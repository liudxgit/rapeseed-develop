#!/usr/bin/env Rscript

## 
setwd("yourfilename")
htseq_counts <- read.table("gene_expre.txt",header = TRUE,sep = '\t')  ## input dataset
dim(htseq_counts)

rownames(htseq_counts) <- htseq_counts$gene
htseq_counts$gene <- NULL
thr<-5000
htseq_counts[htseq_counts>thr]=thr
## 
#htseq_counts <- as.data.frame(htseq_counts)
htseq_counts$mean <- rowMeans(htseq_counts)
htseq_counts_sort <- htseq_counts[order(htseq_counts$mean,decreasing=TRUE),]

htseq_counts_sort_cv <- htseq_counts_sort[1:50000,] ## retain top 5w most expre genes 
htseq_counts_sort_cv$sd <- apply(htseq_counts_sort_cv,1,FUN=sd) ##every line sd
htseq_counts_sort_cv$cv <- htseq_counts_sort_cv$sd/htseq_counts_sort_cv$mean ## cv every line
htseq_counts_sort_cv <- htseq_counts_sort_cv[order(htseq_counts_sort_cv$cv,decreasing=TRUE),] ## order by CV
htseq_counts_sort_cv <- htseq_counts_sort_cv[1:15000,] ## retain top 15,000 genes with high CV
htseq_counts_sort_cv$mean <- NULL ## remove column
htseq_counts_sort_cv$sd <- NULL
htseq_counts_sort_cv$cv <- NULL  ## htseq_counts_sort_cv

htseq_counts_sort_cv <- htseq_counts

library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

library(flashClust)

PYExpr[1:10,1:10]
PYExpr=as.data.frame(t(htseq_counts_sort_cv))
dim(PYExpr)
# data filter
gsg=goodSamplesGenes(PYExpr,verbose=3)
PYExpr=PYExpr[gsg$goodSamples,gsg$goodGenes]
#dim(PYExpr)
# average expression level
meanExpressionByArray=apply(PYExpr,1,mean,na.rm=T)
sizeGrWindow(20,10) ## Rplots1.pdf
par(las=3)
par(mar=c(7,5,5,5))
barplot(meanExpressionByArray,ylab="Mean expression", main="Mean expression across samples",cex.names=0.7)
dev.off()

# sample cluster
sampleTree=flashClust(dist(PYExpr),method="average")
sizeGrWindow(12,9) ## Rplots2.pdf
par(cex=0.8)
par(mar=c(2,6,2,2))
plot(sampleTree,main="Sample clustering",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
dev.off()

# select the power£»
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(PYExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5) # Rplots4.pdf
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower=12  ## according the power
adjacency=adjacency(PYExpr,power=softPower)

# check the network whether it is fit scale-free
ADJ1=abs(cor(PYExpr,use="p"))^12  ## power
k=softConnectivity(datE=PYExpr,power=12)
sizeGrWindow(10,5)
par(mfrow=c(1,2))  ## Rplots5.pdf
hist(k)
scaleFreePlot(k,main="Check scale free topology\n")
dev.off()

#calculate TOM matrix¡¤
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM

#TOM cluster
geneTree=flashClust(as.dist(dissTOM),method="average")
sizeGrWindow(12,9) ##Rplots6.pdf
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",labels=FALSE,hang=0.04)
dev.off()
minModuleSize=30
dynamicMods=cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)
table(dynamicMods)

dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(12,6)  ##Rplots7.pdf
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", main = "Gene dendrogram and module colors",dendroLabels = FALSE, hang = 0.05,addGuide = TRUE, guideHang = 0.05)
dev.off()

#merge the similar module by the expression
MEList=moduleEigengenes(PYExpr,color=dynamicColors)
MEs=MEList$eigengenes
MEDiss=1-cor(MEs)
METree=flashClust(as.dist(MEDiss),method="average")
sizeGrWindow(12,6) ##Rplots8.pdf
plot(METree,main="Clustering of module eigengenes",xlab="",sub="",cex=0.6)
MEDissThres=0.7
dev.off()
merge=mergeCloseModules(PYExpr,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs=merge$newMEs
sizeGrWindow(12,6)  ##Rplots9.pdf
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic") ,dendroLabels = FALSE, hang = 0.05,addGuide = TRUE, guideHang = 0.05)
dev.off()

# rename the module
moduleColors=mergedColors
# use number to represent module
colorOrder=c("grey",standardColors(50))
moduleLabels=match(moduleColors,colorOrder)-1
MEs=mergedMEs
save(MEs,moduleLabels,moduleColors,geneTree,file="BGQ-Long.RData")
write.table(table(moduleColors),file="color_count_308.txt")

#drawing the barplot of module number
uni=unique(moduleColors)
uni=uni[order(uni)]
hehe=data.frame(table(moduleColors),uni)
hehe=hehe[rev(order(hehe$Freq)),]
sizeGrWindow(12,6)  ##Rplots10.pdf
par(mar=c(7,3,2,3))
barplot(hehe[,2],col=hehe[,3],names.arg=hehe[,1],las=3,ylim=c(0,3000))
dev.off()

#select 2k genes and drawing TOM heatmap
nSelect=2000
nGenes=ncol(PYExpr)
nSamples=nrow(PYExpr)
set.seed(10)
select=sample(nGenes,size=nSelect)
selectTOM=dissTOM[select,select]
selectTree=flashClust(as.dist(selectTOM),method="average")
selectColors=moduleColors[select]
sizeGrWindow(12,12)  ## Rplots11.pdf
plotDiss=selectTOM^7
diag(plotDiss)=NA
TOMplot(plotDiss,selectTree,selectColors,main="Network heatmap plot, selected genes")
dev.off()

#export network
for(i in rownames(table(mergedColors))){
  probes <- colnames(PYExpr)
  inModule <- is.finite(match(mergedColors,i))
  modProbes <- probes[inModule]
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modProbes,modProbes)
  cyt <- exportNetworkToCytoscape(modTOM,
                                  edgeFile = paste("edges-",paste(i,collapse = "-"),".sif",sep=""),
                                  nodeFile = paste("nodes-",paste(i,collapse = "-"),".txt",sep=""),
                                  weighted=TRUE,
                                  threshold=0.02,
                                  nodeNames=modProbes,
                                  nodeAttr=mergedColors[inModule])
}