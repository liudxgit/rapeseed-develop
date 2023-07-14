#!/usr/bin/env Rscript

## Fugure S2 Correlation between biological replications
setwd("your pathdir")
TPM<- read.table("RNAseq_TPM.txt",header = TRUE,sep = '\t')
library(ggpubr)

par(mfrow=c(3,3),mar = c(5, 6, 4, 4))
for (name in c("C6","C15","C21")){
  Rep1 <- paste(c(name,"_Rep1"),sep = "",collapse="")
  Rep2 <- paste(c(name,"_Rep2"),sep = "",collapse="")
  Rep3 <- paste(c(name,"_Rep3"),sep = "",collapse="")
  plot(TPM[,Rep1][1:20000]~TPM[,Rep2][1:20000],xlab = Rep2, ylab = Rep1,
       pch=16,cex.axis=2,cex.lab=2)
  abline(lm(TPM[,Rep1][1:20000]~TPM[,Rep2][1:20000]))
  text_y_12 <- max(TPM[,Rep1][1:20000]) - max(TPM[,Rep1][1:20000]/10)
  text_x_12 <- max(TPM[,Rep2][1:20000])/5
  r2_12 <- cor(TPM[,Rep1][1:20000],TPM[,Rep2][1:20000],method = "pearson")
  r2_12 <- round(r2_12,2)
  write_str_12 <- paste(c("R = ",r2_12),sep = "",collapse="")
  text(text_x_12,text_y_12,write_str_12,cex=2)
  
  plot(TPM[,Rep1][1:20000]~TPM[,Rep3][1:20000],xlab = Rep3, ylab = Rep1,
       pch=16,cex.axis=2,cex.lab=2)
  abline(lm(TPM[,Rep1][1:20000]~TPM[,Rep3][1:20000]))
  text_y_13 <- max(TPM[,Rep1][1:20000]) - max(TPM[,Rep1][1:20000]/10)
  text_x_13 <- max(TPM[,Rep3][1:20000])/5
  r2_13 <- cor(TPM[,Rep1][1:20000],TPM[,Rep3][1:20000],method = "pearson")
  r2_13 <- round(r2_13,2)
  write_str_13 <- paste(c("R = ",r2_13),sep = "",collapse="")
  text(text_x_13,text_y_13,write_str_13,cex=2)
  
  plot(TPM[,Rep2][1:20000]~TPM[,Rep3][1:20000],xlab = Rep3, ylab = Rep2,
       pch=16,cex.axis=2,cex.lab=2)
  abline(lm(TPM[,Rep2][1:20000]~TPM[,Rep3][1:20000]))
  text_y_23 <- max(TPM[,Rep2][1:20000]) - max(TPM[,Rep2][1:20000]/10)
  text_x_23 <- max(TPM[,Rep3][1:20000])/5
  r2_23 <- cor(TPM[,Rep2][1:20000],TPM[,Rep3][1:20000],method = "pearson")
  r2_23 <- round(r2_23,2)
  write_str_23 <- paste(c("R = ",r2_23),sep = "",collapse="")
  text(text_x_23,text_y_23,write_str_23,cex=2)
}


## Figure S4
library(phylobase)
library(ggtree)
library(tibble)
library(Biostrings)
library(treeio)
library(ggplot2)
library(ggimage)
library(phyloseq)
library(ggridges)
library(dplyr)
library(tidyr)
setwd("your pathdir")
athtree <- read.tree("AT4G28800_mafft.fa.new.nwk")
tree1 <- as(athtree, "phylo4")
expression=read.table("AT4G28800_expression_tissue.txt",header=T)
#colnames(expression) <- gsub("X","",colnames(expression))
rownames(expression) <- expression$gene
#expression
expression$gene <- NULL
expression[expression>=5000] <-  5000.0
log2_normalization<-function(x){
  return(log(x+1,2))}
expression <- log2_normalization(expression)
expression
tree2 <- phylo4d(tree1, expression, missing.data="warn")
#tree2
athd <- tipData(tree2)
#athd
athd$tip <- rownames(athd)
athdd <- gather(athd, feature, value, -tip)
athcat <- seq(ncol(athd))
names(athcat) <- names(athd)
athdd$cat <- athcat[athdd$feature]
athd2 <- select(athdd, -value) %>% filter(tip == 'AT4G28800')

p <- ggtree(tree2,branch.length="none")  + 
  geom_tiplab(offset = 5.5, size=2,hjust = 0.9) +
  geom_tippoint(size = 2, alpha = .5) +
  xlim_tree(c(0, 1.2))
p

p2 <- facet_plot(p, "Expression",
                 data=athdd, geom=geom_tile,
                 mapping=aes(x = cat, fill = value))
p2
p3 <- facet_plot(p2,"Expression", data=athd2, geom=geom_text,
             mapping=aes(x = cat, y = 0, label = feature),
             size =2 ,angle = 75) +
  scale_fill_viridis_c(na.value = "white",
                       begin = 0.3, end = 1,
                       option = "D") +
  theme(legend.position="right",
        strip.text.x = element_text(size = 8),
        legend.title = element_text(size=8), 
        legend.text=element_text(size = 8)) +
  labs(fill=expression('log'[2]*'(TPM)'))
p3
facet_widths(p3, widths = c(30, 70))

ggsave("AT4G28800_BHLH56.pdf",width = 4,height = 2.2)



## Figure S6  module go analysis
setwd("your pathdir")
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
module_f <- read.table("seed_down_more12_barplot.txt",header = TRUE,sep = '\t')

shorten_names <- function(x, n_word=4, n_char=30){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 30))
  {
    if (nchar(x) > 30) x <- substr(x, 1, 30)
    x <- paste(paste(strsplit(x, " ")
                     [[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
labels=(sapply(levels(module_f$term)[as.numeric(module_f$term)], shorten_names))

term_order=factor(as.integer(rownames(module_f)),labels=labels)
term_order=rev(term_order)

ggplot(module_f, aes(x=term_order, y=number, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=term_order) +
  xlab("GO term") +
  labs(title = "Downregulation") +
  theme(title = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16,colour = "black"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.position = "right")
ggsave("downregulation_go.pdf",width = 9,height = 9)

## Figure S9 gene structure
library(GenomicFeatures)
library(ggplot2)
setwd('your pathdir')
GFF_file <- "zs11.v0.gff3"
txdb <- makeTxDbFromGFF(GFF_file) 
genes_df <- as.data.frame(genes(txdb))
exons_df <- as.data.frame(exons(txdb))
A06_exons=subset(exons_df,seqnames=='scaffoldA06'&start>=3572905&end<=3574580)
A06_exons
ggplot(A06_exons) + 
  geom_segment(data=A06_exons, aes(x=start,xend=end,y=0,yend=0)) +
  geom_rect(aes(xmin=start, xmax=end,ymin=-0.1,ymax=0.1),
            colour="#000000", fill="#000000") +
  ylim(c(-1,1)) + theme_syntenty +
  geom_vline(xintercept = c(3573716,3573742,3573745,3573761,3573895,3573927,3573955))

theme_syntenty <- theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


## Fugure S17 seed-specific gene expression heatmap
library(pheatmap)
setwd("your pathdir")

tissue_all <- read.table("seed_specific_TAU_expre.txt",header = TRUE,sep = '\t')
rownames(tissue_all) <- tissue_all$gene
tissue_all$gene <- NULL
#tissue_all[tissue_all>=5000] <-  5000

# max normalization
max_normalization<-function(x){
 return((x-min(x))/(max(x)-min(x)))}
for (gene in rownames(tissue_all)) {
  tissue_all[gene,] <- max_normalization(tissue_all[gene,])
}                         

pheatmap(tissue_all,cluster_rows = F,cluster_cols = F,show_rownames = F,
         show_colnames = F)
pdf(file = "seed_specific_TAU_expre.pdf",width = 5,height = 8)


## Figure S6 core module mapman
library(ggplot2)
setwd("your pathdir")
mapman <- read.table("mapman_core_module.txt",header = TRUE,sep = '\t')
module1 <- mapman[mapman$module=="module5",]
#module1
module1$Pathway <- factor(module1$Pathway,levels = module1$Pathway)
ggplot(module1,aes(Pvalue,Pathway))+
  geom_point(aes(size=Number,color=-1*log10(FDR)))+
  scale_colour_gradient(low="#89BBB0",high = "#C5AC79")+
  labs(color=expression(-log[10](FDR)),size="Gene",
       x="Pvalue") + module_core
ggsave("module5.pdf",width = 7.4, height = 4)
module_core <- theme_classic() + 
  theme(axis.title.y = element_blank(),
         axis.title.x = element_text(size = 16),
         axis.text = element_text(size = 14,colour = "black"),
         axis.text.x = element_text(size=12,angle=45,hjust = 1),
         legend.text = element_text(size = 14),
         legend.title = element_text(size=16))


