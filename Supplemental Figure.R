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

## demo
#(((BnaA01G0082500ZS:0.11183,BnaC01G0099800ZS:0.10401)0.959:0.04382,((BnaA01G0082600ZS:0.06024,BnaC01G0100200ZS:0.04358)0.982:0.06042,BnaC01G0099700ZS:0.13530):0.03152)0.914:0.19809,AT4G28800:0.19809);

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

## demo
#number	term	type
#32	"photosynthesis, light harvesting"	BP
#23	glycerol ether metabolic process	BP

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

## demo
#gene	A1	A2	A3	A4	A5	A6	B1	B2	B3	B4	B5	B6	B7	B8	B9	B10	B11	B12	B13	B14	B15	B16	B17	B18	B19	B20	B21	B22	B23	B24	B25	B26	C1	C2	C3	C4	C5	C6	C7	C8	C9	C10	C11	C12	C13	C14	C15	C16	C17	C18	C19	C20	C21	C22	C23	C24	D1	D2	F1	F2	F3	F4	G1	G2	G3	G4	G5	G6	G7	G8	G9	G10	G11	G12	G13	G14	G15	G16	G17	G18	G19	G20	G21	G22	G23	J1	J2	J3	106	111	CK1
#BnaA01G0003100ZS	0	0.39	1.07	2	3.37	10.81	21.61	1.41	2.05	2.55	1.49	1.11	0.12	0	0	0	0.02	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.09	0.02	0.02	0	0.01	0	0.02	0.04	0.01	0	0	0.05	0	0.01	0.02	0.03	0	0	0	0	0	0	0	0.02	0	0	0	0	0.01	0.52	0.45	0.28	0.3	0.18	0.15	0.47	0.41	0.5	0.43	0.32	0.29	0.25	0.19	0.17	0.15	0.32	0.34	0.23	0.23	0.2	0.19	0.22	0	0	0	0	0	0.07
#BnaA01G0004800ZS	0	0.03	0.04	0.1	0.21	0.2	1.69	1.18	1.05	0.82	1	0.59	0.75	0.75	0.12	0.37	0.49	0.26	0.25	0.49	0.6	0.39	0.63	0.59	0.69	0.77	0.58	0.59	0.39	0.22	0.61	0.71	0	0	0	0	0.04	0	0	0	0	0	0	0.03	0	0	0	0.03	0	0	0	0	0	0	0.17	0	0.8	0.09	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.03	0	0	0	0	0	0	0	0.09	0.03	0	0.04	0	0.43	0

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

## demo
#Pathway	Number	Pvalue	FDR	module	Number2
#Cell cycle organisation	136	3.07E-13	5.08E-10	module1	204 
#mitosis and meiosis	42	1.66E-08	9.17E-06	module1	63 


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


