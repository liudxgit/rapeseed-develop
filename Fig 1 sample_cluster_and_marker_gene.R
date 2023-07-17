#!/usr/bin/env Rscript

## 26 sample cluster
library(pvclust)
setwd("your pathdir")
seed_tree <- read.table("seed_RNAseq_unrepeat_log1.txt",header = TRUE,sep = '\t') ##expre profile

## demo
#gene	B1	B2	B3	B4	B5	B6	B7	B8	B9	B10	B11	B12	B13	B14	B15	B16	B17	B18	B19	B20	B21	B22	B23	B24	B25	B26
#BnaA01G0001300ZS	0.0247253333333	0.0298143333333	0.0	0.026463	0.0795146666667	0.0	0.024136	0.0	0.0357723333333	0.0938716666667	0.0286716666667	0.0	0.113267333333	0.0555166666667	0.080612	0.087953	0.530192666667	1.11001833333	0.876177666667	0.957424333333	1.85960533333	0.0	0.0	0.0	0.0	0.115738666667
#BnaA01G0003100ZS	21.6110366667	1.41194	2.05249233333	2.55043766667	1.48534866667	1.10822033333	0.117202333333	0.0	0.0	0.0	0.0150486666667	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0

rownames(seed_tree) <- seed_tree$gene
seed_tree$gene <- NULL
fit_tree_seed <- pvclust(seed_tree,method.hclust = "ward.D2",
        method.dist = "euclidean",nboot = 100) ## cluster

library(ggplot2)
require(ggtree)

ggtree(fit_tree_seed$hclust,layout="circular")  + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab() ## using the tree from pvclust

seed1 <- ggtree(fit_tree_seed$hclust) + geom_tiplab(offset = 8, hjust = 0.2)

seed2 <- flip(tree_view = seed1, 38, 39)
seed3 <- flip(tree_view = seed2, 31, 30)
seed4 <- flip(tree_view = seed3, 41, 40)
seed5 <- flip(tree_view = seed4, 33, 32)
seed6 <- flip(tree_view = seed5, 34, 35)
seed6
class(seed6)

seed7 <- rotate(seed6,45)
seed7
seed8 <- rotate(seed7, 48)
seed9 <- rotate(seed8, 42)
seed10 <- rotate(seed9, 51)
seed11 <- rotate(seed10, 46)
seed12 <- rotate(seed11, 40)
seed13 <- rotate(seed12, 49)
seed14 <- rotate(seed13, 32)
seed15 <- rotate(seed14, 44)
seed16 <- rotate(seed15, 47)
seed16
ggsave("seed_specific_tree_final.pdf",width = 5,height = 10)


## marker gene
library(ggplot2)
library(gcookbook)
setwd("your pathdir")
marker <- read.table("marker_gene_expre.txt",header = TRUE,sep = '\t')

## demo
#gene	B1	B2	B3	B4	B5	B6	B7	B8	B9	B10	B11	B12	B13	B14	B15	B16	B17	B18	B19	B20	B21	B22	B23	B24	B25	B26
#BnaA01G0001300ZS	0.0247253333333	0.0298143333333	0.0	0.026463	0.0795146666667	0.0	0.024136	0.0	0.0357723333333	0.0938716666667	0.0286716666667	0.0	0.113267333333	0.0555166666667	0.080612	0.087953	0.530192666667	1.11001833333	0.876177666667	0.957424333333	1.85960533333	0.0	0.0	0.0	0.0	0.115738666667
#BnaA01G0003100ZS	21.6110366667	1.41194	2.05249233333	2.55043766667	1.48534866667	1.10822033333	0.117202333333	0.0	0.0	0.0	0.0150486666667	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0

max_value11 <- max(marker$BnaA07G0095400ZS)
gene11 <- ggplot(marker, aes(x=gene, y=BnaA07G0095400ZS)) + 
  geom_rect(aes(xmin=14, xmax=25, ymin=-Inf, ymax=max_value11),fill='#f5d249',alpha = .5) +
  #geom_rect(aes(xmin=33, xmax=38, ymin=-Inf, ymax=max_value1),fill='#98CA3D',alpha = .5) +
  geom_line(size=0.1) + 
  geom_point(size=0.3, shape=20) + 
  scale_x_continuous(breaks=seq(14, 64, 2)) +
  labs(title="LEC2 BnaA07G0095400ZS") +
  expre

expre <- theme_bw() + theme(axis.title.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.title = element_text(hjust = 0.1,size = 6),
                            axis.text.y = element_text(size=6),
                            axis.text.x = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            axis.line = element_line(colour = "black",size = 0.1),
                            axis.ticks = element_line(size = 0.1))
