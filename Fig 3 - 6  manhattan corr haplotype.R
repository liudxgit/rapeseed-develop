#!/usr/bin/env Rscript

## TWAS manhattan

library(ggplot2)
setwd('D:\\刘东旭\\种子发育转录组\\cageminer')
gwasResults_C18_1 <- read.table("gene_position_annotation.txt", sep='\t',header = T)
#SNP	CHR	BP	P	gene	annotate
#A01_4107093	A01	4107093	0.042352522		
#C09_60267033	C09	60267033	-1.677214362	TT19	yes

# chr length
library(tidyverse)
chr_len <- gwasResults_C18_1 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP))
chr_len
#chr start position
chr_pos <- chr_len %>%
  mutate(total=cumsum(chr_len) - chr_len) %>%
  select(-chr_len)
# SNP position
Snp_pos <- chr_pos %>%
  left_join(gwasResults_C18_1, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)
Snp_pos
#X-axis label position 
X_axis <-  Snp_pos %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
X_axis$center

X_label <- c(paste("A0", 1:9, sep = ""),"A10",paste("C0",1:9,sep = ""))
X_label
library(qqman)
library(ggrepel)

ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #geom_point(aes(color=as.factor(CHR)), alpha=1, size=1) +
  scale_color_manual(values = rep(c("#5a799d", "#87b2d8"), 19 )) +
  scale_x_continuous(label = X_label, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  expand_limits(y=1) + 
  geom_hline(yintercept = c(1.3,-1.3), 
             color = c('grey'),size = 0.3,
             linetype = c("dotted")) +
  geom_hline(yintercept = c(0), 
             color = c('black'),size = 0.3) + 
  xlab("Chromosome") +
  ylab("TWAS (20DAF)\n -log10(FDR)") + 
  geom_point(data=subset(Snp_pos, annotate=="yes"), color="orange", size=0.8) + 
  #geom_point(data=subset(Snp_pos, annotate_a=="yes"), color="#B172B5", size=2) +
  geom_text_repel(aes(label=gene),size=1.8,color="black",direction = "both",
                  segment.colour = "black",max.overlaps = 150) +
  theme_bw() +
  theme(legend.position="none",
    panel.border = element_blank(),
    axis.line = element_line(colour = "black",size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 6,color = "black"),
    axis.title = element_text(size=6),
    axis.text.x = element_text(angle = 45,vjust=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())
ggsave("gene_position_annotation.pdf",width = 4.2,height = 1.6)

# only point
test_oilcontent_20_p <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  geom_point(aes(color=as.factor(CHR)), alpha=1, size=0.5) +
  scale_color_manual(values = rep(c("#5a799d", "#87b2d8"), 19 )) +
  scale_x_continuous(label = X_label, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  expand_limits(y=1) + 
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(test_oilcontent_20_p,filename = "gene_position.png",width = 3.8,height = 1.5)

## correlation of gene expre and phenotype
library(ggplot2)
library(ggpubr)
setwd("your pathdir")
dat <- read.table("seed_specific_expre_phenotype.txt", sep = '\t',header = TRUE)
#pheno	TWAS	expre	gene
#41.77840383	oilcontent	0.223907	BnaA01G0021600ZS
#43.24108343	oilcontent	0.107078	BnaA01G0021600ZS

ggplot(dat, aes(pheno, expre)) +
  geom_point(size=0.5,alpha = 0.7) +
  #scale_fill_manual(values = data$col) + 
  facet_wrap(~gene, ncol = 4, scale = 'free') +
  geom_smooth(method = 'lm') +
  labs(x = '', y = 'Expre') + expre +
  stat_cor(data=dat, method = "pearson") +
  stat_ellipse(level = 0.9,geom = "polygon",alpha=0.2)

ggsave("candidate_gene_expre_phenotype.pdf",width = 6.4, height = 3.2)

expre <- theme_bw() + theme(axis.title.y = element_text(size = 7),
                            axis.title.x = element_blank(),
                            plot.title = element_text(hjust = 0.5,size = 10),
                            axis.text = element_text(size = 7,colour = "black"),
                            axis.text.x = element_text(size=7),
                            panel.grid.major = element_line(),
                            panel.grid.minor = element_line(),
                            legend.position = 'none')

## haplotype olicontent
Haplotype <- theme_classic() + theme(axis.title.y = element_text(size = 12),
                                     axis.title.x = element_text(size = 12),
                                     plot.title = element_text(hjust = 0.5,size = 12),
                                     axis.text = element_text(size = 12,colour = "black"),
                                     axis.text.x = element_text(size=10),
                                     #panel.grid.major = element_line(),
                                     #panel.grid.minor = element_line(),
                                     legend.position = 'none')
library(ggpubr)
setwd("your pathdir")
oilcontent_Hap <- read.table("oilcontent_Hap_pheno_toR.txt", sep = '\t',header = TRUE)
oilcontent_Hap_20DAF <- read.table("oilcontent_Hap_20DAF_toR.txt", sep = '\t',header = TRUE)
oilcontent_Hap_40DAF <- read.table("oilcontent_Hap_40DAF_toR.txt", sep = '\t',header = TRUE)
pique_Hap <- read.table("BnaC08G0351900ZS_pique_Hap2_toR.txt", sep = '\t',header = TRUE)

#BnaC08G0351900ZS TT5 267 91 pheno
BnaC08G0351900ZS_Hap <- oilcontent_Hap[oilcontent_Hap$Gene=="BnaC08G0351900ZS",]
BnaC08G0351900ZS_Hap$Haplotype <- factor(BnaC08G0351900ZS_Hap$Haplotype,
                                         levels = c('TCTCT','CGCTT'),ordered = TRUE)
BnaC08G0351900ZS_Hap_pdf <- ggboxplot(BnaC08G0351900ZS_Hap,x='Haplotype',y='Oilcontent',
                                      fill="#a4a4a4",palette = 'jco',outlier.shape=NA) +
  stat_compare_means(method = 't.test',label='p.format') +
  labs(y = 'Oil content (%)') + Haplotype