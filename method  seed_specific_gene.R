#!/usr/bin/env Rscript

## cageminer method

library(cageminer)
library(SummarizedExperiment)
setwd("your workdir")
twas_sig <- read.table("TWAS_sig_gene.txt",header = TRUE,sep = '\t')

## demo
#Gene
#BnaA04G0015800ZS
#BnaC08G0386300ZS

twas_sig$Gene

guide_acyl <- read.table("acyl_lipid_bna_gene.txt",header = TRUE,sep = '\t')

## demo
#ID
#Bnascaffold0025G0003500ZS
#Bnascaffold0025G0003600ZS

guide_acyl$ID
guide_pheny <- read.table("pheny_bna_gene.txt",header = TRUE,sep = '\t')
guide_pheny$ID

expre_module <- read.table("gene_expre.txt",header = TRUE,row.names = 1,sep = '\t')

expre_module
expre_module2 <- SummarizedExperiment(expre_module)
expre_module2

# Infer GCN
sft <- BioNERO::SFT_fit(expre_module2, net_type = "signed", cor_method = "pearson") ## similar to WGCNA
sft$power

gcn <- BioNERO::exp2gcn(expre_module2, net_type = "signed", cor_method = "pearson",
                        module_merging_threshold = 0.8, SFTpower = sft$power)

# Apply step 2  pheny related candidate genes
candidates_acyl <- mine_step2(expre_module2, gcn = gcn, guides = guide_pheny$ID,
                          twas_sig$Gene)


# Apply step 2  acyl related candidate genes
candidates_pheny <- mine_step2(expre_module2, gcn = gcn, guides = guide_pheny$ID,
                          twas_sig$Gene)








