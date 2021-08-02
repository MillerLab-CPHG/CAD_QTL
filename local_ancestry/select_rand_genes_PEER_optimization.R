######## The purpose of this R script is to select a random  number of genes from one chromosome to run through the eQTL analysis pipeline and 
### determine the optimal number of PEER factors to use in the genome-wide analysis. The expectation is that the number of significant eGenes will 
### increase as PEER factors are added to the model, and then either plateau or start to decrease. The model with the number of PEER factors 
### corresponding to the max number of "discoveries" is the number used for final analyses. File names and paths will need to be changed to match
### your datasets.

library(data.table)
library(dplyr)

###Read in file with list of all genes on chr17 from Ensembl Biomart database:
m <- fread("mart_export.txt", select=c("Gene_stable_ID_version"))
setkey(m, Gene_stable_ID_version)

###Read in file of all genes with an average # of raw reads in TPM file >6:
g <- fread("random_genes_for_select_UVA_coronary.lst", header=F)
setkey(g, V1)

#Assign random numbers to all chr17 genes, then reorder
rand <- sample(1:50000, 2792, replace=F)
m$rand <- rand
setorder(m, rand)

#Select random middle 60% of genes using ENSG names (2792 total) to narrow down from:
nrow(m)
2792-(2792/5) #filter bottom 20% of genes (randomly)
2792-(4*2792/5) # filter top 20% of genes (randomly)
m2 <- m[560:2234, ] # select middle 60% of genes (randomly), 1674 total to filter based on minimum expression
setkey(m2, Gene_stable_ID_version)

#Combine to select only genes that are both in the random chr17 subset and expressed at detectable leves in our study population:
both <- m2[g]
setorder(both, rand)

#Unique genes in the random dataset and meet expression criteria
new <- unique(both[!is.na(rand), c(1)])
nrow(new)

write.table(new, "Random_genes_for_PEER_optimization.lst", row.names=F, quote=F, sep="\t")

###Read in full-genome expression files to select genes of interest and remake bed file for PEER factor # optimization.
tpm <- fread("TPM_UVAcoronary_May.gct")
raw <- fread("Raw_reads_UVAcoronary_May.txt")

t1 <- tpm[Name %in% new$Gene_stable_ID_version, ]
head(t1[, c(1:5)])
write.table(t1, "UVA_coronary_PEERopt_TPM.txt", sep="\t", row.names=F, quote=F)

r1 <- raw[Name %in% new$Gene_stable_ID_version, ]
head(r1[, c(1:5)])
write.table(r1, "UVA_coronary_PEERopt_raw.txt", sep="\t", row.names=F, quote=F)

