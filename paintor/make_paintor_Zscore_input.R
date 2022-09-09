##### This file was used to get all SNPs from GWAS summary statistics to run with Paintor
##### NOTE: didn't include Japanese studies' summary statistics because the LD is so different
##### If somehow all your GWAS sum stats are magically formatted the same, can use a function to do this 
##### across published studies before reducing

library(readr)
library(dplyr)
library(data.table)
library(purrr)

chrom="1"

setwd('/path/to/main/directory')

qtl <- fread(paste0("QTL/folder/QTL_res_chr", chrom, "_date.tsv.gz"))
qtl <- qtl[!is.na(pval.meta)]
qtl <- qtl[, c(1:3, 15)]
setnames(qtl, "stat.meta", "UVAQTL_Z")
setkey(qtl, variant)

snps <- fread(paste0("/scratch/ch2um/chr",chrom,"_mixQTL_nominal_for_paintor.txt.gz"), select=c(1:3))
snps <- unique(snps)
setkey(snps, variant)

a <- snps[qtl, nomatch=0]
vars <- as.vector(unique(a$variant))

b <- fread("/path/to/GWAS/results/study1.tsv.gz", select=c(1,2,6))
setnames(b, c("rsid", "P", "Effect"))
b <- b[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
b[(Effect<0 & Z<0) | (Effect>0 & Z>0), GCST011365_Z:=Z]
b[(Effect>0 & Z<0) | (Effect<0 & Z>0), GCST011365_Z:=(-1)*Z]
b[Effect==0, GCST011365_Z:=0]
b <- b[, c(1,5)]

d <- fread("/path/to/GWAS/results/study2.txt.gz", select=c(1,2,6))
setnames(d, c("rsid", "P", "Effect"))
d <- d[, Z:=qnorm(P)] %>% filter (rsid %in% vars)
d[(Effect<0 & Z<0) | (Effect>0 & Z>0), GCST005194_Z:=Z]
d[(Effect>0 & Z<0) | (Effect<0 & Z>0), GCST005194_Z:=(-1)*Z]
d[Effect==0, GCST005194_Z:=0]
d <- d[, c(1,5)]

e <- fread("/path/to/GWAS/results/study3.txt.gz", select=c(1,14,16))
setnames(e, c("rsid", "P", "Effect"))
e <- e[, Z:=qnorm(P)] %>% filter (rsid %in% vars)
e[(Effect<0 & Z<0) | (Effect>0 & Z>0), HDL_Z:=Z]
e[(Effect>0 & Z<0) | (Effect<0 & Z>0), HDL_Z:=(-1)*Z]
e[Effect==0, HDL_Z:=0]
e <- e[,c(1,5)]

f <- fread("/path/to/GWAS/results/study4.txt.gz", select=c(1,14,16))
setnames(f, c("rsid", "P", "Effect"))
f <- f[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
f[(Effect<0 & Z<0) | (Effect>0 & Z>0), LDL_Z:=Z]
f[(Effect>0 & Z<0) | (Effect<0 & Z>0), LDL_Z:=(-1)*Z]
f[Effect==0, LDL_Z:=0]
f <- f[,c(1,5)]

g <- fread("/path/to/GWAS/results/study5.txt.gv", select=c(1,14,16))
setnames(g, c("rsid", "P", "Effect"))
g <- g[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
g[(Effect<0 & Z<0) | (Effect>0 & Z>0), TC_Z:=Z]
g[(Effect>0 & Z<0) | (Effect<0 & Z>0), TC_Z:=(-1)*Z]
g[Effect==0, TC_Z:=0]
g <- g[,c(1,5)]

h <- fread("/path/to/GWAS/results/study6.txt.gz", select=c(1,14,16))
setnames(h, c("rsid", "P", "Effect"))
h <- h[, Z:=qnorm(P)] %>% filter (rsid %in% vars)
h[(Effect<0 & Z<0) | (Effect>0 & Z>0), logTG_Z:=Z]
h[(Effect>0 & Z<0) | (Effect<0 & Z>0), logTG_Z:=(-1)*Z]
h[Effect==0, logTG_Z:=0]
h <- h[,c(1,5)]

i <- fread("/path/to/GWAS/results/study7.txt.gz", select=c(1,2,6))
setnames(i, c("rsid", "P", "Effect"))
i <- i[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
i[(Effect<0 & Z<0) | (Effect>0 & Z>0), PP_Z:=Z]
i[(Effect>0 & Z<0) | (Effect<0 & Z>0), PP_Z:=(-1)*Z]
i[Effect==0, PP_Z:=0]
i <- i[,c(1,5)]

j <- fread("/path/to/GWAS/results/study8.txt.gz", select=c(1,2,6))
setnames(j, c("rsid", "P", "Effect"))
j <- j[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
j[(Effect<0 & Z<0) | (Effect>0 & Z>0), SBP_Z:=Z]
j[(Effect>0 & Z<0) | (Effect<0 & Z>0), SBP_Z:=(-1)*Z]
j[Effect==0, SBP_Z:=0]
j <- j[,c(1,5)]

k <- fread("/path/to/GWAS/results/study9.txt.gz", select=c(1,2,6))
setnames(k, c("rsid", "P", "Effect"))
k <- k[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
k[(Effect<0 & Z<0) | (Effect>0 & Z>0), DBP_Z:=Z]
k[(Effect>0 & Z<0) | (Effect<0 & Z>0), DBP_Z:=(-1)*Z]
k[Effect==0, DBP_Z:=0]
k <- k[,c(1,5)]

n <- fread("/path/to/GWAS/results/study10.txt.gz", select=c(20,19,8))
setnames(n, c("rsid", "P", "Effect"))
n <- n[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
n[(Effect<0 & Z<0) | (Effect>0 & Z>0), IMT_Z:=Z]
n[(Effect>0 & Z<0) | (Effect<0 & Z>0), IMT_Z:=(-1)*Z]
n[Effect==0, IMT_Z:=0]
n <- n[,c(1,5)]

o <- fread("/path/to/GWAS/results/study11.txt.gz", select=c(18,21,6))
setnames(o, c("rsid", "P", "Effect"))
o <- o[, Z:=qnorm(P)] %>% filter(rsid %in% vars)
o[(Effect<0 & Z<0) | (Effect>0 & Z>0), Plaque_Z:=Z]
o[(Effect>0 & Z<0) | (Effect<0 & Z>0), Plaque_Z:=(-1)*Z]
o[Effect==0, Plaque_Z:=0]
o <- o[,c(1,5)]

studies <- list(a,b,d,e,f,g,h,i,j,k,l,n,o) #c,m ended up not being used, hence excluded from this script
setnames(a, "variant", "rsid")
z <- reduce(studies, full_join, by="rsid")

setorder(z, gene, pos)
setnames(z, c("rsid", "chr", "pos"), c("RSID", "CHR", "POS"))
z <- z[, c(2,3,1,6:18,4,5)]

write_delim(z, file=gzfile(paste0("/path/name/Date_eqtl_gwas_forpaintor_chr", chrom,".txt.gz")), escape="none")
