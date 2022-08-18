##### This file was used to get all SNPs from GWAS summary statistics to run with Paintor
##### NOTE: didn't include Japanese summary statistics because the LD is so different

library(readr)
library(dplyr)
library(data.table)
library(purrr)

chrom="1"

setwd('/project/cphg-millerlab/chani/QTL_res')

uva <- fread(paste0("mixQTL/mixQTL_res_chr", chrom, "_Aug2_all.tsv.gz"))
uva <- uva[!is.na(pval.meta)]
uva <- uva[, c(1:3, 15)]
setnames(uva, "stat.meta", "UVAQTL_Z")
setkey(uva, variant)

snps <- fread(paste0("/project/cphg-millerlab/chani/QTL_input/mixqtl/Aug11_converted_chr",chrom,"_variant_annotation.txt.gz"), select=c(1:3))
setkey(snps, variant)

a <- snps[uva, nomatch=0]
vars <- as.vector(unique(a$variant))

b <- fread("../external_datasets/MI_GCST011365_buildGRCh37_forcoloc.tsv.gz", select=c(1,2,6))
setnames(b, c("rsid", "P", "Effect"))
b[, Z:=qnorm(P)]
b[(Effect<0 & Z<0) | (Effect>0 & Z>0), GCST011365_Z:=Z]
b[(Effect>0 & Z<0) | (Effect<0 & Z>0), GCST011365_Z:=(-1)*Z]
b <- b[rsid %in% vars, c(1,5)]

c <- fread("../external_datasets/GCST011364_buildGRCh37_forcoloc.tsv.gz", select=c(1,2,6))
setnames(c, c("rsid", "P", "Effect"))
c[, Z:=qnorm(P)]
c[(Effect<0 & Z<0) | (Effect>0 & Z>0), GCST011364_Z:=Z]
c[(Effect>0 & Z<0) | (Effect<0 & Z>0), GCST011364_Z:=(-1)*Z]
c <- c[rsid %in% vars, c(1,5)]

d <- fread("../external_datasets/GCST005194_forcoloc.txt.gz", select=c(1,2,6))
setnames(d, c("rsid", "P", "Effect"))
d[, Z:=qnorm(P)]
d[(Effect<0 & Z<0) | (Effect>0 & Z>0), GCST005194_Z:=Z]
d[(Effect>0 & Z<0) | (Effect<0 & Z>0), GCST005194_Z:=(-1)*Z]
d <- d[rsid %in% vars, c(1,5)]

e <- fread("/project/cphg-millerlab/reference_datasets/Lipid_GWAS_PMID_34887591/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz", select=c(1,14,16))
setnames(e, c("rsid", "P", "Effect"))
e[, Z:=qnorm(P)]
e[(Effect<0 & Z<0) | (Effect>0 & Z>0), HDL_Z:=Z]
e[(Effect>0 & Z<0) | (Effect<0 & Z>0), HDL_Z:=(-1)*Z]
e <- e[rsid %in% vars,c(1,5)]

f <- fread("/project/cphg-millerlab/reference_datasets/Lipid_GWAS_PMID_34887591/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz", select=c(1,14,16))
setnames(f, c("rsid", "P", "Effect"))
f[, Z:=qnorm(P)]
f[(Effect<0 & Z<0) | (Effect>0 & Z>0), LDL_Z:=Z]
f[(Effect>0 & Z<0) | (Effect<0 & Z>0), LDL_Z:=(-1)*Z]
f <- f[rsid %in% vars,c(1,5)]

g <- fread("/project/cphg-millerlab/reference_datasets/Lipid_GWAS_PMID_34887591/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_TC_INV_ALL_with_N_1.gz", select=c(1,14,16))
setnames(g, c("rsid", "P", "Effect"))
g[, Z:=qnorm(P)]
g[(Effect<0 & Z<0) | (Effect>0 & Z>0), TC_Z:=Z]
g[(Effect>0 & Z<0) | (Effect<0 & Z>0), TC_Z:=(-1)*Z]
g <- g[rsid %in% vars,c(1,5)]

h <- fread("/project/cphg-millerlab/reference_datasets/Lipid_GWAS_PMID_34887591/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz", select=c(1,14,16))
setnames(h, c("rsid", "P", "Effect"))
h[, Z:=qnorm(P)]
h[(Effect<0 & Z<0) | (Effect>0 & Z>0), logTG_Z:=Z]
h[(Effect>0 & Z<0) | (Effect<0 & Z>0), logTG_Z:=(-1)*Z]
h <- h[rsid %in% vars,c(1,5)]

i <- fread("../external_datasets/Evangelou_PP_for_coloc.txt.gz", select=c(1,2,6))
setnames(i, c("rsid", "P", "Effect"))
i[, Z:=qnorm(P)]
i[(Effect<0 & Z<0) | (Effect>0 & Z>0), PP_Z:=Z]
i[(Effect>0 & Z<0) | (Effect<0 & Z>0), PP_Z:=(-1)*Z]
i <- i[rsid %in% vars,c(1,5)]

j <- fread("../external_datasets/Evangelou_SBP_for_coloc.txt.gz", select=c(1,2,6))
setnames(j, c("rsid", "P", "Effect"))
j[, Z:=qnorm(P)]
j[(Effect<0 & Z<0) | (Effect>0 & Z>0), SBP_Z:=Z]
j[(Effect>0 & Z<0) | (Effect<0 & Z>0), SBP_Z:=(-1)*Z]
j <- j[rsid %in% vars,c(1,5)]

k <- fread("../external_datasets/Evangelou_DBP_for_coloc.txt.gz", select=c(1,2,6))
setnames(k, c("rsid", "P", "Effect"))
k[, Z:=qnorm(P)]
k[(Effect<0 & Z<0) | (Effect>0 & Z>0), DBP_Z:=Z]
k[(Effect>0 & Z<0) | (Effect<0 & Z>0), DBP_Z:=(-1)*Z]
k <- k[rsid %in% a$variant,c(1,5)]

l <- fread("../external_datasets/1000G_CAC_metaanalysis_EA_AA_sumstats.txt.gz", select=c(1,2,6))
setnames(l, c("rsid", "P", "Effect"))
l[, Z:=qnorm(P)]
l[(Effect<0 & Z<0) | (Effect>0 & Z>0), CAC_1KG_Z:=Z]
l[(Effect>0 & Z<0) | (Effect<0 & Z>0), CAC_1KG_Z:=(-1)*Z]
l <- l[rsid %in% vars,c(1,5)]

m <- fread("../external_datasets/Topmed_CAC_GWAS.txt.gz", select=c(1,2,6))
setnames(m, c("rsid", "P", "Effect"))
m[, Z:=qnorm(P)]
m[(Effect<0 & Z<0) | (Effect>0 & Z>0), Topmed_CAC_Z:=Z]
m[(Effect>0 & Z<0) | (Effect<0 & Z>0), Topmed_CAC_Z:=(-1)*Z]
m <- m[rsid %in% vars,c(1,5)]

n <- fread("/project/cphg-millerlab/reference_datasets/CHARGE_IMT_plaque/IMT.EA.META.MAF1.HetDF4_jun.txt", select=c(1,19,8))
setnames(n, c("rsid", "P", "Effect"))
n[, Z:=qnorm(P)]
n[(Effect<0 & Z<0) | (Effect>0 & Z>0), IMT_Z:=Z]
n[(Effect>0 & Z<0) | (Effect<0 & Z>0), IMT_Z:=(-1)*Z]
n <- n[rsid %in% vars,c(1,5)]

o <- fread("/project/cphg-millerlab/reference_datasets/CHARGE_IMT_plaque/Plaque_meta_032218.txt", select=c(1,19,8))
setnames(o, c("rsid", "P", "Effect"))
o[, Z:=qnorm(P)]
o[(Effect<0 & Z<0) | (Effect>0 & Z>0), Plaque_Z:=Z]
o[(Effect>0 & Z<0) | (Effect<0 & Z>0), Plaque_Z:=(-1)*Z]
o <- o[rsid %in% vars,c(1,5)]

studies <- list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
setnames(a, "variant", "rsid")
z <- reduce(studies, full_join, by="rsid")

setorder(z, gene, position)
setnames(z, c("rsid", "chromosome", "position"), c("RSID", "CHR", "POS"))
z <- z[, c(2,3,1,6:20,4,5)]

write_delim(z, file=gzfile(paste0("/scratch/ch2um/Aug22_mixqtl_gwas_forpaintor_chr", chrom,".txt.gz")), escape="none")
