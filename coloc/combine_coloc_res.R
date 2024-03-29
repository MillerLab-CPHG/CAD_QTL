#### R script to combine the coloc results from GWAS used in this study

library(data.table)
library(purrr)
library(dplyr)

a <- fread("vdh_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(a, c("vdh_hit2", "vdh_hit1", "vdh_nsnps", "vdh_PPH3", "vdh_PPH4", "vdh_best1","vdh_best2", "vdh_best4", "gene_name"))
b <- fread("GCST011365_mixQTL_all_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(b, c("Hartiala_hit2", "Hartiala_hit1", "Hartiala_nsnps", "Hartiala_PPH3", "Hartiala_PPH4", "Hartiala_best1", "Hartiala_best2", "Hartiala_best4", "gene_name"))
c <- fread("Koyama_CAD_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(c, c("Koyama_hit2", "Koyama_hit1", "Koyama_nsnps", "Koyama_PPH3", "Koyama_PPH4", "Koyama_best1", "Koyama_best2", "Koyama_best4", "gene_name"))
d <- fread("BBJ_Riken_CAD_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(d, c("Matsunaga_hit2", "Matsunaga_hit1", "Matsunaga_nsnps", "Matsunaga_PPH3", "Matsunaga_PPH4", "Matsunaga_best1", "Matsunaga_best2", "Matsunaga_best4", "gene_name"))
e <- fread("Evangelou_DBP_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(e, c("DBP_hit2", "DBP_hit1", "DBP_nsnps", "DBP_PPH3", "DBP_PPH4", "DBP_best1","DBP_best2", "DBP_best4", "gene_name"))
f <- fread("Evangelou_SBP_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(f, c("SBP_hit2", "SBP_hit1", "SBP_nsnps", "SBP_PPH3", "SBP_PPH4", "SBP_best1","SBP_best2", "SBP_best4", "gene_name"))
g <- fread("Evangelou_PP_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(g, c("PP_hit2", "PP_hit1", "PP_nsnps", "PP_PPH3", "PP_PPH4", "PP_best1","PP_best2", "PP_best4", "gene_name"))
h <- fread("TC_metaanalysis_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(h, c("TC_hit2", "TC_hit1", "TC_nsnps", "TC_PPH3", "TC_PPH4", "TC_best1", "TC_best2", "TC_best4", "gene_name"))
j <- fread("HDL_metaanalysis_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(j, c("HDL_hit2", "HDL_hit1", "HDL_nsnps", "HDL_PPH3", "HDL_PPH4", "HDL_best1", "HDL_best2", "HDL_best4", "gene_name"))
i <- fread("LDL_metaanalysis_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(i, c("LDL_hit2", "LDL_hit1", "LDL_nsnps", "LDL_PPH3", "LDL_PPH4", "LDL_best1", "LDL_best2", "LDL_best4", "gene_name"))
k <- fread("logTG_metaanalysis_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(k, c("logTG_hit2", "logTG_hit1", "logTG_nsnps", "logTG_PPH3", "logTG_PPH4", "logTG_best1", "logTG_best2", "logTG_best4", "gene_name"))
l <- fread("CHARGE_plaque_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(l, c("CHARGE_plaque_hit2", "CHARGE_plaque_hit1", "CHARGE_plaque_nsnps", "CHARGE_plaque_PPH3", "CHARGE_plaque_PPH4", "CHARGE_plaque_best1", "CHARGE_plaque_best2", "CHARGE_plaque_best4", "gene_name"))
m <- fread("CHARGE_IMT_meta_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(m, c("CHARGE_IMT_hit2", "CHARGE_IMT_hit1", "CHARGE_IMT_nsnps", "CHARGE_IMT_PPH3", "CHARGE_IMT_PPH4", "CHARGE_IMT_best1", "CHARGE_IMT_best2", "CHARGE_IMT_best4", "gene_name"))
n <- fread("1000G_CAC_EA_AA_meta_coloc_June2022.txt", select=c(1:3,7:11,14))
setnames(n, c("CAC_1KG_hit2", "CAC_1KG_hit1", "CAC_1KG_nsnps", "CAC_1KG_PPH3", "CAC_1KG_PPH4", "CAC_1KG_best1", "CAC_1KG_best2", "CAC_1KG_best4", "gene_name"))
o <- fread("Topmed_CAC_coloc_Aug2022.txt", select=c(1:3,7:11,14))
setnames(o, c("Topmed_CAC_hit2", "Topmed_CAC_hit1", "Topmed_CAC_nsnps", "Topmed_CAC_PPH3", "Topmed_CAC_PPH4", "Topmed_CAC_best1", "Topmed_CAC_best2", "Topmed_CAC_best4", "gene_name"))

genes <- list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
coloc <- reduce(genes, full_join, by = "gene_name")
coloc[, coloc_evidence:=0]
coloc[Topmed_CAC_PPH4>0.8 | CAC_1KG_PPH4>0.8 | CHARGE_IMT_PPH4>0.8 | CHARGE_plaque_PPH4>0.8 | logTG_PPH4>0.8 | LDL_PPH4>0.8 | HDL_PPH4>0.8 | TC_PPH4>0.8 | PP_PPH4>0.8 | SBP_PPH4>0.8 | DBP_PPH4>0.8 | Matsunaga_PPH4>0.8 | Koyama_PPH4>0.8 | Hartiala_PPH4>0.8 | vdh_PPH4>0.8, coloc_evidence:=1]
coloc[, indy_evidence:=0]
coloc[Topmed_CAC_PPH3>0.8 | CAC_1KG_PPH3>0.8 | CHARGE_IMT_PPH3>0.8 | CHARGE_plaque_PPH3>0.8 | logTG_PPH3>0.8 | LDL_PPH3>0.8 | HDL_PPH3>0.8 | TC_PPH3>0.8 | PP_PPH3>0.8 | SBP_PPH3>0.8 | DBP_PPH3>0.8 | Matsunaga_PPH3>0.8 | Koyama_PPH3>0.8 | Hartiala_PPH3>0.8 | vdh_PPH3>0.8, indy_evidence:=1]

write.table(coloc, "UVA_QTL_coloc_res_Aug2022.txt", row.names=F, quote=F, sep="\t")
