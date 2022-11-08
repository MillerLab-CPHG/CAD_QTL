### UVA_coronary_ALL_COLOC  
library(coloc)
library(data.table)
library(dplyr)
library(tidyverse)
setwd('/project/cphg-millerlab/chani/QTL_res/coloc_res')

## Binary traits use type="cc" and s for dataset1
study="case_control_1"

## Continuous traits use type="quant" and N for dataset1
#study="quant_trait_1"

sqtl <- fread("/path/to/sQTL_res/Tissue_type_for_coloc_date.txt.gz")
sqtl <- sqtl[stderr>0 & !(snp=="."),] 
sqtl[, ALT_AF:=1-ref_freq]
sqtl[, N:=100000] #This is specific to number of study participants

gwas <- fread("/path/to/published/GWAS/PMID_sumstats_forcoloc.tsv.gz")

setnames(gwas, c("snp", "p_value", "ref", "alt", "REF_AF", "beta", "standard_error"))
gwas <- gwas[snp %in% sqtl$snp & REF_AF>0 & REF_AF<1,]

        gene_list <- sqtl[!duplicated(splice_site), c("splice_site","chr","position")]
        genes <- sqtl[, c("splice_site", "snp", "beta", "N", "ALT_AF", "position", "stderr")]

dat_all<-data.table()

for (t in 1:length(gene_list$splice_site))        {

        uva_dat <- genes[splice_site==gene_list[t, 1]]
        uva_dat <- uva_dat[!duplicated(snp)]
        gwas_dat <- gwas[snp %in% uva_dat$snp,]
        gwas_dat <- gwas_dat[!duplicated(snp)]

        if(!(is.null(gwas_dat))) {
        print(genes[t, splice_site])
        }
        if(nrow(gwas_dat)==0)   {
        print("Next gene"); next
        }

        uva_dat <- uva_dat[snp %in% gwas_dat$snp,]
        setorder(gwas_dat, snp); setorder(uva_dat, snp)
        gwas_dat[, position:=uva_dat$position]

        dat2<-coloc.signals(
### Case proportions for "s" value below (for binary traits only): 
       dataset1=list(snp=gwas_dat$snp, type="cc", MAF=gwas_dat$REF_AF, beta=gwas_dat$beta, varbeta=gwas_dat$standard_error^2, position=gwas_dat$position, s=0.297),
### Sample sizes for quant traits for "N" value below (for continuous/quantitative traits only): 
#       dataset1=list(snp=gwas_dat$snp, type="quant", MAF=gwas_dat$REF_AF, N=35000, beta=gwas_dat$beta, varbeta=gwas_dat$standard_error^2, position=gwas_dat$position),
        dataset2=list(snp=uva_dat$snp, type="quant", MAF=uva_dat$ALT_AF, N=uva_dat$N, beta=uva_dat$beta, varbeta=uva_dat$stderr^2, position=uva_dat$position),
        method = c("single"), mode = c("iterative"),
        p1 = 1e-04, p2 = 1e-04, p12= 1e-05,
        maxhits = 3)

  dat2<-as.data.table(dat2$summary)
  dat2[, splice_site:=gene_list[t, 1]]
  dat_all<-rbind(dat_all, dat2)
}

setnames(dat_all, c(paste0(study,"_hit2"), paste0(study,"_hit1"), paste0(study,"_nsnps"), paste0(study,"_PPH0"), paste0(study,"_PPH1"), paste0(study,"_PPH2"), paste0(study,"_PPH3"), paste0(study,"_PPH4"), paste0(study, "_best1"), paste0(study, "_best2"), paste0(study, "_best4"), paste0(study, "_hit1_margz"), paste0(study, "_hit2_margz"), "splice_site"))

write.table(dat_all, paste0("sqtl_coloc/",study,"_sQTL_coloc_date.txt"), sep="\t", quote=F, row.names=F)
