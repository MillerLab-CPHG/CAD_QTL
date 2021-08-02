### The purpose of this script is generating PCs for ~171 VCF files using SNPRelate BioConductor package
### Code adapted from https://github.com/zhengxwen/SNPRelate

### If package not installed:
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#       BiocManager::install("SNPRelate")
#       BiocManager::install("gdsfmt")
#install.packages("data.table", type = "source", repos = "https://Rdatatable.gitlab.io/data.table")

margs<-commandArgs(T)
vcf_file=margs[1]
out_pc_file=margs[2]

library(SNPRelate)
library(gdsfmt)
library(data.table)

##VCF file name to read in for GDS filetype conversion: 
vcf.fn <- vcf_file

## Convert to GDS using:
snpgdsVCF2GDS(vcf.fn, paste0("/path/you/store/stuff/gds_",vcf_file,"_", Sys.Date()), method="biallelic.only", ignore.chr.prefix = "chr", snpfirstdim=F, verbose=T)

## Check out the fancy file you've made:
### NOTE: this file is HUGE. If you want to keep it, definitely zip it; otherwise remove after you're done
snpgdsSummary(paste0("/path/you/store/stuff/gds_",vcf_file,"_", Sys.Date()))

##Get rid of all those totally lame SNPs in LD with each other because that will definitely mess up your PCs in a sample size this small. Start by opening the file:
genofile <- snpgdsOpen(paste0("/path/you/store/stuff/gds_",vcf_file,"_", Sys.Date()))

##Next, hit up some sweet LD pruning using whatever parameters you like. 
snpsforpca <- snpgdsLDpruning(genofile, slide.max.bp=500000, remove.monosnp=T, ld.threshold=0.2, autosome.only=T, method="corr")
names(snpsforpca)
head(snpsforpca$chr1)
pcasnp.id <- unlist(unname(snpsforpca))

##Then run the actual PCA:
pca <- snpgdsPCA(genofile, snp.id=pcasnp.id, num.thread=4)

### Add in population IDs
sample.id <- as.vector(read.gdsn(index.gdsn(genofile, "sample.id")))

pcaout <- as.data.frame(pca$eigenvect)
pcaout$sampleid <- sample.id

write.table(pcaout, paste0("/path/where/output/file/belongs/PCs_",vcf_file,".txt"), sep="\t", row.names=T, quote=F)

