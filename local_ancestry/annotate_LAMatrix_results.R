### The purpose of this script is to select the minimum p-value SNP for each gene from LAMatrix results
### adjusting for local ancestry regions for one of five groups: AFR, AMR, EAS, EUR, or SAS.
### Arguments for each numbered chromosome and file path to by-gene results are read in from a shell script.

library(data.table)
library(stringr)
library(dplyr)
library(R.utils)

margs=commandArgs(T)
chromosome=margs[1]
chrompath=margs[2]
freqfile=margs[3]
outpath=margs[4]

### Read in allele frequencies from our genotype data
freq <- fread(freqfile)
freq <- as_tibble(freq[chr==chromosome, ])
print(head(freq, n=3))

### Merge results of all SNPs for each gene and 
multimerge = function(mergefiles) {
        filenames=list.files(path=mergefiles, full.names=T)
	      datalist=combined_files <- bind_rows(lapply(filenames, fread))
        }

reslist <- as.data.table(multimerge(chrompath))
dim(reslist)
setnames(reslist, c("beta", "tstat", "pval", "geneid", "position", "rsid"))
setorder(reslist, geneid, pval)
print(head(reslist, n=3))
reslist[, chr:=as.integer(chromosome)]

### Combine frequency table and results
combo <- left_join(reslist, freq, by=c("chr", "position")) %>%
	filter(!is.na(alt_freq))
combo$rsid.y <- NULL
print(paste0("combo rows:", nrow(combo)))

############### 
### Write out all variants with effect size -5 < beta < 5 for each gene tested
### NOTE: this is a BIG file, but including all the rows is necessary for generating accurate
### FDRs/Benjamini-Hochberg adjusted p-value for lead variants to identify significant eGenes.
############### 
write.table(combo, paste0(outpath, "/Oct28_filtered_LAmatrix_eqtl_res_chr",chromosome,".txt"), quote=F, row.names=F, sep="\t")
print("filtered file p<0.1 written")


### Select lead SNP for each gene
leadsnps <- combo %>%
		group_by(geneid) %>%
		filter(beta < abs(5)) %>%
		slice_min(.reslist, .preserve=F, order_by=pval, with_ties=F) %>%
		filter(!is.na(position))

print(head(leadsnps, n=3))

############### 
### Write out the variant with effect size -5 < beta < 5 and the lowest p-value for each gene tested
############### 
write.table(leadsnps, paste0(outpath, "/Oct28_LAmatrix_eqtl_lead_snps_chr",chromosome,".txt"), quote=F, row.names=F, sep="\t")
print("lead SNP file written")



