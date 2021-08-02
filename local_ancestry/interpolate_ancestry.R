###### The purpose of this R script is to assign RFmix-inferred ancestry to each chromosomal position corresponding to a SNP
### in our VCF files. Position-specific ancestry estimates will then be incorporated as "local ancestry" covariates in eQTL analysis.

library(data.table)
library(stringr)
library(dplyr)

margs<-commandArgs(T)

chrom=margs[1]
individual=margs[2]

subset=paste0("/path/to/your/VCF/UVA_b37_chr",chrom,"_filtered.vcf.gz")
mapfile=paste0("/path/to/your/rfmixresults/chr",chrom,"_UVA_local_QTL.msp.tsv")

vcf <- fread(subset)
setnames(vcf, c("#CHROM", "POS", "ID"), c("chr", "pos", "rsid"))

indiv <- individual
print(indiv)
indiv1=paste0(indiv, "_0")
indiv2=paste0(indiv, "_1")

bed <- fread(mapfile)
myvect <- c("#chm", "spos", "epos", indiv1, indiv2)
newbed <- bed[, ..myvect]
setnames(newbed, c(indiv1, indiv2), c("V4", "V5"))

# Screening out individuals who have 100% ancestry on both haplotypes and writing that ancestry to output file
# Our population didn't have anyone reported at 100% African, Native American/Amerindian, or South Asian ancestry
# but those can be added if relevant
a <- nrow(newbed)
#af1 <- nrow(newbed[V4==0, ])
#af2 <- nrow(newbed[V5==0, ])
ea1 <- nrow(newbed[V4==2, ])
ea2 <- nrow(newbed[V5==2, ])
eu1 <- nrow(newbed[V4==3, ])
eu2 <- nrow(newbed[V5==3, ])

        if (a==ea1 & a==ea2) {
                newvcf <- as.data.table(vcf[, c(1:3)])
                newvcf[, paste0(indiv, "_AFR"):=0]
                newvcf[, paste0(indiv, "_AMR"):=0]
                newvcf[, paste0(indiv, "_EAS"):=2]
                newvcf[, paste0(indiv, "_EUR"):=0]
                newvcf[, paste0(indiv, "_SAS"):=0]
                newvcf[, paste0(indiv, "_UNK"):=0]
                write.table(newvcf, paste0("/project/cphg-millerlab/chani/QTL_input/local_anc/chr",chrom,"/",indiv, "_chr", chrom,"_hg19.maf1.localcov.tsv"), row.names=F, quote=F, sep="\t", col.names=F)

                print("This person is all East Asian ancestry on both haplotypes")
                q("no")
                }

                if (a==eu1 & a==eu2) {
                newvcf <- as.data.table(vcf[, c(1:3)])
                newvcf[, paste0(indiv, "_AFR"):=0]
                newvcf[, paste0(indiv, "_AMR"):=0]
                newvcf[, paste0(indiv, "_EAS"):=0]
                newvcf[, paste0(indiv, "_EUR"):=2]
                newvcf[, paste0(indiv, "_SAS"):=0]
                newvcf[, paste0(indiv, "_UNK"):=0]
                write.table(newvcf, paste0("/project/cphg-millerlab/chani/QTL_input/local_anc/chr",chrom,"/",indiv, "_chr", chrom,"_hg19.maf1.localcov.tsv"), row.names=F, quote=F, sep="\t", col.names=F)

                print("This person is all European ancestry on both haplotypes")
                q("no")
                }


########### Moving on for individuals with any admixiture
newvcf <- as.data.table(c("chr"))
newvcf$V2 <- "pos"; newvcf$V3 <- "rsid"; newvcf$V4 <- paste0(indiv, "_AFR"); newvcf$V5 <- paste0(indiv, "_AMR")
newvcf$V6 <- paste0(indiv, "_EAS"); newvcf$V7 <- paste0(indiv, "_EUR"); newvcf$V8 <- paste0(indiv, "_SAS"); newvcf$V9 <- paste0(indiv, "_UNK")
setnames(newvcf, c("chr", "pos", "rsid",paste0(indiv, "_AFR"), paste0(indiv, "_AMR"),paste0(indiv, "_EAS"),paste0(indiv, "_EUR"),paste0(indiv, "_SAS"),paste0(indiv, "_UNK")))
write.table(newvcf, paste0("/project/cphg-millerlab/chani/QTL_input/local_anc/chr",chrom,"/",indiv, "_chr", chrom,"_hg19.maf1.localcov.tsv"), row.names=F, quote=F, sep="\t", col.names=F)


################ Loop for each segment of relevant chromosome in the RFmix output file ##################
        for (s in seq(1, nrow(newbed), by=1)) {
        spos <- as.numeric(newbed[s, c("spos")])
        epos <- as.numeric(newbed[s, c("epos")])
        print(spos)

################ Loop for each SNP in the VCF file within the relevant chromosomal segment ###################
        for (i in seq_along(1:nrow(vcf)))       {
                newvcf <- as.data.table(vcf[i, c(1:3)])
                if(newvcf$pos<spos) next
                if(newvcf$pos>epos) {
                                print("Next chunk on chromosome ", chrom)
                                break
                                }
                newvcf[, ind_1:=newbed[s, c(4)]]
                newvcf[, ind_2:=newbed[s, c(5)]]

                newvcf[, UNK_1 := ifelse(ind_1!=0 & ind_1!=1 & ind_1!=2 & ind_1!=3 & ind_1!=4, 1, 0)]
                newvcf[, UNK_2 := ifelse(ind_2!=0 & ind_2!=1 & ind_2!=2 & ind_2!=3 & ind_2!=4, 1, 0)]
                newvcf[, AFR_1 := ifelse(ind_1==0, 1, 0)]; newvcf[, AFR_2 := ifelse(ind_2==0, 1, 0)]
                newvcf[, AMR_1 := ifelse(ind_1==1, 1, 0)]; newvcf[, AMR_2 := ifelse(ind_2==1, 1, 0)]
                newvcf[, EAS_1 := ifelse(ind_1==2, 1, 0)]; newvcf[, EAS_2 := ifelse(ind_2==2, 1, 0)]
                newvcf[, EUR_1 := ifelse(ind_1==3, 1, 0)]; newvcf[, EUR_2 := ifelse(ind_2==3, 1, 0)]
                newvcf[, SAS_1 := ifelse(ind_1==4, 1, 0)]; newvcf[, SAS_2 := ifelse(ind_2==4, 1, 0)]

                newvcf[, paste0(indiv, "_AFR") := AFR_1 + AFR_2]
                newvcf[, paste0(indiv, "_AMR") := AMR_1 + AMR_2]
                newvcf[, paste0(indiv, "_EAS") := EAS_1 + EAS_2]
                newvcf[, paste0(indiv, "_EUR") := EUR_1 + EUR_2]
                newvcf[, paste0(indiv, "_SAS") := SAS_1 + SAS_2]
                newvcf[, paste0(indiv, "_UNK") := UNK_1 + UNK_2]

        newvcf <- newvcf[, c(1:3,18:23)]
        print(newvcf)

        write.table(newvcf, paste0("/path/to/your/QTLinput/datasets/chr",chrom,"/",indiv, "_chr", chrom,"_hg19.maf1.localcov.tsv"), row.names=F, quote=F, sep="\t", append=T, col.names=F)
        }
}

