##### The purpose of this file is to merge ENCODE annotations with the QTL and GWAS Zscores to make 
##### gene-specific paintor input files

library(data.table)
library(readr)
library(dplyr)

chromosome=1

annot <- fread("/path/to/folder/coronary_artery_ENCODE_ABCpredict.txt")
annot <- annot[chr==paste0("chr", chromosome)]

chrom <- fread(paste0("/path/to/folder/qtl_gwas_forpaintor_chr",chromosome,".txt.gz"), select=c(1:3))
chrom <- chrom[!duplicated(RSID)]
setorder(chrom, POS)

        for (s in seq(1, nrow(annot), by=1)) {
        spos <- as.numeric(annot[s, c("start")])
        epos <- as.numeric(annot[s, c("end")])
        print(spos)

        for (i in seq_along(1:nrow(chrom)))       {
                newchrom <- as.data.table(chrom[i,])
                if(newchrom$POS<spos) next
                if(newchrom$POS>epos) {
                                print("Next locus")
                                break
                                }
                newchrom[, ABC:=ifelse(POS>spos & POS<epos, annot[s, c("ABC.Score")], 0)]
                newchrom[, ABC_gene:=ifelse(POS>spos & POS<epos, annot[s, c("TargetGene")], 0)]

                write_delim(newchrom, file=gzfile(paste0("ENCODE_ABC_chr",chromosome,"_paintor.txt.gz"), append=T, escape="none")

                }
        }
