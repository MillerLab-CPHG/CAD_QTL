##### 
# The purpose of this script is to combine my results files with annotations from a couple sources 
# (in this case ENCODE and ABC files from our own data) into files that can be subsetted for PAINTOR input. 
# The order is somewhat relevant because ABC annotations are gene specific and the output file has to have SNPs 
# in the exact same order as the QTL results file because there's no SNP column.
#####

library(data.table)
library(purrr)

chrom="1"

a <- fread(paste0("ENCODE_annot_chr",chrom,"_paintor.txt.gz"))
setnames(a, c("CHR", "POS", "RSID", "CTCF"))
b <- fread(paste0("ENCODE_h3k4me3_chr",chrom,"_paintor.txt.gz"))
setnames(b, c("CHR", "POS", "RSID", "H3K4ME3"))
c <- fread(paste0("ENCODE_h3k27ac_chr",chrom,"_paintor.txt.gz"))
setnames(c, c("CHR", "POS", "RSID", "H3K27AC"))

chr <- fread(paste0("/path/to/nominal/results/qtl_gwas_paintor_chr",chrom,".txt.gz", select=c(1:3,10))
chr[, n:=1:nrow(chr)]

abc <- fread(paste0("ENCODE_ABC_chr",chrom,"_paintor.txt.gz")
setnames(abc, c("CHR", "POS", "RSID", "ABC_score", "gene_name"))

anns <- list(chr, a, b, c, abc)
total <- reduce(anns, full_join, by=c("CHR", "POS", "RSID"))

total <- total[!is.na(n),]
total[is.na(CTCF), CTCF:=0]
total[is.na(H3K4ME3), H3K4ME3:=0]
total[is.na(H3K27AC), H3K27AC:=0]
total[is.na(ABC_score), ABC_score:=0]
setorder(total, n)
total$n <- NULL

head(total)
write.table(total, paste0("Paintor_annotations_chr",chrom,".txt"), row.names=F, quote=F, sep="\t")
