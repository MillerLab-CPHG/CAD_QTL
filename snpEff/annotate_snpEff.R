library(data.table)
library(dplyr)
library(stringr)

##### Purpose of this script is annotating SNPEff results from logfile #####
snpeff <- fread("qtl_snpeff_annotations.txt", select=c(1:3,8)) ## read in log results
setnames(snpeff, c("chr", "pos", "rsid", "annotation"))
setkey(snpeff, rsid)
eqtl_sig <- fread("eQTL_results_file_lead_variants.txt", select=c("gene_name", "rsid", "rsid_pos", "TSS_dist", "mix_bh")) ##read in QTL results (prepared separately)
setkey(eqtl_sig, rsid)
d <- eqtl_sig[snpeff, nomatch=0]

gv32_t <- fread("gencode_v32_transcript_annotation.txt") ##Need to use isoforms since snpEff annotates by transcript not gene
setkey(gv32_t, gene_name); setkey(d, gene_name)
e <- gv32_t[d, nomatch=0]

transcripts <- e$transcript
b <- as.data.table("gene_id")

## Clunky but functional loop to grab the right combo of variant-gene based on your eQTL results
for (i in seq_along(1:length(transcripts))) {
trans_test=e[i]
splits <- as.data.table(strsplit(trans_test[, annotation], "\\|\\|\\|\\|\\|\\|")) ##Split string of results for each SNP with "||||||" separation
splits[, var_annot:=str_match(splits[,V1], "\\|.+?variant")]
splits$var_annot <- substr(splits$var_annot, 2, 50) 
splits[, enst:=str_match(splits[,V1], "ENST0000.+?\\.[0-9]")]
splits$V1 <- NULL; splits <- splits[enst %like% "ENST"]
ensts <- as.vector(splits$enst)
trans_test[, realsnp:=ifelse(transcript %in% ensts, 1, 0)] ##Label SNPs for which the annotation matches your QTL input gene
trans_test$annotation <- NULL
setkey(trans_test, transcript); setkey(splits, enst)
a <- trans_test[splits]
b <- bind_rows(b, a)
}

c <- b[realsnp==1, c(2,3,5:11,4,12:14,18)] ##Exclude SNPs which were annotated to a different gene (i.e., not a QTL for that gene)
write.table(c, "snpeff_annotated_results_date_year.txt", row.names=F, quote=F, sep="\t")
