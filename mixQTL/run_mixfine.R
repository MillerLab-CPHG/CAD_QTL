##### The purpose of this script is to run mixfine, a fine-mapping functionality from the mixQTL package
##### All elements of this package and relevant descriptions can be found here: https://github.com/hakyimlab/mixqtl

##### NOTE: the window for this script is set to +/-500kb from the _TSS_, not from the gene body which varies a lot
##### Additional note: this script generates output files per gene but for genes with a single causal variant writes
##### to a single file (""), remember to include this when concatenating results.

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mixqtl))
suppressPackageStartupMessages(library(Rfast))

###############################################################################
## Loading expression and variant annotation file just to sort by chromosome
variant_annotation <- fread(paste0("/path/to/mixqtl/input/converted_",chrom,"_variant_annotation_date.txt.gz"))
siggenes <- fread(paste0("/path/to/mixqtl/results/mixqtl_egene_ensgs_",chrom,"date.txt"), header=F)
mixres <- fread(paste0("/path/to/mixqtl/results/",chrom,"_mixQTL_nominal_for_finemapping.txt.gz"), select=c(3:6))
setnames(mixres, c("variant", "gene", "gene_name", "pval_meta"))
mixres <- mixres[gene %in% siggenes$V1 & !is.na(pval_meta),]
variant_annotation <- variant_annotation[variant %in% mixres$variant,]
mixres <- NULL

expression_annotation <- fread("/path/to/mixqtl/input/short_v32_gencode_mixqtl_input.txt")
expression_annotation <- expression_annotation[gene_id %in% siggenes$V1,]
siggenes <- NULL

# Data Loading & Prep
cat("Loading library size\n")
lib_size <- read_tsv("/path/to/mixqtl/input/RNAseq_library_size.txt") %>% select(indiv, lib_size)

cat("Loading haplotypes\n")
geno1 <- read_tsv(paste0("/path/to/mixqtl/input/converted_",chrom,"_date.txt.gz"))
geno2 <- read_tsv(paste0("/path/to/mixqtl/input/converted_",chrom,"_hap2_date.txt.gz"))

cat("Loading expression\n")
y1 <- read_tsv("/path/to/mixqtl/input/haps1_expression_short_date.tsv.gz")
y2 <- read_tsv("/path/to/mixqtl/input/haps2_expression_short_date.tsv.gz")
ytotal <- read_tsv("/path/to/mixqtl/input/total_raw_exp_cts_date.txt")

cat("Loading covariates\n")
cov_offset <- read_tsv("/path/to/mixqtl/input/QTL_covariate_file_date.txt")

cat("trimming individuals\n")
individuals <- colnames(geno1)[-1]
individuals <- individuals[individuals %in% colnames(cov_offset)]
#individuals <- individuals[individuals %in% colnames(geno2)]
individuals <- individuals[individuals %in% colnames(ytotal)]
individuals <- individuals[individuals %in% colnames(y1)]
#individuals <- individuals[individuals %in% colnames(y2)]
individuals <- individuals[individuals %in% lib_size$indiv]
cat(glue::glue("Kept {length(individuals)} individuals\n"))

cat("Preparing data\n")
cov_df <- as.data.frame(cov_offset[, c("ID", individuals)])
print(nrow(cov_df))
library_size <- lib_size[lib_size$indiv %in% individuals,] #Chani added this
lib_size <- as.data.table(lib_size)

haplotype_1 <- geno1[,c("variant", individuals)]
print(paste0("Hap 1 samples: ", ncol(haplotype_1)-1))
haplotype_2 <- geno2[,c("variant", individuals)]
print(paste0("Hap 2 samples: ", ncol(haplotype_2)-1))
exp_total <- ytotal[,c("gene_list", individuals)]
print(paste0("Total expression samples: ", ncol(exp_total)-1))
expression_count_1 <- y1[,c("gene_list", individuals)]
print(paste0("Exp count one samples: ", ncol(expression_count_1)-1))
expression_count_2 <- y2[,c("gene_list", individuals)]
print(paste0("Exp count two samples: ", ncol(expression_count_2)-1))

transpose_data <- function(d) {
  genes <- d$gene_list
  d <- d %>% select(-gene_list)
  td_ <- transpose(d)
  colnames(td_) <- genes
  td_$SUBJID <- individuals
  td_ <- td_ %>%
        relocate(SUBJID)
  return(td_)
}

etrc <- transpose_data(exp_total)
print(paste0("total expression samples: ", nrow(etrc)))
print(etrc[1:5,1:5])
eac1 <- transpose_data(expression_count_1)
print(paste0("expression hap1 samples: ", nrow(eac1)))
eac2 <- transpose_data(expression_count_2)
print(paste0("expression hap2 samples: ", nrow(eac2)))

get_haplotypes <- function(variants, haplotype) {
  g <- variants %>% select(variant) %>% inner_join(haplotype, by="variant")
  ids <- g$variant
  g <- g %>% select(-variant) %>% as.matrix %>% t
  colnames(g) <- ids
  return(g)
}

###############################################################################
# Processing
genes <- expression_annotation %>% filter(chromosome == chrom) %>% .$gene_id
n_genes <- length(genes)
print(n_genes)
first_write <- TRUE

############################

for (i in 1:n_genes) {
  gene_ <- genes[i]
        print(gene_)

  tryCatch({
    gene <- expression_annotation %>% filter(gene_id==gene_)
    print(gene)
    window_start <- max(0, (gene$start - 500000))
    window_end <- min(max(variant_annotation$position), gene$end + 500000)
    print(paste0("Start position: ", window_start, ",   End position: ", window_end))

        variants <- variant_annotation %>% filter(window_start <= position, position <= window_end)
        genotypes_1 <- get_haplotypes(variants, haplotype_1)
        genotypes_2 <- get_haplotypes(variants, haplotype_2)
        print(genotypes_2[1:2,1:5])

        if (is.null(etrc[, c(gene_)])) {
        a <- paste0(gene_, ": This gene had a problem with input data")
                write.table(a, "/path/to/mixqtl/results/mixfine/reruns.txt",
                append=T, row.names=F, quote=F, sep="\t", col.names=F)
                next
                }
##### Generate covariate offset based on what was used in mixQTL analyses
        covariate_offset = mixqtl::regress_against_covariate(trc = etrc[[gene$gene_id]],
                        covariates = cov_df, lib_size = lib_size$lib_size)
#####Run mixfine adjusting for offset covariates
        r <- mixfine(geno1 = genotypes_1, geno2 = genotypes_2,
                        #y1 = eac1[[gene$ensg]], y2 = eac2[[gene$ensg]], ytotal = etrc[[gene$ensg]],  
                        y1 = eac1[[gene$gene_id]], y2 = eac2[[gene$gene_id]], ytotal = etrc[[gene$gene_id]],
                        lib_size = lib_size$lib_size, cov_offset = covariate_offset)

#Cancels job for nonconverging genes and writes out gene name to reruns file
        if (is.null(r$sets$cs)) {
                a <- paste0(gene_, ": This gene did not converge")
                write.table(a, "/path/to/mixqtl/results/mixfine/reruns.txt",
                append=T, row.names=F, quote=F, sep="\t", col.names=F)
                next
                }

        r_ <- data.table(gene = gene$ensg, gene_name = gene$gene_name, variant = colnames(genotypes_1), stringsAsFactors=FALSE)
        r_[, rownumber:=(1:nrow(r_))]

        pips <- as.data.frame(r$pip)
        pips$rsid <- rownames(pips)
        setnames(pips, c("pip", "rsid"))
        pips <- as.data.table(pips)
        print(head(pips))

        vars <- r$sets$cs
#Workaround for genes which have only one signal with one 100% causal variant (these need to be added back in later)
        if (lengths(vars)<=1) {
                new <- r_[vars$L1,]
                pips <- pips[rsid %in% new$variant]
                r2 <- data.table(gene = gene$ensg, gene_name = gene$gene_name, r$sets$purity, coverage = r$sets$coverage,
                CS_num = r$sets$cs_index, bayes = summary(r)$cs$cs_log10bf, variant = new$variant, pip = pips$pip,
                variants[variant==pips$rsid, c(2:5)], stringsAsFactors=F)
                write.table(r2, "/path/to/mixqtl/results/mixfine/mixfine_summary_res_June2022.txt",
                row.names=F, quote=F, sep="\t", col.names=F, append=T)
        next
        }
        vars.m <- as.data.table(reshape2::melt(vars))
        vars.m[, CS_num:=as.numeric(substr(L1, 2, 3))]
        vars.m$L1 <- NULL
        setnames(vars.m, c("rownumber", "CS_num"))
print(head(vars.m))

        setkey(r_, rownumber)
        setkey(vars.m, rownumber)
        new <- r_[vars.m]
print(head(new))

        pips <- pips[rsid %in% new$variant]
        setkey(pips, rsid); setkey(new, variant)
        newout <- new[pips, nomatch=0]
print(head(newout)) #Now write an output file of all variants for each gene to use for rgional plotting
        write.table(newout, paste0("/path/to/mixqtl/results/mixfine/", gene_, "_for_plots.txt"), row.names=F, quote=F, sep="\t")

        setorder(newout, -pip)
        newout <-newout[!duplicated(CS_num)]
print(newout)
        setkey(newout, variant)
        setkey(variants, variant)
        new <- newout[variants, nomatch=0]
        setkey(new, CS_num, gene, gene_name)

        r2 <- data.table(gene = gene$ensg, gene_name = gene$gene_name, r$sets$purity,
        coverage = r$sets$coverage, CS_num = r$sets$cs_index, bayes = summary(r)$cs$cs_log10bf, stringsAsFactors=F)
        setkey(r2, CS_num, gene, gene_name)
        r3 <- r2[new]
        r3$rownumber <- NULL
#Append summary results and lead SNP for each independent signal per locus onto one output file
        write.table(r3, "/path/to/mixqtl/results/mixfine/mixfine_summary_res_Aug2022.txt", row.names=F, quote=F, sep="\t", col.names=F, append=T)
        },
            error=function(cond) {
            message(paste("stopped running on this gene:", gene_))
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        })
}




