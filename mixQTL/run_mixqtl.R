#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mixqtl))

chrom="chr22"

parser <- ArgumentParser(description='Process some integers')
parser$add_argument("-trc_cutoff", type="integer", default=20)
parser$add_argument("-asc_cutoff", type="integer", default=5)
parser$add_argument("-weight_cap", type="double", default=100)
parser$add_argument("-asc_cap", type="integer", default=5000)
parser$add_argument('-library_size')
parser$add_argument('-variant_annotation')
parser$add_argument('-haplotype_1')
parser$add_argument('-haplotype_2')
parser$add_argument('-covariates')
parser$add_argument('-expression_annotation')
parser$add_argument('-expression_total_count')
parser$add_argument('-expression_count_1')
parser$add_argument('-expression_count_2')
parser$add_argument('-window', type="integer", default=1000000)
parser$add_argument('-output')
args <- parser$parse_args()

cat(glue::glue("Running with trc_cutoff: {args$trc_cutoff}"), "\n")
cat(glue::glue("Running with asc_cutoff: {args$asc_cutoff}"), "\n")
cat(glue::glue("Running with weight_cap: {args$weight_cap}"), "\n")
cat(glue::glue("Running with asc_cap: {args$asc_cap}"), "\n")

###############################################################################
# Data Loading & Prep
cat("Loading library size\n")
library_size <- read_tsv(args$library_size) %>% select(indiv, lib_size)

cat("Loading variant annotation\n")
variant_annotation <- read_tsv(args$variant_annotation)

cat("Loading haplotypes\n")
haplotype_1 <- read_tsv(args$haplotype_1)
haplotype_2 <- read_tsv(args$haplotype_2)

cat("Loading expression annotation\n")
expression_annotation <- read_tsv(args$expression_annotation)
expression_annotation <- expression_annotation %>% filter(chromosome==chrom) #Chani added this

cat("Loading expression\n")
expression_total_count <- read_tsv(args$expression_total_count)
expression_count_1 <- read_tsv(args$expression_count_1)
expression_count_2 <- read_tsv(args$expression_count_2)

if(!is.null(args$covariates)) {
  cat("Loading covariates\n")
  covariates <- read_tsv(args$covariates)
}

cat("trimming individuals\n")
individuals <- colnames(haplotype_1)[-1]
#if(!is.null(args$covariates)) {
   individuals <- individuals[individuals %in% colnames(covariates)]
#  individuals <- individuals[individuals %in% covariates$SUBJID]
#}
individuals <- individuals[individuals %in% colnames(haplotype_2)]
individuals <- individuals[individuals %in% colnames(expression_total_count)]
individuals <- individuals[individuals %in% colnames(expression_count_1)]
individuals <- individuals[individuals %in% colnames(expression_count_2)]
individuals <- individuals[individuals %in% library_size$indiv]
cat(glue::glue("Kept {length(individuals)} individuals\n"))

cat("Preparing data\n")
if(!is.null(args$covariates)) {
  covariates <- covariates[, c("ID", individuals)]
}
print(nrow(covariates))
library_size <- library_size[library_size$indiv %in% individuals,] #Chani added this
haplotype_1 <- haplotype_1[,c("variant", individuals)]
print(paste0("Hap 1 samples: ", ncol(haplotype_1)))
haplotype_2 <- haplotype_2[,c("variant", individuals)]
print(paste0("Hap 2 samples: ", ncol(haplotype_2)))
expression_total_count <- expression_total_count[,c("gene_list", individuals)]
print(paste0("Total expression samples: ", ncol(expression_total_count)))
expression_count_1 <- expression_count_1[,c("gene_list", individuals)]
print(paste0("Exp count one samples: ", ncol(expression_count_1)))
expression_count_2 <- expression_count_2[,c("gene_list", individuals)]
print(paste0("Exp count two samples: ", ncol(expression_count_2)))

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

etrc <- transpose_data(expression_total_count)
print(paste0("total expression samples: ", nrow(etrc)))
eac1 <- transpose_data(expression_count_1)
print(paste0("expression hap1 samples: ", nrow(eac1)))
eac2 <- transpose_data(expression_count_2)
print(paste0("expression hap2 samples: ", nrow(eac2)))

covariates <- data.frame(covariates)

get_haplotypes <- function(variants, haplotype) {
  g <- variants %>% select(variant) %>% inner_join(haplotype, by="variant")
  ids <- g$variant
  g <- g %>% select(-variant) %>% as.matrix %>% t
  colnames(g) <- ids
  return(g)
}

###############################################################################
# Processing
cat("Processing\n")
expression_annotation <- expression_annotation %>% filter(chromosome==chrom & (gene_id %in% expression_total_count$gene_list | gene_id %in% expression_count_1$gene_list))

#genes <- c(colnames(etrc), colnames(eac1)) %>% unique
genes <- expression_annotation %>% filter(chromosome == chrom) %>% .$gene_id
n_genes <- length(genes)
print(n_genes)
first_write <- TRUE

########covariates$ID <- NULL

for (i in 1:n_genes) {
#for (i in 40:50) {
  gene_ <- genes[i]
        print(gene_)

  tryCatch({
    gene <- expression_annotation %>% filter(gene_id == gene_)
    print(gene)
    window_start <- max(0, (gene$start - args$window))
    window_end <- (gene$end + args$window)
        print(paste0("Start position: ", window_start, ",   End position: ", window_end))
    variants <- variant_annotation %>% filter(chromosome == gene$chromosome, window_start <= position, position <= window_end)
if (nrow(variants) == 0) {
      cat("No available variants\n")
      next;
      } else {
      cat(paste0("Number of variants being tested: ", nrow(variants), "\n"))
    }

    indiv_offset <- regress_against_covariate(etrc[[gene$gene_id]], library_size$lib_size, covariates)
    genotypes_1 <- get_haplotypes(variants, haplotype_1)
    genotypes_2 <- get_haplotypes(variants, haplotype_2)

    r <- mixqtl(genotypes_1, genotypes_2, eac1[[gene$gene_id]], eac2[[gene$gene_id]], etrc[[gene$gene_id]],
                library_size$lib_size, indiv_offset, trc_cutoff = args$trc_cutoff, asc_cutoff = args$asc_cutoff,
                weight_cap = args$weight_cap, asc_cap = args$asc_cap)
    r_ <- data.frame(gene = gene$gene_id, gene_name = gene$gene_name, variant = colnames(genotypes_1), stringsAsFactors=FALSE)

    r <- cbind(r_, r)

    gzf = gzfile(args$output, ifelse(first_write,"w", "a"));
    write.table(r, gzf, append=!first_write, col.names = first_write, row.names = FALSE, sep="\t", quote = FALSE)
    close(gzf)
    first_write <- FALSE
  },
  error = function(e){
    cat("Error: ", gene$gene_id, "\n")
    message(e)
    return(NA)
  })
}
