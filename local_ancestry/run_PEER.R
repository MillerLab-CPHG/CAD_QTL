###### The purpose of this script is to run probabilistic estimate of expression residuals (PEER, PMID: )
### on a computing cluster with arguments supplied by a shell script. Demographic covariates such as age
### and sex were included in the linear regression model and hence were not adjusted for here.

# Author of original script: Francois Aguet

library(peer, quietly=TRUE)  # https://github.com/PMBio/peer
library(argparser, quietly=TRUE)
library(data.table)

WriteTable <- function(data, filename, index.name) {
    datafile <- file(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))
    writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
    write.table(data, datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

p <- arg_parser("Run PEER factor estimation")
p <- add_argument(p, "expr_file", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "n", help="Number of hidden confounders to estimate")
p <- add_argument(p, "--covariates", help="Observed covariates")
p <- add_argument(p, "--alphaprior_a", help="", default=0.001)
p <- add_argument(p, "--alphaprior_b", help="", default=0.01)
p <- add_argument(p, "--epsprior_a", help="", default=0.1)
p <- add_argument(p, "--epsprior_b", help="", default=10)
p <- add_argument(p, "--max_iter", help="", default=200)
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv <- parse_args(p)

cat("PEER: loading expression data ... ")

df <- fread(argv$expr_file, header=T)
        print(head(df[1:5, 1:10]))

nrows <- nrow(df)
        print(nrows)

id <- as.vector(df[, 4])
        print(head(id))
        row.names(df) <- id$pid

df <- df[, 7:ncol(df)]
        print(head(df[, 1:10]))

#M <- t(as.matrix(df))
M <- t(as.matrix(df))
        print(dim(M))
        print(head(M[, 1:5]))

cat("done.\n")


# run PEER
cat(paste0("PEER: estimating hidden confounders (", argv$n, ")\n"))
model = PEER()
invisible(PEER_setNk(model, argv$n))
invisible(PEER_setPhenoMean(model, M))
dim(PEER_getPhenoMean(model))
invisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))
invisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))
invisible(PEER_setNmax_iterations(model, argv$max_iter))
if (!is.null(argv$covariates) && !is.na(argv$covariates)) {
    covar.df <- read.table(argv$covariates, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    covar.df <- sapply(covar.df, as.numeric)
    cat(paste0("  * including ", dim(covar.df)[2], " covariates", "\n"))
    invisible(PEER_setCovariates(model, as.matrix(covar.df)))  # samples x covariates
}
time <- system.time(PEER_update(model))

X <- PEER_getX(model)  # samples x PEER factors
A <- PEER_getAlpha(model)  # PEER factors x 1
R <- t(PEER_getResiduals(model))  # genes x samples

# add relevant row/column names
c <- paste0("InferredCov",1:ncol(X))
rownames(X) <- rownames(M)
colnames(X) <- c
rownames(A) <- c
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(M)
colnames(R) <- rownames(M)

# write results
cat("PEER: writing results ... ")
WriteTable(t(X), file.path(argv$output_dir, paste0(argv$prefix, ".PEER_covariates.txt")), "ID")  # format(X, digits=6)
WriteTable(A, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_alpha.txt")), "ID")
WriteTable(R, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_residuals.txt")), "ID")
cat("done.\n")
