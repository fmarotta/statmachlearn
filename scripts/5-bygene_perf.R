# For each gene, perform a 5-fold cross-validation and save the R^2 in
# order to evaluate the performance of the model. Do a nested
# cross-validation to tune the parameters.
#
# Federico Marotta (federico.marotta@edu.unito.it),
# Oct 2019

library(cvtools)
library(data.table)
library(docopt)
library(fplyr)

doc <- "Models Performance

Usage:
    bygene_perf.R --tba=<file> --expr=<file> --genes=<file> --samples=<file> --model=<str> --autoparams=<str> [--nested_folds=<int>] [--threads=<int>] [--reshape=<str>] --out_prefix=<prefix>
    bygene_perf.R (-h | --help)

Options:
    --tba=<file>            File with the predictors.
    --expr=<file>           File with the response.
    --genes=<file>          List of genes on which to run the models.
    --samples=<file>        List of individuals on which to run the models.
    --model=<str>           Name of the model; must match a function in the models.R script.
    --autoparams=<str>      Name of the autoparams; must match a function in the models.R script.
    --nested_folds=<int>    Number of folds in the inner CV loop. [default: 10]
    --threads=<int>         Number of cores to use. [default: 5]
    --reshape=<str>         Method to reshape the TBA. [default: reshape_esum]
    --out_prefix=<prefix>   Path of the outputs.
    -h --help               Show this help and exit.

"
argv <- docopt(doc)

source("scripts/models.R")

# Read the training data sets
tba.file <- argv$tba
expression <- fread(argv$expr, key = c("GENE", "IID"))
samples <- readLines(argv$samples)
genes <- readLines(argv$genes)

# Filter for the samples and the genes
expression <- expression[IID %in% samples]
expression <- expression[GENE %in% genes]

# Define the threads and the number of folds for cross-validation
n.folds <- 5
n.nested_folds <- as.integer(argv$nested_folds)
threads <- as.integer(argv$threads)

# Define the reshaping methods
reshape_sum <- function(tba) {
    tba <- tba[, .(TBA = sum(TBA)), by = c("GENE", "MODEL", "IID")]
    dcast(tba, GENE + IID ~ MODEL, value.var = "TBA")
}
reshape_esum <- function(tba) {
    tba <- tba[, .(TBA = log(sum(exp(TBA)))), by = c("GENE", "MODEL", "IID")]
    dcast(tba, GENE + IID ~ MODEL, value.var = "TBA")
}
reshape_stagger <- function(tba) {
    tba <- dcast(tba, GENE + IID + MODEL ~ REGION, value.var = "TBA")
    dcast(tba, GENE + IID ~ MODEL, value.var = names(tba)[grep("^ENSR", names(tba))])
}
reshape_tba <- match.fun(argv$reshape)

# Five-Fold Nested Cross-Validation
cv_function <- function(tba) {
    tba <- reshape_tba(tba)
    d.train <- merge(expression, tba)
    if (nrow(d.train) == 0)
        return(NULL)

    # Remove columns of zero variance
    w <- sapply(d.train[, 4:ncol(d.train)], function(reg) var(reg) == 0)
    d.train[, which(c(F, F, F, w)) := NULL]
    if (ncol(d.train) <= 3)
        return(NULL)

    cat("Considering gene ", d.train$GENE[1], "\n")
    set.seed(2019)

    cv.out <- NestedCV(d.train,
                       match.fun(argv$model),
                       params.auto = match.fun(argv$autoparams)(d.train),
                       n.outer.folds = n.folds,
                       n.inner.folds = n.nested_folds,
                       threads = threads,
                       verbose = 2)

    # Return everything
    cv.out
}

l <- flply(tba.file, FUN = cv_function, sep = "\t", header = TRUE, col.names = c("GENE", "REGION", "MODEL", "IID", "TBA"))
saveRDS(l, paste0(argv$out, ".perf_full.Rds"))

inner.cv <- lapply(l, function(x) x$inner.cv)
saveRDS(inner.cv, paste0(argv$out, ".perf_inner.Rds"))

table <- lapply(l, function(x) x$table)
fivefold.cv.table <- do.call("rbind", table)
n_params <- which(names(fivefold.cv.table) == "data.name") - 1
fwrite(fivefold.cv.table[, -(1:n_params)], paste0(argv$out, ".perf_merg.tsv"),
       quote = FALSE, sep = "\t")

aggr <- lapply(l, function(x) x$aggr)
fivefold.cv.aggr_table <- do.call("rbind", aggr)
fwrite(fivefold.cv.aggr_table[, -(1:n_params)], paste0(argv$out, ".perf_aggr.tsv"),
       quote = FALSE, sep = "\t")

