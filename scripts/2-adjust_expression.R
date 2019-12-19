# adjust_expression.R
#
# Fit a linear model EXPRESSION ~ COVARIATES and obtain the residuals.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Oct 2019

library(data.table)
library(docopt)

doc <- "Adjust Expression

Usage:
    adjust_expression.R --qnorm_std=<file> --covariates=<file> --adj_expr=<file>

Options:
    --qnorm_std=<file>	    Path to the standardized expression.
    --covariates=<file>	    Path to the covariates.
    --adj_expr=<file>	    Path to the output file.
    -h --help		    Show this help and exit.

"
argv <- docopt(doc)

Transpose_table <- function(dt, id) {
    head <- c(id, dt[[1]])
    ret <- as.data.table(cbind(names(dt[, -1]), t(dt[, -1])))
    names(ret) <- head
    ret
}

# Read the input files
expr <- Transpose_table(fread(argv$qnorm_std), "IID")
covariates <- fread(argv$covariates)

# Transpose the covariates
covnames <- covariates$ID
samplenames <- names(covariates[, -1])
covariates <- as.data.table(t(covariates[, -1]))
names(covariates) <- covnames
covariates$IID <- samplenames

# Set the keys
setkeyv(expr, "IID")
setkeyv(covariates, "IID")

# For each gene, fit the linear model and return the residuals
cat("Fitting linear models... ")
d <- merge(expr, covariates, by = "IID")
genes <- names(expr[, -1])
residuals <- t(sapply(genes, function(g) {
    data <- d[, c(g, covnames), with = FALSE]
    names(data)[1] <- "EXPRESSION"
    lm.fit <- lm(EXPRESSION ~ ., data)
    lm.fit$residuals
}))
cat("Done.\n")

# Create a data table from the residuals
r <- as.data.table(residuals)
names(r) <- d$IID
r <- cbind(GENE = genes, r)

# Melt the data table and print it
mr <- melt(r,
           id.vars = "GENE",
           # measure.vars = samples,
           value.name = "EXPRESSION",
           variable.name = "IID",
           variable.factor = FALSE)[order(GENE, IID), ]
fwrite(mr, argv$adj_expr,
       quote = FALSE, sep = "\t")

