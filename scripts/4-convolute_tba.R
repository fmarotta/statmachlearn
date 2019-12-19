# Apply an input function to combine [TF] and TBA in biologically
# meaningful (or not) ways.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Nov 2019

library(biomaRt)
library(fplyr)
library(data.table)
library(docopt)

doc <- "Convolute TBA

Usage:
    convolute_tba.R --tba=<file> --tf=<file> --expr=<file> --tba_tf=<file> --input_function=<str> --threads=<int>

Options:
    --tba=<file>            Path to the TBA file.
    --tf=<file>             Path to HOCOMOCO's file.
    --expr=<file>           Path to the expression file.
    --tba_tf=<file>         Path to the output file.
    --input_function=<str>  Name of the input function.
    --threads=<int>         Number of threads.
    -h --help               Show this help and exit.

"
argv <- docopt(doc)

# Define the input functions
infunc.plus <- function(tf.tba, tf.expr) { # use with std expression
    tf.tba + tf.expr
}
infunc.times <- function(tf.tba, tf.expr) {
    tf.tba * tf.expr
}
infunc.etimes <- function(tf.tba, tf.expr) { # use this if your TBA is log-transformed
   exp(tf.tba) * tf.expr
}
infunc.hill <- function(tf.tba, tf.expr) {
    tf.expr / (1/tf.tba + tf.expr)
}
infunc.ehill <- function(tf.tba, tf.expr) { # use this if your TBA is log-transformed
    tf.expr / (1/exp(tf.tba) + tf.expr)
}

# Assign the appropriate function to infunc(), which therefore behaves sort of
# like a generic function. All the possible functions that
# infunc() can assume must be vectorised.
infunc <- match.fun(paste("infunc", argv$input_function, sep = "."))


# Define a function to remove the version (the part after the dot) from the ENSG
strip_ensg_version <- function(ensg) {
    strsplit(as.character(ensg), ".", fixed = TRUE)[[1]][1]
}
strip_ensg_version <- Vectorize(strip_ensg_version, USE.NAMES = FALSE)


# Read the expression and the TF files
expr <- fread(argv$expr, col.names = c("ensembl_gene_id", "IID", "EXPRESSION"), key = c("ensembl_gene_id", "IID"))[IID != ""]
tf <- fread(argv$tf, select = c("Model", "EntrezGene"), col.names = c("MODEL", "entrezgene_id"), key = "MODEL")
expr$ensembl_gene_id <- strip_ensg_version(expr$ensembl_gene_id)
# Some entrezgene_id are of the form "entrez1; entrez2". With the following code
# we expand all these instances into two different lines. Afterwards, the
# expression of the MODEL will be given by the average of all the genes
# associated to that MODEL.
entrez <- strsplit(tf$entrezgene_id, "; ")
tf <- data.table(MODEL = rep(tf$MODEL, sapply(entrez, length)), entrezgene_id = unlist(entrez))
tf$entrezgene_id <- as.integer(tf$entrezgene_id)

# Obtain a map to convert entrezgene_id to ENSG
ensembl <- useEnsembl(biomart = "ensembl", GRCh = "37", dataset = "hsapiens_gene_ensembl", host = "http://grch37.ensembl.org" )
map <- setDT(getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "entrezgene_id",
                   values = tf$entrezgene_id,
                   mart = ensembl))

# Add the expression to the map
map_expr <- merge(map, expr, by = "ensembl_gene_id")[, !"ensembl_gene_id"]
tf_expr <- merge(tf, map_expr, by = "entrezgene_id", allow.cartesian = TRUE)[, !"entrezgene_id"]
tf_expr <- tf_expr[, .(EXPRESSION = mean(EXPRESSION)), by = c("MODEL", "IID")]
setkeyv(tf_expr, "MODEL")

# Combine TBA and TF expression
ffply(argv$tba, argv$tba_tf, parallel = argv$threads, FUN = function(tba, by) {
    m <- merge(tba, tf_expr, by = c("MODEL", "IID"))
    m[order(REGION, MODEL, IID), .(REGION, MODEL, IID, TBA_TF = infunc(TBA, EXPRESSION))]
})
