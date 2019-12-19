# Convert vcf_rider's output for the TBA in a matrix form. Moreover,
# exclude the regions whose indels are too large.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Oct,Dec 2019

library(data.table)
library(fplyr)
library(docopt)

doc <- "Aggregate TBA

Usage:
    aggregate_tba.R --snps=<file> --sorted_tba=<file> --tba=<file> --fraction=<float> --threads=<int>

Options:
    --snps=<file>           Path to the SNP file returned by vcf_rider.
    --sorted_tba=<file>     Path to the sorted TBA file.
    --tba=<file>            Path to the output file.
    --fraction=<float>      Exclude regions whose indels are larger than frac.
    --threads=<int>         The number of threads to use.
    -h --help               Show this help and exit.

"
argv <- docopt(doc)

snps <- fread(argv$snps, header = FALSE, col.names = c("REGION", "START", "END", "OV_MUT"))
fraction <- as.numeric(argv$fraction)

# Function to obtain the genes whose indels are too large
exclude_from_indel_length <- function(snps, fraction) {
    snps$CUTOFF <- (snps$END - snps$START) * fraction
    indels <- strsplit(snps$OV_MUT, ",")

    snps$MAX_INS <- sapply(indels, function(i) {
        i <- strsplit(i, "_")
        sum(sapply(i, function(j) {
            if (j[3] == "true" && j[2] > 0)
                as.numeric(j[2])
            else
                0
        }))
    })
    snps$MAX_DEL <- sapply(indels, function(i) {
        i <- strsplit(i, "_")
        sum(sapply(i, function(j) {
            if (j[3] == "true" && j[2] < 0)
                -as.numeric(j[2])
            else
                0
        }))
    })

    # snps[MAX_INS > CUTOFF | MAX_DEL > CUTOFF, c("REGION", "START", "END")]
    snps[MAX_INS > CUTOFF | MAX_DEL > CUTOFF]$REGION
}

# Find the regions to be excluded due to length problems
invalid.regions <- exclude_from_indel_length(snps, fraction)

# Aggregate the TBA
ffply(argv$sorted_tba, argv$tba, parallel = argv$threads,
      header = FALSE, col.names = c("GENE", "REGION", "START", "END", "MODEL", "IID", "ALLELE", "TBA"),
      FUN = function(d, by) {
    if (is.null(d) || nrow(d) == 0)
        return(NULL)
    # Exclude the problematic regions
    d <- d[!REGION %in% invalid.regions]
    if (nrow(d) == 0)
        return(NULL)

    # Sum the TBAs for each region
    d <- d[, .(TBA = sum(as.numeric(TBA))), by = c("REGION", "MODEL", "IID")]
    d$TBA <- log(d$TBA)

    return(d)
})

# Print the warnings so that I can check if everything is OK
warnings()
