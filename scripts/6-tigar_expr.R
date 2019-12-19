# Prepare an expression file satisfying TIGAR's requirements.
# Federico Marotta (federico.marotta@edu.unito.it)
# Dec 2019

library(data.table)
library(docopt)

"Tigar Expression

Usage:
    tigar_expr.R --annot=<file> --expr=<file> --tigexpr=<file>
    tigar_expr.R (-h | --help)

Options:
    --annot=<file>      File with the annotated contigs.
    --expr=<file>       Expression values.
    --tigexpr=<file>    Output file.
    -h --help           Print this help and exit.

" -> doc

argv <- docopt(doc)


annot <- fread(argv$annot,
               skip = 6,
               header = F,
               select = c(1, 3:5, 9),
               col.names = c("CHROM", "Type", "GeneStart", "GeneEnd", "Info"))
annot <- annot[Type == "gene"]
annot <- annot[CHROM %in% 1:22]

info <- strsplit(annot$Info, "\"", fixed = TRUE)
ensg <- sapply(info, function(i) i[2])
name <- sapply(info, function(i) i[10])

annot$Info <- NULL
annot$Type <- NULL
annot <- cbind(annot, TargetID = ensg, GeneName = name)

expr <- fread(argv$expr)
expr <- dcast(expr, GENE ~ IID, value.var = "EXPRESSION")

tigexpr <- merge(annot, expr, by.x = "TargetID", by.y = "GENE")
perm <- c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", names(tigexpr)[-(1:5)])

fwrite(tigexpr[, ..perm], argv$tigexpr, quote = FALSE, sep = "\t")
