start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("11norm_convert.R: Are normalizations consistent over time (Conversions)?\n")

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata[["sampleid"]] <- rownames(metadata)
## Uses these genes for quick tests
test_genes <- c("FBgn0000014","FBgn0000008","FBgn0000017","FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## Make sure that my invocation of cpm() is the same as edgeR's.
pasilla_convert <- sm(normalize_expt(pasilla_expt, convert="cpm"))
expected <- edgeR::cpm(pasilla_expt[["expressionset"]])
actual <- Biobase::exprs(pasilla_convert[["expressionset"]])
test_that("cpm conversions are equivalent?", {
    expect_equal(expected, actual)
})

## Check that the different ways of calling rpkm() are identical
pasilla_convert <- convert_counts(pasilla_expt, convert="rpkm")
pasilla_norm <- sm(normalize_expt(pasilla_expt, convert="rpkm"))
expected <- pasilla_convert[["count_table"]]
actual <- Biobase::exprs(pasilla_norm[["expressionset"]])
test_that("calling convert_counts and normalize_expt are equivalent?", {
    expect_equal(expected, actual)
})

## Similarly check that edgeR's rpkm() comes out the same
fdata_lengths <- as.vector(Biobase::fData(pasilla_expt[["expressionset"]])[["length"]])
names(fdata_lengths) <- rownames(Biobase::fData(pasilla_expt[["expressionset"]]))
expected <- edgeR::rpkm(Biobase::exprs(pasilla_expt[["expressionset"]]), gene.length=fdata_lengths)
actual <- Biobase::exprs(pasilla_norm[["expressionset"]])
test_that("rpkm conversions are equivalent?", {
    expect_equal(expected, actual)
})

## I have a modification of rpkm(), cp_seq_m(), which should give some expected results.
tt <- sm(require.auto("BSgenome.Dmelanogaster.UCSC.dm6"))
tt <- sm(library("BSgenome.Dmelanogaster.UCSC.dm6"))
pasilla_convert <- sm(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6))
expected <- c(0.03443909, 0.46719586, 22.38432168, 54.32421464, 0.03908639)
actual <- as.numeric(Biobase::exprs(pasilla_convert[["expressionset"]])[test_genes, 1])
test_that("cp_seq_m works for TA?", {
    expect_equal(expected, actual)
})

## Repeat cp_seq_m() for ATG
pasilla_convert <- sm(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6, pattern="ATG"))
expected <- c(0.04536343, 0.51893853, 27.76677691, 46.94320722, 0.05237078)
actual <- as.numeric(Biobase::exprs(pasilla_convert[["expressionset"]])[test_genes, c("untreated1")])
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 11norm_convert.R in ", elapsed, " seconds."))
