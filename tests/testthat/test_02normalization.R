library(testthat)
library(hpgltools)
context("Are normalizations consistent over time?")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

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

## First make sure the pasilla_expt still has the stuff we expect
expected <- "This is an expt class."
actual <- pasilla_expt[["title"]]
test_that("Pasilla title?", {
    expect_equal(expected, actual)
})

## Ensure that the beginning count table library sizes are identical.
expected <- colSums(counts)
actual <- pasilla_expt[["original_libsize"]]
names(expected) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
test_that("Pasilla libsize?", {
    expect_equal(expected, actual)
})

## Check a few arbitrary counts to make sure they are maintained.
expected <- counts[test_genes, "untreated1"]
testing_counts <- Biobase::exprs(pasilla_expt[["expressionset"]])
actual <- as.numeric(testing_counts[test_genes, "untreated1"])
test_that("Pasilla count tables? (untreated1)", {
    expect_equal(expected, actual)
})

## Check that all samples agree for 1 gene.
test_gene <- "FBgn0062565"
expected <- as.numeric(counts[test_gene, ])
actual <- as.numeric(Biobase::exprs(pasilla_expt[["expressionset"]])[test_gene, ])
expected <- c(4, 7, 3, 3, 9, 10, 9)
test_that("Pasilla count tables? (gene FBgn0063565)", {
    expect_equal(expected, actual)
})

## Ensure that normalize_expt does not mess up the data when called without arguments (this wasn't true once)
unmolested <- s_p(normalize_expt(pasilla_expt))[["result"]]
expected <- as.matrix(Biobase::exprs(pasilla_expt[["expressionset"]]))  ## I accidently changed this to potentially return a data.frame
actual <- Biobase::exprs(unmolested[["expressionset"]])
test_that("Pasilla (un)normalized counts?", {
    expect_equal(expected, actual)
})

## Make sure that the pData information is maintained through normalization
expected <- Biobase::pData(pasilla_expt[["expressionset"]])
actual <- Biobase::pData(unmolested[["expressionset"]])
test_that("Pasilla (un)normalized pdata?", {
    expect_equal(expected, actual)
})

## Also ensure that the library sizes (which are very important for limma) are not messed up.
expected <- pasilla_expt[["libsize"]]
actual <- unmolested[["libsize"]]
test_that("Pasilla (un)normalized libsize?", {
    expect_equal(expected, actual)
})

## Make sure that my invocation of cpm() is the same as edgeR's.
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cpm"))[["result"]]
expected <- edgeR::cpm(pasilla_expt[["expressionset"]])
actual <- Biobase::exprs(pasilla_convert[["expressionset"]])
test_that("cpm conversions are equivalent?", {
    expect_equal(expected, actual)
})

## Check that the different ways of calling rpkm() are identical
pasilla_convert <- convert_counts(pasilla_expt, convert="rpkm")
pasilla_norm <- s_p(normalize_expt(pasilla_expt, convert="rpkm"))[["result"]]
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
tt <- s_p(require.auto("BSgenome.Dmelanogaster.UCSC.dm6"))
tt <- s_p(library("BSgenome.Dmelanogaster.UCSC.dm6"))
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_convert[["expressionset"]]))[, 1])
expected <- c(0.46719586, 0.03443909, 22.38432168, 54.32421464, 0.03908639, 106.58455139)
test_that("cp_seq_m works for TA?", {
    expect_equal(expected, actual)
})

## Repeat cp_seq_m() for ATG
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6, pattern="ATG"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_convert[["expressionset"]]))[, c("untreated1")])
expected <- c(0.51893853, 0.04536343, 27.76677691, 46.94320722, 0.05237078, 99.09109542)
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected, actual)
})

## Test normalizations -- I should change this to be automatically generated for expected
expected <- as.numeric(c(5.857143, 91.500000, 4400.000000, 543.785714, 10.714286))
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="quant"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- as.numeric(actual_df[test_genes, c("untreated1")])
test_that("quant normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Similar test for size-factor normalization
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(6702.317616, 7.906784, 3.514126, 222.268496, 524.483368, 193.276953)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("size-factor normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector)
})

## Check another size-factor normalization
expected <- c(4.392658, 80.824908, 4097.471407, 512.183926, 8.785316)
names(expected) <- test_genes
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf2"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("size-factor2 normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Oh I never noticed before that this is a log, too
expected <- c(5.488150, 7.082043, 12.021996, 9.160395, 5.707992)
names(expected) <- test_genes
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="vsd"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("vsd normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qsmooth")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink_median")

expected <- c(4.927997, 91.830657, 4765.366532, 613.466245, 9.342734)
names(expected) <- test_genes
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="tmm"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("tmm normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.902692, 90.336774, 4803.090308, 608.726226, 9.488822)
names(expected) <- test_genes
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="upperquartile"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("upperquartile normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.927854, 91.079703, 4840.296148, 615.582521, 9.205998)
names(expected) <- test_genes
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="rle"))[["result"]]
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("RLE normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## Test transformations
expected <- c(2.584963, 6.539159, 12.187661, 9.189825, 3.459432)
names(expected) <- test_genes
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log2"))[["result"]]
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("log2 transformation gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(0.7781513, 1.9684829, 3.6688516, 2.7664128, 1.0413927)
names(expected) <- test_genes
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log10"))[["result"]]
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("log10 transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(0.7781513, 1.9684829, 3.6688516, 2.7664128, 1.0413927)
names(expected) <- test_genes
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log"))[["result"]]
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("loge transformation gives expected values (why log10!?)?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## Test filter
expected <- c(7526, 7)
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cbcb"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("cbcb filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="pofa"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("pofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="kofa"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("kofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cv"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("cv filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(9784, 7)
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="simple"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("simple filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## Test batch
expected <- c(3.333333, 64.500000, 3040.166667, 383.916667, 7.083333)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limma"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME this is obviously broken
expected <- c(1.60048774, 0.02101530, -0.07524254, 0.15555548, 0.49697157)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limmaresid"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma-residuals batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(3.095141, 60.542865, 3032.546240, 355.354483, 6.666536)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combatmod"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("combatmod from cbcbSEQ batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.875168, 87.749156, 4418.793105, 557.089237, 9.641209)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="sva"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("sva batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO Figure out what is up with these
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat")  ## broken
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior") ## takes forever
expected <- c(2.267660, 71.646546, 3481.434409, 410.347808, 6.658058)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combat_scale"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("combat_scale gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO Figure this guy out
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior_scale") ## takes forever
expected <- c(4.610009, 82.109047, 4099.039071, 519.407501, 9.116170)
names(expected) <- test_genes
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="svaseq"))[["result"]]
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("svaseq gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME I broke ruvg
## pasilla_batch <- normalize_expt(pasilla_expt, batch="ruvg") ## broken for the moment
