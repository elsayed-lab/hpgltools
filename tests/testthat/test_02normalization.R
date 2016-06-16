library(testthat)
library(hpgltools)

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

context("Are normalizations consistent over time?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## First make sure the pasilla_expt still has the stuff we expect
actual <- pasilla_expt[["title"]]
expected <- "This is an expt class."
test_that("Pasilla title?", {
    expect_equal(expected, actual)
})

actual <- pasilla_expt[["original_libsize"]]
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(expected) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
test_that("Pasilla libsize?", {
    expect_equal(expected, actual)
})

actual <- as.numeric(head(Biobase::exprs(pasilla_expt[["expressionset"]]))[["untreated1"]])
expected <- c(7629, 9, 4, 253, 597, 220)
test_that("Pasilla count tables? (untreated1)", {
    expect_equal(expected, actual)
})

actual <- as.numeric(head(Biobase::exprs(pasilla_expt[["expressionset"]]))["FBgn0062565", ])
expected <- c(4, 7, 3, 3, 9, 10, 9)
test_that("Pasilla count tables? (gene FBgn0063565)", {
    expect_equal(expected, actual)
})

unmolested <- s_p(normalize_expt(pasilla_expt))[["result"]]
actual <- Biobase::exprs(unmolested[["expressionset"]])
expected <- as.matrix(Biobase::exprs(pasilla_expt[["expressionset"]]))  ## I accidently changed this to potentially return a data.frame
test_that("Pasilla (un)normalized counts?", {
    expect_equal(expected, actual)
})

actual <- Biobase::pData(unmolested[["expressionset"]])
expected <- Biobase::pData(pasilla_expt[["expressionset"]])
test_that("Pasilla (un)normalized pdata?", {
    expect_equal(expected, actual)
})

actual <- unmolested[["libsize"]]
expected <- pasilla_expt[["libsize"]]
test_that("Pasilla (un)normalized libsize?", {
    expect_equal(expected, actual)
})

## First test conversions alone
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cpm"))[["result"]]
actual <- Biobase::exprs(pasilla_convert[["expressionset"]])
expected <- edgeR::cpm(pasilla_expt[["expressionset"]])
test_that("cpm conversions are equivalent?", {
    expect_equal(expected, actual)
})

pasilla_convert <- convert_counts(pasilla_expt, convert="rpkm")
pasilla_norm <- s_p(normalize_expt(pasilla_expt, convert="rpkm"))[["result"]]
test_that("calling convert_counts and normalize_expt are equivalent?", {
    expect_equal(pasilla_convert[["count_table"]], Biobase::exprs(pasilla_norm[["expressionset"]]))
})

fdata_lengths <- as.vector(Biobase::fData(pasilla_expt[["expressionset"]])[["length"]])
names(fdata_lengths) <- rownames(Biobase::fData(pasilla_expt[["expressionset"]]))
actual <- Biobase::exprs(pasilla_norm[["expressionset"]])
expected <- edgeR::rpkm(x=exprs(pasilla_expt[["expressionset"]]), gene.length=fdata_lengths)
test_that("rpkm conversions are equivalent?", {
    expect_equal(expected, actual)
})

tt <- s_p(require.auto("BSgenome.Dmelanogaster.UCSC.dm6"))
tt <- s_p(library("BSgenome.Dmelanogaster.UCSC.dm6"))
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_convert[["expressionset"]]))[, 1])
expected <- c(0.46719586, 0.03443909, 22.38432168, 54.32421464, 0.03908639, 106.58455139)
test_that("cp_seq_m works for TA?", {
    expect_equal(expected, actual)
})

pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6, pattern="ATG"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_convert[["expressionset"]]))[, c("untreated1")])
expected <- c(0.51893853, 0.04536343, 27.76677691, 46.94320722, 0.05237078, 99.09109542)
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected, actual)
})

## Test normalizations
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="quant"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- as.vector(c(7225.000000, 10.142857, 4.714286, 237.285714, 558.714286, 207.285714))
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")

test_that("quant normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(6702.317616, 7.906784, 3.514126, 222.268496, 524.483368, 193.276953)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("size-factor normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf2"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(6702.317616, 7.906784, 3.514126, 222.268496, 524.483368, 193.276953)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("size-factor2 normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector)
})

## Oh I never noticed before that this is a log, too
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="vsd"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(12.723632, 5.669751, 5.431568, 8.132187, 9.191164, 7.972554)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("vsd normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## FIXME
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qsmooth")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink_median")

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="tmm"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(7614.957402, 8.870394, 4.086935, 266.221201, 557.761246, 221.757090)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("tmm normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="upperquartile"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(7614.957402, 8.870394, 4.086935, 266.221201, 557.761246, 221.757090)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("upperquartile normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="rle"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_norm[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(7552.685400, 8.870137, 4.151197, 267.139584, 549.598099, 223.089510)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("RLE normalization gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## Test transformations
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log2"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_trans[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(12.897467, 3.321928, 2.321928, 7.988685, 9.224002, 7.787903)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("log2 transformation gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log10"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_trans[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(3.882525, 1.000000, 0.698970, 2.404834, 2.776701, 2.344392)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("log10 transformation gives expected values (why log10!?)?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_trans[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(8.939843, 2.302585, 1.609438, 5.537334, 6.393591, 5.398163)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("loge transformation gives expected values (why log10!?)?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## Test filter
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cbcb"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
expected <- c(7526, 7)
test_that("cbcb filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="pofa"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
expected <- c(10153, 7)
test_that("pofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="kofa"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
expected <- c(10153, 7)
test_that("kofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cv"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
expected <- c(10153, 7)
test_that("cv filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="simple"))[["result"]]
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
expected <- c(9784, 7)
test_that("simple filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## Test batch
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limma"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(5271.833333, 4.458333, 3.791667, 230.000000, 413.625000, 160.875000)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("limma batch gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## FIXME this is obviously broken
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limmaresid"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(0.05337675, -0.36662233, -0.02590011, 0.54896800, 0.05749937, 0.08030627)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("limma-residuals batch gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combatmod"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(5228.115819, 3.807374, 3.761282, 226.661066, 387.082964, 152.333373)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("combatmod from cbcbSEQ batch gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="sva"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(7301.355553, 8.185490, 3.877108, 250.278232, 570.622325, 211.156206)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("sva batch gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## TODO Figure out what is up with these
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat")  ## broken
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior") ## takes forever

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combat_scale"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(5650.664883, 6.702011, 3.283283, 229.060736, 447.552075, 166.215787)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("combat_scale gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## TODO Figure this guy out
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior_scale") ## takes forever

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="svaseq"))[["result"]]
actual_df <- head(Biobase::exprs(pasilla_batch[["expressionset"]]))
actual_vector <- actual_df[, c("untreated1")]
expected_vector <- c(6836.560113, 7.163647, 3.750438, 246.863925, 534.647045, 198.986139)
names(expected_vector) <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("svaseq gives expected values?", {
    expect_equal(expected_vector, actual_vector, tolerance=0.0001)
})

## FIXME I broke ruvg
## pasilla_batch <- normalize_expt(pasilla_expt, batch="ruvg") ## broken for the moment
