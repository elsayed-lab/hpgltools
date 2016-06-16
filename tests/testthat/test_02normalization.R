library(testthat)
library(hpgltools)

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

actual <- as.numeric(head(Biobase::exprs(pasilla_expt[["expressionset"]][, 1])))
expected <- c(92, 5, 4664, 583, 10, 1446)
test_that("Pasilla count tables?", {
    expect_equal(expected, actual)
})

unmolested <- s_p(normalize_expt(pasilla_expt))[["result"]]
actual <- Biobase::exprs(unmolested[["expressionset"]])
expected <- Biobase::exprs(pasilla_expt[["expressionset"]])
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

pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="rpkm"))[["result"]]
lengths <- as.vector(Biobase::fData(pasilla_expt[["expressionset"]])[["length"]])
actual <- Biobase::exprs(pasilla_convert[["expressionset"]])
expected <- edgeR::rpkm(x=pasilla_expt[["expressionset"]], gene.length=lengths)
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
actual <- as.numeric(head(Biobase::exprs(pasilla_convert[["expressionset"]]))[, 1])
expected <- c(0.51893853, 0.04536343, 27.76677691, 46.94320722, 0.05237078, 99.09109542)
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected, actual)
})

## Test normalizations
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="quant"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(91.500000, 5.857143, 4400.000000, 543.785714, 10.714286, 1352.428571)
test_that("quant normalization gives expected values?", {
    expect_equal(expected, actual)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(80.824908, 4.392658, 4097.471407, 512.183926, 8.785316, 1270.356701)
test_that("size-factor normalization gives expected values?", {
    expect_equal(expected, actual)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf2"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(80.824908, 4.392658, 4097.471407, 512.183926, 8.785316, 1270.356701)
test_that("size-factor2 normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Oh I never noticed before that this is a log, too
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="vsd"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(7.082043, 5.488150, 12.021996, 9.160395, 5.707992, 10.378651)
test_that("vsd normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qsmooth")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink_median")

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="tmm"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(91.830657, 4.927997, 4765.366532, 613.466245, 9.342734, 1457.548873)
test_that("tmm normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="upperquartile"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(90.336774, 4.902692, 4803.090308, 608.726226, 9.488822, 1462.362203)
test_that("upperquartile normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="rle"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_norm[["expressionset"]]))[, 1])
expected <- c(91.079703, 4.927854, 4840.296148, 615.582521, 9.205998, 1466.306507)
test_that("RLE normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## Test transformations
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log2"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_trans[["expressionset"]]))[, 1])
expected <- c(6.539159, 2.584963, 12.187661, 9.189825, 3.459432, 10.498849)
test_that("log2 transformation gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log10"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_trans[["expressionset"]]))[, 1])
expected <- c(1.9684829, 0.7781513, 3.6688516, 2.7664128, 1.0413927, 3.1604685)
test_that("log10 transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_trans[["expressionset"]]))[, 1])
expected <- c(4.532599, 1.791759, 8.447843, 6.369901, 2.397895, 7.277248)
test_that("loge transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance=0.0001)
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
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(64.500000, 3.333333, 3040.166667, 383.916667, 7.083333, 979.000000)
test_that("limma batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME this is obviously broken
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limmaresid"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(0.02101530, 1.60048774, -0.07524254, 0.15555548, 0.49697157, 0.21066246)
test_that("limma-residuals batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combatmod"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(60.542865, 3.095141, 3032.546240, 355.354483, 6.666536, 912.345141)
test_that("combatmod from cbcbSEQ batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="sva"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(87.749155, 4.875168, 4418.793068, 557.089233, 9.641209, 1387.892080)
test_that("sva batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO Figure out what is up with these
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat")  ## broken
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior") ## takes forever

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combat_scale"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(71.646617, 2.267661, 3481.438297, 410.347768, 6.658058, 1031.946620)
test_that("combat_scale gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO Figure this guy out
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior_scale") ## takes forever

pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="svaseq"))[["result"]]
actual <- as.numeric(head(Biobase::exprs(pasilla_batch[["expressionset"]]))[, 1])
expected <- c(82.109044, 4.610009, 4099.038912, 519.407483, 9.116170, 1301.057485)
test_that("svaseq gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME I broke ruvg
## pasilla_batch <- normalize_expt(pasilla_expt, batch="ruvg") ## broken for the moment
