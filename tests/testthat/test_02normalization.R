library(testthat)
library(hpgltools)

context("Are normalizations consistent over time?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## First make sure the pasilla_expt still has the stuff we expect
expected_title <- "This is an expt class."
actual_title <- pasilla_expt$title
expected_olibsize <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(expected_olibsize) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
actual_olibsize <- pasilla_expt$original_libsize
expected_expt <- c(92, 5, 4664, 583, 10, 1446)
actual_expt <- as.numeric(head(Biobase::exprs(pasilla_expt$expressionset[, 1])))
test_that("Does the expt contain everything we expect?", {
    expect_equal(expected_title, actual_title)
    expect_equal(expected_olibsize, actual_olibsize)
    expect_equal(expected_expt, actual_expt)
})

unmolested <- s_p(normalize_expt(pasilla_expt))$result

expected_expt <- Biobase::exprs(pasilla_expt$expressionset)
actual_expt <- Biobase::exprs(unmolested$expressionset)
expected_pexpt <- Biobase::pData(pasilla_expt$expressionset)
actual_pexpt <- Biobase::pData(unmolested$expressionset)
expected_libsize <- pasilla_expt$libsize
actual_libsize <- unmolested$libsize
test_that("If we use normalize_expt with no arguments, nothing should happen, yesno?", {
    expect_equal(expected_expt, actual_expt)
    expect_equal(expected_pexpt, actual_pexpt)
    expect_equal(expected_libsize, actual_libsize)
})

## First test conversions alone
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cpm"))$result
hpgl_cpm <- Biobase::exprs(pasilla_convert$expressionset)
edger_cpm <- edgeR::cpm(pasilla_expt$expressionset)
test_that("cpm conversions are equivalent?", {
    expect_equal(edger_cpm, hpgl_cpm)
})

pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="rpkm"))$result
lengths <- as.vector(Biobase::fData(pasilla_expt$expressionset)$length)
hpgl_rpkm <- Biobase::exprs(pasilla_convert$expressionset)
edger_rpkm <- edgeR::rpkm(x=pasilla_expt$expressionset, gene.length=lengths)
test_that("rpkm conversions are equivalent?", {
    expect_equal(edger_rpkm, hpgl_rpkm)
})

tt <- require.auto("BSgenome.Dmelanogaster.UCSC.dm6")
tt <- s_p(library("BSgenome.Dmelanogaster.UCSC.dm6"))
pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6))$result
expected_tas <- c(0.46719586, 0.03443909, 22.38432168, 54.32421464, 0.03908639, 106.58455139)
actual_tas <- as.numeric(head(Biobase::exprs(pasilla_convert$expressionset))[, 1])
test_that("cp_seq_m works for TA?", {
    expect_equal(expected_tas, actual_tas)
})

pasilla_convert <- s_p(normalize_expt(pasilla_expt, convert="cp_seq_m", genome=BSgenome.Dmelanogaster.UCSC.dm6, pattern="ATG"))$result
expected_atgs <- as.numeric(head(Biobase::exprs(pasilla_convert$expressionset))[,1])
actual_atgs <- c(0.51893853, 0.04536343, 27.76677691, 46.94320722, 0.05237078, 99.09109542)
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected_atgs, actual_atgs)
})

## Test normalizations
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="quant"))$result
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf"))$result
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="sf2"))$result
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="vsd"))$result
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qsmooth")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink_median")
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="tmm"))$result
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="upperquartile"))$result
pasilla_norm <- s_p(normalize_expt(pasilla_expt, norm="rle"))$result

## Test transformations
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log2"))$result
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log10"))$result
pasilla_trans <- s_p(normalize_expt(pasilla_expt, transform="log"))$result

## Test filter
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cbcb"))$result
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="pofa"))$result
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="kofa"))$result
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="cv"))$result
pasilla_filter <- s_p(normalize_expt(pasilla_expt, filter="simple"))$result

## Test batch
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limma"))$result
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="limmaresid"))$result
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combatmod"))$result
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="sva"))$result
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat")  ## broken
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior") ## takes forever
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="combat_scale"))$result
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior_scale") ## takes forever
pasilla_batch <- s_p(normalize_expt(pasilla_expt, batch="svaseq"))$result
## pasilla_batch <- normalize_expt(pasilla_expt, batch="ruvg") ## broken for the moment
