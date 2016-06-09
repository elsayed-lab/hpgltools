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

unmolested <- normalize_expt(pasilla_expt)

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
cpm_expt <- s_p(normalize_expt(pasilla_expt, convert="cpm"))$result
hpgl_cpm <- Biobase::exprs(cpm_expt$expressionset)
edger_cpm <- edgeR::cpm(pasilla_expt$expressionset)

rpkm_expt <- normalize_expt(pasilla_expt, convert="rpkm")
lengths <- as.vector(Biobase::fData(pasilla_expt$expressionset)$length)
hpgl_rpkm <- Biobase::exprs(rpkm_expt$expressionset)
edger_rpkm <- edgeR::rpkm(x=pasilla_expt$expressionset, gene.length=lengths)

test_that("cpm/rpkm conversions are equivalent?", {
    expect_equal(edger_cpm, hpgl_cpm)
    expect_equal(edger_rpkm, hpgl_rpkm)
})
