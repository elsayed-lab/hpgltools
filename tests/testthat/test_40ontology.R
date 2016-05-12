library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

expected_term_result <- "adenyl ribonucleotide binding"
actual_term_result <- as.character(goterm())
expected_syn_result <- "mitochondrial inheritance"
actual_syn_result <- as.character(gosyn())
expected_sec_result <- c("GO:0000141", "GO:0030482")
actual_sec_result <- as.character(gosec())
expected_def_result <- "An assembly of actin filaments that are on the same axis but may be oriented with the same or opposite polarities and may be packed with different levels of tightness."
actual_def_result <- as.character(godef())
expected_ont_result <- c("CC","CC")
actual_ont_result <- as.character(goont())
expected_level_result <- c("3", "3")
actual_level_result <- as.character(golevel())

test_that("Are GO.db functions working?", {
    expect_equal(expected_term_result, actual_term_result)
    expect_equal(expected_syn_result, actual_syn_result)
    expect_equal(expected_sec_result, actual_sec_result)
    expect_equal(expected_def_result, actual_def_result)
    expect_equal(expected_ont_result, actual_ont_result)
    expect_equal(expected_level_result, actual_level_result)
})


