library(testthat)
library(hpgltools)

context("Does biomart function?")

if (!identical(Sys.getenv("TRAVIS"), "true")) {

    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sp(get_sig_genes(table, column="untreated")$up_genes)$result
    dmel_annotations <- sp(get_biomart_annotations(species="dmelanogaster"))$result
    dmel_go <- sp(get_biomart_ontologies(species="dmelanogaster"))$result

    expected_lengths <- c(1776, 819, 2361, NA, 633, 1164)
    actual_lengths <- head(dmel_annotations$length)
    test_that("Did the gene lengths come out?", {
        expect_equal(expected_lengths, actual_lengths)
    })

    expected_ids <- c("FBgn0041711", "FBgn0041711", "FBgn0041711", "FBgn0041711", "FBgn0032283", "FBgn0042110")
    actual_ids <- head(dmel_go$ID)
    expected_go <- c("GO:0005576", "GO:0048067", "GO:0016853", "GO:0042438", "", "GO:0016772")
    actual_go <- head(dmel_go$GO)
    test_that("Did the ontologies come out?", {
        expect_equal(expected_ids, actual_ids)
        expect_equal(expected_go, actual_go)
    })

}
