start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("40ann_biomart.R: Does biomart function?\n")

limma <- new.env()
load("de_limma.rda", envir=limma)
table <- limma$hpgl_table
sig_genes <- sm(get_sig_genes(table, column="untreated")$up_genes)
dmel_annotations <- sm(get_biomart_annotations(species="dmelanogaster"))
dmel_go <- sm(get_biomart_ontologies(species="dmelanogaster"))

expected_lengths <- c(1776, 819, 2361, NA, 633, 1164)
actual_lengths <- head(dmel_annotations$length)
test_that("Did the gene lengths come out?", {
    expect_equal(expected_lengths, actual_lengths)
})

expected_ids <- c("FBgn0041711", "FBgn0041711", "FBgn0041711", "FBgn0041711", "FBgn0032283", "FBgn0042110")
actual_ids <- head(dmel_go$ID)
expected_go <- c("GO:0005576", "GO:0048067", "GO:0016853", "GO:0042438", "", "GO:0016772")
actual_go <- head(dmel_go$GO)
test_that("Did the ontologies come out (ids)?", {
    expect_equal(expected_ids, actual_ids)
})
test_that("Did the ontologies come out (go)?", {
    expect_equal(expected_go, actual_go)
})


if (!identical(Sys.getenv("TRAVIS"), "true")) {
    test_genes <- head(rownames(sig_genes))
    linkage_test <- biomart_orthologs(test_genes, first_species="dmelanogaster",
                                      second_species="mmusculus",
                                      first_attributes=c("ensembl_gene_id"),
                                      second_attributes=c("ensembl_gene_id","hgnc_symbol"))
    linked_genes <- linkage_test$linked_genes
    expected_linkage <- c("ENSMUSG00000025815", "ENSMUSG00000033006",
                          "ENSMUSG00000024176", "ENSMUSG00000000567")
    actual_linkage <- linked_genes[["mmusculus"]]
    test_that("Can I link some melanogaster and mouse genes?", {
        expect_equal(expected_linkage, actual_linkage)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 40ann_biomart.R in ", elapsed,  " seconds."))
