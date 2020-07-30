start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("40ann_biomartgenbank.R: Does biomart function?
  1234567\n")

limma <- new.env()
load("de_limma.rda", envir=limma)
limma_result <- limma[["hpgl_limma"]]
table <- limma_result[["all_tables"]][[1]]
sig_genes <- get_sig_genes(table, column="logFC")[["up_genes"]]
dmel_annotations <- load_biomart_annotations(species="dmelanogaster",
                                             host="useast.ensembl.org",
                                             overwrite=TRUE)
dmel_go <- load_biomart_go(species="dmelanogaster",
                           ## Yay archive.ensembl is back!
                           ## host="useast.ensembl.org",
                           overwrite=TRUE)

expected_lengths <- c(1776, 819, 2361, NA, 633, 1164)
actual_lengths <- head(dmel_annotations[["annotation"]][["cds_length"]])
## 01
test_that("Did the gene lengths come out?", {
    expect_equal(expected_lengths, actual_lengths)
})

## I am not sure why, but these tests are failing when I do not run them interactively...
## But when I run them myself, no problems...
expected_ids <- c(
    "FBgn0041711", "FBgn0041711", "FBgn0041711",
    "FBgn0041711", "FBgn0032283", "FBgn0042110")
actual_ids <- head(dmel_go[["go"]][["ID"]])
expected_go <- c(
    "GO:0005576", "GO:0048067", "GO:0016853",
    "GO:0042438", "", "GO:0016772")
actual_go <- head(dmel_go[["go"]][["GO"]])
## 0203
test_that("Did the ontologies come out?", {
  expect_equal(expected_ids, actual_ids)
  expect_equal(expected_go, actual_go)
})

test_genes <- head(rownames(sig_genes))
linkage_test <- load_biomart_orthologs(test_genes,
                                       ## host="useast.ensembl.org",
                                       first_species="dmelanogaster",
                                       second_species="mmusculus",
                                       attributes="ensembl_gene_id")
linked_genes <- linkage_test[["subset_linked_genes"]]
## Hard-set the order of the genes to avoid test failures for silly reasons.
actual_linkage <- sort(linked_genes[["mmusculus"]])
expected_linkage <- c(
    "ENSMUSG00000000567", "ENSMUSG00000024176",
    "ENSMUSG00000025815", "ENSMUSG00000033006")
## 04
test_that("Can I link some melanogaster and mouse genes?", {
    expect_equal(expected_linkage, actual_linkage)
})

gbk2txdb_test <- load_genbank_annotations(accession="AE009949")
expected <- 1895017
actual <- GenomicRanges::width(gbk2txdb_test[["seq"]])
## 05
test_that("The genbank txdb S.pyogenes genome's size is correct?", {
    expect_equal(expected, actual)
})

expected <- c("spyM18_0001", "spyM18_0002", "spyM18_0004",
              "spyM18_0005", "spyM18_0007", "spyM18_0008")
actual <- head(as.data.frame(gbk2txdb_test[["cds"]])[["locus_tag"]])
## 06
test_that("The first few genbank S.pyogenes cds spyIDs downloaded?", {
    expect_equal(expected, actual)
})

gene_annotations <- gbk_annotations(gbk2txdb_test[["txdb"]])
expected <- c(1356, 1137, 1116, 570, 3504, 273)
actual <- head(GenomicRanges::width(gene_annotations))
## 07
test_that("I can extract the gene lengths of the first few S.pyogenes genes?", {
    expect_equal(expected, actual)
})

grabbed_accession <- download_gbk()
actual <- grabbed_accession[["strings"]][[1]][1]
expected <- "LOCUS       AE009949             1895017 bp    DNA     circular BCT 31-JAN-2014"
## 0809
test_that("I downloaded an accession using download_gbk()?", {
    expect_equal(expected, actual)
    expect_true(file.exists("AE009949.gb"))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 40ann_biomartgenbank.R in ", elapsed,  " seconds."))
