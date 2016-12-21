start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("40ann_biomartgenbank.R: Does biomart function?\n")

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

query_many <- translate_ids_querymany(actual_ids, species="dmelanogaster", from=NULL, fields=NULL)
expected <- c("yellow-e", "yellow-e", "yellow-e", "yellow-e", "CG7296", "CG18765")
actual <- query_many[["symbol"]]
test_that("Does queryMany return sensible outputs?", {
    expect_equal(expected, actual)
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

    gbk2txdb_test <- sm(gbk2txdb())
    expected <- 1895017
    actual <- GenomicRanges::width(gbk2txdb_test[["seq"]])
    test_that("The genbank txdb S.pyogenes genome's size is correct?", {
        expect_equal(expected, actual)
    })
    expected <- c("spyM18_0001", "spyM18_0002", "spyM18_0004",
                  "spyM18_0005", "spyM18_0007", "spyM18_0008")
    actual <- head(as.data.frame(gbk2txdb_test$cds)$locus_tag)
    test_that("The first few genbank S.pyogenes cds spyIDs downloaded?", {
        expect_equal(expected, actual)
    })

    gene_annotations <- gbk_annotations(gbk2txdb_test$txdb)
    expected <- c(1356, 1137, 1116, 570, 3504, 273)
    actual <- head(GenomicRanges::width(gene_annotations))
    test_that("I can extract the gene lengths of the first few S.pyogenes genes?", {
        expect_equal(expected, actual)
    })

    grabbed_accession <- sm(download_gbk())
    actual <- grabbed_accession[["strings"]][[1]]
    expected <- "LOCUS       AE009949             1895017 bp    DNA     circular BCT 31-JAN-2014"
    test_that("I downloaded an accession using download_gbk()?", {
        expect_equal(expected, actual)
    })
    test_that("The file exists?", {
        expect_true(file.exists("AE009949.gb"))
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 40ann_biomartgenbank.R in ", elapsed,  " seconds."))
