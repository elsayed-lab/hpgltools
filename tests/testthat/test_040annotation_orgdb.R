start_time <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("040annotation_orgdb.R")
## 2017-12, exported functions in annotation_orgdb:
## load_host_annotations()(likely deprecated)
## load_orgdb_annotations(), load_orgdb_go(),
## take_orgdb_from_ah(), map_orgdb_ids()

## Taken by doing head(keys(Homo.sapiens, keytype = "ENSEMBL"))
test_genes <- c("ENSG00000121410", "ENSG00000175899", "ENSG00000256069",
                "ENSG00000171428", "ENSG00000156006", "ENSG00000196136")

## load_orgdb_annotations()
testing <- load_orgdb_annotations()
expected <- 320000
actual <- nrow(testing[["genes"]])
## 0102
test_that("Do we get the expected amount of orgdb gene data?", {
  expect_gt(actual, expected)
})

start_positions <- testing[["genes"]][["ensembl"]]
test_ids <- head(unique(start_positions), n = 500)
guess <- guess_orgdb_keytype(ids = test_ids, orgdb = NULL)
test_that("Can we guess appropriate keytypes from gene IDs?", {
  expect_equal(guess, "ENSEMBL")
})

expected <- 82000
actual <- nrow(as.data.frame(testing[["transcripts"]]))
## 0304
test_that("Do we get the expected amount of orgdb transcript data?", {
  expect_gt(actual, expected)
})

## load_orgdb_go()
## Interesting, querying homo sapiens reminds me that we need to be more careful about which
## evidences we accept, as this table is astonishingly redundant.
testing <- load_orgdb_go(gene_ids = test_genes)
## Another function on which I get different answers on different hosts.
expected <- 8100
actual <- nrow(testing)
## 0506
test_that("Do we get the expected amount of orgdb GO data?", {
  expect_gt(actual, expected)
})

## take_from_ah()
testing <- orgdb_from_ah()
expected <- 26
actual <- length(AnnotationDbi::keytypes(testing))
## 07
test_that("Do we get the expected keytypes from an ah orgdb?", {
  expect_equal(expected, actual)
})

## map_orgdb_ids()
testing <- map_orgdb_ids(orgdb = testing, keytype = "entrezid")
expected <- c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419",
              "ENSG00000000457", "ENSG00000000460", "ENSG00000000938")
actual <- head(sort(testing[["ensembl"]]))
## 08
test_that("Do we get the expected ensembl IDs from entrez IDs?", {
  expect_equal(expected, actual)
})

end_time <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end_time - start_time))
message("\nFinished 040annotation_orgdb.R in ", elapsed,  " seconds.")
