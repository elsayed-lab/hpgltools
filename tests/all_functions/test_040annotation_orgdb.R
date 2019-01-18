start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("040annotation_orgdb.R
  12345678\n")
## 2017-12, exported functions in annotation_orgdb:
## load_host_annotations()(likely deprecated)
## load_orgdb_annotations(), load_orgdb_go(),
## take_orgdb_from_ah(), map_orgdb_ids()

## Taken by doing head(keys(Homo.sapiens, keytype="ENSEMBL"))
test_genes <- c("ENSG00000121410", "ENSG00000175899", "ENSG00000256069",
                "ENSG00000171428", "ENSG00000156006", "ENSG00000196136")

## load_orgdb_annotations()
testing <- load_orgdb_annotations()
## hmm, I get different answers here on different hosts.
##expected <- c(327565, 6)
expected <- c(334845, 6)
actual <- dim(testing[["genes"]])
## 0102
test_that("Do we get the expected amount of orgdb gene data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})
expected <- c(82960, 7)
actual <- dim(as.data.frame(testing[["transcripts"]]))
## 0304
test_that("Do we get the expected amount of orgdb transcript data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## Wow, even querying only 1 gene seems to take too long, I must be doing something wrong.
## testing <- load_orgdb_annotations(fields="all", gene_ids=test_genes[1])

## load_orgdb_go()
## Interesting, querying homo sapiens reminds me that we need to be more careful about which
## evidences we accept, as this table is astonishingly redundant.
testing <- load_orgdb_go(gene_ids=test_genes)
##expected <- c(13627, 10)
## Another function on which I get different answers on different hosts.
expected <- c(13062, 10)
actual <- dim(testing)
## 0506
test_that("Do we get the expected amount of orgdb GO data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## Looks like a recent glibc change is causing this to segfault.
## take_from_ah()
testing <- orgdb_from_ah()
expected <- 26
actual <- length(AnnotationDbi::keytypes(testing))
## 07
test_that("Do we get the expected keytypes from an ah orgdb?", {
  expect_equal(expected, actual)
})
## map_orgdb_ids()
testing <- map_orgdb_ids(orgdb=testing, keytype="entrezid")
expected <- c("ENSG00000121410", "ENSG00000175899", "ENSG00000256069",
              "ENSG00000171428", "ENSG00000156006", NA)
actual <- head(testing[["ensembl"]])
## 08
test_that("Do we get the expected ensembl IDs from entrez IDs?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 040annotation_orgdb.R in ", elapsed,  " seconds."))
