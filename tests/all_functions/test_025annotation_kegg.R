start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("025annotation_kegg.R
  1234567\n")
## 2017-12, exported functions in annotation_kegg:
## load_kegg_annotations(), map_kegg_to_ensembl()
## There are some functions in ontology_kegg which probably should be moved here.

## load_kegg_annotations()
test_kegg <- load_kegg_annotations()
actual <- head(test_kegg[["GID"]])
expected <- c("b0001", "b0002", "b0003", "b0004", "b0005", "b0006")
## 01
test_that("Do we get the expected KEGG GIDs?", {
  expect_equal(expected, actual)
})
kegg_ids <- head(test_kegg[["kegg_geneid"]])

actual <- head(test_kegg[["GID"]])
expected <- c("b0001", "b0002", "b0003", "b0004", "b0005", "b0006")
## 02
test_that("Do we get the expected ncbi gene IDs?", {
  expect_equal(expected, actual)
})

actual <- head(test_kegg[["ncbi_proteinid"]])
expected <- c("NP_414542", "NP_414543", "NP_414544",
              "NP_414545", "NP_414546", "NP_414547")
## 03
test_that("Do we get the expected ncbi protein IDs?", {
  expect_equal(expected, actual)
})

actual <- head(test_kegg[["uniprotid"]])
expected <- c("P0AD86, A0A387D2P5", "P00561, A0A387CVZ1", "P00547, A0A387D0L4",
              "P00934, A0A387CXN6", "P75616, A0A387D2A3", "P0A8I3, A0A387DCJ9")
## 04
test_that("Do we get the expected uniprot IDs?", {
  expect_equal(expected, actual)
})

actual <- head(test_kegg[["pathways"]])
expected <- c(
  "", "eco00260, eco00261, eco00270, eco00300, eco01100, eco01110, eco01120, eco01130, eco01230",
  "eco00260, eco01100, eco01110, eco01120, eco01230",
  "eco00260, eco00750, eco01100, eco01110, eco01120, eco01230", "", "")
## 05
test_that("Do we get the expected KEGG pathways?", {
  expect_equal(expected, actual)
})

## map_kegg_ensembl()
mapped <- map_kegg_dbs(kegg_ids)
expected <- c(6, 7)
actual <- dim(mapped)
## 0607
test_that("Do we get the expected db mapping size?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 025annotation_kegg.R in ", elapsed,  " seconds."))
