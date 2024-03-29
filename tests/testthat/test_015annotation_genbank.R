start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("015annotation_genbank.R")
## 2017-12, exported functions in annotation_genbank:
##   load_genbank_annotations() gbk_annotations() download_gbk()

## load_genbank_annotations(), the others are called inside it.
## Something changed in genbankr, this now fails.
testing <- load_genbank_annotations()
other_data <- summary(testing[["others"]])
## 01
test_that("Do we get some granges?", {
  expect_equal("GRanges object with 79 ranges and 5 metadata columns", other_data)
})

exon_data <- summary(testing[["exons"]])
## 02
test_that("Do we get some granges?", {
  expect_equal("GRanges object with 1845 ranges and 12 metadata columns", exon_data)
})

cds_data <- summary(testing[["cds"]])
## 03
test_that("Do we get some granges?", {
  expect_equal("GRanges object with 1845 ranges and 12 metadata columns", cds_data)
})

genes_data <- summary(testing[["genes"]])
## 04
test_that("Do we get some granges?", {
  expect_equal("GRanges object with 1857 ranges and 5 metadata columns", genes_data)
})

txdb_data <- load_orgdb_annotations(testing[["txdb"]], keytype = "TXID",
                                    fields = c("CDSID", "CDSNAME", "EXONID", "EXONNAME",
                                             "GENEID", "TXID", "TXNAME"))
txdb_transcripts <- as.data.frame(txdb_data[["transcripts"]])
## 05
test_that("Do we get some transcript sequence names?", {
  expect_equal(as.character(head(txdb_transcripts[["seqnames"]])),
               c("MGAS8232", "MGAS8232", "MGAS8232", "MGAS8232", "MGAS8232", "MGAS8232"))
})

## 06
test_that("Do we get some transcript starts?", {
  expect_equal(as.numeric(head(txdb_transcripts[["start"]])),
               c(202,1709, 3447, 4632, 5204, 8869))
})

## 07
test_that("Do we get some transcript ends?", {
  expect_equal(as.numeric(head(txdb_transcripts[["end"]])),
               c(1557, 2845, 4562, 5201, 8707, 9141))
})

## 08
test_that("Do we get some transcript widths?", {
  expect_equal(as.numeric(head(txdb_transcripts[["width"]])),
               c(1356, 1137, 1116, 570, 3504, 273))
})

## 09
test_that("Do we get some transcript txnames?", {
  expect_equal(as.character(head(txdb_transcripts[["tx_name"]])),
               c("spyM18_0001.1", "spyM18_0002.1", "spyM18_0004.1",
                 "spyM18_0005.1", "spyM18_0007.1", "spyM18_0008.1"))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 015annotation_genbank.R in ", elapsed,  " seconds.")
