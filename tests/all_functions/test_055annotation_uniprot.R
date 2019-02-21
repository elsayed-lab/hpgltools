start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("055annotation_uniprot.R:
  1234\n")
## 2017-12, exported functions in annotation_gff:
## load_uniprot_annotations() download_uniprot_proteome()
website_up <- TRUE
## download_uniprot_proteome()
testing <- try(download_uniprot_proteome(
  species="Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)",
  first=TRUE))
if (class(testing)[1] == "try-error") {
  website_up <- FALSE
}

if (isTRUE(website_up)) {
  expected <- "UP000001584.txt.gz"
  actual <- as.character(testing[["filename"]])
  ## 0102
  test_that("Did we get the correct uniprot text file?", {
    expect_equal(expected, actual)
    expect_true(file.exists(expected))
  })

  ## load_uniprot_annotations()
  testing <- load_uniprot_annotations(file=expected)
  expected <- c(3993, 76)
  actual <- dim(testing)
  ## 0304
  test_that("Did we get the correct uniprot annotations?", {
    expect_equal(expected[1], actual[1])
    expect_equal(expected[2], actual[2])
  })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 055annotation_uniprot.R in ", elapsed,  " seconds."))
