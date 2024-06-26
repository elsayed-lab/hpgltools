start <- as.POSIXlt(Sys.time())
context("050annotation_txt.R")
## 2017-12, exported functions in annotation_txt:
## load_trinotate_annotations() load_trinotate_go()

tmp <- system.file("share/sb/trinotate_head.csv.xz", package = "hpgldata")
testing <- load_trinotate_annotations(trinotate = tmp)
## Moved rownames to a column for tibble, ergo 1 more column
## I recently changed this to drop columns which have no information, losing 3.
expected <- c(49999, 31)
actual <- dim(testing)
test_that("Do we get expected trinotate data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## FIXME: The creation of this go table is not quite right I think.
## I added some options to make this flexible over different versions of trinotate
testing <- load_trinotate_go(trinotate = tmp, blast2go_column = "gene_ontology_blast",
                             pfam_column = "gene_ontology_pfam")
## In addition, it now collapses isoforms of one putative gene and so
## drops the number of rows significantly.
## expecteed <- c(136152, 5)
expected <- c(13871, 5)
actual <- dim(testing[["go_data"]])
test_that("Do we get expected trinotate GO data?", {
  expect_equal(expected[1], actual[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 050annotation_txt.R in ", elapsed,  " seconds.")
