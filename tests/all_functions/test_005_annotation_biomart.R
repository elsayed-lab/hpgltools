start <- as.POSIXlt(Sys.time())
context("005annotation_biomart.R\n")
## 2017-12, exported functions in annotation_biomart:
##   load_biomart_annotations(), load_biomart_go(), load_biomart_orthologs()

## load_biomart_annotations()
testing <- load_biomart_annotations()
gene_ids <- head(rownames(testing[["annotation"]]))
data <- testing[["annotation"]]
expected <- c(197995, 12)
actual <- dim(data)
test_that("Do we receive expected output from load_biomart_annotations()?", {
  expect_equal(expected, actual)
})

## load_biomart_go()
testing <- load_biomart_go()
data <- testing[["go"]]
expected <- c(318558, 2)
actual <- dim(data)
test_that("Do we receive expected output from load_biomart_go()?", {
  expect_equal(expected, actual)
})

## load_biomart_orthologs()
## Oh yeah, I moved the default biomart to hg38/90 or 89.  Thus these are unlikely
## to be correct.  I should just query the number of orthologs found.
testing <- load_biomart_orthologs(gene_ids=gene_ids)
data <- testing[["all_linked_genes"]]
actual <- dim(data)
expected <- c(26635, 3)
test_that("Do we get expected orthologs from load_biomart_orthologs()?", {
  expect_equal(expected[1], actual[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 005annotation_biomart.R in ", elapsed,  " seconds."))
