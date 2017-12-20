start <- as.POSIXlt(Sys.time())
context("005annotation_biomart.R\n")
## 2017-12, exported functions in annotation_biomart:
##   load_biomart_annotations(), load_biomart_go(), load_biomart_orthologs()

## load_biomart_annotations()
testing <- sm(load_biomart_annotations())
gene_ids <- head(rownames(testing[["annotation"]]))
data <- testing[["annotation"]]
expected <- c(197995, 12)
actual <- dim(data)
test_that("Do we receive expected output from load_biomart_annotations()?", {
  expect_equal(expected, actual)
})

## load_biomart_go()
testing <- sm(load_biomart_go())
data <- testing[["go"]]
expected <- c(318558, 2)
actual <- dim(data)
test_that("Do we receive expected output from load_biomart_go()?", {
  expect_equal(expected, actual)
})

## load_biomart_orthologs()
testing <- sm(load_biomart_orthologs(gene_ids=gene_ids))
data <- head(testing[["all_gene_list"]])
hs_actual <- data[["hsapiens"]]
hs_expected <- c("ENSG00000135070", "ENSG00000113108", "ENSG00000212907",
                 "ENSG00000249115", "ENSG00000213918", "ENSG00000198888")
mm_actual <- data[["mmusculus"]]
mm_expected <- c("ENSMUSG00000078139", "ENSMUSG00000006050", "ENSMUSG00000065947",
                 "ENSMUSG00000078762", "ENSMUSG00000005980", "ENSMUSG00000064341")
test_that("Do we get expected orthologs from load_biomart_orthologs()?", {
  expect_equal(hs_expected, hs_actual)
  expect_equal(mm_expected, mm_actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 005annotation_biomart.R in ", elapsed,  " seconds."))
