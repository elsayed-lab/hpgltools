start <- as.POSIXlt(Sys.time())
context("050annotation_txt.R:\n")
## 2017-12, exported functions in annotation_txt:
## load_trinotate_annotations() load_trinotate_go()

tmp <- system.file("sb/trinotate.csv.xz", package="hpgltools")
testing <- sm(load_trinotate_annotations(trinotate=tmp))
expected <- c(136152, 33)
actual <- dim(testing)
test_that("Do we get expected trinotate data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## FIXME: The creation of this go table is not quite right I think.
testing <- sm(load_trinotate_go(trinotate=tmp))
expecteed <- c(136152, 5)
actual <- dim(testing$go_data)
test_that("Do we get expected trinotate GO data?", {
  expect_equal(expected[1], actual[1])
##  expect_equal(expected[2], actual[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 050annotation_txt.R in ", elapsed,  " seconds."))
