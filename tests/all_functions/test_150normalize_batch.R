start <- as.POSIXlt(Sys.time())
context("150normalize_batch.R:\n")
## 2017-12, exported functions in annotation_gff:

pombe_expt <- make_pombe_expt()

testing <- normalize_expt(pombe_expt, filter=TRUE, batch="limma")
test_counts <- exprs(testing)

expected <- c(27.21933, 33.08824, 147.85699, 19.53751, 89.46974,
              206.14546, 39.54991, 117.83365, 139.79517, 20.77142)
actual <- as.numeric(test_counts[1:10, 1])
test_that("limma batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 150normalize_batch.R in ", elapsed,  " seconds."))
