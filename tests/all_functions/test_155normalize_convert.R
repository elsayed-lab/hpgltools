start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("155normalize_convert.R:
  12\n")
## 2018-04, exported functions in normalize_convert:
## convert_counts(), divide_seq(), hpgl_rpkm()

pombe_expt <- make_pombe_expt()

## I will mostly access these functions from normalize_expt, so first invoke them from there.
testing <- normalize_expt(pombe_expt, convert="cpm")
test_counts <- exprs(testing)

expected <- c(0.51068105, 1.46820803, 0.00000000, 0.06383513, 2.36189987,
              0.00000000, 0.00000000, 0.12767026, 0.00000000, 0.00000000)
actual <- as.numeric(test_counts[1:10, 1])
test_that("cpm modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <- normalize_expt(pombe_expt, convert="rpkm", column="cds_length")
test_counts <- exprs(testing)
expected <- c(0.09019446, NA, NA, 0.17159982, 2.72422131,
              0.00000000, 0.00000000, 0.15144752, 0.00000000, 0.00000000)
actual <- as.numeric(test_counts[1:10, 1])
test_that("rpkm modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

## Since I do not seem to have a usable fission BSgenome, I would ask the reader to
## look to ../travis/test_11norm_convert.R for the test of convert_counts()
## using pasilla_expt and the method 'cp_seq_m'.

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 155normalize_convert.R in ", elapsed,  " seconds."))
