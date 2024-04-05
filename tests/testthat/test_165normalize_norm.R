start <- as.POSIXlt(Sys.time())
context("165normalize_norm.R")
## 2017-12, exported functions in annotation_gff:

pombe_expt <- make_pombe_expt(annotation = FALSE)

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "quant")
test_counts <- exprs(testing)
expected <- c(14.27778, 23.72222, 109.58333, 11.50000, 61.94444,
              129.55556, 32.97222, 72.55556, 107.02778, 13.47222)
actual <- as.numeric(test_counts[1:10, 1])
## I am loosening this tolerance due to preprocessCore problems 20181210
test_that("Does quantile normalization return expected values?", {
  expect_equal(expected, actual, tolerance = 1.0)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "sf2")
test_counts <- exprs(testing)
expected <- c(14.54869, 23.40442, 98.04553, 12.01848, 57.56221,
              116.38953, 30.99504, 66.41794, 95.51532, 13.91614)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does size-factor normalization (DESeq2) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "sf")
test_counts <- exprs(testing)
expected <- c(14.54869, 23.40442, 98.04553, 12.01848, 57.56221,
              116.38953, 30.99504, 66.41794, 95.51532, 13.91614)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does size-factor normalization (DESeq) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.01)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "vsd")
test_counts <- exprs(testing)
expected <- c(6.298605, 6.508338, 7.462539, 6.226385, 7.045104,
              7.613672, 6.654783, 7.149607, 7.440215, 6.281225)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does vsd (DESeq2) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.01)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "tmm")
test_counts <- exprs(testing)
expected <- c(19.13119, 38.70352, 143.26251, 22.09686, 111.67109,
              266.38818, 52.58351, 135.67889, 194.09264, 24.84304)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does tmm (edgeR) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.1)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "upperquartile")
test_counts <- exprs(testing)
expected <- c(19.48496, 39.75666, 147.30954, 22.23285, 111.32201,
              260.75273, 50.90794, 130.31228, 185.40237, 24.39166)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does tmm (edgeR) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, filter = TRUE, norm = "rle")
test_counts <- exprs(testing)
expected <- c(18.81332, 38.83165, 143.44134, 22.54498, 113.47495,
              270.35551, 51.77898, 137.22436, 196.26130, 24.53660)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does rle (edgeR) return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 165normalize_norm.R in ", elapsed,  " seconds.")
