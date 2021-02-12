start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("175normalize_transform.R:
  12345\n")
## 2017-12, exported functions in normalize_transform:
## transform_counts(), this function has switches for: rounding, voom, weighted voom, log2, log10

pombe_expt <- make_pombe_expt()
test_genes <- c("SPAC212.11", "SPAC212.09c", "SPNCRNA.70", "SPAC212.12", "SPAC212.04c",
                "SPAC212.03", "SPAC212.02", "SPAC212.01c", "SPAC977.03", "SPAC977.04")
testing <- normalize_expt(pombe_expt, transform = "log2")
test_counts <- exprs(testing)
expected <- c(3.169925, 4.584963, 0.000000, 1.000000, 5.247928,
              0.000000, 0.000000, 1.584963, 0.000000, 0.000000)
names(expected) <- test_genes
actual <- test_counts[test_genes, 1]
test_that("Does log2 transformation return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, transform = "log10")
test_counts <- exprs(testing)
expected <- c(0.9542425, 1.3802112, 0.0000000, 0.3010300, 1.5797836,
              0.0000000, 0.0000000, 0.4771213, 0.0000000, 0.0000000)
names(expected) <- test_genes
actual <- test_counts[test_genes, 1]
test_that("Does log10 transformation return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, transform = "log")
test_counts <- exprs(testing)
expected <- c(2.1972246, 3.1780538, 0.0000000, 0.6931472, 3.6375862,
              0.0000000, 0.0000000, 1.0986123, 0.0000000, 0.0000000)
names(expected) <- test_genes
actual <- test_counts[test_genes, 1]
test_that("Does loge transformation return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, transform = "voom")
test_counts <- exprs(testing)
expected <- c(-0.8820428, 0.5850832, -4.9695057, -3.3845431, 1.2593130,
              -4.9695057, -4.9695057, -2.6475776, -4.9695057, -4.9695057)
names(expected) <- test_genes
actual <- test_counts[test_genes, 1]
test_that("Does voom transformation return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

testing <- normalize_expt(pombe_expt, transform = "voomweight")
test_counts <- exprs(testing)
expected <- c(-0.8820428, 0.5850832, -4.9695057, -3.3845431, 1.2593130,
              -4.9695057, -4.9695057, -2.6475776, -4.9695057, -4.9695057)
names(expected) <- test_genes
actual <- test_counts[test_genes, 1]
test_that("Does weighted voom transformation return expected values?", {
  expect_equal(expected, actual, tolerance = 0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 175normalize_transform.R in ", elapsed,  " seconds."))
