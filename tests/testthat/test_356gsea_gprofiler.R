start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("56gsea_gprofiler.R: Do gProfiler searches work?
  12345\n")

sig_file <- "351_gsea_siggenes.rda"
if (file.exists(sig_file)) {
  load(sig_file)
} else {
  stop("The significance file.")
}

##gprofiler_result <- simple_gprofiler(z_sig_genes, species = "dmelanogaster", first_col = "untreated")
gprofiler_result <- simple_gprofiler(z_sig_genes, species = "dmelanogaster", first_col = "logFC")

expected <- c(2.44e-07, 4.08e-07, 3.69e-06, 1.62e-05, 1.62e-05, 2.13e-05)
actual <- head(sort(gprofiler_result[["go"]][["p.value"]]))
test_that("Does gprofiler return expected values?", {
    expect_equal(expected, actual, tolerance = 0.1)
})

expected <- c(1.62e-05, 1.46e-03, 4.42e-03, 4.80e-03, 7.72e-03, 9.59e-03)
actual <- head(sort(gprofiler_result[["pvalue_plots"]][["mfp_plot_over"]][["data"]][["pvalue"]]))
test_that("Does gprofiler return expected values? (mfpplot_data)", {
    expect_equal(expected, actual, tolerance = 0.01)
})

## 20181003: gprofiler updated their searches, they now return 2 more hits.
expected <- c(0.00151, 0.00226, 0.00275, 0.01050, 0.01790, 0.01850)
actual <- head(sort(gprofiler_result[["pvalue_plots"]][["bpp_plot_over"]][["data"]][["pvalue"]]))
test_that("Does gprofiler return expected values? (bppplot_data)", {
    expect_equal(expected, actual, tolerance = 0.03)
})

expected <- NULL
actual <- head(gprofiler_result[["pvalue_plots"]][["cp_plot_over"]][["data"]][["pvalue"]])
test_that("Does gprofiler return expected values? (ccpplot_data)", {
    expect_equal(expected, actual, tolerance = 0.001)
})

actual <- write_gprofiler_data(gprofiler_result, excel = "test_gprofiler.xlsx")
expected <- 1
test_that("Can we write a gprofiler xlsx file?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 56gsea_gprofiler.R in ", elapsed,  " seconds.")
