start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("56gsea_gprofiler.R: Do gProfiler searches work?\n")

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.
## Run your tests here
load("gsea_siggenes.rda")

##gprofiler_result <- simple_gprofiler(z_sig_genes, species="dmelanogaster", first_col="untreated")
gprofiler_result <- sm(simple_gprofiler(z_sig_genes, species="dmelanogaster", first_col="logFC"))

expected <- c(2.44e-07, 4.08e-07, 3.69e-06, 1.62e-05, 1.62e-05, 2.13e-05)
actual <- head(sort(gprofiler_result[["go"]][["p.value"]]))
test_that("Does gprofiler return expected values?", {
    expect_equal(expected, actual, tolerance=0.1)
})

expected <- c(1.62e-05, 1.46e-03, 4.42e-03, 4.80e-03, 7.72e-03, 9.59e-03)
actual <- head(sort(gprofiler_result[["pvalue_plots"]][["mfp_plot_over"]][["data"]][["pvalue"]]))
test_that("Does gprofiler return expected values? (mfpplot_data)", {
    expect_equal(expected, actual, tolerance=0.01)
})

## When I run this in an interactive session it works.
## But when I run it with make test,  no!  WTF!?
expected <- c(0.00101, 0.01380, 0.02210)
actual <- head(sort(gprofiler_result[["pvalue_plots"]][["bpp_plot_over"]][["data"]][["pvalue"]]))
test_that("Does gprofiler return expected values? (bppplot_data)", {
    expect_equal(expected, actual, tolerance=0.03)
})

expected <- NULL
actual <- head(gprofiler_result[["pvalue_plots"]][["cp_plot_over"]][["data"]][["pvalue"]])
test_that("Does gprofiler return expected values? (ccpplot_data)", {
    expect_equal(expected, actual, tolerance=0.001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 56gsea_gprofiler.R in ", elapsed,  " seconds."))
tt <- clear_session()
