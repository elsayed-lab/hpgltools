library(testthat)
library(hpgltools)
context("Do gProfiler searches work?")

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ## Run your tests here
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sm(get_sig_genes(table, column="untreated")$up_genes)

    gprofiler_result <- sm(simple_gprofiler(sig_genes, species="dmelanogaster", first_col="untreated"))

    expected_go <- c(2.75e-03, 6.20e-04, 1.98e-06, 1.47e-04, 5.24e-04, 3.11e-05)
    actual_go <- head(gprofiler_result[["go"]][["p.value"]])
    expected_mfplot_data <- c(1.54e-06, 7.20e-06, 1.30e-05, 1.55e-05, 1.64e-05, 7.70e-05)
    actual_mfplot_data <- head(gprofiler_result$plots$mfp_plot_over$data$pvalue)
    expected_bpplot_data <- c(0.00062, 0.00275)
    actual_bpplot_data <- head(gprofiler_result$plots$bpp_plot_over$data$pvalue)
    expected_ccplot_data <- NULL
    actual_ccplot_data <- head(gprofiler_result$plots$cp_plot_over$data$pvalue)

    test_that("Does gprofiler return expected values?", {
        expect_equal(expected_go, actual_go, tolerance=0.1)
    })

    test_that("Does gprofiler return expected values? (mfpplot_data)", {
        expect_equal(expected_mfplot_data, actual_mfplot_data, tolerance=0.001)
    })

    test_that("Does gprofiler return expected values? (bppplot_data)", {
        expect_equal(expected_bpplot_data, actual_bpplot_data, tolerance=0.002)
    })

    test_that("Does gprofiler return expected values? (ccpplot_data)", {
        expect_equal(expected_ccplot_data, actual_ccplot_data, tolerance=0.001)
    })
}
