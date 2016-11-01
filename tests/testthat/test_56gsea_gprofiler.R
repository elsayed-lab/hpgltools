library(testthat)
library(hpgltools)
context("56gsea_gprofiler.R: Do gProfiler searches work?\n")

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

    expected <- c(6.57e-03, 1.62e-05, 2.80e-02, 8.38e-05, 1.37e-04, 1.10e-04)
    actual <- head(gprofiler_result[["go"]][["p.value"]])
    test_that("Does gprofiler return expected values?", {
        expect_equal(expected, actual, tolerance=0.1)
    })

    expected <- c(0.0264, 0.0214, 0.0187, 0.0183, 0.0471, 0.0145)
    actual <- head(gprofiler_result$plots$mfp_plot_over$data$pvalue)
    test_that("Does gprofiler return expected values? (mfpplot_data)", {
        expect_equal(expected, actual, tolerance=0.001)
    })

    expected <- c(0.01120, 0.02110, 0.00186, 0.02640, 0.01160, 0.00222)
    actual <- head(gprofiler_result$plots$bpp_plot_over$data$pvalue)
    test_that("Does gprofiler return expected values? (bppplot_data)", {
        expect_equal(expected, actual, tolerance=0.002)
    })

    expected <- NULL
    actual <- head(gprofiler_result$plots$cp_plot_over$data$pvalue)
    test_that("Does gprofiler return expected values? (ccpplot_data)", {
        expect_equal(expected, actual, tolerance=0.001)
    })
}

message("\nFinished 56gsea_gprofiler.R")
