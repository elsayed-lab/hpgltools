start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("55gsea_gostats.R: Does GOstats work?\n")

load("gsea_siggenes.rda")

gst_result <- sm(simple_gostats(fcp_sig_genes, gff_df=dmel,
                                go_db=dmel_ontologies,
                                gff_type="protein_coding"))

## There is some run-to-run variability in these ontology searches
expected <- 70
actual <- nrow(gst_result[["tables"]][["mf_over_enriched"]])
test_that("Are the GOstats interesting results as expected? (MF)", {
    expect_gt(actual, expected)
})

expected <- 125
actual <- nrow(gst_result[["tables"]][["bp_over_enriched"]])
test_that("Are the GOstats interesting results as expected? (BP)", {
    expect_gt(actual, expected)
})

expected <- "gg"
actual <- class(gst_result[["pvalue_plots"]][["mfp_plot_over"]])[[1]]
test_that("Are the GOstats pvalue plots generated?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 55gsea_gostats.R in ", elapsed,  " seconds."))
