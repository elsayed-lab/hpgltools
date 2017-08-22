start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("55gsea_gostats.R: Does GOstats work?\n")

load("gsea_siggenes.rda")

gst_result <- sm(simple_gostats(fcp_sig_genes, gff_df=dmel, goids_df=dmel_ontologies,
                                gff_type="protein_coding"))
## There is some run-to-run variability in these ontology searches
expected <- c("GO:0000146", "GO:0000295", "GO:0001871",
              "GO:0003824", "GO:0003974", "GO:0003978")
actual <- head(sort(gst_result[["mf_over_enriched"]][["GOMFID"]]))
test_that("Are the GOstats interesting results as expected? (MF)", {
    expect_equal(expected, actual)
})

expected <- c("GO:0000422", "GO:0001508", "GO:0001676",
              "GO:0002118", "GO:0002121", "GO:0002218")
actual <- head(sort(gst_result[["bp_over_enriched"]][["GOBPID"]]))
test_that("Are the GOstats interesting results as expected? (BP)", {
    expect_equal(expected, actual)
})

expected <- c("GO:0005576", "GO:0005578", "GO:0005604",
              "GO:0005637", "GO:0005639", "GO:0005811")
actual <- head(sort(gst_result[["cc_over_enriched"]][["GOCCID"]]))
test_that("Are the GOstats interesting results as expected? (CC)", {
    expect_equal(expected, actual)
})

expected <- "gg"
actual <- class(gst_result[["pvalue_plots"]][["mfp_plot_over"]])[[1]]
test_that("Are the GOstats pvalue plots generated?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 55gsea_gostats.R in ", elapsed,  " seconds."))
tt <- clear_session()
