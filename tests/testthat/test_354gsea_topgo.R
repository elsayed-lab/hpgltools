start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("54gsea_topgo.R: Does topGO work?\n")

sig_file <- "351_gsea_siggenes.rda"
if (file.exists(sig_file)) {
  load(sig_file)
} else {
  stop("The significance file.")
}

tp_result <- simple_topgo(fcp_sig_genes, go_db=dmel_ontologies, overwrite=TRUE,
                          excel="topgo.xlsx", pval_column="adj.P.Val")
test_that("Did we get an excel output?", {
  expect_true(file.exists("topgo.xlsx"))
})

## There is some run-to-run variability in these searches.
expected <- c("GO:0000146", "GO:0000295", "GO:0000981",
              "GO:0003824", "GO:0003974", "GO:0003978")
actual <- head(sort(tp_result[["tables"]][["mf_interesting"]][["GO.ID"]]))
test_that("Are the topGO interesting results expected (MF GOIDs)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0000422", "GO:0001508", "GO:0001676",
              "GO:0002118", "GO:0002121", "GO:0002218")
actual <- head(sort(tp_result[["tables"]][["bp_interesting"]][["GO.ID"]]))
test_that("Are the topGO interesting results expected (BP GOIDs)?", {
    expect_equal(expected, actual)
})

##expected <- c("GO:0005576", "GO:0005578", "GO:0005604",
##              "GO:0005637", "GO:0005639", "GO:0005811")
expected <- c("GO:0005604", "GO:0005637", "GO:0005639",
              "GO:0005811", "GO:0005859", "GO:0005967")
actual <- head(sort(tp_result[["tables"]][["cc_interesting"]][["GO.ID"]]))
test_that("Are the topGO interesting results expected (CC GOIDs)?", {
    expect_equal(expected, actual)
})

tp_trees <- sm(topgo_trees(tp_result))
expected <- c("GO:0003674", "GO:0003824", "GO:0003974",
              "GO:0004467", "GO:0004553", "GO:0004556")
actual <-  head(sort(tp_trees[["mf_fisher_nodes"]][["dag"]]@nodes))
test_that("Are the topGO interesting results as expected? (MF trees)?", {
    expect_equal(expected, actual)
})

## I reordered these plots, these values are no longer valid.
##expected <- c("GO:0003974", "GO:0004467", "GO:0004556",
##              "GO:0004742", "GO:0004793", "GO:0005044")
expected <- c("GO:0000146", "GO:0003974", "GO:0003978",
              "GO:0004467", "GO:0004556", "GO:0004591")
actual <- head(sort(tp_result[["pvalue_plots"]][["mfp_plot_over"]][["data"]][["GO.ID"]]))
test_that("Are the topGO interesting results as expected? (MF pval)?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 54gsea_topgo.R in ", elapsed,  " seconds."))
