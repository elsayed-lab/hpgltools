start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("353gsea_clusterp.R: Does clusterProfiler work?\n")

sig_file <- "351_gsea_siggenes.rda"
if (file.exists(sig_file)) {
  load(sig_file)
} else {
  stop("The significance file.")
}

dmel_cp <- simple_clusterprofiler(
    sig_genes = z_sig_genes,
    de_table = table,
    do_david = FALSE,
    orgdb = "org.Dm.eg.db")

expected <- 90
actual <- nrow(dmel_cp[["enrich_go"]][["MF_all"]])
test_that("Does the set of MF_all have the expected IDs?", {
    expect_gt(actual, expected)
})

actual <- length(dmel_cp[["enrich_go"]][["MF_all"]][["ID"]])
test_that("Do we get a similar number of MF ids?", {
  expect_gt(actual, 85)
})

expected <- c("GO:0006022", "GO:0006030", "GO:0006040",
              "GO:0006208", "GO:0006812", "GO:0006814")
actual <- head(sort(dmel_cp[["enrich_go"]][["BP_all"]][["ID"]]))
test_that("Does the set of BP_all have the expected IDs?", {
    expect_equal(expected, actual)
})

actual <- length(dmel_cp[["enrich_go"]][["BP_all"]][["ID"]])
test_that("Do we get a similar number of BP ids?", {
  expect_gt(actual, 10)
})

expected <- c("GO:0005887", "GO:0031012", "GO:0031226")
actual <- head(sort(dmel_cp[["enrich_go"]][["CC_all"]][["ID"]]))
test_that("Does the set of CC_all have the expected IDs?", {
    expect_equal(expected, actual)
})

## Gather some scores
expected <- c(9.448405e-07, 9.448405e-07, 5.461374e-06,
              5.461374e-06, 5.253210e-04, 1.071284e-03)
actual <- head(sort(dmel_cp[["enrich_go"]][["MF_sig"]][["p.adjust"]]))
test_that("Does the set of MF_sig have the expected p.adjusts?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

expected <- c(1.318702e-07, 1.318702e-07, 2.451734e-04)
actual <- head(sort(dmel_cp[["enrich_go"]][["CC_sig"]][["p.adjust"]]))
test_that("Does the set of CC_sig have the expected p.adjusts?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

expected <- 8
actual <- nrow(dmel_cp[["kegg_data"]][["kegg_sig"]])
test_that("Did cp pick up consistent KEGG categories?", {
    expect_gt(actual, expected)
})

expected <- 500
actual <- nrow(dmel_cp[["group_go"]][["BP"]])
test_that("Does cp find consistent group_go categories in biological processes?", {
    expect_gt(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 53gsea_clusterp.R in ", elapsed,  " seconds."))
