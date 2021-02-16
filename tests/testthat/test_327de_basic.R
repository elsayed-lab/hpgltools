start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("327de_basic: Does the basic differential expression analysis work?\n")

pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]

norm_expt <- normalize_expt(pasilla_expt, transform = "log2",
                            norm = "quant", filter = TRUE,
                            convert = "cbcbcpm")

hpgl_pasilla_basic <- basic_pairwise(pasilla_expt)
hpgl_norm_basic <- basic_pairwise(norm_expt)

expected <- hpgl_pasilla_basic[["all_tables"]][[1]][["logFC"]]
actual <- hpgl_norm_basic[["all_tables"]][[1]][["logFC"]]
test_that("Does a non-normalized basic run equal a normalized basic run?", {
  expect_equal(expected, actual, tolerance = 0.1)
})

## Another casuality of normalize.quantiles vs normalize.quantiles.robust
##expected_medians_treated <- c(-3.008634, -3.008632, -2.487802, -2.487802, -1.750838, -1.509827)
expected_medians_treated <- c(-2.318990, -2.297095, -2.123484, -2.069237, -1.384643, -1.184830)
actual_medians_treated <- head(sort(hpgl_norm_basic[["medians"]][["treated"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_treated, actual_medians_treated, tolerance = 0.001)
})

##expected_medians_untreated <- c(-1.9991561, -1.7320136, -1.2739814,
##                                -1.0475990, -0.9636259, -0.7944623)
expected_medians_untreated <- c(-1.9202939, -1.2937864, -1.2228098,
                                -1.1269285, -1.0597740, -0.9812334)
actual_medians_untreated <- head(sort(hpgl_norm_basic[["medians"]][["untreated"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_untreated, actual_medians_untreated, tolerance = 0.01)
})

expected_logfc <- c(-3.307, -3.136, -3.061, -2.997, -2.961, -2.948)
actual_logfc <- head(sort(hpgl_norm_basic[["all_tables"]][["untreated_vs_treated"]][["logFC"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_logfc, actual_logfc, tolerance = 0.1)
})

basic_written <- write_basic(hpgl_norm_basic, excel = "basic_test.xlsx")
test_that("Is it possible to write the results of a basic analysis?", {
    expect_true(file.exists("basic_test.xlsx"))
})

hpgl_basic <- sm(basic_pairwise(pasilla_expt))
save(list = ls(), file = "327_de_basic.rda")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 327de_basic.R in ", elapsed,  " seconds.")
