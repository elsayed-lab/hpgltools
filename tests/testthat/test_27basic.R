library(testthat)
library(hpgltools)

context("Does the basic differential expression analysis work?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

norm_expt <- s_p(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm"))$result

hpgl_pas_basic <- s_p(basic_pairwise(pasilla_expt))$result
hpgl_basic <- s_p(basic_pairwise(norm_expt))$result

test_that("Does a non-normalized basic run equal a normalized basic run?", {
    expect_equal(hpgl_pas_basic$all_tables$untreated_vs_treated, hpgl_basic$all_tables$untreated_vs_treated)
})

expected_medians_treated <- c(2.8356541, -4.0744928, 8.3510620, 5.0276022, -0.9340304, 6.2167638)
actual_medians_treated <- head(hpgl_basic[["medians"]][["treated"]])
expected_medians_untreated <- c(2.919061, -4.074470, 8.513323, 5.099943, -1.109891, 6.266102)
actual_medians_untreated <- head(hpgl_basic[["medians"]][["untreated"]])
expected_logfc <- c(8.341e-02, 2.234e-05, 1.623e-01, 7.234e-02, -1.759e-01, 4.934e-02)
actual_logfc <- head(hpgl_basic[["all_tables"]][["untreated_vs_treated"]][["logFC"]])

test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_treated, actual_medians_treated, tolerance=0.001)
    expect_equal(expected_medians_untreated, actual_medians_untreated, tolerance=0.001)
    expect_equal(expected_logfc, actual_logfc, tolerance=0.001)
})

save(list=ls(), file="de_basic.rda")
