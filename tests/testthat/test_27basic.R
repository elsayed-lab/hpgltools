library(testthat)
library(hpgltools)
context("Does the basic differential expression analysis work?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

norm_expt <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm"))

hpgl_pas_basic <- sm(basic_pairwise(pasilla_expt))
hpgl_basic <- sm(basic_pairwise(norm_expt))

test_that("Does a non-normalized basic run equal a normalized basic run?", {
    expect_equal(hpgl_pas_basic[["all_tables"]][["untreated_vs_treated"]],
                 hpgl_basic[["all_tables"]][["untreated_vs_treated"]])
})


expected_medians_treated <- c(-4.074512, -4.074512, -4.074512, -4.074512, -4.074512, -4.074512)
actual_medians_treated <- head(sort(hpgl_basic[["medians"]][["treated"]]))
expected_medians_untreated <- c(-4.074491, -4.074491, -4.074491, -4.074491, -4.074491, -4.074491)
actual_medians_untreated <- head(sort(hpgl_basic[["medians"]][["untreated"]]))
expected_logfc <- c(-3.835, -3.645, -3.471, -3.328, -3.072, -3.007)
actual_logfc <- head(sort(hpgl_basic[["all_tables"]][["untreated_vs_treated"]][["logFC"]]))

test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_treated, actual_medians_treated, tolerance=0.001)
    expect_equal(expected_medians_untreated, actual_medians_untreated, tolerance=0.001)
    expect_equal(expected_logfc, actual_logfc, tolerance=0.001)
})

save(list=ls(), file="de_basic.rda")
