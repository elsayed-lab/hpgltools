library(testthat)
library(hpgltools)
context("27de_basic: Does the basic differential expression analysis work?\n")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

norm_expt <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm", filter=TRUE))

hpgl_pas_basic <- basic_pairwise(pasilla_expt)
hpgl_basic <- basic_pairwise(norm_expt)
head(hpgl_pas_basic$input_data)
head(hpgl_basic$input_data)

test_that("Does a non-normalized basic run equal a normalized basic run?", {
    expect_equal(hpgl_pas_basic[["all_tables"]][["untreated_vs_treated"]],
                 hpgl_basic[["all_tables"]][["untreated_vs_treated"]])
})

expected_medians_treated <- c(-4.074444, -4.074444, -4.074444, -4.074444, -4.074444, -4.074444)
actual_medians_treated <- head(sort(hpgl_basic[["medians"]][["treated"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_treated, actual_medians_treated, tolerance=0.001)
}

expected_medians_untreated <- c(-3.850703, -3.850703, -3.850703, -3.850703, -3.850703, -3.850703)
actual_medians_untreated <- head(sort(hpgl_basic[["medians"]][["untreated"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_untreated, actual_medians_untreated, tolerance=0.001)
})

expected_logfc <- c(-3.645, -3.611, -3.471, -3.328, -3.079, -3.007)
actual_logfc <- head(sort(hpgl_basic[["all_tables"]][["untreated_vs_treated"]][["logFC"]]))
test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_logfc, actual_logfc, tolerance=0.001)
})

save(list=ls(), file="de_basic.rda")

message("\nFinished 27de_basic.R")
