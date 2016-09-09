library(testthat)
library(hpgltools)

context("Is it possible to graph the various metrics with hpgltools?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
## Uses these genes for quick tests
test_genes <- c("FBgn0000014","FBgn0000008","FBgn0000017","FBgn0000018", "FBgn0000024")

## What graphs can we make!?

libsize_plot <- plot_libsize(pasilla_expt)
actual <- libsize_plot[["data"]][["sum"]]
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
test_that("The libsize plot is as expected?", {
    expect_equal(expected, actual)
})

nonzero_plot <- plot_nonzero(pasilla_expt)
actual <- nonzero_plot[["data"]][["nonzero_genes"]]
expected <- c(9863, 10074, 9730, 9786, 10087, 9798, 9797)
test_that("The non-zero genes is as expected?", {
    expect_equal(expected, actual)
})

## These tests have also been affected by the changed order of expressionsets.
density_plot <- sm(plot_density(pasilla_expt))
expected <- c(92, 5, 4664, 583, 10, 1446)
actual <- head(density_plot[["data"]][["counts"]])
test_that("Density plot data is as expected?", {
    expect_equal(expected, actual)
})

hist_plot <- sm(plot_histogram(data.frame(Biobase::exprs(pasilla_expt[["expressionset"]]))))
actual <- head(hist_plot[["data"]][["values"]])
## The values of expected have not changed
test_that("Histogram data is as expected?", {
    expect_equal(expected, actual)
})

box_plot <- sm(plot_boxplot(pasilla_expt))
actual <- head(box_plot[["data"]][["value"]])
## The values of expected have not changed
test_that("Box plot data is as expected?", {
    expect_equal(expected, actual)
})

## Ahh yes I changed the cbcb_filter options to match those from the cbcbSEQ vignette.
norm <- sm(normalize_expt(pasilla_expt, transform="log2", convert="cbcbcpm", norm="quant", filter=TRUE))
expected <- "recordedplot"  ## for all the heatmaps

corheat_plot <- plot_corheat(norm)
actual <- class(corheat_plot)
test_that("corheat is a recorded plot?", {
    expect_equal(expected, actual)
})

disheat_plot <- plot_disheat(norm)
actual <- class(disheat_plot)
test_that("disheat is a recorded plot?", {
    expect_equal(expected, actual)
})

sampleheat_plot <- plot_sample_heatmap(norm)
actual <- class(sampleheat_plot)
test_that("sampleheat is a recorded plot?", {
    expect_equal(expected, actual)
})

smc_plot <- sm(plot_sm(norm, method="pearson"))
actual <- head(smc_plot[["data"]][["sm"]])
expected <- c(0.9779550, 0.9822023, 0.9798571, 0.9816597, 0.9783356, 0.9805947)
test_that("Is the normalized smc data expected?", {
    expect_equal(expected, actual, tolerance=0.004)
})

smd_plot <- sm(plot_sm(norm, method="euclidean"))
actual <- head(smd_plot[["data"]][["sm"]])
expected <- c(72.26, 64.90, 68.90, 65.68, 71.64, 67.65)
test_that("Is the normalized smd data expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pca_stuff <- plot_pca(norm)
pca_plot <- pca_stuff[["plot"]]
pca_pca <- head(pca_stuff[["pca"]])

actual <- pca_plot[["data"]][["PC1"]]
expected <- c(-0.484593397, -0.464554956, 0.008466696, -0.030043569, -0.070616473, 0.496169040, 0.545172660)
test_that("Is the pca data as expected for PC1?", {
    expect_equal(expected, actual, tolerance=0.001)
})

actual <- as.numeric(head(pca_stuff[["pca"]][["v"]][, 1]))
expected <- c(-0.484593397, -0.464554956, 0.008466696, -0.030043569, -0.070616473, 0.496169040)
test_that("Is the SVD 'v' element expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

actual <- pca_stuff[["res"]][[1]]
expected <- c(29.82, 25.10, 15.84, 10.73, 9.71, 8.80)
test_that("Is the pca residual table as expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

actual_pca_var <- pca_stuff[["variance"]]
expected_pca_var <- c(29.82, 25.10, 15.84, 10.73, 9.71, 8.80)
test_that("Is the variance list as expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

actual <- pca_stuff[["table"]][["PC2"]]
expected <- c(0.06119973, 0.04514489, -0.53934612, -0.44395229, 0.69872704, 0.12172316, 0.05650359)
test_that("Is the PCA PC2 as expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})
