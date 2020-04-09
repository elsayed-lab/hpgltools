start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("03graph_metrics.R: Is it possible to graph the various metrics with hpgltools?
  1234567890123456789\n")

pasilla <- new.env()
load("pasilla.rda", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")

## What graphs can we make!?
libsize_plot <- plot_libsize(pasilla_expt)
actual <- libsize_plot[["table"]][["sum"]]
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
## 01
test_that("The libsize plot is as expected?", {
    expect_equal(expected, actual)
})

nonzero_plot <- plot_nonzero(pasilla_expt)
actual <- nonzero_plot[["table"]][["nonzero_genes"]]
expected <- c(9863, 10074, 9730, 9786, 10087, 9798, 9797)
## 02
test_that("The non-zero genes is as expected?", {
    expect_equal(expected, actual)
})

## These tests have also been affected by the changed order of expressionsets.
density <- sm(plot_density(pasilla_expt))
density_plot <- density[["plot"]]
density_table <- density[["table"]]
expected <- c(92, 5, 4664, 583, 10, 1446)
actual <- head(density_table[["counts"]])
## 03
test_that("Density plot data is as expected?", {
    expect_equal(expected, actual)
})

hist_plot <- sm(plot_histogram(data.frame(exprs(pasilla_expt))))
actual <- head(hist_plot[["data"]][["values"]])
## The values of expected have not changed
## 04
test_that("Histogram data is as expected?", {
    expect_equal(expected, actual)
})

box_plot <- sm(plot_boxplot(pasilla_expt))
actual <- head(box_plot[["data"]][["value"]])
## 05
test_that("Box plot data is as expected?", {
    expect_equal(expected, actual, tolerance=1)
})

## Ahh yes I changed the cbcb_filter options to match those from the cbcbSEQ vignette.
## Note that the filtering has changed slightly, and this affects the results.
norm <- sm(normalize_expt(pasilla_expt, transform="log2", convert="cbcbcpm",
                          norm="quant", filter=TRUE))
expected <- "recordedplot"  ## for all the heatmaps

corheat_plot <- plot_corheat(norm)
actual <- class(corheat_plot[["plot"]])
## 06
test_that("corheat is a recorded plot?", {
    expect_equal(expected, actual)
})

disheat_plot <- plot_disheat(norm)
actual <- class(disheat_plot[["plot"]])
## 07
test_that("disheat is a recorded plot?", {
    expect_equal(expected, actual)
})

sampleheat_plot <- plot_sample_heatmap(norm)
actual <- class(sampleheat_plot)
## 08
test_that("sampleheat is a recorded plot?", {
    expect_equal(expected, actual)
})

smc_plot <- sm(plot_sm(norm, method="pearson"))
actual <- head(smc_plot[["plot"]][["data"]][["sm"]])
expected <- c(0.9759981, 0.9824316, 0.9759981, 0.9821373, 0.9784851, 0.9786376)
## 09
test_that("Is the normalized smc data expected?", {
    expect_equal(expected, actual, tolerance=0.004)
})

smd_plot <- sm(plot_sm(norm, method="euclidean"))
actual <- head(smd_plot[["plot"]][["data"]][["sm"]])
## 201812 Changed due to peculiarities in normalization methods.
## 201907 I changed the normalization back, so the values need to return.
expected <- c(42.43941, 36.43116, 42.43941, 36.60569, 40.01228, 40.04465)
## 10
test_that("Is the normalized smd data expected?", {
    expect_equal(expected, actual, tolerance=0.01)
})

pca_stuff <- plot_pca(norm)
pca_plot <- pca_stuff[["plot"]]
pca_pca <- head(pca_stuff[["pca"]])

actual <- pca_plot[["data"]][["PC1"]]
## Ibid
expected <- c(-0.3595704, -0.4034877, -0.2721046, -0.2432659, 0.2850077, 0.4995150, 0.4939059)
## 11
test_that("Is the pca data as expected for PC1?", {
    expect_equal(expected, actual, tolerance=0.01)
})

actual <- as.numeric(head(pca_stuff[["result"]][["v"]][, 1]))
## Ibid
expected <- c(-0.3595704, -0.4034877, -0.2721046, -0.2432659, 0.2850077, 0.4995150)
## 12
test_that("Is the SVD 'v' element expected?", {
    expect_equal(expected, actual, tolerance=0.01)
})

actual <- pca_stuff[["residual_df"]][[1]]
expected <- c(42.44, 31.26, 13.09, 5.84, 4.13, 3.24)
## 13
test_that("Is the pca residual table as expected?", {
    expect_equal(expected, actual, tolerance=0.01)
})

actual <- pca_stuff[["prop_var"]]
expected <- c(42.44, 31.26, 13.09, 5.84, 4.13, 3.24)
## 14
test_that("Is the variance list as expected?", {
    expect_equal(expected, actual, tolerance=0.01)
})

actual <- pca_stuff[["table"]][["PC2"]]
expected <- c(0.3014798, 0.2752759, -0.4578880, -0.3894024, 0.6349581, -0.1471152, -0.2173081)
## 15
test_that("Is the PCA PC2 as expected?", {
    expect_equal(expected, actual, tolerance=0.01)
})

tsne_stuff <- plot_tsne(norm, seed=1)
tsne_stuff$plot
actual <- tsne_stuff[["table"]][["Factor1"]]
expected <- c(117.7296, 117.7064, -287.8868, -287.8774, 105.9297, 114.4243, 119.9742)
## These values seem to have changed in the new version of Rtsne.
## 16
test_that("Is the tsne data as expected for Comp1?", {
    expect_equal(expected, actual, tolerance=0.1)
})

actual <- as.numeric(head(tsne_stuff[["result"]][["Y"]][, 2]))
expected <- c(-165.0350, -170.4171, -164.4600, -169.8417, 223.7109, 223.1911)
## These also changed.
## 17
test_that("Is the tsne second component data expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

actual <- head(tsne_stuff[["residual_df"]][["condition_rsquared"]])
expected <- c(0.2910028, 0.9998866, 0.2910028, 0.9998866, 0.2910028, 0.9998866)
## 18
test_that("Is the tsne r-squared by condition as expected?", {
    expect_equal(expected, actual, tolerance=0.001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 03graph_metrics.R in ", elapsed, " seconds."))
tt <- try(clear_session())
