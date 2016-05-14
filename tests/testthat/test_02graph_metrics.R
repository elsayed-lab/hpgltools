library(testthat)
library(hpgltools)

context("Is it possible to graph the various metrics with hpgltools?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## What graphs can we make!?
libsize_plot <- plot_libsize(pasilla_expt)
nonzero_plot <- plot_nonzero(pasilla_expt)
density_plot <- suppressMessages(plot_density(pasilla_expt))
hist_plot <- suppressMessages(plot_histogram(data.frame(Biobase::exprs(pasilla_expt$expressionset))))
box_plot <- suppressMessages(plot_boxplot(pasilla_expt))

expected_libsize_data <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
actual_libsize_data <- libsize_plot$data$sum
expected_nonzero_data <- c(9863, 10074, 9730, 9786, 10087, 9798, 9797)
actual_nonzero_data <- nonzero_plot$data$nonzero_genes
expected_density_data <- c(92, 5, 4664, 583, 10, 1446)
actual_density_data <- head(density_plot$data$counts)
expected_hist_data <- c(92, 5, 4664, 583, 10, 1446)
actual_hist_data <- head(hist_plot$data$values)
expected_box_data <- c(92, 5, 4664, 583, 10, 1446)
actual_box_data <- head(box_plot$data$value)
test_that("Do the pre-normalized plots work", {
    expect_equal(expected_libsize_data, actual_libsize_data)
    expect_equal(expected_nonzero_data, actual_nonzero_data)
    expect_equal(expected_density_data, actual_density_data)
    expect_equal(expected_hist_data, actual_hist_data)
    expect_equal(expected_box_data, actual_box_data)
})

norm <- suppressMessages(normalize_expt(pasilla_expt, transform="log2", convert="cpm", norm="quant", filter_low=TRUE))

corheat_plot <- plot_corheat(norm)
disheat_plot <- plot_disheat(norm)
sampleheat_plot <- plot_sample_heatmap(norm)
smc_plot <- suppressMessages(plot_sm(norm, method="pearson"))
smd_plot <- suppressMessages(plot_sm(norm, method="euclidean"))
pca_stuff <- plot_pca(norm)
pca_plot <- pca_stuff$plot
pca_pca <- head(pca_stuff$pca)

expected_corheat_class <- expected_disheat_class <- expected_sampleheat_class <- "recordedplot"
actual_corheat_class <- class(corheat_plot)
actual_disheat_class <- class(disheat_plot)
actual_sampleheat_class <- class(sampleheat_plot)
expected_smc_data <- c(0.9779550, 0.9822023, 0.9798571, 0.9816597, 0.9783356, 0.9805947)
actual_smc_data <- head(smc_plot$data$sm)
expected_smd_data <- c(72.25638, 64.89870, 68.90467, 65.67845, 71.63792, 67.65131)
actual_smd_data <- head(smd_plot$data$sm)
expected_pca_data <- c(-0.484593397, -0.464554956, 0.008466696, -0.030043569, -0.070616473, 0.496169040, 0.545172660)
actual_pca_data <- pca_plot$data$PC1
expected_pca_v <- c(-0.484593397, -0.464554956, 0.008466696, -0.030043569, -0.070616473, 0.496169040)
actual_pca_v <- as.numeric(head(pca_stuff$pca$v[, 1]))
expected_pca_res <- c(29.82, 25.10, 15.84, 10.73, 9.71, 8.80)
actual_pca_res <- pca_stuff$res[[1]]
expected_pca_var <- c(29.82, 25.10, 15.84, 10.73, 9.71, 8.80)
actual_pca_var <- pca_stuff$variance
expected_pca_tab <- c(0.06119973, 0.04514489, -0.53934612, -0.44395229, 0.69872704, 0.12172316, 0.05650359)
actual_pca_tab <- pca_stuff$table$PC2
test_that("Do the normalized plots work", {
    expect_equal(expected_corheat_class, actual_corheat_class)
    expect_equal(expected_disheat_class, actual_disheat_class)
    expect_equal(expected_sampleheat_class, actual_sampleheat_class)
    expect_equal(expected_smc_data, actual_smc_data, tolerance=0.001)
    expect_equal(expected_smd_data, actual_smd_data, tolerance=0.001)
    expect_equal(expected_pca_data, actual_pca_data, tolerance=0.001)
    expect_equal(expected_pca_v, actual_pca_v, tolerance=0.001)
    expect_equal(expected_pca_res, actual_pca_res, tolerance=0.001)
    expect_equal(expected_pca_var, actual_pca_var, tolerance=0.001)
    expect_equal(expected_pca_tab, actual_pca_tab, tolerance=0.001)
})
