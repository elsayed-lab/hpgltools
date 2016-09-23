library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Do surrogate estimators provide expected outputs?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

pasilla_svasup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_supervised"))
expected <- c(0.1184778, 0.5678904, 0.4166829, 0.2938793, 0.5250849, 0.2844593, 0.2164316)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_unsupervised"))
expected <- c(0.1185214, 0.5617694, 0.4119131, 0.2947645, 0.5325897, 0.2850143, 0.2211886)
actual <- abs(as.numeric(pasilla_svaunsup[["model_adjust"]]))
test_that("Have the sva unsupervised model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_pca <- sm(get_model_adjust(pasilla_expt, estimate_type="pca"))
actual <- as.numeric(pasilla_pca[["model_adjust"]])
expected <- c(-0.1089195, -0.5708557, 0.4132932, 0.2933180, -0.5271456, 0.2846236, 0.2156861)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvsup <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_supervised"))
actual <- as.numeric(pasilla_ruvsup[["model_adjust"]])
expected <- c(-0.1085012, -0.5712475, 0.4130910, 0.2928087, -0.5269890, 0.2849590, 0.2158789)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvresid <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_residuals"))
actual <- as.numeric(pasilla_ruvresid[["model_adjust"]])
expected <- c(-0.2029818, -0.3981149, 0.2824784, 0.2863400, -0.6407893, 0.3101006, 0.3629670)
test_that("Have the ruv resid model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvemp <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_empirical"))
actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
expected <- c(-0.1409532, -0.5704923, 0.3987880, 0.2900914, -0.5148203, 0.3007386, 0.2366479)
test_that("Have the ruv resid empirical adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_surrogates <- sm(compare_surrogate_estimates(pasilla_expt))
actual <- as.numeric(as.character(pasilla_surrogates[["pca_adjust"]][["model_adjust"]]))
expected <- as.numeric(c("-0.108919503928491", "-0.570855719628606",
                         "0.413293224481248", "0.293317999343559",
                         "-0.527145619820086", "0.284623554215443",
                         "0.215686065336933"))
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(expected, actual, tolerance=0.0001)
})
