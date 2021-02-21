start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
context("18norm_surrogates.R: Do surrogate estimators provide expected outputs?\n")

## Some changes were made to the surrogate detectors in 2017-01/2016-12.  In at least one instance
## an error crept in.  In addition, in this time frame I added a hook to the batch_counts() function
## which allows it to call on get_model_adjust() when a batch adjustment method is actually in it.
## The result is a more flexible batch method, but sadly one which has/had at least one error.

load("pasilla_df.rda")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]

pasilla_svasup <- all_adjusters(pasilla_expt,
                                estimate_type = "sva_supervised",
                                surrogates = 2)
expected <- c(0.11443442, 0.56917862, 0.41547755, 0.29366937, 0.52595030,
              0.28371316, 0.21670326, 0.84677563, 0.07331834, 0.05389192,
              0.14758535, 0.32210777, 0.25490952, 0.29013342)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

pasilla_svasup <- all_adjusters(pasilla_expt,
                                estimate_type = "sva_supervised",
                                surrogates = "be")
expected <- c(0.1197830, 0.5674614, -0.4171903, -0.2939268,
              0.5248040, -0.2844763, -0.2164550)
actual <- as.numeric(pasilla_svasup[["model_adjust"]])
test_that("Have the sva supervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance = 0.000001)
})

pasilla_svaunsup <- all_adjusters(pasilla_expt,
                                  estimate_type = "sva_unsupervised",
                                  surrogates = 2)
expected <- c(0.11724026, 0.56231486, -0.41018053, -0.29286547,
              0.53309354, -0.28644237, -0.22316028, 0.85562428,
              -0.32330553, -0.37725111, -0.01347759, -0.14237170,
              0.01786325, -0.01708160)
actual <- as.numeric(pasilla_svaunsup[["model_adjust"]])
test_that("Have the sva unsupervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance = 0.001)
})

pasilla_svaunsup <- all_adjusters(pasilla_expt,
                                  estimate_type = "sva_unsupervised",
                                  surrogates = "be")
expected <- c(0.1185214, 0.5617694, -0.4119131, -0.2947645, 0.5325897, -0.2850143, -0.2211886)
actual <- as.numeric(pasilla_svaunsup[["model_adjust"]])
test_that("Have the sva unsupervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance = 0.001)
})

pasilla_pca <- all_adjusters(pasilla_expt,
                             estimate_type = "pca",
                             surrogates = 2)
actual <- as.numeric(pasilla_pca[["model_adjust"]])
expected <- c(-0.1089195, -0.5708557, 0.4132932, 0.2933180, -0.5271456,
              0.2846236, 0.2156861, 0.4473820, 0.2704362, 0.2892021,
              0.2884778, -0.4160025, -0.4443341, -0.4351615)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

pasilla_ruvsup <- all_adjusters(pasilla_expt,
                                estimate_type = "ruv_supervised",
                                surrogates = 2)
actual <- as.numeric(pasilla_ruvsup[["model_adjust"]])
expected <- c(-0.1089663, -0.5707737, 0.4133079, 0.2933209, -0.5272058, 0.2846454, 0.2156716)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

pasilla_ruvresid <- all_adjusters(pasilla_expt,
                                  estimate_type = "ruv_residuals",
                                  surrogates = 2)
actual <- as.numeric(pasilla_ruvresid[["model_adjust"]])
expected <- c(-0.19913107, -0.40282892, 0.28770617, 0.29350932,
              -0.64006827, 0.30664759, 0.35416518, 0.79201575,
              -0.44226077, -0.40448415, 0.06206771, -0.08226267,
              0.03103248, 0.04389165)
test_that("Have the ruv resid model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

pasilla_ruvemp <- all_adjusters(pasilla_expt,
                                estimate_type = "ruv_empirical",
                                surrogates = 2)
actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
expected <- c(-0.14104080, -0.57046187, 0.39876727, 0.29005796, -0.51481268,
              0.30077535, 0.23671477, 0.83924154, -0.36282048, -0.38458306,
              -0.03511887, -0.11437361, 0.03760568, 0.02004879)
test_that("Have the ruv resid empirical adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance = 0.001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 18norm_surrogates.R in ", elapsed,  " seconds.")
