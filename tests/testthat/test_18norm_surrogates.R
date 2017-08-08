start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("18norm_surrogates.R: Do surrogate estimators provide expected outputs?\n")

## Some changes were made to the surrogate detectors in 2017-01/2016-12.  In at least one instance
## an error crept in.  In addition, in this time frame I added a hook to the batch_counts() function
## which allows it to call on get_model_adjust() when a batch adjustment method is actually in it.
## The result is a more flexible batch method, but sadly one which has/had at least one error.
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

pasilla_svasup <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="sva_supervised", surrogates="leek"))
expected <- c(0.32445, 0.35891, 0.30689, 0.31796, 0.59513,
              0.28930, 0.36436, 0.77704, 0.42547, 0.37468,
              0.07181, 0.22753, 0.12013, 0.05870)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance=0.0001)
})

pasilla_svasup <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="sva_supervised",
                                      surrogates="be"))
expected <- c(0.32445150, 0.35891445, -0.30688842, -0.31795736,
              0.59513446, -0.28929712, -0.36435751, 0.77703632,
              -0.42547368, -0.37467660, 0.07180690, -0.22752995,
              0.12013255, 0.05870446)
actual <- as.numeric(pasilla_svasup[["model_adjust"]])
test_that("Have the sva supervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt,
                                        estimate_type="sva_unsupervised",
                                        surrogates="leek"))
expected <- c(0.32956, 0.34928, -0.30688, -0.31297, 0.59854,
              -0.29098, -0.36655, 0.76588, -0.43026, -0.38806,
              0.07680, -0.22376, 0.13028, 0.06912)
actual <- as.numeric(pasilla_svaunsup[["model_adjust"]])
test_that("Have the sva unsupervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt,
                                        estimate_type="sva_unsupervised",
                                        surrogates="be"))
expected <- c(0.32955716, 0.34927818, -0.30687904, -0.31297048, 0.59853774, -0.29097615,
              -0.36654742, 0.76587766, -0.43025965, -0.38806172, 0.07679965, -0.22375754,
              0.13027875, 0.06912285)
actual <- as.numeric(pasilla_svaunsup[["model_adjust"]])
test_that("Have the sva unsupervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_pca <- sm(get_model_adjust(pasilla_expt,
                                   estimate_type="pca",
                                   surrogates="leek"))
actual <- as.numeric(pasilla_pca[["model_adjust"]])
expected <- c(-0.3319, -0.3221, -0.3462, -0.3086, 0.4400,
              0.4356, 0.4331, 0.3195, 0.3393, -0.3383,
              -0.3148, 0.6132, -0.2815, -0.3375)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvsup <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="ruv_supervised",
                                      surrogates="leek"))
actual <- as.numeric(pasilla_ruvsup[["model_adjust"]])
expected <- c(-0.32613293, -0.33890442, 0.33582425, 0.30807688, -0.60864913, 0.28081003,
              0.34897533, 0.73915563, -0.46097713, -0.40320483, 0.04569659, -0.19624851,
              0.13867460, 0.13690364)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvresid <- sm(get_model_adjust(pasilla_expt,
                                        estimate_type="ruv_residuals",
                                        surrogates="leek"))
actual <- as.numeric(pasilla_ruvresid[["model_adjust"]])
expected <- c(-0.18528, -0.40475, 0.28166, 0.29323, -0.64607,
              0.31415, 0.34704, 0.78648, -0.44130, -0.41741,
              0.08139, -0.06528, 0.02565, 0.03047)
test_that("Have the ruv resid model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvemp <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="ruv_empirical",
                                      surrogates="leek"))
actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
expected <- c(-0.29016, -0.33367, 0.30928, 0.31323, -0.63656,
              0.28689, 0.35098, 0.76982, -0.38843, -0.42542,
              0.04191, -0.22051, 0.09823, 0.12440)
test_that("Have the ruv resid empirical adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 18norm_surrogates.R in ", elapsed,  " seconds."))
tt <- clear_session()
