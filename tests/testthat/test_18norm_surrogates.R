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
expected <- c(0.3334780, 0.3537427, 0.3167739, 0.3146217, 0.5925343, 0.2851664, 0.3631930)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svasup <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="sva_supervised",
                                      surrogates="be"))
expected <- c(0.32445150, 0.35891445, -0.30688842, -0.31795736, 0.59513446, -0.28929712,
              -0.36435751, 0.77703632, -0.42547368, -0.37467660, 0.07180690, -0.22752995,
              0.12013255, 0.05870446)
actual <- as.numeric(pasilla_svasup[["model_adjust"]])
test_that("Have the sva supervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt,
                                        estimate_type="sva_unsupervised",
                                        surrogates="leek"))
expected <- c(0.3339817, 0.3453492, -0.3132165, -0.3104993, 0.5983120, -0.2889621, -0.3649650)
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
expected <- c(-0.3318802, -0.3221039, -0.3461874, -0.3086267, 0.4400211, 0.4356463, 0.4331308)
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
expected <- c(-0.1852807, -0.4047459, 0.2816629, 0.2932314, -0.6460658, 0.3141537, 0.3470443)
test_that("Have the ruv resid model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvemp <- sm(get_model_adjust(pasilla_expt,
                                      estimate_type="ruv_empirical",
                                      surrogates="leek"))
actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
expected <- c(-0.2901575, -0.3336673, 0.3092839, 0.3132296, -0.6365617, 0.2868927, 0.3509803)
test_that("Have the ruv resid empirical adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_surrogates <- sm(compare_surrogate_estimates(pasilla_expt, surrogates="leek"))
actual <- as.numeric(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])
expected <- c(-0.3318802, -0.3221039, -0.3461874, -0.3086267, 0.4400211, 0.4356463, 0.4331308)
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 18norm_surrogates.R in ", elapsed,  " seconds."))
tt <- clear_session()
