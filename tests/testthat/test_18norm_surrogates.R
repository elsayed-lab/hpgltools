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
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    pasilla <- new.env()
    load("pasilla.Rdata", envir=pasilla)
    pasilla_expt <- pasilla[["expt"]]

    pasilla_svasup <- sm(get_model_adjust(pasilla_expt,
                                          estimate_type="sva_supervised",
                                          surrogates="leek"))
    expected <- c(0.32445150, 0.35891445, 0.30688842, 0.31795736, 0.59513446, 0.28929712,
                  0.36435751, 0.77703632, 0.42547368, 0.37467660, 0.07180690, 0.22752995,
                  0.12013255, 0.05870446)
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
    expected <- c(0.32955716, 0.34927818, -0.30687904, -0.31297048, 0.59853774, -0.29097615,
                  -0.36654742, 0.76587766, -0.43025965, -0.38806172, 0.07679965, -0.22375754,
                  0.13027875, 0.06912285)
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
    expected <- c(-0.3318802, -0.3221039, -0.3461874, -0.3086267, 0.4400211, 0.4356463,
                  0.4331308, 0.3195302, 0.3393188, -0.3383166, -0.3148410, 0.6132450,
                  -0.2814536, -0.3374828)
    test_that("Have the pca model adjustments stayed the same?", {
        expect_equal(expected, actual, tolerance=0.001)
    })

    pasilla_ruvsup <- sm(get_model_adjust(pasilla_expt,
                                          estimate_type="ruv_supervised",
                                          surrogates="leek"))
    actual <- as.numeric(pasilla_ruvsup[["model_adjust"]])
    ## expected <- c(-0.31886894, -0.34235209, 0.33123193, 0.31618389, -0.61161597, 0.28287038,
    ##               0.34255080, 0.74525582, -0.45037457, -0.40959131, 0.07230244, -0.19867584,
    ##               0.13129157, 0.10979189)
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
    expected <- c(-0.18528065, -0.40474586, 0.28166288, 0.29323140, -0.64606575, 0.31415365,
                  0.34704432, 0.78648114, -0.44129992, -0.41740825, 0.08138963, -0.06528100,
                  0.02564705, 0.03047136)
    test_that("Have the ruv resid model adjustments stayed the same?", {
        expect_equal(expected, actual, tolerance=0.001)
    })

    pasilla_ruvemp <- sm(get_model_adjust(pasilla_expt,
                                          estimate_type="ruv_empirical",
                                          surrogates="leek"))
    actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
    expected <- c(-0.29015754, -0.33366726, 0.30928395, 0.31322962, -0.63656174, 0.28689269,
                  0.35098027, 0.76982415, -0.38843346, -0.42542331, 0.04191394, -0.22050955,
                  0.09823111, 0.12439712)
    test_that("Have the ruv resid empirical adjustments stayed the same?", {
        expect_equal(expected, actual, tolerance=0.001)
    })

    pasilla_surrogates <- sm(compare_surrogate_estimates(pasilla_expt, surrogates="leek"))
    actual <- as.numeric(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])
    expected <- c(-0.3318802, -0.3221039, -0.3461874, -0.3086267, 0.4400211, 0.4356463,
                  0.4331308, 0.3195302, 0.3393188, -0.3383166, -0.3148410, 0.6132450,
                  -0.2814536, -0.3374828)
    test_that("Does the compare_surrogate stuff work?", {
        expect_equal(expected, actual, tolerance=0.0001)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 18norm_surrogates.R in ", elapsed,  " seconds."))
