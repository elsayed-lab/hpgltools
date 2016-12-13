start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("18norm_surrogates.R: Do surrogate estimators provide expected outputs?\n")

##if (!identical(Sys.getenv("TRAVIS"), "true")) {
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

pasilla_svasup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_supervised",  surrogates="leek"))
expected <- c(0.32640488, 0.35743450, 0.29314922, 0.32047629, 0.59410327,
              0.29067607, 0.37364107, 0.78321891, 0.40191262, 0.37760230,
              0.06839857, 0.24567374, 0.12108913, 0.05248206, 0.14380313,
              0.60234550, 0.64364278, 0.24317391, 0.31576815, 0.05750947, 0.20018519)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svasup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_supervised",  surrogates="be"))
expected <- c(0.32445150, 0.35891445, 0.30688842, 0.31795736, 0.59513446, 0.28929712,
              0.36435751, 0.77703632, 0.42547368, 0.37467660, 0.07180690, 0.22752995,
              0.12013255, 0.05870446)
actual <- abs(as.numeric(pasilla_svasup[["model_adjust"]]))
test_that("Have the sva supervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance=0.000001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_unsupervised", surrogates="leek"))
expected <- c(0.33434443, 0.35128648, 0.29953952, 0.31788830, 0.59331046,
              0.29064867, 0.37086488, 0.77202802, 0.39881500, 0.39351192,
              0.07152666, 0.25150094, 0.12880931, 0.07146387, 0.15922272,
              0.61066712, 0.63487349, 0.24013176, 0.31141362, 0.05123417, 0.20347678)
actual <- abs(as.numeric(pasilla_svaunsup[["model_adjust"]]))
test_that("Have the sva unsupervised model adjustments stayed the same? (leek estimation)", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_svaunsup <- sm(get_model_adjust(pasilla_expt, estimate_type="sva_unsupervised", surrogates="be"))
expected <- c(0.32955716, 0.34927818, 0.30687904, 0.31297048, 0.59853774, 0.29097615,
              0.36654742, 0.76587766, 0.43025965, 0.38806172, 0.07679965, 0.22375754,
              0.13027875, 0.06912285)
actual <- abs(as.numeric(pasilla_svaunsup[["model_adjust"]]))
test_that("Have the sva unsupervised model adjustments stayed the same? (be estimation)", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_pca <- sm(get_model_adjust(pasilla_expt, estimate_type="pca", surrogates="leek"))
actual <- as.numeric(pasilla_pca[["model_adjust"]])
expected <- c(-0.33188016, -0.32210388, -0.34618744, -0.30862665, 0.44002111,
              0.43564626, 0.43313076, 0.31953018, 0.33931882, -0.33831657,
              -0.31484103, 0.61324502, -0.28145361, -0.33748281, -0.75746814,
              0.44315486, 0.39164587, -0.07983911, 0.21689621, -0.11309356, -0.10129612)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvsup <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_supervised", surrogates="leek"))
actual <- as.numeric(pasilla_ruvsup[["model_adjust"]])
expected <- c(-0.31886894, -0.34235209, 0.33123193, 0.31618389, -0.61161597, 0.28287038,
              0.34255080, 0.74525582, -0.45037457, -0.40959131, 0.07230244, -0.19867584,
              0.13129157, 0.10979189, -0.36632996, -0.32087498, -0.34475205, -0.27288995,
              0.47024987, 0.40996248, 0.42463458)
test_that("Have the pca model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvresid <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_residuals", surrogates="leek"))
actual <- as.numeric(pasilla_ruvresid[["model_adjust"]])
expected <- c(-0.18528065, -0.40474586, 0.28166288, 0.29323140, -0.64606575, 0.31415365,
              0.34704432, 0.78648114, -0.44129992, -0.41740825, 0.08138963, -0.06528100,
              0.02564705, 0.03047136, -0.23243097, 0.32359335, -0.64436300, 0.57123122,
              -0.10347752, -0.16391176, 0.24935868)
test_that("Have the ruv resid model adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_ruvemp <- sm(get_model_adjust(pasilla_expt, estimate_type="ruv_empirical", surrogates="leek"))
actual <- as.numeric(pasilla_ruvemp[["model_adjust"]])
expected <- c(-0.290157538, -0.333667257, 0.309283947, 0.313229622, -0.636561735,
              0.286892694, 0.350980267, 0.769824150, -0.388433458, -0.425423308,
              0.041913942, -0.220509554, 0.098231109, 0.124397118, -0.271071141,
              0.273712196, -0.661336852, 0.630444290, -0.005197405, -0.073245559,
              0.106694471)

test_that("Have the ruv resid empirical adjustments stayed the same?", {
    expect_equal(expected, actual, tolerance=0.001)
})

pasilla_surrogates <- sm(compare_surrogate_estimates(pasilla_expt, surrogates="leek"))
actual <- as.numeric(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])
expected <- c(-0.33188016, -0.32210388, -0.34618744, -0.30862665, 0.44002111, 0.43564626,
              0.43313076, 0.31953018, 0.33931882, -0.33831657, -0.31484103, 0.61324502,
              -0.28145361, -0.33748281, -0.75746814, 0.44315486, 0.39164587, -0.07983911,
              0.21689621, -0.11309356, -0.10129612)
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 18norm_surrogates.R in ", elapsed,  " seconds."))
##}
