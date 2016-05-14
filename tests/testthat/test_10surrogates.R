library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Do surrogate estimators provide expected outputs?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

pasilla_surrogates <- suppressMessages(compare_surrogate_estimates(pasilla_expt))
pca_char_adjust <- as.numeric(c("-0.108919503928491", "-0.570855719628606",
                                "0.413293224481248", "0.293317999343559",
                                "-0.527145619820086", "0.284623554215443",
                                "0.215686065336933"))
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(pca_char_adjust,
                 as.numeric(as.character(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])),
                 tolerance <- 0.0001)
})
