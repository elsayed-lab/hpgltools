start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("19varpart.R: Does variancePartition work as expected?\n")

## Some changes were made to the surrogate detectors in 2017-01/2016-12.  In at least one instance
## an error crept in.  In addition, in this time frame I added a hook to the batch_counts() function
## which allows it to call on get_model_adjust() when a batch adjustment method is actually in it.
## The result is a more flexible batch method, but sadly one which has/had at least one error.
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    pasilla <- new.env()
    load("pasilla.Rdata", envir=pasilla)
    pasilla_expt <- pasilla[["expt"]]

    message("Friendly warning, this is slow.")
    pasilla_varpart <- sm(varpart(pasilla_expt, predictor=NULL, factors=c("condition", "batch")))

    ## Grab the model and see if it survived.
    expected <- "(1 | condition) + (1 | batch)"
    actual <- as.character(pasilla_varpart[["model_used"]])[[2]]
    test_that("Does my varpart function return a sane model?", {
        expect_equal(expected, actual)
    })

    ## See if the variance by condition is as expected
    expected <- c(0.8451941, 0.8436158, 0.8364982, 0.8361348, 0.8240863, 0.8174928)
    actual <- head(pasilla_varpart[["sorted_df"]])[["condition"]]
    test_that("Do we get expected values of variance by condition?", {
        expect_equal(expected, actual, tolerance=0.001)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 19varpart.R in ", elapsed,  " seconds."))
tt <- clear_session()
