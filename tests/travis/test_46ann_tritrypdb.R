start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("46ann_tritrypdb.R: Is it possible to extract TriTrypDB data?\n")

eupath_config <- get_eupath_config()
expected <- 28
actual <- eupath_config["lmajor_friedlin", ][["db_version"]]  ## Pull the db version for lmajor.
test_that("Do we maintain the lmajor data for creating organismdbi objects?", {
    expect_equal(expected, actual)
})

if (!identical(Sys.getenv("TRAVIS"), "true")) {
    message("Testing make_organismdbi() takes _foever_!  skipping it.")
    ##difficult <- sm(make_organismdbi(id="tcruzi_nonesmer"))
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 46ann_tritrypdb.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
