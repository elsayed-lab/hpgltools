start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("41ann_organdb.R: Is it possible to create OrganismDbi instances using TriTrypDB data?\n")

eupath_config <- get_eupath_config()

expected <- 27
actual <- eupath_config["lmajor_friedlin", ][["db_version"]]  ## Pull the db version for lmajor.
test_that("Do we maintain the lmajor data for creating organismdbi objects?", {
    expect_equal(expected, actual)
})

message("Testing organismdbi stuff takes _foever_!  skipping it.")
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ##difficult <- sm(make_organismdbi(id="tcruzi_nonesmer"))
}

message("\nFinished 41ann_organdb.R")
end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 41ann_organdb.R in ", elapsed,  " seconds."))
