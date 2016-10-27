library(testthat)
library(hpgltools)
context("41ann_organdb.R: Is it possible to create OrganismDbi instances using TriTrypDB data?\n")

eupath_config <- get_eupath_config()

expected <- 27
actual <- eupath_config["lmajor_friedlin", ][["db_version"]]  ## Pull the db version for lmajor.
expect_equal(expected, actual)

if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ##difficult <- sm(make_organismdbi(id="tcruzi_nonesmer"))
}

message("\nFinished 41ann_organdb.R")
