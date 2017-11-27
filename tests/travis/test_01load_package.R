start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("01load_package.R: Does hpgltools load properly?\n")

installedp <- utils::vignette(package="hpgltools")
available_vignettes <- as.character(as.data.frame(installedp[["results"]])[["Item"]])

## In theory, there should currently be 5 vignettes installed, a-d and a reference manual.
expected <- 5
actual <- length(available_vignettes)
test_that("Did the vignettes install?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 01load_package.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
