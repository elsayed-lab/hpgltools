start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("01load_package.R: Does hpgltools load properly?")

installedp <- utils::vignette(package="hpgltools")
available_vignettes <- as.character(as.data.frame(installedp[["results"]])[["Item"]])
expected <- 5
actual <- length(available_vignettes)
test_that("Did the vignettes install?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 01load_data.R in ", elapsed,  " seconds."))
