start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("10ann_eupathdb.R: Is it possible to extract EuPathDB data?\n")

eupath_metadata <- download_eupathdb_metadata()
expected <- c(281, 19)
actual <- dim(eupath_metadata)
test_that("Is the eupathdb metadata the expected size?", {
  expect_equal(expected, actual)
})

tritryp_metadata <- download_eupathdb_metadata(webservice="tritrypdb")
expected <- c(41, 19)
actual <- dim(tritryp_metadata)
test_that("Is the eupathdb metadata the expected size?", {
  expect_equal(expected, actual)
})

## You know what, making all of these will take days, just pull a random 5
chosen <- sample(1:nrow(eupath_metadata), 5)
eupath_names <- eupath_metadata[["Species"]]
for (index in chosen) {
  species <- eupath_names[index]
  eupath_test <- make_eupath_organismdbi(species)
  test_that("Did the orgdb get installed?", {
    expect_true(eupath_test[["orgdb_name"]] %in% installed.packages())
  })
  test_that("Did the txdb get installed?", {
    expect_true(eupath_test[["txdb_name"]] %in% installed.packages())
  })
  test_that("Did the organismdbi get installed?", {
    expect_true(eupath_test[["organdb_name"]] %in% installed.packages())
  })
  bsgenome_test <- make_eupath_bsgenome(species)
  expect_true(bsgenome_test)

  dm28c_test <- make_eupath_bsgenome(
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 46ann_tritrypdb.R in ", elapsed,  " seconds."))
