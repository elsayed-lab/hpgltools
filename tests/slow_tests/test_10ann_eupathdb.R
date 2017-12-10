start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("10ann_eupathdb.R: Is it possible to extract EuPathDB data?\n")

eupath_metadata <- download_eupathdb_metadata()
expected <- c(285, 19)
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

## You know what, making all of these will take days, just pull a random 2
ok <- RNGkind()
RNGkind("Super")#matches  "Super-Duper"
set.seed(.Random.seed)
## Apparently this is insufficient, I still get the same 2...
chosen <- sample(1:nrow(eupath_metadata), 2)
eupath_names <- eupath_metadata[["Species"]]
for (index in chosen) {
  species <- eupath_names[index]
  eupath_test <- try(make_eupath_organismdbi(species=species, metadata=eupath_metadata))
  if (class(eupath_test) != "try-error") {
    test_that("Did the organismdbi get installed?", {
      expect_true(eupath_test[["organdb_name"]] %in% installed.packages())
    })
  }
  bsgenome_test <- make_eupath_bsgenome(species, reinstall=TRUE)
  test_that("Did the bsgenome get installed?", {
    expect_true(bsgenome_test[["bsgenome_name"]] %in% installed.packages())
  })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 46ann_tritrypdb.R in ", elapsed,  " seconds."))
