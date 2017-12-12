start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("10ann_eupathdb.R: Is it possible to extract EuPathDB data?\n")
## 2017-12, exported functions in annotation_eupathdb:
##  make_eupath_bsgenome(), make_eupath_organismdbi() download_eupath_metadata(),
##  make_eupath_orgdb(), make_eupath_txdb(), get_eupath_gff_table(),
##  get_eupath_gene_types(), get_eupath_go_term_table(), get_eupath_pathway_table(),
##  get_eupath_interpro_table(), get_eupath_ortholog_table() get_eupath_text()
##  Most of these are implicitly tested via make_eupath_organismdbi().

eupath_metadata <- sm(download_eupath_metadata())
expected <- c(285, 19)
actual <- dim(eupath_metadata)
test_that("Is the eupathdb metadata the expected size?", {
  expect_equal(expected, actual)
})

tritryp_metadata <- sm(download_eupath_metadata(webservice="tritrypdb"))
expected <- c(41, 19)
actual <- dim(tritryp_metadata)
test_that("Is the eupathdb metadata the expected size?", {
  expect_equal(expected, actual)
})

## You know what, making all of these will take days, just pull a random 2
randomOrg <- random::randomNumbers(n=3, min=1, max=nrow(eupath_metadata), col=1)
chosen <- randomOrg[, 1]
eupath_names <- eupath_metadata[["Species"]]
for (i in 1:length(chosen)) {
  index <- chosen[i]
  species <- eupath_names[index]
  message(paste0("Going to attempt making packages for: ", species))
  eupath_test <- try(sm(make_eupath_organismdbi(species=species, metadata=eupath_metadata)))
  if (class(eupath_test) != "try-error") {
    test_that("Did the organismdbi get installed?", {
      expect_true(eupath_test[["organdb_name"]] %in% installed.packages())
    })
  } else {
    message(paste0("Creation of orgdb/etc failed for ", species))
  }
  bsgenome_test <- make_eupath_bsgenome(species, reinstall=TRUE)
  test_that("Did the bsgenome get installed?", {
    expect_true(bsgenome_test[["bsgenome_name"]] %in% installed.packages())
  })
  text_test <- sm(get_eupath_text(species, metadata=eupath_metadata))
  test_that("Do we get interesting text data?", {
    expect_gt(nrow(text_test), 1000)
    expect_gt(ncol(text_test), 10)
  })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 10annotation_eupathdb.R in ", elapsed,  " seconds."))
