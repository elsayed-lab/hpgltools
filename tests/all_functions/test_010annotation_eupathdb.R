start <- as.POSIXlt(Sys.time())
context("010annotation_eupathdb.R\n")
## 2018-02: Hey, there is a new eupathdb release!  Some stuff has changed!
## 2017-12, exported functions in annotation_eupathdb:
##  make_eupath_bsgenome(), make_eupath_organismdbi() download_eupath_metadata(),
##  make_eupath_orgdb(), make_eupath_txdb(), get_eupath_gff_table(),
##  get_eupath_gene_types(), get_eupath_go_term_table(), get_eupath_pathway_table(),
##  get_eupath_interpro_table(), get_eupath_ortholog_table() get_eupath_text()
##  Most of these are implicitly tested via make_eupath_organismdbi().

testing <- sm(download_eupath_metadata())
expected <- c(300, 19)
actual <- dim(testing)
test_that("Is the eupathdb metadata the expected size?", {
  expect_equal(expected, actual)
})

## You know what, making all of these will take days, just pull a random 2
random_org <- try(random::randomNumbers(n=1, min=1, max=nrow(testing), col=1), silent=TRUE)
chosen <- 1
if (class(random_org) == "try-error") {
  chosen <- sample(1:nrow(testing), 1)
} else {
  chosen <- random_org[, 1]
}
eupath_names <- testing[["Species"]]
species <- eupath_names[chosen]
message(paste0("\nGoing to attempt making packages for: ", species))
eupath_test <- make_eupath_organismdbi(species=species, metadata=testing, reinstall=TRUE)
if (class(eupath_test) != "try-error") {
  test_that("Did the organismdbi get installed?", {
    expect_true(eupath_test[["organdb_name"]] %in% installed.packages())
  })
} else {
  message(paste0("Creation of orgdb/etc failed for ", species))
}
bsgenome_test <- sm(make_eupath_bsgenome(species, reinstall=TRUE))
test_that("Did the bsgenome get installed?", {
  expect_true(bsgenome_test[["bsgenome_name"]] %in% installed.packages())
})
text_test <- sm(post_eupath_table(species=species, metadata=testing))
test_that("Do we get interesting text data?", {
  expect_gt(nrow(text_test), 1000)
  expect_gt(ncol(text_test), 10)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 010annotation_eupathdb.R in ", elapsed,  " seconds."))
