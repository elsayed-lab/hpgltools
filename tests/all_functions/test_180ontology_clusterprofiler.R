start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("180ontology_clusterprofiler.R:\n")
## 2017-12, exported functions in ontology_cluster_profiler: simple_clusterprofiler

## hmm I think I should split that up into separate functions for the various things it can do.
cb_sig <- environment()
load(file="test_065_significant.rda", envir=cb_sig)
##ups <- cb_sig[["cb_sig"]][["ups"]][[1]]
##all <- cb_sig[["test_condbatch"]][["all_tables"]][[1]]
## It looks like I messed up the save.
ups <- cb_sig[["limma"]][["ups"]][[1]]
all <- cb_sig[["limma"]][["ma_plots"]][[1]][["data"]]

## Gather the pombe annotation data.
tmp <- sm(library(AnnotationHub))
ah <- sm(AnnotationHub())
orgdbs <- sm(query(ah, "OrgDb"))
sc_orgdb <- sm(query(ah, c("OrgDB", "Saccharomyces"))) ##   AH49589 | org.Sc.sgd.db.sqlite
pombe <- sc_orgdb[[3]]

cp_test <- sm(simple_clusterprofiler(ups, de_table=all, orgdb=pombe))
test_that("Did clusterprofiler provide the expected number of entries?", {
  actual <- nrow(cp_test[["group_go"]][["MF"]])
  expected <- 155
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["group_go"]][["BP"]])
  expected <- 571
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["group_go"]][["CC"]])
  expected <- 745
  expect_equal(expected, actual, tolerance=2)

  actual <- nrow(cp_test[["enrich_go"]][["MF_all"]])
  expected <- 13
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["enrich_go"]][["BP_all"]])
  expected <- 8
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["enrich_go"]][["CC_all"]])
  expected <- 3
  expect_equal(expected, actual, tolerance=2)
})

test_that("Do we get some plots?", {
  expected <- "gg"
  actual <- class(cp_test[["plots"]][["ggo_mf_bar"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ggo_bp_bar"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ggo_cc_bar"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_all_mf"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_all_bp"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_all_cc"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_sig_mf"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_sig_bp"]])[1]
  expect_equal(expected, actual)
  actual <- class(cp_test[["plots"]][["ego_sig_cc"]])[1]
  expect_equal(expected, actual)
})


end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_clusterprofiler.R in ", elapsed,  " seconds."))
