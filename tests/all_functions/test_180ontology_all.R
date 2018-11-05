start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("180gene_ontology_enrichment.R:\n")
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
tmp <- library(AnnotationHub)
ah <- AnnotationHub()
orgdbs <- query(ah, "OrgDb")
sc_orgdb <- query(ah, c("OrgDB", "Saccharomyces"))
## AH49589 | org.Sc.sgd.db.sqlite
pombe <- sc_orgdb[[3]]

pombe_expt <- make_pombe_expt()
pombe_lengths <- fData(pombe_expt)[, c("ensembl_gene_id", "cds_length")]
colnames(pombe_lengths) <- c("ID", "length")

cp_test <- simple_clusterprofiler(ups, de_table=all, orgdb=pombe)
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

pombe_go <- load_biomart_go(species="spombe", host="fungi.ensembl.org")[["go"]]
go_test <- simple_goseq(ups, go_db=pombe_go, length_db=pombe_lengths)

actual <- dim(go_test[["bp_interesting"]])
expected <- c(3, 6)
test_that("Does goseq provide a few biological processes?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_clusterprofiler.R in ", elapsed,  " seconds."))
