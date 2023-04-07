start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("180gene_ontology_enrichment.R:
           1         2         3         4
  1234567890123456789012345678901234567890\n")
## 2021-02, Removing tests for explicit ontology categories, because they keep changing.
## 2017-12, exported functions in ontology_cluster_profiler: simple_clusterprofiler

## hmm I think I should split that up into separate functions for the various things it can do.
cb_sig <- environment()
load(file = "test_065_significant.rda", envir = cb_sig)
cb_combined <- environment()
load(file = "test_065_combined.rda", envir = cb_combined)

ups <- cb_sig[["deseq"]][["ups"]][[2]]
combined <- cb_combined[["test_condbatch_combined"]]
table <- combined[["data"]][["somesig"]]

## 2020-07 It looks like AnnotationHub does not currently have pombe...
## yet another victim to 2020.

## Gather the pombe annotation data.
## pombe_orgdb <- orgdb_from_ah(species = "^Schizosaccharomyces pombe$")
## Here is a fallback for if/when annotationhub is not working:
if (! "EuPathDB" %in% installed.packages()) {
  devtools::install_github("abelew/EuPathDB", force = TRUE)
}
fungidb_metadata <- sm(EuPathDB::download_eupath_metadata(webservice = "fungidb",
                                                          eu_version = "46"))
pombe_entry <- EuPathDB::get_eupath_entry(species = "pombe", metadata = fungidb_metadata)
pkgnames <- EuPathDB::get_eupath_pkgnames(entry = pombe_entry)
if (! pkgnames[["orgdb"]] %in% installed.packages()) {
  pombe_org <- sm(EuPathDB::make_eupath_orgdb(entry = pombe_entry))
}
pombe_orgdb <- pkgnames[["orgdb"]]

pombe_expt <- make_pombe_expt()
pombe_lengths <- fData(pombe_expt)[, c("ensembl_gene_id", "cds_length")]
colnames(pombe_lengths) <- c("ID", "length")

pombe_go <- load_biomart_go(species = "spombe", host = "fungi.ensembl.org",
                            overwrite = TRUE)[["go"]]

## Note that I default to using entrez IDs, but the eupathdb does not,
## so change the orgdb_to argument.
cp_test <- simple_clusterprofiler(ups, de_table = table, orgdb = pombe_orgdb, orgdb_to = "GID")
test_that("Did clusterprofiler provide the expected number of entries (MF group)?", {
  actual <- nrow(cp_test[["group_go"]][["MF"]])
  expected <- 155
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Did clusterprofiler provide the expected number of entries (BP group)?", {
  actual <- nrow(cp_test[["group_go"]][["BP"]])
  expected <- 571
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Did clusterprofiler provide the expected number of entries (CC group)?", {
  actual <- nrow(cp_test[["group_go"]][["CC"]])
  expected <- 745
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Did clusterprofiler provide the expected number of entries (MF enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["MF_all"]])
  expected <- 13
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Did clusterprofiler provide the expected number of entries (BP enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["BP_all"]])
  expected <- 160
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Did clusterprofiler provide the expected number of entries (CC enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["CC_all"]])
  expected <- 11
  expect_equal(expected, actual, tolerance = 2)
})

test_that("Do we get some plots?", {
  ## 07 - 15
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

cp_written <- write_cp_data(cp_test, excel = "test_cp_write.xlsx")
test_that("Were we able to write clusterprofiler data?", {
  expect_true(file.exists("test_cp_write.xlsx"))
})
if (file.exists("test_cp_write.xlsx")) {
  removed <- file.remove("test_cp_write.xlsx")
  removed <- unlink("test_cp_write", recursive = TRUE)
}

## } ## End checking for github actions
go_test <- simple_goseq(ups, go_db = pombe_go, length_db = pombe_lengths)
## As of 202205, there are 3 new elements in the result list;
## mf_enrich, bp_enrich, and cc_enrich.  These should make it possible to strip out
## all of my ontology plotting functions in favor of clusterprofiler/enrichplot's.

actual <- nrow(go_test[["bp_interesting"]])
expected <- 90
## 16 and 17
test_that("Does goseq provide a few biological processes?", {
  expect_gt(actual, expected)
})

## 18
## only 1 mf interesting category
expected <- c(8.481665e-06, 2.841514e-04, 4.305725e-04,
              2.514813e-03, 3.659096e-03, 3.659096e-03)
actual <- head(go_test[["mf_interesting"]][["over_represented_pvalue"]])
test_that("Did goseq give the expected mf_interesting?", {
  expect_equal(expected, actual, tolerance = 0.01)
})

## 19
expected <- c(4.322235e-12, 3.097237e-06, 3.451860e-05,
              1.204558e-04, 4.305725e-04, 4.308406e-04)
actual <- head(go_test[["bp_interesting"]][["over_represented_pvalue"]])
test_that("Did goseq give the expected bp_interesting?", {
  expect_equal(expected, actual, tolerance = 0.01)
})

p_expected <- c(3.649345e-12, 3.097698e-06, 8.664056e-06,
                3.363297e-05, 1.433993e-04, 2.763290e-04)
p_actual <- head(go_test[["all_data"]][["over_represented_pvalue"]])
q_expected <- c(2.069544e-08, 8.783524e-03, 1.637795e-02,
                4.768314e-02, 1.626435e-01, 2.174207e-01)
q_actual <- head(go_test[["all_data"]][["qvalue"]])
cat_expected <- c("GO:0008150", "GO:0055114", "GO:0016491",
                  "GO:0003674", "GO:0010844", "GO:1904765")
cat_actual <- go_test[["all_data"]][["category"]]
## 202122
test_that("Did the table of all results include the expected material?", {
  expect_equal(p_expected, p_actual, tolerance = 0.001)
  expect_equal(q_expected, q_actual, tolerance = 0.2)
  expect_equal(5, sum(cat_expected %in% cat_actual))
})

## Test if the goseq enrichResult elements are plottable.
test_enrichplot <- enrichplot::dotplot(go_test[["mf_enrich"]])
expected <- "gg"
test_that("We can enrichplot::dotplot the new goseq outputs?", {
  expect_true(expected %in% class(test_enrichplot))
})

test_enrichplot <- enrichplot::cnetplot(go_test[["mf_enrich"]])
test_that("We can enrichplot::cnetplot the new goseq outputs?", {
  expect_true(expected %in% class(test_enrichplot))
})

go_written <- write_goseq_data(go_test, excel = "test_go_write.xlsx")
test_that("Were we able to write goseq data?", {
  expect_true(file.exists("test_go_write.xlsx"))
})
if (file.exists("test_go_write.xlsx")) {
  removed <- file.remove("test_go_write.xlsx")
  removed <- unlink("test_go_write", recursive = TRUE)
}

top_test <- simple_topgo(ups, go_db = pombe_go, overwrite = TRUE)
cat_expected <- c("GO:0016491", "GO:0016614", "GO:0016616",
                  "GO:0004032", "GO:0008106", "GO:0010844")
cat_actual <- rownames(top_test[["tables"]][["mf_subset"]])
test_that("Do we get expected catalogs from topgo?", {
  expect_equal(6, sum(cat_expected %in% cat_actual))
})

top_written <- write_topgo_data(top_test, excel = "test_topgo_write.xlsx")
test_that("Were we able to write topgo data?", {
  expect_true(file.exists("test_topgo_write.xlsx"))
})
if (file.exists("test_topgo_write.xlsx")) {
  removed <- file.remove("test_topgo_write.xlsx")
  removed <- unlink("test_topgo_write", recursive = TRUE)
}

## I think it would not be difficult for me to add a little logic to make gostats smarter
## with respect to how it finds the correct annotations.
annot <- fData(pombe_expt)
colnames(annot) <- c("txid", "txid2", "ID", "description", "type", "width",
                     "chromosome", "strand", "start", "end")
gos_test <- simple_gostats(ups, go_db = pombe_go, gff_df = annot, gff_type = "protein_coding")
cat_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["GOMFID"]])
cat_expected <- c("GO:0016491", "GO:0016614", "GO:0016616",
                  "GO:0008106", "GO:0004033", "GO:0004032")
test_that("Do we get expected stuff from gostats? (cat)", {
  expect_equal(6, sum(cat_expected %in% cat_actual))
})

gos_written <- write_gostats_data(gos_test, excel = "test_gostats_write.xlsx")
test_that("Were we able to write gostats data?", {
  expect_true(file.exists("test_gostats_write.xlsx"))
})
if (file.exists("test_gostats_write.xlsx")) {
  removed <- file.remove("test_gostats_write.xlsx")
  removed <- unlink("test_gostats_write", recursive = TRUE)
}

gprof_test <- simple_gprofiler(sig_genes = ups, species = "spombe")
gprof_table <- gprof_test[["GO"]]
actual_dim <- nrow(gprof_table)
## expected_dim <- c(67, 14)
expected_dim <- 30
test_that("Does gprofiler provide some expected tables?", {
  expect_gt(actual_dim, expected_dim)
})

actual_go <- head(sort(gprof_table[["term_id"]]))
##expected_go <- c("GO:0001678", "GO:0006884", "GO:0007186",
##                 "GO:0007187", "GO:0007188", "GO:0007189")
expected_go <- c("GO:0000407", "GO:0004032", "GO:0004033",
                 "GO:0005946", "GO:0005975", "GO:0005984")
test_that("Does gprofiler give some expected GO categories?", {
  expect_gt(sum(actual_go %in% expected_go), 1)
})

gp_written <- write_gprofiler_data(gprof_test, excel = "test_gp_write.xlsx")
test_that("Were we able to write gprofiler data?", {
  expect_true(file.exists("test_gp_write.xlsx"))
})
if (file.exists("test_gp_write.xlsx")) {
  removed <- file.remove("test_gp_write.xlsx")
  removed <- unlink("test_gp_write", recursive = TRUE)
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 180ontology_all.R in ", elapsed,  " seconds.")
