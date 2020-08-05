start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("180gene_ontology_enrichment.R:
  123456789012345678901234567890123456789\n")
## 2017-12, exported functions in ontology_cluster_profiler: simple_clusterprofiler

## hmm I think I should split that up into separate functions for the various things it can do.
cb_sig <- environment()
load(file="test_065_significant.rda", envir=cb_sig)
cb_combined <- environment()
load(file="test_065_combined.rda", envir=cb_combined)

##ups <- cb_sig[["cb_sig"]][["ups"]][[1]]
##all <- cb_sig[["test_condbatch"]][["all_tables"]][[1]]
## It looks like I messed up the save.
ups <- cb_sig[["limma"]][["ups"]][[1]]
combined <- cb_combined[["test_condbatch_combined"]]
table <- combined[["data"]][["data"]]

## 2020-07 It looks like AnnotationHub does not currently have pombe...
## yet another victim to 2020.

## Gather the pombe annotation data.
## tmp <- library(AnnotationHub)
## ah <- AnnotationHub()
##  orgdbs <- AnnotationHub::query(ah, "OrgDb")
## sc_orgdb <- query(ah, c("OrgDB", "Saccharomyces"))
## AH67545 | org.Sc.sgd.db.sqlite3
## sc_orgdb
## pombe <- sc_orgdb[[1]]
## pombe <- orgdb_from_ah(species="^Schizosaccharomyces pombe$")

pombe_expt <- make_pombe_expt()
pombe_lengths <- fData(pombe_expt)[, c("ensembl_gene_id", "cds_length")]
colnames(pombe_lengths) <- c("ID", "length")

pombe_go <- load_biomart_go(species="spombe", host="fungi.ensembl.org")[["go"]]

## Since annotationhub seems to be having difficulty, I will use the eupathdb package
## that Keith wrote and I picked up.
fungidb_metadata <- EuPathDB::download_eupath_metadata(webservice="fungidb")
pombe_entry <- EuPathDB::get_eupath_entry(species="pombe", metadata=fungidb_metadata)
pombe_org <- sm(EuPathDB::make_eupath_orgdb(entry=pombe_entry))
pkgnames <- EuPathDB::get_eupath_pkgnames(pombe_entry)
pombe_orgdb <- pkgnames[["orgdb"]]

## Note that I default to using entrez IDs, but the eupathdb does not, so change the orgdb_to argument.
cp_test <- simple_clusterprofiler(ups, de_table=table, orgdb=pombe_orgdb, orgdb_to="GID")
test_that("Did clusterprofiler provide the expected number of entries (MF group)?", {
  actual <- nrow(cp_test[["group_go"]][["MF"]])
  expected <- 155
  expect_equal(expected, actual, tolerance=2)
})

test_that("Did clusterprofiler provide the expected number of entries (BP group)?", {
  actual <- nrow(cp_test[["group_go"]][["BP"]])
  expected <- 571
  expect_equal(expected, actual, tolerance=2)
})

test_that("Did clusterprofiler provide the expected number of entries (CC group)?", {
  actual <- nrow(cp_test[["group_go"]][["CC"]])
  expected <- 745
  expect_equal(expected, actual, tolerance=2)
})

test_that("Did clusterprofiler provide the expected number of entries (MF enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["MF_all"]])
  expected <- 13
  expect_equal(expected, actual, tolerance=2)
})

test_that("Did clusterprofiler provide the expected number of entries (BP enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["BP_all"]])
  expected <- 160
  expect_equal(expected, actual, tolerance=2)
})

test_that("Did clusterprofiler provide the expected number of entries (CC enriched)?", {
  actual <- nrow(cp_test[["enrich_go"]][["CC_all"]])
  expected <- 3
  expect_equal(expected, actual, tolerance=2)
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

## I made a change to how I process goseq data which leaves is significantly less restrictive.
## and therefore requires one to come back and decide what to drop.
go_test <- simple_goseq(ups, go_db=pombe_go, length_db=pombe_lengths)

actual <- nrow(go_test[["bp_interesting"]])
expected <- 97
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
  expect_equal(expected, actual, tolerance=0.01)
})

## 19
expected <- c(4.322235e-12, 3.097237e-06, 3.451860e-05,
              1.204558e-04, 4.305725e-04, 4.308406e-04)
actual <- head(go_test[["bp_interesting"]][["over_represented_pvalue"]])
test_that("Did goseq give the expected bp_interesting?", {
  expect_equal(expected, actual, tolerance=0.01)
})

p_expected <- c(3.649345e-12, 3.097698e-06, 8.664056e-06,
                3.363297e-05, 1.433993e-04, 2.763290e-04)
p_actual <- head(go_test[["all_data"]][["over_represented_pvalue"]])
q_expected <- c(2.069544e-08, 8.783524e-03, 1.637795e-02,
                4.768314e-02, 1.626435e-01, 2.174207e-01)
q_actual <- head(go_test[["all_data"]][["qvalue"]])
cat_expected <- c("GO:0008150", "GO:0055114", "GO:0016491",
                  "GO:0003674", "GO:0010844", "GO:1904765")
cat_actual <- head(go_test[["all_data"]][["category"]])
## 202122
test_that("Did the table of all results include the expected material?", {
  expect_equal(p_expected, p_actual, tolerance=0.001)
  expect_equal(q_expected, q_actual, tolerance=0.2)
  expect_equal(cat_expected, cat_actual, tolerance=0.001)
})

top_test <- simple_topgo(ups, go_db=pombe_go, overwrite=TRUE)
cat_expected <- c("GO:0016491", "GO:0016614", "GO:0016616",
                  "GO:0004032", "GO:0008106", "GO:0010844")
cat_actual <- rownames(head(top_test[["tables"]][["mf_subset"]]))
test_that("Do we get expected catalogs from topgo?", {
  expect_equal(cat_expected, cat_actual)
})

annot_expected <- c(289, 66, 62, 6, 7, 2)
annot_actual <- head(top_test[["tables"]][["mf_subset"]][["Annotated"]])
test_that("Do we get expected annotations from topgo?", {
  expect_equal(annot_expected, annot_actual)
})

sig_actual <- head(top_test[["tables"]][["mf_subset"]][["Significant"]])
sig_expected <- c(21, 9, 8, 3, 3, 2)
test_that("Do we get expected significances from topgo?", {
  expect_equal(sig_expected, sig_actual)
})

exp_actual <- head(top_test[["tables"]][["mf_subset"]][["Expected"]])
exp_expected <- c(7.60, 1.71, 1.61, 0.15, 0.18, 0.05)
test_that("Do we get expected MF values from topgo?", {
  expect_equal(exp_expected, exp_actual, tolerance=0.2)
})

fi_actual <- head(top_test[["tables"]][["mf_subset"]][["fisher"]])
fi_expected <- c(1.6e-05, 4.3e-05, 1.8e-04, 3.1e-04, 5.3e-04, 6.5e-04)
test_that("Do we get expected fisher values from topgo?", {
  expect_equal(fi_expected, fi_actual, tolerance=0.0001)
})

ks_actual <- head(top_test[["tables"]][["mf_subset"]][["KS"]])
ks_expected <- c(0.1134, 0.1728, 0.2361, 0.0538, 0.0825, 0.0208)
test_that("Do we get expected KS values from topgo?", {
  expect_equal(ks_expected, ks_actual, tolerance=0.05)
})

el_actual <- head(top_test[["tables"]][["mf_subset"]][["EL"]])
el_expected <- c(0.344, 0.167, 0.196, 0.054, 0.083, 0.021)
test_that("Do we get expected EL values from topgo?", {
  expect_equal(el_expected, el_actual, tolerance=0.01)
})

we_actual <- head(top_test[["tables"]][["mf_subset"]][["weight"]])
we_expected <- c(0.81453, 1.00000, 0.01174, 0.00031, 1.00000, 0.00065)
test_that("Do we get expected weight values from topgo?", {
  expect_equal(we_expected, we_actual, tolerance=0.05)
})

## I think it would not be difficult for me to add a little logic to make gostats smarter
## with respect to how it finds the correct annotations.
annot <- fData(pombe_expt)
colnames(annot) <- c("txid", "txid2", "ID", "description", "type", "width",
                     "chromosome", "strand", "start", "end")
gos_test <- simple_gostats(ups, go_db=pombe_go, gff_df=annot, gff_type="protein_coding")
cat_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["GOMFID"]])
cat_expected <- c("GO:0016491", "GO:0016616", "GO:0016614",
                  "GO:0004032", "GO:0008106", "GO:0010844")
test_that("Do we get expected stuff from gostats? (cat)", {
  expect_equal(cat_expected, cat_actual)
})

test_that("Do we get expected stuff from gostats? (annot)", {
  expect_equal(annot_expected, annot_actual)
})

p_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Pvalue"]])
p_expected <- c(2.650842e-06, 5.213081e-05, 2.302690e-04,
                2.691926e-04, 4.627229e-04, 5.922791e-04)
test_that("Do we get expected stuff from gostats? (p)", {
  expect_equal(p_expected, p_actual, tolerance=0.001)
})

odd_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["OddsRatio"]])
odd_expected <- c(4.177321, 7.763441, 7.262411, 40.938776, 30.696429, Inf)
test_that("Do we get expected stuff from gostats? (odd)", {
  expect_equal(odd_expected, odd_actual, tolerance=0.06)
})

exp_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["ExpCount"]])
exp_expected <- c(5.64382896, 1.27599611, 1.34961127, 0.14723032, 0.17176871, 0.04907677)
test_that("Do we get expected stuff from gostats? (exp)", {
  expect_equal(exp_expected, exp_actual, tolerance=0.01)
})

count_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Count"]])
count_expected <- c(19, 8, 8, 3, 3, 2)
test_that("Do we get expected stuff from gostats? (count)", {
  expect_equal(count_expected, count_actual)
})

size_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Size"]])
size_expected <- c(229, 52, 55, 6, 7, 2)
test_that("Do we get expected stuff from gostats? (size)", {
  expect_equal(size_expected, size_actual)
})

gprof_test <- simple_gprofiler(sig_genes=ups, species="spombe")
gprof_table <- gprof_test[["go"]]
actual_dim <- dim(gprof_table)
expected_dim <- c(35, 14)
test_that("Does gprofiler provide some expected tables?", {
  expect_equal(actual_dim, expected_dim)
})

actual_go <- head(sort(gprof_table[["term.id"]]))
expected_go <- c("GO:0001678", "GO:0006884", "GO:0007186",
                 "GO:0007187", "GO:0007188", "GO:0007189")
test_that("Does gprofiler give some expected GO categories?", {
  expect_equal(actual_go, expected_go)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_all.R in ", elapsed,  " seconds."))
