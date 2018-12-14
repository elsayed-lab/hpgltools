start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("180gene_ontology_enrichment.R:
  12345678901234567\n")
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
##orgdbs <- AnnotationHub::query(ah, "OrgDb")
sc_orgdb <- query(ah, c("OrgDB", "Saccharomyces pombe"))
## AH67545 | org.Sc.sgd.db.sqlite3
sc_orgdb
pombe <- sc_orgdb[[2]]
##pombe <- orgdb_from_ah(species="^Schizosaccharomyces pombe$")

pombe_expt <- make_pombe_expt()
pombe_lengths <- fData(pombe_expt)[, c("ensembl_gene_id", "cds_length")]
colnames(pombe_lengths) <- c("ID", "length")

cp_test <- simple_clusterprofiler(ups, de_table=all, orgdb=pombe)
test_that("Did clusterprofiler provide the expected number of entries?", {
  ## 010203
  actual <- nrow(cp_test[["group_go"]][["MF"]])
  expected <- 155
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["group_go"]][["BP"]])
  expected <- 571
  expect_equal(expected, actual, tolerance=2)
  actual <- nrow(cp_test[["group_go"]][["CC"]])
  expected <- 745
  expect_equal(expected, actual, tolerance=2)
  ## 040506
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

pombe_go <- load_biomart_go(species="spombe", host="fungi.ensembl.org")[["go"]]
go_test <- simple_goseq(ups, go_db=pombe_go, length_db=pombe_lengths)

actual <- dim(go_test[["bp_interesting"]])
expected <- c(3, 6)
## 16 and 17
test_that("Does goseq provide a few biological processes?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

## 18
## only 1 mf interesting category
expected <- c(8.6640560513167e-06)
actual <- go_test[["mf_interesting"]][["over_represented_pvalue"]]
test_that("Did goseq give the expected mf_interesting?", {
  expect_equal(expected, actual, tolerance=0.01)
})

## 19
expected <- c(3.649345e-12, 3.097698e-06, 3.363297e-05)
actual <- go_test[["bp_interesting"]][["over_represented_pvalue"]]
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
                  "GO:0110034", "GO:0010619", "GO:0003674")
cat_actual <- head(go_test[["all_data"]][["category"]])
## 202122
test_that("Did the table of all results include the expected material?", {
  expect_equal(p_expected, p_actual, tolerance=0.001)
  expect_equal(q_expected, q_actual, tolerance=0.03)
  expect_equal(cat_expected, cat_actual, tolerance=0.001)
})

top_test <- simple_topgo(ups, go_db=pombe_go)
cat_expected <- c("GO:0016491", "GO:0016614", "GO:0016616",
                  "GO:0004032", "GO:0008106", "GO:0010844")
cat_actual <- rownames(head(top_test[["tables"]][["mf_subset"]]))
annot_expected <- c(297, 67, 63, 6, 7, 2)
annot_actual <- head(top_test[["tables"]][["mf_subset"]][["Annotated"]])
sig_actual <- head(top_test[["tables"]][["mf_subset"]][["Significant"]])
sig_expected <- c(21, 9, 8, 3, 3, 2)
exp_actual <- head(top_test[["tables"]][["mf_subset"]][["Expected"]])
exp_expected <- c(7.60, 1.71, 1.61, 0.15, 0.18, 0.05)
fi_actual <- head(top_test[["tables"]][["mf_subset"]][["fisher"]])
fi_expected <- c(1.6e-05, 4.3e-05, 1.8e-04, 3.1e-04, 5.3e-04, 6.5e-04)
ks_actual <- head(top_test[["tables"]][["mf_subset"]][["KS"]])
ks_expected <- c(0.1134, 0.1728, 0.2361, 0.0538, 0.0825, 0.0208)
el_actual <- head(top_test[["tables"]][["mf_subset"]][["EL"]])
el_expected <- c(0.426, 0.305, 0.405, 0.054, 0.083, 0.021)
we_actual <- head(top_test[["tables"]][["mf_subset"]][["weight"]])
we_expected <- c(0.84363, 0.48465, 0.35583, 0.00031, 1.00000, 0.00065)
q_actual <- head(top_test[["tables"]][["mf_subset"]][["qvalue"]])
q_expected <- c(0.04032759, 0.05407051, 0.14596579, 0.19310957, 0.26530034, 0.27036311)
test_that("Do we get expected stuff from topgo?", {
  expect_equal(cat_expected, cat_actual)
  expect_equal(annot_expected, annot_actual)
  expect_equal(sig_expected, sig_actual)
  expect_equal(exp_expected, exp_actual)
  expect_equal(fi_expected, fi_actual)
  expect_equal(ks_expected, ks_actual)
  expect_equal(el_expected, el_actual)
  expect_equal(we_expected, we_actual)
  expect_equal(q_expected, q_actual, tolerance=0.01)
})

## I think it would not be difficult for me to add a little logic to make gostats smarter
## with respect to how it finds the correct annotations.
annot <- fData(pombe_expt)
colnames(annot) <- c("txid", "txid2", "ID", "description", "type", "width",
                     "chromosome", "strand", "start", "end")
gos_test <- simple_gostats(ups, go_db=pombe_go, gff_df=annot, gff_type="protein_coding")
cat_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["GOMFID"]])
cat_expected <- c("GO:0016491", "GO:0016614", "GO:0016616",
                  "GO:0004032", "GO:0008106", "GO:0010844")
p_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Pvalue"]])
p_expected <- c(2.650842e-06, 5.213081e-05, 2.302690e-04,
                2.691926e-04, 4.627229e-04, 5.922791e-04)
odd_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["OddsRatio"]])
odd_expected <- c(4.090278, 7.134409, 6.592908, 41.081633, 30.803571, Inf)
exp_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["ExpCount"]])
exp_expected <- c(5.74697337, 1.36949153, 1.27167070,
                  0.14673123, 0.17118644, 0.04891041)
count_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Count"]])
count_expected <- c(19, 8, 7, 3, 3, 2)
size_actual <- head(gos_test[["tables"]][["mf_over_enriched"]][["Size"]])
size_expected <- c(235, 56, 52, 6, 7, 2)
test_that("Do we get expected stuff from gostats?", {
  expect_equal(cat_expected, cat_actual)
  expect_equal(annot_expected, annot_actual)
  expect_equal(p_expected, p_actual)
  expect_equal(odd_expected, odd_actual, tolerance=0.001)
  expect_equal(exp_expected, exp_actual)
  expect_equal(count_expected, count_actual)
  expect_equal(size_expected, size_actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_clusterprofiler.R in ", elapsed,  " seconds."))
