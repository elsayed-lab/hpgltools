start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("245gsva.R:
  123456\n")

## In order to test my GSVA functions, I am going to need a rda with some data...
## I am thinking to copy some of my human macrophage data.

## I am not sure what I want to test with this either... hmmm
hs_envir <- environment()
hs_file <- system.file("share/hs_expt.rda", package = "hpgltools")
load(file = hs_file, envir = hs_envir)
hs_expt <- hs_envir[["expt"]]
hs_expt <- subset_expt(hs_expt, subset="condition=='sh'|condition=='chr'")

hs_annot <- load_biomart_annotations()[["annotation"]]
rownames(hs_annot) <- make.names(hs_annot[["ensembl_gene_id"]], unique = TRUE)
drop <- grepl(pattern = "\\.\\d+$", x = rownames(hs_annot))
hs_annot <- hs_annot[!drop, ]

fData(hs_expt[["expressionset"]]) <- hs_annot

hs_filt <- normalize_expt(hs_expt, filter = "cv")
annotation(hs_filt[["expressionset"]]) <- "org.Hs.eg.db"
gsva_result <- simple_gsva(hs_filt, cores = 1)

## This test passes on my computer, but fails on the github infrastructure...
actual <- dim(exprs(gsva_result[["gsva"]]))
test_that("Do we get an expected gsva result?", {
  expect_equal(actual[1], 2970)
  expect_equal(actual[2], 19)
})

gsva_sig <- get_sig_gsva_categories(gsva_result, excel = NULL, model_batch = FALSE)
test_that("Can we acquire significant gsva scores?", {
  expect_equal("gg", class(gsva_sig[["score_pca"]])[1])
  expect_equal("recordedplot", class(gsva_sig[["score_plot"]])[1])
})

## Test making geneset Collections.
cb_sig <- environment()
load(file = "test_065_significant.rda", envir = cb_sig)
ups <- cb_sig[["deseq"]][["ups"]][["somesig"]]
downs <- cb_sig[["deseq"]][["downs"]][["somesig"]]
sig_gsc <- make_gsc_from_ids(first_ids = rownames(ups), second_ids = rownames(downs),
                             researcher_name = "Idunno",
                             current_id = NULL, required_id = NULL,
                             study_name = "fission", category_name = "30vs0")
test_that("We can make gene set collections from DE outputs?", {
  expect_equal(3, length(names(sig_gsc)))
  expect_gt(length(GSEABase::geneIds(sig_gsc[[1]])), 400)
})

## The following xcell call throws a warning because there is a gene with no variance.
## I do not particularly care, so I will suppress it.
xcell_result <- suppressWarnings(simple_xcell(expt = hs_filt, column = "cds_length", cores = 1))
test_that("We get some expected results from xCell?", {
  expect_equal("recordedplot", class(xcell_result[["heatmap"]])[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 180ontology_all.R in ", elapsed,  " seconds.")
