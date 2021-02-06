start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("245gsva.R:
")

## In order to test my GSVA functions, I am going to need a rda with some data...
## I am thinking to copy some of my human macrophage data.

## I am not sure what I want to test with this either... hmmm
hs_envir <- environment()
hs_file <- system.file("hs_expt.rda", package="hpgltools")
load(file=hs_file, envir=hs_envir)
hs_expt <- hs_envir[["expt"]]

hs_annot <- load_biomart_annotations()[["annotation"]]
rownames(hs_annot) <- make.names(hs_annot[["ensembl_gene_id"]], unique=TRUE)
drop <- grepl(pattern="\\.\\d+$", x=rownames(hs_annot))
hs_annot <- hs_annot[!drop, ]

fData(hs_expt[["expressionset"]]) <- hs_annot

hs_filt <- normalize_expt(hs_expt, filter="cv")
annotation(hs_filt[["expressionset"]]) <- "org.Hs.eg.db"
gsva_result <- sm(simple_gsva(hs_filt))

actual <- head(as.numeric(exprs(gsva_result[["gsva"]])["WINTER_HYPOXIA_UP", ]))
expected <- c(0.2811093, 0.2914057, 0.2990333, 0.2999755, 0.3044319, 0.2981790)
test_that("Do we get an expected gsva result?", {
  expect_equal(actual, expected, tolerance=0.001)
})

gsva_expt <- gsva_result[["expt"]]
gsva_dis <- plot_sample_heatmap(gsva_expt)
test_that("Can we plot a gsva result?", {
  expect_equal("recordedplot", class(gsva_dis))
})

xcell_result <- simple_xcell(expt=hs_filt, column="cds_length")
test_that("We get some expected results from xCell?", {
  expect_equal("recordedplot", class(xcell_result[["heatmap"]])[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_all.R in ", elapsed,  " seconds."))
