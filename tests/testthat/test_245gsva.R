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

hs_filt <- normalize_expt(hs_expt, filter="cv")
annotation(hs_filt[["expressionset"]]) <- "org.Hs.eg.db"
gsva_result <- sm(simple_gsva(hs_filt))

gsva_expt <- gsva_result[["expt"]]
gsva_dis <- plot_sample_heatmap(gsva_expt)

##gsva_pca <- plot_pca_genes(gsva_expt, pc_method="tsne", theta=0.9,
##                           iterations=10000, perplexity=50)
##expected <- "gg"
##actual <- class(gsva_pca[["plot"]])[1]
##test_that("Do we get a pca plot?", {
##  expect_equal(expected, actual)
##})

reactome_subset <- grepl(x=rownames(gsva_expt$expressionset), pattern="^REACTOME")
reactome_gsva <- gsva_expt$expressionset[reactome_subset, ]
tt <- heatmap.3(exprs(reactome_gsva), cexRow=0.1, cexCol=0.5, trace="none")

types <- c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils")
expressionset <- hs_expt$expressionset
colnames(pData(expressionset))[34:38] <- types

##gsva_intersections <- intersect_signatures(gsva_expt)

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 180ontology_all.R in ", elapsed,  " seconds."))
