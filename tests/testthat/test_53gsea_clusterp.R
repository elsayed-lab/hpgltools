start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("53gsea_clusterp.R: Does clusterProfiler work?\n")

## Load the set of limma results and pull the significantly 'up' genes.
limma <- new.env()
load("de_limma.rda", envir=limma)
table <- limma$hpgl_table
sig_genes <- sm(get_sig_genes(table, column="untreated")$up_genes)

## This information should be available through the annotation tests and passed here.
## Use biomart's result to get the gene lengths etc.
dmel_annotations <- sm(get_biomart_annotations(species="dmelanogaster"))
## And ontology cateogies.
dmel_ontologies <- sm(get_biomart_ontologies(species="dmelanogaster"))
dmel_cp <- sm(simple_clusterprofiler(sig_genes, table, orgdb="org.Dm.eg.db", fc_column="untreated"))

expected <- c("GO:0006952", "GO:0098542", "GO:0006260",
              "GO:0043207", "GO:0051707", "GO:0009607")
actual <- head(dmel_cp$enrich_go$BP_all$ID)
test_that("Does the set of BP_all have the expected IDs?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0018022", "GO:0002812", "GO:0002813",
              "GO:0007568", "GO:0048813", "GO:0007409")
actual <- tail(dmel_cp$enrich_go$BP_sig$ID)
test_that("Does the set of BP_sig have the expected IDs?", {
    expect_equal(expected, actual)
})

## Every time I run this stupid tests it switches between 1 of 2 sets of ontologies,  WTF!?
##expected <- c("GO:1903076", "GO:1903146", "GO:1903532",
##              "GO:1903729", "GO:1904375", "GO:1904950")
##actual <- tail(dmel_cp$gse_go$BP_all$ID)
##test_that("Does the set of BP_gsea_all have the expected IDs?", {
##    expect_equal(expected, actual)
##})

expected <- c("dme00240", "dme00564", "dme00280")
actual <- head(dmel_cp$kegg_data$kegg_sig$ID)
test_that("Does the set of KEGG data have the expected IDs?", {
    expect_equal(expected, actual)
})

expected <- "gg"
actual <- class(dmel_cp$plots$dot_all_bp)[[1]]
test_that("Did we get a valid dotplot for BP data?", {
    expect_equal(expected, actual)
})

actual <- class(dmel_cp$plots$ego_all_bp)[[1]]
test_that("Did we get a valid barplot for BP data?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 53gsea_clusterp.R in ", elapsed,  " seconds."))
