library(testthat)
library(hpgltools)
context("53gsea_clusterp.R: Does clusterProfiler work?\n")

if (!identical(Sys.getenv("TRAVIS"), "true")) {

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

    expected <- c("GO:0055085", "GO:0006811", "GO:0006820", "GO:0030001", "GO:0006030", "GO:0006812")
    actual <- head(dmel_cp$enrich_go$BP_all$ID)
    test_that("Does the set of BP_all have the expected IDs?", {
        expect_equal(expected, actual)
    })
    expected <- c("GO:0055085", "GO:0006811")
    actual <- dmel_cp$enrich_go$BP_sig$ID
    test_that("Does the set of BP_sig have the expected IDs?", {
        expect_equal(expected, actual)
    })
    expected <- c("GO:0016052", "GO:0019991", "GO:0034329", "GO:0035151", "GO:0035159", "GO:0044724")
    actual <- head(dmel_cp$gse_go$BP_all$ID)
    test_that("Does the set of BP_gsea_all have the expected IDs?", {
        expect_equal(expected, actual)
    })
    expected <- c("dme00983", "dme00500", "dme00860", "dme00903", "dme00053", "dme00830")
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
}

message("\nFinished 53gsea_clusterp.R")
