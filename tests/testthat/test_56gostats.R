library(testthat)
library(hpgltools)

context("Does GOstats work?")

if (!identical(Sys.getenv("TRAVIS"), "true")) {

    ## Load the set of limma results and pull the significantly 'up' genes.
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- s_p(get_sig_genes(table, column="untreated")$up_genes)$result

    ## Use biomart's result to get the gene lengths etc.
    dmel_annotations <- s_p(get_biomart_annotations(species="dmelanogaster"))$result
    ## And ontology cateogies.
    dmel_ontologies <- s_p(get_biomart_ontologies(species="dmelanogaster"))$result

    ## Get the annotations ready to be recast as a gff file.
    dmel_annotations$strand <- ifelse(dmel_annotations$strand == "1", "+", "-")
    ## Then make them into a granges object
    dmel_granges <- GenomicRanges::makeGRangesFromDataFrame(dmel_annotations, keep.extra.columns=TRUE)
    ## I got a weird error when the column was Type and not type, I suspect though that this line is not needed.
    dmel_granges$type <- dmel_annotations$Type
    ## Recast the data frame first as a List of GRanges
    dmel <- as.data.frame(dmel_granges)
    dmel$ID <- dmel$geneID

    gst_result <- s_p(simple_gostats(sig_genes, gff_df=dmel, goids_df=dmel_ontologies, gff_type="protein_coding"))$result
    ## There is some run-to-run variability in these ontology searches
    expected_gst_mf <- c("GO:0004252", "GO:0008236", "GO:0017171")
    actual_gst_mf <- head(gst_result$mf_over_enriched$GOMFID, n=3)
    expected_gst_bp <- c("GO:0006811", "GO:0055085", "GO:0055114")
    actual_gst_bp <- head(gst_result$bp_over_enriched$GOBPID, n=3)
    expected_gst_cc <- c("GO:0016021", "GO:0031224", "GO:0016020")
    actual_gst_cc <- head(gst_result$cc_over_enriched$GOCCID, n=3)
    test_that("Are the GOstats interesting results as expected?", {
        expect_equal(expected_gst_mf, actual_gst_mf)
        expect_equal(expected_gst_bp, actual_gst_bp)
        expect_equal(expected_gst_cc, actual_gst_cc)
    })
}
