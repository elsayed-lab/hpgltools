library(testthat)
library(hpgltools)

context("Does GOstats work?")

if (!identical(Sys.getenv("TRAVIS"), "true")) {

    ## Load the set of limma results and pull the significantly 'up' genes.
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sp(get_sig_genes(table, column="untreated")$up_genes)$result

    ## Use biomart's result to get the gene lengths etc.
    dmel_annotations <- sp(get_biomart_annotations(species="dmelanogaster"))$result
    ## And ontology cateogies.
    dmel_ontologies <- sp(get_biomart_ontologies(species="dmelanogaster"))$result

    ## Get the annotations ready to be recast as a gff file.
    dmel_annotations$strand <- ifelse(dmel_annotations$strand == "1", "+", "-")
    ## Then make them into a granges object
    dmel_granges <- GenomicRanges::makeGRangesFromDataFrame(dmel_annotations, keep.extra.columns=TRUE)
    ## I got a weird error when the column was Type and not type, I suspect though that this line is not needed.
    dmel_granges$type <- dmel_annotations$Type
    ## Recast the data frame first as a List of GRanges
    dmel <- as.data.frame(rtracklayer::GenomicData(dmel_granges))
    dmel$ID <- dmel$geneID

    gst_result <- sp(simple_gostats(sig_genes, gff_df=dmel, goids_df=dmel_ontologies, gff_type="protein_coding"))$result
    ##gst_trees <- sp(gostats_trees(gst_result))$result
    expected_gst_mf <- c("GO:0004252", "GO:0008236", "GO:0017171", "GO:0008509", "GO:0004175", "GO:0022857")
    actual_gst_mf <- head(gst_result$mf_over_enriched$GOMFID)
    expected_gst_bp <- c("GO:0006811", "GO:0055085", "GO:0055114", "GO:0030001", "GO:0006814", "GO:0006820")
    actual_gst_bp <- head(gst_result$bp_over_enriched$GOBPID)
    expected_gst_cc <- c("GO:0016021", "GO:0031224", "GO:0016020", "GO:0044425", "GO:0005576", "GO:0044421")
    actual_gst_cc <- head(gst_result$cc_over_enriched$GOCCID)
    ##expected_tree_mf_nodes <- c("GO:0001071", "GO:0003674", "GO:0003700", "GO:0003824", "GO:0004175", "GO:0004222")
    ##actual_tree_mf_nodes <-  head(gst_trees$mf_fisher_nodes$dag@nodes)
    ##expected_gst_mfp <- c("GO:0004252", "GO:0008236", "GO:0017171", "GO:0008509", "GO:0004175", "GO:0022857")
    ##actual_gst_mfp <- head(gst_result$pvalue_plots$mfp_plot_over$data$GO.ID)
    test_that("Are the GOstats interesting results as expected?", {
        expect_equal(expected_gst_mf, actual_gst_mf)
        expect_equal(expected_gst_bp, actual_gst_bp)
        expect_equal(expected_gst_cc, actual_gst_cc)
    ##    expect_equal(expected_tree_mf_nodes, actual_tree_mf_nodes)
    ##    expect_equal(expected_gst_mfp, actual_gst_mfp)
    })
}
