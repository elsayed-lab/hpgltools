library(testthat)
library(hpgltools)

context("Does clusterProfiler work?")

if (!identical(Sys.getenv("TRAVIS"), "true")) {

    ## Load the set of limma results and pull the significantly 'up' genes.
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- suppressMessages(get_sig_genes(table, column="untreated")$up_genes)

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
    dmel <- rtracklayer::GenomicData(dmel_granges)
    ## Set the gff filename
    gff_file="dmel.gff"
    ## And write the entries as a gff file.  This gff file may be used by clusterprofiler, topgo, and gostats.
    dmel_gff <- rtracklayer::export(object=dmel, con=gff_file)

    message("Testing cluster profiler, this will probably take a while.")
    message("Since the last time I played with clusterprofiler, it looks like some pretty big changes have happened.")
    message("It is probably time to re-do my entire clusterprofiler analysis.\n")
    cp_result <- sp(simple_clusterprofiler(sig_genes, gff=gff_file, goids=dmel_ontologies))
    cp_result <- cp_result$result

    expected_cp_mf <- c("GO:0004252", "GO:0008509", "GO:0008028", "GO:0046943", "GO:0005342", "GO:0008236")
    actual_cp_mf <- head(cp_result$mf_interesting$ID)
    expected_cp_bp <- c("GO:0006116", "GO:0006734", "GO:0032881", "GO:0070873", "GO:0006206", "GO:0006814")
    actual_cp_bp <- head(cp_result$bp_interesting$ID)
    expected_cp_cc <- c("GO:0016021", "GO:0031224", "GO:0000421", "GO:0045254", "GO:0044665", "GO:0009897")
    actual_cp_cc <- head(cp_result$cc_interesting$ID)
    expected_cp_mfp <- c(0.0001466297, 0.0003302495, 0.0003874293, 0.0004893901, 0.0007106647, 0.0010191967)
    actual_cp_mfp <- head(cp_result$pvalue_plots$mfp_plot_over$data$pvalue)

    test_that("Are the clusterprofiler interesting results as expected?", {
        expect_equal(expected_cp_mf, actual_cp_mf)
        expect_equal(expected_cp_bp, actual_cp_bp)
        expect_equal(expected_cp_cc, actual_cp_cc)
        expect_equal(expected_cp_mfp, actual_cp_mfp, tolerance=0.0001)
    })
}
