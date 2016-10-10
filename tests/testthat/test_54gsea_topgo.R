library(testthat)
library(hpgltools)

context("54gsea_topgo.R: Does topGO work?\n")

if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ## Load the set of limma results and pull the significantly 'up' genes.
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sm(get_sig_genes(table, column="untreated")$up_genes)

    ## Use biomart's result to get the gene lengths etc.
    dmel_annotations <- sm(get_biomart_annotations(species="dmelanogaster"))
    ## And ontology cateogies.
    dmel_ontologies <- sm(get_biomart_ontologies(species="dmelanogaster"))

    ## Get the annotations ready to be recast as a gff file.
    dmel_annotations$strand <- ifelse(dmel_annotations$strand == "1", "+", "-")
    ## Then make them into a granges object
    dmel_granges <- GenomicRanges::makeGRangesFromDataFrame(dmel_annotations, keep.extra.columns=TRUE)
    ## I got a weird error when the column was Type and not type, I suspect though that this line is not needed.
    dmel_granges$type <- dmel_annotations$Type
    ## Recast the data frame first as a List of GRanges
    ## dmel <- rtracklayer::GenomicData(dmel_granges)
    ## Set the gff filename
    gff_file="dmel.gff"
    ## And write the entries as a gff file.  This gff file may be used by clusterprofiler, topgo, and gostats.
    dmel_gff <- sm(rtracklayer::export(object=dmel_granges, con=gff_file))

    tp_result <- sm(simple_topgo(sig_genes, gff=gff_file, goids_df=dmel_ontologies))

    test_parallel <- FALSE
    if (isTRUE(test_parallel)) {
        tp_par_result <- sm(simple_topgo(sig_genes, gff=gff_file, goids_df=dmel_ontologies, parallel=TRUE))
        test_that("Are parallel and serial topGO runs equivalent?", {
            expect_equal(tp_result, tp_par_result)
        })
        rm(tp_par_result)
    }

    ## There is some run-to-run variability in these searches.
    expected_tp_mf <- c("GO:0004252", "GO:0008236", "GO:0017171")
    actual_tp_mf <- head(tp_result$tables$mf_interesting$GO.ID, n=3)
    test_that("Are the topGO interesting results expected (MF GOIDs)?", {
        expect_equal(expected_tp_mf, actual_tp_mf)
    })
    expected_tp_bp <- c("GO:0006811", "GO:0055085", "GO:0055114")
    actual_tp_bp <- head(tp_result$tables$bp_interesting$GO.ID, n=3)
    test_that("Are the topGO interesting results expected (BP GOIDs)?", {
        expect_equal(expected_tp_bp, actual_tp_bp)
    })
    expected_tp_cc <- c("GO:0016021", "GO:0031224", "GO:0016020")
    actual_tp_cc <- head(tp_result$tables$cc_interesting$GO.ID, n=3)
    test_that("Are the topGO interesting results expected (CC GOIDs)?", {
        expect_equal(expected_tp_cc, actual_tp_cc)
    })


    tp_trees <- sm(topgo_trees(tp_result))
    expected_tree_mf_nodes <- c("GO:0001071", "GO:0003674", "GO:0003700")
    actual_tree_mf_nodes <-  head(tp_trees$mf_fisher_nodes$dag@nodes, n=3)
    test_that("Are the topGO interesting results as expected? (MF trees)?", {
        expect_equal(expected_tree_mf_nodes, actual_tree_mf_nodes)
    })

    expected_tp_mfp <- c("GO:0004252", "GO:0008236", "GO:0017171")
    actual_tp_mfp <- head(tp_result$pvalue_plots$mfp_plot_over$data$GO.ID, n=3)
    test_that("Are the topGO interesting results as expected? (MF pval)?", {
        expect_equal(expected_tp_mfp, actual_tp_mfp)
    })
}

message("\nFinished 54gsea_topgo.R")
