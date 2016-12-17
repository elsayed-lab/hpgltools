start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("54gsea_topgo.R: Does topGO work?\n")

if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ## Load the set of limma results and pull the significantly 'up' genes.
    load("gsea_siggenes.rda")

    ## Set the gff filename
    gff_file="dmel.gff"
    ## And write the entries as a gff file.  This gff file may be used by clusterprofiler, topgo, and gostats.
    dmel_gff <- sm(rtracklayer::export(object=dmel_granges, con=gff_file))

    tp_result <- sm(simple_topgo(fcp_sig_genes, gff=gff_file, goids_df=dmel_ontologies))

    test_parallel <- FALSE
    if (isTRUE(test_parallel)) {
        tp_par_result <- sm(simple_topgo(fcp_sig_genes, gff=gff_file, goids_df=dmel_ontologies, parallel=TRUE))
        test_that("Are parallel and serial topGO runs equivalent?", {
            expect_equal(tp_result, tp_par_result)
        })
        rm(tp_par_result)
    }

    ## There is some run-to-run variability in these searches.
    expected <- c("GO:0000146", "GO:0000295", "GO:0001871",
                  "GO:0003824", "GO:0003974", "GO:0003978")
    actual <- head(sort(tp_result$tables$mf_interesting$GO.ID))
    test_that("Are the topGO interesting results expected (MF GOIDs)?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0000422", "GO:0001508", "GO:0001676",
                  "GO:0001708", "GO:0002118", "GO:0002121")
    actual <- head(sort(tp_result$tables$bp_interesting$GO.ID))
    test_that("Are the topGO interesting results expected (BP GOIDs)?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0005576", "GO:0005578", "GO:0005604",
                  "GO:0005637", "GO:0005639", "GO:0005811")
    actual <- head(sort(tp_result$tables$cc_interesting$GO.ID))
    test_that("Are the topGO interesting results expected (CC GOIDs)?", {
        expect_equal(expected, actual)
    })

    tp_trees <- sm(topgo_trees(tp_result))
    expected <- c("GO:0003674", "GO:0003824", "GO:0003974",
                  "GO:0004175", "GO:0004252", "GO:0004467")
    actual <-  head(sort(tp_trees$mf_fisher_nodes$dag@nodes))
    test_that("Are the topGO interesting results as expected? (MF trees)?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0003974", "GO:0004252", "GO:0004467",
                  "GO:0004556", "GO:0004742", "GO:0005044")
    actual <- head(sort(tp_result$pvalue_plots$mfp_plot_over$data$GO.ID))
    test_that("Are the topGO interesting results as expected? (MF pval)?", {
        expect_equal(expected, actual)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 54gsea_topgo.R in ", elapsed,  " seconds."))
