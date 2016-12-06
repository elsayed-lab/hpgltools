start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("55gsea_gostats.R: Does GOstats work?\n")

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
    dmel <- as.data.frame(dmel_granges)
    dmel$ID <- dmel$geneID

    gst_result <- sm(simple_gostats(sig_genes, gff_df=dmel, goids_df=dmel_ontologies, gff_type="protein_coding"))
    ## There is some run-to-run variability in these ontology searches
    expected <- c("GO:0003824", "GO:0019840", "GO:0005044")
    actual <- head(gst_result$mf_over_enriched$GOMFID, n=3)
    test_that("Are the GOstats interesting results as expected? (MF)", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0044699", "GO:0044710", "GO:0050896")
    actual <- head(gst_result$bp_over_enriched$GOBPID, n=3)
    test_that("Are the GOstats interesting results as expected? (BP)", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0000421", "GO:0044665", "GO:0005776")
    actual <- head(gst_result$cc_over_enriched$GOCCID, n=3)
    test_that("Are the GOstats interesting results as expected? (CC)", {
        expect_equal(expected, actual)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 55gsea_gostats.R in ", elapsed,  " seconds."))
