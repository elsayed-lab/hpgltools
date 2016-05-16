library(testthat)
library(hpgltools)

context("Does goseq work?")

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sp(get_sig_genes(table, column="untreated")$up_genes)$result

    dmel_annotations <- sp(get_biomart_annotations(species="dmelanogaster"))$result
    dmel_ontologies <- sp(get_biomart_ontologies(species="dmelanogaster"))$result
    dmel_lengths <- dmel_annotations[, c("geneID", "length")]
    colnames(dmel_lengths) <- c("ID","width")

    goseq_result <- sp(simple_goseq(sig_genes, lengths=dmel_lengths, goids_df=dmel_ontologies))$result

    expected_goseq_mf <- c("GO:0004252", "GO:0003824", "GO:0043169", "GO:0020037", "GO:0016918", "GO:0008146")
    actual_goseq_mf <- head(rownames(goseq_result$mf_interesting))
    expected_goseq_bp <- c("GO:0006508", "GO:0055085", "GO:0055114", "GO:0048565", "GO:0008152", "GO:0006030")
    actual_goseq_bp <- head(rownames(goseq_result$bp_interesting))
    expected_goseq_cc <- c("GO:0005615", "GO:0016021", "GO:0005576", "GO:0016020", "GO:0043231", "GO:0031012")
    actual_goseq_cc <- head(rownames(goseq_result$cc_interesting))
    expected_goseq_mfp <- c(0.2332155, 0.1883853, 0.5333333, 0.2279412, 0.4375000, 0.2129630)
    actual_goseq_mfp <- head(goseq_result$pvalue_plots$mfp_plot_over$data$score)
    expected_goseq_bpp <- c(0.1866913, 0.2135922, 0.1853547, 0.1791444, 0.2075472, 0.2777778)
    actual_goseq_bpp <- head(goseq_result$pvalue_plots$bpp_plot_over$data$score)
    expected_goseq_ccp <- c(0.1822785, 0.1783088, 0.1606498, 0.1739675, 0.2471910, 0.1843972)
    actual_goseq_ccp <- head(goseq_result$pvalue_plots$ccp_plot_over$data$score)

    test_that("Are the goseq interesting results as expected?", {
        expect_equal(expected_goseq_mf, actual_goseq_mf)
        expect_equal(expected_goseq_bp, actual_goseq_bp)
        expect_equal(expected_goseq_cc, actual_goseq_cc)
        expect_equal(expected_goseq_mfp, actual_goseq_mfp, tolerance=0.0001)
        expect_equal(expected_goseq_bpp, actual_goseq_bpp, tolerance=0.0001)
        expect_equal(expected_goseq_ccp, actual_goseq_ccp, tolerance=0.0001)
    })

}
