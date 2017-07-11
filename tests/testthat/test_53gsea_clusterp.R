start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("53gsea_clusterp.R: Does clusterProfiler work?\n")

load("gsea_siggenes.rda")

if (!identical(Sys.getenv("TRAVIS"), "true")) {
    dmel_cp <- sm(simple_clusterprofiler(z_sig_genes, table,
                                         orgdb="org.Dm.eg.db",
                                         fc_column="untreated"))

    expected <- c("GO:0000975", "GO:0000981", "GO:0000982",
                  "GO:0001067", "GO:0001071", "GO:0001871")
    actual <- head(sort(dmel_cp$enrich_go$MF_all$ID))
    test_that("Does the set of MF_all have the expected IDs?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0000082", "GO:0000122", "GO:0000186",
                  "GO:0001558", "GO:0001676", "GO:0002118")
    actual <- head(sort(dmel_cp$enrich_go$BP_all$ID))
    test_that("Does the set of BP_all have the expected IDs?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0005811", "GO:0030894", "GO:0034708",
                  "GO:0035097", "GO:0043596", "GO:0043601")
    actual <- head(sort(dmel_cp$enrich_go$CC_all$ID))
    test_that("Does the set of CC_all have the expected IDs?", {
        expect_equal(expected, actual)
    })

    ## Gather some scores
    expected <- c(0.01822725, 0.01822725, 0.01822725, 0.01822725, 0.01822725, 0.01822725)
    actual <- head(sort(dmel_cp$enrich_go$MF_sig$p.adjust))
    test_that("Does the set of MF_sig have the expected p.adjusts?", {
        expect_equal(expected, actual, tolerance=0.001)
    })

    expected <- c(1.546513e-06, 1.546513e-06, 1.546513e-06,
                  9.483361e-06, 9.483361e-06, 9.772667e-06)
    actual <- head(sort(dmel_cp$enrich_go$BP_sig$p.adjust))
    test_that("Does the set of BP_sig have the expected p.adjusts?", {
        expect_equal(expected, actual, tolerance=0.00001)
    })

    expected <- c(0.01689551)
    actual <- head(sort(dmel_cp$enrich_go$CC_sig$p.adjust))
    test_that("Does the set of CC_sig have the expected p.adjusts?", {
        expect_equal(expected, actual, tolerance=0.0001)
    })

    expected <- c("dme00240", "dme00410", "dme00564", "dme00280",
                  "dme00640", "dme00650", "dme00230", "dme04140")
    actual <- rownames(dmel_cp$kegg_data$kegg_sig)
    test_that("Did cp pick up consistent KEGG categories?", {
        expect_equal(expected, actual)
    })

    expected <- c("GO:0000728", "GO:0000734", "GO:0000742",
                  "GO:0000749", "GO:0000751", "GO:0000754")
    actual <- head(sort(rownames(dmel_cp$group_go$BP)))
    test_that("Does cp find consistent group_go categories in biological processes?", {
        expect_equal(expected, actual)
    })
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 53gsea_clusterp.R in ", elapsed,  " seconds."))
tt <- clear_session()
