start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("57gsea_kegg.R: Do KEGGREST and pathview work?\n")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("de_limma.rda", envir=limma)

test_orgn <- kegg_get_orgn("Drosophila melanogaster", short=FALSE)

actual <- as.character(test_orgn[["orgid"]])
expected <- c("dme", "wol")
test_that("Is it possible to look up a kegg species ID?", {
    expect_equal(expected, actual)
})

## Make a map of the weird flybase IDs FBgn to the also weird Cg ids.
dm_orgdb <- sm(choose_orgdb("drosophila_melanogaster"))
mapping <- sm(orgdb_idmap(dm_orgdb, mapto=c("ENSEMBL","ENTREZID","FLYBASE","FLYBASECG","GENENAME")))
expected <- c("FBgn0040373", "FBgn0040372", "FBgn0261446", "FBgn0000316", "FBgn0005427", "FBgn0040370")
actual <- head(mapping$flybase)
test_that("Did orgdb give useful ID mappings? (FBgn IDs)", {
    expect_equal(expected, actual)
})

expected <- c("30970", "30971", "30972", "30973", "30975", "30976")
actual <- head(mapping$entrezid)
test_that("Did orgdb give useful ID mappings? (entrez)", {
    expect_equal(expected, actual)
})

limma_result <- limma$hpgl_limma
all_genes <- limma_result[["all_tables"]][["untreated_vs_treated"]]
all_genes <- merge(x=all_genes, y=mapping, by.x="row.names", by.y="flybase", all.x=TRUE)
sig_up <- sm(get_sig_genes(all_genes, z=2)$up_genes)
all_ids <- paste0("Dmel_", all_genes[["flybasecg"]])
sig_ids <- paste0("Dmel_", sig_up[["flybasecg"]])

## Note, I split the result of this into percent_nodes and percent_edges
pct_citrate <- sm(pct_kegg_diff(all_ids, sig_ids, organism="dme", pathway="00500"))
expected <- 6.1
actual <- pct_citrate$percent_nodes
test_that("Can we extract the percent differentially expressed genes in one pathway?", {
    expect_equal(expected, actual, tolerance=0.1)
})

pathways <- c("00010", "00020", "00030", "00040","nonexistent", "00051")
all_percentages <- sm(pct_all_kegg(all_ids, sig_ids, pathways=pathways, organism="dme"))
expected <- c(3.704, 4.651, 4.167, 5.769, NA, 3.448)
actual <- all_percentages[["percent_nodes"]]
test_that("Can we extract the percent differentially expressed genes from multiple pathways?", {
    expect_equal(expected, actual, tolerance=0.1)
})

## Try testing out pathview
mel_id <- kegg_get_orgn("melanogaster")
rownames(sig_up) <- make.names(sig_up[["flybasecg"]], unique=TRUE)
funkytown <- sm(hpgl_pathview(sig_up, fc_column="logFC", species="dme",
                              from_list=c("CG"), to_list=c("Dmel_CG")))

expected <- c(22, 34, 5, 3, 12, 106)
actual <- head(funkytown$total_mapped_nodes)
test_that("Did pathview work? (total mapped nodes)", {
    expect_equal(expected, actual, tolerance=0.1)
})

expected <- c(17, 8, 5, 3, 9, 100)
actual <- head(funkytown$unique_mapped_nodes)
test_that("Did pathview work? (unique mapped nodes)", {
    expect_equal(expected, actual, tolerance=0.1)
})

unlink("kegg_pathways", recursive=TRUE)

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 57gsea_kegg.R in ", elapsed,  " seconds."))
