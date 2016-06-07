library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
sm <- suppressMessages
sw <- suppressWarnings
ss <- function(code) { suppressMessages(suppressWarnings(code)) }
context("Do KEGGREST and pathview work?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("de_limma.rda", envir=limma)

test_orgn <- kegg_get_orgn("Drosophila melanogaster", short=FALSE)
actual_id <- as.character(test_orgn[["orgid"]])
expected_id <- c("dme", "wol")

test_that("Is it possible to look up a kegg species ID?", {
    expect_equal(expected_id, actual_id)
})

## Make a map of the weird flybase IDs FBgn to the also weird Cg ids.
dm_orgdb <- s_p(choose_orgdb("drosophila_melanogaster"))$result
mapping <- ss(orgdb_idmap(dm_orgdb, mapto=c("ENSEMBL","ENTREZID","FLYBASE","FLYBASECG","GENENAME")))
expected_flybase_mapping <- c("FBgn0040373", "FBgn0040372", "FBgn0261446", "FBgn0000316", "FBgn0005427", "FBgn0040370")
actual_flybase_mapping <- head(mapping$flybase)
expected_entrezid_mapping <- c("30970", "30971", "30972", "30973", "30975", "30976")
actual_entrezid_mapping <- head(mapping$entrezid)
test_that("Did orgdb give useful ID mappings?", {
    expect_equal(expected_entrezid_mapping, actual_entrezid_mapping)
    expect_equal(expected_flybase_mapping, actual_flybase_mapping)
})

limma_result <- limma$hpgl_limma
all_genes <- limma_result$all_tables[[3]]
all_genes <- merge(x=all_genes, y=mapping, by.x="row.names", by.y="flybase", all.x=TRUE)
sig_up <- s_p(get_sig_genes(all_genes, z=2)$up_genes)$result
all_ids <- paste0("Dmel_", all_genes[["flybasecg"]])
sig_ids <- paste0("Dmel_", sig_up[["flybasecg"]])

pct_citrate <- s_p(sm(pct_kegg_diff(all_ids, sig_ids, organism="dme")))$result
test_that("Can we extract the percent differentially expressed genes in one pathway?", {
    expect_equal(19.67, pct_citrate$percent, tolerance=0.1)
})

pathways <- c("00010", "00020", "00030", "00040","nonexistent", "00051")
all_percentages <- s_p(sm(pct_all_kegg(all_ids, sig_ids, pathways=pathways, organism="dme")))$result

expected_percentages <- c(7.547, 9.302, 0.000, 13.040, NA, 3.704)
actual_percentages <- all_percentages[["percent"]]
test_that("Can we extract the percent differentially expressed genes from multiple pathways?", {
    expect_equal(expected_percentages, actual_percentages, tolerance=0.1)
})

unlink("kegg_pathways", recursive=TRUE)
