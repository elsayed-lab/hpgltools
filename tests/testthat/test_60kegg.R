library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

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
require.auto("org.Dm.eg.db")
suppressMessages(library(org.Dm.eg.db))

x <- org.Dm.egFLYBASE
mapped_genes <- mappedkeys(org.Dm.egFLYBASECG)
flybase_ids <- as.data.frame(org.Dm.egFLYBASE[mapped_genes])
cg_ids <- as.data.frame(org.Dm.egFLYBASECG[mapped_genes])
mapped_ids <- merge(flybase_ids, cg_ids, by="gene_id")

limma_result <- limma$hpgl_limma
all_genes <- limma_result$all_tables[[3]]
all_genes <- merge(x=all_genes, y=mapped_ids, by.x="row.names", by.y="flybase_id", all.x=TRUE)
sig_up <- suppressMessages(get_sig_genes(all_genes, z=2)$up_genes)
all_ids <- paste0("Dmel_", all_genes[["flybase_cg_id"]])
sig_ids <- paste0("Dmel_", sig_up[["flybase_cg_id"]])


pct_citrate <- suppressMessages(pct_kegg_diff(all_ids, sig_ids, organism="dme"))
test_that("Can we extract the percent differentially expressed genes in one pathway?", {
    expect_equal(19.67, pct_citrate$percent, tolerance=0.1)
})

pathways <- c("00010", "00020", "00030", "00040","nonexistent", "00051")
all_percentages <- suppressMessages(pct_all_kegg(all_ids, sig_ids, pathways=pathways, organism="dme"))

expected_percentages <- c(7.547, 9.302, 0.000, 13.040, NA, 3.704)
actual_percentages <- all_percentages[["percent"]]
test_that("Can we extract the percent differentially expressed genes from multiple pathways?", {
    expect_equal(expected_percentages, actual_percentages, tolerance=0.1)
})

unlink("kegg_pathways", recursive=TRUE)
