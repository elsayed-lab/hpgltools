library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Does a small bacterial RNAseq experiment load?")

mgas_data <- new.env()
cdm_data <- system.file("cdm_expt.rda", package="hpgltools")
load(cdm_data, envir=mgas_data)
rm(cdm_data)

mgas_expt <- create_expt(count_dataframe=mgas_data$cdm_counts,
                         meta_dataframe=mgas_data$cdm_metadata,
                         gene_info=mgas_data$gene_info)
rm(mgas_data)

expected_gene_names <- c("dnaA", "dnaN", "M5005_Spy_0003", "ychF", "pth", "trcF")
actual_gene_names <- head(Biobase::fData(mgas_expt$expressionset)$Name)
test_that("Did the gene information load?", {
    expect_equal(expected_gene_names, actual_gene_names)
})