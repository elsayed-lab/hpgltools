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

mgas_norm <- s_p(normalize_expt(mgas_expt, transform="log2", norm="quant", convert="cpm", filter=TRUE, batch="combat_scale", low_to_zero=TRUE))$result
test_that("Are the expt notes and state maintained?", {
    expect_match(object=mgas_norm$notes, regexp="log2\\(combat_scale\\(cpm\\(quant\\(filter\\(data\\)\\)\\)\\)\\)")
    expect_equal("cbcb", mgas_norm$state$filter)
    expect_equal("quant", mgas_norm$state$normalization)
    expect_equal("cpm", mgas_norm$state$conversion)
    expect_equal("combat_scale", mgas_norm$state$batch)
    expect_equal("log2", mgas_norm$state$transform)
})

forced_pairwise <- all_pairwise(mgas_expt)

mgas_data <- hpgltools::gbk2txdb()
actual_width <- GenomicRanges::width(mgas_data$seq)
expected_width = 1895017
actual_exons <- as.data.frame(mgas_data$exons)
expected_num_exons <- 1845
actual_num_exons <- nrow(actual_exons)
expected_gene_names <- c("dnaA", "dnaN", NA, "pth", "trcF", NA)
actual_gene_names <- head(actual_exons$gene)

test_that("Can I extract the chromosome sequence from a genbank file?", {
    expect_equal(expected_width, actual_width)
    expect_equal(expected_num_exons, actual_num_exons)
    expect_equal(expected_gene_names, actual_gene_names)
})

actual_microbe_ids <- as.character(hpgltools::get_microbesonline_ids("pyogenes MGAS5005"))
expected_microbe_ids <- c("293653", "Streptococcus pyogenes MGAS5005")
test_that("Can I get data from microbesonline?", {
    expect_equal(expected_microbe_ids, actual_microbe_ids)
})

mgas_df <- hpgltools::get_microbesonline_annotation(expected_microbe_ids[[1]])
actual_mgas_names <- as.character(head(mgas_df$name))
expected_mgas_names <- c("dnaA","dnaN","SPy0004","SPy0006","pth","trcF")
test_that("Did the mgas annotations download?", {
    expect_equal(expected_mgas_names, actual_mgas_names)
})

circos_test <- circos_prefix()
circos_kary <- circos_karyotype("mgas", length=actual_width)
##circos_plus <- circos_plus_minus(table, circos_test)

