start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("185sequence_queries.R:
  123456789012\n")

a909_fasta <- system.file("share/gbs_tnseq/sagalactiae_a909.fasta", package = "hpgltools")
number_atg <- count_nmer(a909_fasta)
expected <- 37134
actual <- as.numeric(number_atg)

test_that("Can we count patterns in sequence data?", {
  expect_equal(expected, actual)
})

gene_info <- sm(load_biomart_annotations(species = "dmelanogaster"))[["annotation"]]
info_idx <- gene_info[["gene_biotype"]] == "protein_coding"
gene_info <- gene_info[info_idx, ]
rownames(gene_info) <- make.names(gene_info[["ensembl_gene_id"]], unique = TRUE)
bsg <- "BSgenome.Dmelanogaster.UCSC.dm6"
bsg_installed <- please_install(bsg)
lib_result <- sm(requireNamespace(bsg))
att_result <- sm(try(attachNamespace(bsg), silent = TRUE))
dm_bsg <- get0(bsg)
gene_info[["chromosome_name"]] <- paste0("chr", gene_info[["chromosome_name"]])

dm_utrs <- gather_utrs_padding(dm_bsg, annot_df = gene_info, name_column = "ensembl_gene_id",
                               chr_column = "chromosome_name", start_column = "start_position",
                               end_column = "end_position", strand_column = "strand",
                               type_column = "gene_biotype", gene_type = "protein_coding", padding = 100)
expected <- sum(gene_info[["strand"]] == "+" & gene_info[["gene_biotype"]] == "protein_coding")
test_that("Do we get some sequences from the plus strand genes?", {
  expect_equal(nrow(dm_utrs[["fiveprime_plus_table"]]), expected)
  expect_equal(nrow(dm_utrs[["threeprime_plus_table"]]), expected)
  expect_equal(nrow(dm_utrs[["cds_plus_table"]]), expected)
})

expected <- sum(gene_info[["strand"]] == "-" & gene_info[["gene_biotype"]] == "protein_coding")
test_that("Do we get some sequences from the minus strand genes?", {
  expect_equal(nrow(dm_utrs[["fiveprime_minus_table"]]), expected)
  expect_equal(nrow(dm_utrs[["threeprime_minus_table"]]), expected)
  expect_equal(nrow(dm_utrs[["cds_minus_table"]]), expected)
})



end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 185sequence_queries.R in ", elapsed,  " seconds.")
