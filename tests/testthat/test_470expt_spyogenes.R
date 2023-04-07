start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("470expt_spyogenes.R: Does a small bacterial RNAseq experiment load?")

mgas_data <- new.env()
cdm_data <- system.file("share/cdm_expt.rda", package = "hpgltools")
load(cdm_data, envir = mgas_data)
rm(cdm_data)

mgas_expt <- sm(create_expt(count_dataframe = mgas_data[["cdm_counts"]],
                            metadata = mgas_data[["cdm_metadata"]],
                            gene_info = mgas_data[["gene_info"]]))

expected <- c("dnaA", "dnaN", "M5005_Spy_0003", "ychF", "pth", "trcF")
actual <- head(fData(mgas_expt)[["Name"]])
## 01
test_that("Did the gene information load?", {
  expect_equal(expected, actual)
})

mgas_norm <- normalize_expt(mgas_expt, transform = "log2",
                            convert = "cbcbcpm", filter = TRUE)
## 02
test_that("Is the filter state maintained?", {
  expect_equal("cbcb", mgas_norm[["state"]][["filter"]])
})
## 03
test_that("Is the normalization state maintained?", {
  expect_equal("raw", mgas_norm[["state"]][["normalization"]])
})
## 04
test_that("Is the conversion state maintained?", {
  expect_equal("cbcbcpm", mgas_norm[["state"]][["conversion"]])
})
## 06
test_that("Is the transformation state maintained?", {
  expect_equal("log2", mgas_norm[["state"]][["transform"]])
})

mgas_norm <- normalize_expt(mgas_norm, batch = "combat_scale")
## 05
test_that("Is the batch state maintained?", {
  expect_equal("combat_scale", mgas_norm[["state"]][["batch"]])
})

mgas_pairwise <- all_pairwise(mgas_expt)
expected <- 0.64
actual <- min(mgas_pairwise[["comparison"]][["comp"]])
## 07
test_that("Do we get reasonably high similarities among the various DE tools?", {
  expect_gt(actual, expected)
})

mgas_combined <- combine_de_tables(mgas_pairwise, excel = FALSE)
mgas_sig <- extract_significant_genes(mgas_combined, excel = FALSE)
expected <- 150
actual <- nrow(mgas_sig[["deseq"]][["ups"]][["wtllcf_vs_mga1llcf"]])
## 08
test_that("Do we find some significant genes in the mga/wt fructose analysis?", {
  expect_gt(actual, expected)
})

mgas_data <- load_genbank_annotations(accession = "AE009949")
expected <- 1895017
actual <- GenomicRanges::width(mgas_data[["seq"]])  ## This fails on travis?
actual_width <- actual
## 09
test_that("Can I extract the chromosome sequence from a genbank file? (widths)", {
  expect_equal(expected, actual)
})

expected <- c(1845, 17)
actual <- dim(as.data.frame(mgas_data[["exons"]]))
## 10
test_that("Can I extract the chromosome sequence from a genbank file? (exons)", {
  expect_equal(expected, actual)
})

expected <- c("dnaA", "dnaN", NA, "pth", "trcF", NA)
actual <- head(as.data.frame(mgas_data[["genes"]])[["gene"]])
## 11
test_that("Can I extract the chromosome sequence from a genbank file? (gene names)", {
  expect_equal(expected, actual)
})

taxon <- "293653"
mgas_df <- load_microbesonline_annotations(id = taxon)
mgas_df[["sysName"]] <- gsub(pattern = "Spy_", replacement = "Spy", x = mgas_df[["sysName"]])
expected <- c("dnaA","dnaN","M5005_Spy_0003","M5005_Spy_0004","pth","trcF")
actual <- as.character(head(mgas_df[["name"]]))
## 12
test_that("Did the mgas annotations download?", {
  expect_equal(expected, actual)
})

mgas_go <- load_microbesonline_go(taxon)
colnames(mgas_go) <- c("ID", "GO")
mgas_go <- unique(mgas_go)
expected <- c(4161, 2)
actual <- dim(mgas_go)
## 13
test_that("Do we get expected gene ontology information?", {
  expect_equal(expected, actual)
})

circos_annot_df <- as.data.frame(mgas_df)
circos_annot_df <- circos_annot_df[, c("start", "stop", "strand", "COGFun")]
circos_annot_df[["chromosome"]] <- "chr1"
rownames(circos_annot_df) <- make.names(gsub(x = mgas_df[["sysName"]],
                                             pattern = "Spy_", replacement = "Spy"),
                                        unique = TRUE)

## There is no way circos will work on travis, lets be realistic.
##if (identical(Sys.getenv("HAS_CIRCOS"), "true")) {
## Plot the coefficients of latelog glucose
glucose_table <- mgas_pairwise[["limma"]][["identity_tables"]][["mga1llcg"]]
wtvmga_glucose <- mgas_pairwise[["limma"]][["all_tables"]][["wtllcg_vs_mga1llcg"]]
relevant_widths <- merge(glucose_table, mgas_df, by.x = "row.names",
                         by.y = "sysName", all.x = TRUE)
## Since genbankr died, get the gene lengths from microbesonline
relevant_widths <- suppressWarnings(as.numeric(relevant_widths[["width"]]))
na_widths <- is.na(relevant_widths)
relevant_widths[na_widths] <- 0

na_cog <- is.na(circos_annot_df[["COGFun"]])
circos_annot_df[na_cog, "COGFun"] <- "Z"

circos_test <- circos_prefix(circos_annot_df, name = "mgas",
                             chr_column = "chromosome", stop_column = "stop")
lengths <- 1835600
names(lengths) <- "chr1"
circos_kary <- circos_karyotype(circos_test, lengths = lengths)
circos_plus <- circos_plus_minus(circos_test)
circos_hist_ll_cg <- circos_hist(circos_test, glucose_table,
                                 colname = "logFC", outer = circos_plus)
circos_tile_wtmga <- circos_tile(circos_test, wtvmga_glucose,
                                 colname = "logFC", outer = circos_hist_ll_cg)
circos_suffix(circos_test)
circos_made <- circos_make(circos_test, target = "mgas")

expected <- "circos/mgas.svg"
test_that("Did circos run?", {
  expect_true(file.exists(expected))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = (as.numeric(end - start)))
message("\nFinished 70expt_spyogenes.R in ", elapsed,  " seconds.")
