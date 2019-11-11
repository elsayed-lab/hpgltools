start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("70expt_spyogenes.R: Does a small bacterial RNAseq experiment load?
  1234567890123456789012345678901\n")

mgas_data <- new.env()
cdm_data <- system.file("cdm_expt.rda", package="hpgltools")
load(cdm_data, envir=mgas_data)
rm(cdm_data)

mgas_expt <- sm(create_expt(count_dataframe=mgas_data[["cdm_counts"]],
                            metadata=mgas_data[["cdm_metadata"]],
                            gene_info=mgas_data[["gene_info"]]))

expected <- c("dnaA", "dnaN", "M5005_Spy_0003", "ychF", "pth", "trcF")
actual <- head(fData(mgas_expt)[["Name"]])
## 01
test_that("Did the gene information load?", {
    expect_equal(expected, actual)
})

mgas_norm <- normalize_expt(mgas_expt, transform="log2",
                            convert="cbcbcpm", filter=TRUE)
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

mgas_norm <- normalize_expt(mgas_norm, batch="combat_scale")
## 05
test_that("Is the batch state maintained?", {
    expect_equal("combat_scale", mgas_norm[["state"]][["batch"]])
})

mgas_pairwise <- sm(all_pairwise(mgas_expt, parallel=FALSE))

expected <- 0.64
actual <- min(mgas_pairwise[["comparison"]][["comp"]])
## 07
test_that("Do we get reasonably high similarities among the various DE tools?", {
  expect_gt(actual, expected)
})

mgas_combined <- sm(combine_de_tables(mgas_pairwise, excel=FALSE))
mgas_sig <- sm(extract_significant_genes(mgas_combined, excel=FALSE))
expected <- 209
actual <- nrow(mgas_sig[["deseq"]][["ups"]][["wt_ll_cf_vs_mga1_ll_cf"]])
## 08
test_that("Do we find some significant genes in the mga/wt fructose analysis?", {
  expect_equal(expected, actual)
})

mgas_data <- sm(load_genbank_annotations(accession="AE009949"))
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
mgas_df <- sm(load_microbesonline_annotations(id=taxon))
mgas_df[["sysName"]] <- gsub(pattern="Spy_", replacement="Spy", x=mgas_df[["sysName"]])
expected <- c("dnaA","dnaN","M5005_Spy_0003","M5005_Spy_0004","pth","trcF")
actual <- as.character(head(mgas_df[["name"]]))
## 12
test_that("Did the mgas annotations download?", {
  expect_equal(expected, actual)
})

mgas_go <- sm(load_microbesonline_go(taxon))
colnames(mgas_go) <- c("ID", "GO")
mgas_go <- unique(mgas_go)
expected <- c(4161, 2)
actual <- dim(mgas_go)
## 13
test_that("Do we get expected gene ontology information?", {
  expect_equal(expected, actual)
})

mgas_df <- fData(mgas_expt)
mgas_df <- mgas_df[, c("start", "end", "width", "strand", "gene")]
colnames(mgas_df) <- c("start", "stop", "width", "strand", "COGFun")
mgas_df[["start"]] <- suppressWarnings(as.numeric(mgas_df[["start"]]))
mgas_df[["stop"]] <- suppressWarnings(as.numeric(mgas_df[["stop"]]))
na_entries <- is.na(mgas_df[["start"]])
mgas_df <- mgas_df[!na_entries, ]

## There is no way circos will work on travis, lets be realistic.
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  ## Plot the coefficients of latelog glucose
  glucose_table <- mgas_pairwise[["limma"]][["identity_tables"]][["mga1_ll_cg"]]
  wtvmga_glucose <- mgas_pairwise[["limma"]][["all_tables"]][["wt_ll_cg_vs_mga1_ll_cg"]]
  relevant_widths <- merge(glucose_table, mgas_df, by="row.names", all.x=TRUE)
  ## Since genbankr died, get the gene lengths from microbesonline
  relevant_widths <- suppressWarnings(as.numeric(relevant_widths[["width"]]))
  na_widths <- is.na(relevant_widths)
  relevant_widths[na_widths] <- 0

  circos_test <- circos_prefix()
  circos_kary <- circos_karyotype(name="mgas", length=1899877)
  circos_plus <- circos_plus_minus(mgas_df, cfgout=circos_test)
  circos_hist_ll_cg <- circos_hist(glucose_table, mgas_df, circos_test, outer=circos_plus)
  circos_heat_ll_cf <- circos_heatmap(glucose_table, mgas_df, circos_test, outer=circos_hist_ll_cg)
  circos_tile_wtmga <- circos_tile(wtvmga_glucose, mgas_df, circos_test, outer=circos_heat_ll_cf)
  circos_suffix(cfgout=circos_test)
  ## circos_made <- circos_make(target="mgas")
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 70expt_spyogenes.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
