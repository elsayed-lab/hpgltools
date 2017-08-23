start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("70expt_spyogenes.R: Does a small bacterial RNAseq experiment load?\n")

mgas_data <- new.env()
cdm_data <- system.file("cdm_expt.rda", package="hpgltools")
load(cdm_data, envir=mgas_data)
rm(cdm_data)

mgas_expt <- sm(create_expt(count_dataframe=mgas_data[["cdm_counts"]],
                            metadata=mgas_data[["cdm_metadata"]],
                            gene_info=mgas_data[["gene_info"]]))

expected <- c("dnaA", "dnaN", "M5005_Spy_0003", "ychF", "pth", "trcF")
actual <- head(fData(mgas_expt)[["Name"]])
test_that("Did the gene information load?", {
    expect_equal(expected, actual)
})

mgas_norm <- sm(normalize_expt(mgas_expt, transform="log2", norm="quant",
                               convert="cbcbcpm", filter=TRUE,
                               batch="combat_scale", low_to_zero=TRUE))

test_that("Is the filter state maintained?", {
    expect_equal("hpgl", mgas_norm[["state"]][["filter"]])
})

test_that("Is the normalization state maintained?", {
    expect_equal("quant", mgas_norm[["state"]][["normalization"]])
})

test_that("Is the conversion state maintained?", {
    expect_equal("cbcbcpm", mgas_norm[["state"]][["conversion"]])
})

test_that("Is the batch state maintained?", {
    expect_equal("combat_scale", mgas_norm[["state"]][["batch"]])
})

test_that("Is the transformation state maintained?", {
    expect_equal("log2", mgas_norm[["state"]][["transform"]])
})

mgas_pairwise <- sm(all_pairwise(mgas_expt, parallel=FALSE))

mgas_data <- sm(gbk2txdb(accession="AE009949"))
expected <- 1895017
actual <- GenomicRanges::width(mgas_data[["seq"]])  ## This fails on travis?
actual_width <- actual
test_that("Can I extract the chromosome sequence from a genbank file? (widths)", {
    expect_equal(expected, actual)
})

expected <- c(1845, 17)
actual <- dim(as.data.frame(mgas_data[["exons"]]))
test_that("Can I extract the chromosome sequence from a genbank file? (exons)", {
    expect_equal(expected, actual)
})

expected <- c("dnaA", "dnaN", NA, "pth", "trcF", NA)
actual <- head(as.data.frame(mgas_data[["genes"]])[["gene"]])
test_that("Can I extract the chromosome sequence from a genbank file? (gene names)", {
    expect_equal(expected, actual)
})

expected <- c("293653", "Streptococcus pyogenes MGAS5005")
actual <- sm(as.character(get_microbesonline_ids("pyogenes MGAS5005")))
test_that("Can I get data from microbesonline?", {
    expect_equal(expected, actual)
})

mgas_df <- sm(load_microbesonline_annotations(expected[[1]])[[1]])
mgas_df[["sysName"]] <- gsub(pattern="Spy_", replacement="Spy", x=mgas_df[["sysName"]])
rownames(mgas_df) <- make.names(mgas_df[["sysName"]], unique=TRUE)

expected <- c("dnaA","dnaN","M5005_Spy_0003","M5005_Spy_0004","pth","trcF")
actual <- as.character(head(mgas_df[["name"]]))
test_that("Did the mgas annotations download?", {
    expect_equal(expected, actual)
})

## Plot the coefficients of latelog glucose
glucose_table <- mgas_pairwise[["limma"]][["identity_tables"]][["mga1_ll_cg"]]
wtvmga_glucose <- mgas_pairwise[["limma"]][["all_tables"]][["wt_ll_cg_vs_mga1_ll_cg"]]

## There is no way circos will work on travis, lets be realistic.
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  circos_test <- sm(circos_prefix())
  circos_kary <- sm(circos_karyotype("mgas", length=actual_width))
  circos_plus <- sm(circos_plus_minus(mgas_df, circos_test))
  circos_hist_ll_cg <- sm(circos_hist(glucose_table, mgas_df, circos_test, outer=circos_plus))
  circos_heat_ll_cf <- sm(circos_heatmap(glucose_table, mgas_df, circos_test, outer=circos_hist_ll_cg))
  circos_tile_wtmga <- sm(circos_tile(wtvmga_glucose, mgas_df, circos_test, outer=circos_heat_ll_cf))
  circos_suffix(cfgout=circos_test)
  circos_made <- sm(circos_make(target="mgas"))
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 70expt_spyogenes.R in ", elapsed,  " seconds."))
tt <- clear_session()
