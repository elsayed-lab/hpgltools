start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("29de_shared.R: Do the combined differential expression searches work?\n")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

deseq <- new.env()
load("de_deseq.rda", envir=deseq)
edger <- new.env()
load("de_edger.rda", envir=edger)
limma <- new.env()
load("de_limma.rda", envir=limma)
basic <- new.env()
load("de_basic.rda", envir=basic)

## The following lines should not be needed any longer.
normalized_expt <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm", filter=TRUE))
hpgl_result <- sm(all_pairwise(normalized_expt, model_batch=TRUE, which_voom="hpgl", parallel=FALSE, edger_method="short"))

expected <- deseq$hpgl_deseq$all_tables[["untreated_vs_treated"]]
actual <- hpgl_result$deseq$all_tables[["untreated_vs_treated"]]
test_that("Do we get similar results to previous DE runs: (DESeq2)?", {
    expect_equal(expected, actual)
})

expected <- edger$hpgl_edger$all_tables[["untreated_vs_treated"]]
actual <- hpgl_result$edger$all_tables[["untreated_vs_treated"]]
test_that("Do we get similar results to previous DE runs: (edgeR)?", {
    expect_equal(expected, actual)
})

expected <- limma$hpgl_limma$all_tables[["untreated_vs_treated"]]
actual <- hpgl_result$limma$all_tables[["untreated_vs_treated"]]
test_that("Do we get similar results to previous DE runs: (limma)?", {
    expect_equal(expected, actual)
})

expected <- basic$hpgl_pasilla_basic$all_tables[["untreated_vs_treated"]]
actual <- hpgl_result$basic$all_tables[["untreated_vs_treated"]]
test_that("Do we get similar results to previous DE runs: (basic)?", {
    expect_equal(expected, actual)
})

le <- hpgl_result$comparison$comp[[1]]
ld <- hpgl_result$comparison$comp[[2]]
ed <- hpgl_result$comparison$comp[[3]]
lb <- hpgl_result$comparison$comp[[4]]
eb <- hpgl_result$comparison$comp[[5]]
db <- hpgl_result$comparison$comp[[6]]
test_that("Are the comparisons between DE tools sufficiently similar? (limma/edger)", {
    expect_gt(le, 0.95)
})
test_that("Are the comparisons between DE tools sufficiently similar? (limma/deseq)", {
    expect_gt(ld, 0.80)
})
test_that("Are the comparisons between DE tools sufficiently similar? (edger/deseq)", {
    expect_gt(ed, 0.70)
})
test_that("Are the comparisons between DE tools sufficiently similar? (limma/basic)", {
    expect_gt(lb, 0.88)
})
test_that("Are the comparisons between DE tools sufficiently similar? (edger/basic)", {
    expect_gt(eb, 0.83)
})
test_that("Are the comparisons between DE tools sufficiently similar? (deseq/basic)", {
    expect_gt(db, 0.68)
})

combined_table <- sm(combine_de_tables(hpgl_result, excel=FALSE))
expected <- c(7526, 42)
actual <- dim(combined_table$data[[1]])
test_that("Has the untreated/treated combined table been filled in?", {
    expect_equal(expected, actual)
})

sig_tables <- sm(extract_significant_genes(combined_table, according_to="all", excel=FALSE))
expected <- 118
actual <- nrow(sig_tables[["limma"]][["ups"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 94
actual <- nrow(sig_tables[["limma"]][["downs"]][[1]])
test_that("Are the limma significant downs expected?", {
    expect_equal(expected, actual)
})

expected <- 106
actual <- nrow(sig_tables[["edger"]][["ups"]][[1]])
test_that("Are the edger significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 109
actual <- nrow(sig_tables[["edger"]][["downs"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 66
actual <- nrow(sig_tables[["deseq"]][["ups"]][[1]])
test_that("Are the deseq significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 51
actual <- nrow(sig_tables[["deseq"]][["downs"]][[1]])
test_that("Are the deseq significant downs expected?", {
    expect_equal(expected, actual)
})

expected <- 45
actual <- nrow(sig_tables[["basic"]][["ups"]][[1]])
test_that("Are the basic significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 31
actual <- nrow(sig_tables[["basic"]][["downs"]][[1]])
test_that("Are the basic significant downs expected?", {
    expect_equal(expected, actual)
})

funkytown <- plot_num_siggenes(combined_table$data[[1]])
expected <- c(5.35194, 5.29788, 5.24382, 5.18976, 5.13570, 5.08164)
actual <- as.numeric(head(funkytown$up_data[[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected, actual, tolerance=0.02)
})

expected <- c(-3.22938, -3.19676, -3.16414, -3.13152, -3.09890, -3.06628)
actual <- as.numeric(head(funkytown$down_data[[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected, actual, tolerance=0.02)
})

## It is strange, when using an interactive session, the following tests complete without problem
## However, when I use make test they fail with "File does not exist."

## Ensure that the excel table printer is printing excel tables
test_keepers <- list("treatment" = c("treated","untreated"))
combined_excel <- sm(combine_de_tables(hpgl_result, excel="test_excel.xlsx", keepers=test_keepers))
test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("test_excel.xlsx"))
})

## We previously checked that we can successfully combine tables, let us now ensure that plots get created etc.
## Check that there are some venn plots in the excel workbook:
expected <- "recordedplot"
actual <- class(combined_excel[["venns"]][["treatment"]][["up_noweight"]])
test_that("Are venn plots getting generated for the excel sheets?", {
    expect_equal(expected, actual)
})

expected <- "gg"
actual <- class(combined_excel$limma_plots$treatment$scatter)[[1]]
test_that("Do we get a pretty limma scatter plot?", {
    expect_equal(expected, actual)
})

actual <- class(combined_excel$deseq_plots$treatment$scatter)[[1]]
test_that("Do we get a pretty deseq scatter plot?", {
    expect_equal(expected, actual)
})

actual <- class(combined_excel$edger_plots$treatment$scatter)[[1]]
test_that("Do we get a pretty edger scatter plot?", {
    expect_equal(expected, actual)
})

expected <- c("transcriptid", "geneid", "description", "type",
              "length", "chromosome", "strand", "start",
              "end", "limma_logfc", "deseq_logfc", "edger_logfc",
              "limma_adjp", "deseq_adjp", "edger_adjp", "limma_ave",
              "limma_t", "limma_p", "limma_b", "limma_q",
              "deseq_basemean", "deseq_lfcse", "deseq_stat", "deseq_p",
              "deseq_q", "edger_logcpm", "edger_lr", "edger_p",
              "edger_q", "basic_nummed", "basic_denmed", "basic_numvar",
              "basic_denvar", "basic_logfc", "basic_t", "basic_p",
              "basic_adjp", "fc_meta", "fc_var", "fc_varbymed",
              "p_meta", "p_var")
actual <- colnames(combined_excel$data$treatment)
test_that("Do we get expected columns from the excel sheet?", {
    expect_equal(expected, actual)
})

## Test that we can extract the significant genes and get pretty graphs
significant_excel <- sm(extract_significant_genes(combined_excel, excel="test_excel_sig.xlsx"))
test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("test_excel_sig.xlsx"))
})

## How many significant up genes did limma find?
actual <- dim(significant_excel$limma$ups$treatment)
expected <- c(94, 42)
test_that("Is the number of significant genes as expected? (limma)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel$deseq$ups$treatment)
expected <- c(51, 42)
test_that("Is the number of significant genes as expected? (deseq)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel$edger$ups$treatment)
expected <- c(109, 42)
test_that("Is the number of significant genes as expected? (edger)", {
    expect_equal(expected, actual)
})

actual <- class(significant_excel$sig_bar_plots$limma)[[1]]
expected <- "gg"
test_that("Are the significance bar plots generated? (limma)",  {
    expect_equal(expected, actual)
})

## Check to make sure that if we specify a direction for the comparison, that it is maintained.
forward_keepers <- list("treatment" = c("treated","untreated"))
reverse_keepers <- list("treatment" = c("untreated","treated"))
reverse_combined_excel <- sm(combine_de_tables(hpgl_result, keepers=reverse_keepers, excel=FALSE))
forward_combined_excel <- sm(combine_de_tables(hpgl_result, keepers=forward_keepers, excel=FALSE))
forward_fold_changes <- forward_combined_excel$data$treatment$limma_logfc
expected <- sort(forward_fold_changes)
actual <- sort(reverse_combined_excel$data$treatment$limma_logfc * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$edger_logfc)
actual <- sort(reverse_combined_excel$data$treatment$edger_logfc * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$deseq_logfc)
actual <- sort(reverse_combined_excel$data$treatment$deseq_logfc * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (deseq)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$basic_logfc)
actual <- sort(reverse_combined_excel$data$treatment$basic_logfc * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (basic)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$limma_adjp)
actual <- sort(reverse_combined_excel$data$treatment$limma_adjp)
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$edger_adjp)
actual <- sort(reverse_combined_excel$data$treatment$edger_adjp)
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel$data$treatment$deseq_adjp)
actual <- sort(reverse_combined_excel$data$treatment$deseq_adjp)
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (deseq)", {
    expect_equal(expected, actual)
})

## Make sure that MA plots from combined tables are putting the logFCs in the right direction
forward_plot <- extract_de_ma(forward_combined_excel, type="limma")
reverse_plot <- extract_de_ma(reverse_combined_excel, type="limma")
expected <- sort(forward_plot$df$logfc)
actual <- sort(reverse_plot$df$logfc * -1)
test_that("Plotting an MA plot from a combined DE table provides logFCs in the correct orientation?", {
    expect_equal(expected, actual)
})

tt <- file.remove("test_excel.xlsx")
tt <- file.remove("test_excel_sig.xlsx")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 29de_shared.R in ", elapsed,  " seconds."))
