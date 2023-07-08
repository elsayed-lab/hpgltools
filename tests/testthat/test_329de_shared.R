start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("329de_shared.R: Do the combined differential expression searches work?")

pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]
deseq <- new.env()
load("324_de_deseq.rda", envir = deseq)
edger <- new.env()
load("326_de_edger.rda", envir = edger)
limma <- new.env()
load("320_de_limma.rda", envir = limma)
basic <- new.env()
load("327_de_basic.rda", envir = basic)

## The following lines should not be needed any longer.
normalized_expt <- normalize_expt(pasilla_expt, transform = "log2", norm = "quant",
                                  convert = "cbcbcpm", filter = "cbcb", thresh = 1)

## Interestingly, doParallel does not work when run from packrat.
test_keepers <- list("treatment" = c("treated", "untreated"))
hpgl_all <- all_pairwise(pasilla_expt, filter = TRUE)
hpgl_tables <- combine_de_tables(hpgl_all, keepers = test_keepers,
                                 excel = "excel_test.xlsx")

combined_excel <- hpgl_tables
test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("excel_test.xlsx"))
})
removed <- file.remove("excel_test.xlsx")

hpgl_sva_result <- all_pairwise(pasilla_expt, model_batch = "sva", which_voom = "limma",
                                limma_method = "robust", edger_method = "long",
                                edger_test = "qlr", filter = TRUE)

expected <- deseq[["hpgl_deseq"]][["all_tables"]][["untreated_vs_treated"]]
actual <- hpgl_all[["deseq"]][["all_tables"]][["untreated_vs_treated"]]
shared <- rownames(actual) %in% rownames(expected)
actual <- actual[shared, ]
shared <- rownames(expected) %in% rownames(actual)
expected <- expected[shared, ]
table_order <- rownames(expected)
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (DESeq2)?", {
    expect_equal(expected, actual, tolerance = 0.05)
})

expected <- edger[["hpgl_edger"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["edger"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (edgeR)?", {
    expect_equal(expected, actual, tolerance = 0.08)
})

expected <- limma[["hpgl_limma"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["limma"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (limma)?", {
    expect_equal(expected, actual, tolerance = 0.07)
})

expected <- basic[["hpgl_basic"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["basic"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (basic)?", {
    expect_equal(expected, actual)
})

le <- hpgl_all[["comparison"]][["comp"]]["limma_vs_edger", ]
ld <- hpgl_all[["comparison"]][["comp"]]["limma_vs_deseq", ]
ed <- hpgl_all[["comparison"]][["comp"]]["deseq_vs_edger", ]
lb <- hpgl_all[["comparison"]][["comp"]]["limma_vs_basic", ]
eb <- hpgl_all[["comparison"]][["comp"]]["edger_vs_basic", ]
db <- hpgl_all[["comparison"]][["comp"]]["deseq_vs_basic", ]
test_that("Are the comparisons between DE tools sufficiently similar? (limma/edger)", {
    expect_gt(le, 0.95)
})
test_that("Are the comparisons between DE tools sufficiently similar? (limma/deseq)", {
    expect_gt(ld, 0.95)
})
test_that("Are the comparisons between DE tools sufficiently similar? (edger/deseq)", {
    expect_gt(ed, 0.90)
})
test_that("Are the comparisons between DE tools sufficiently similar? (limma/basic)", {
    expect_gt(lb, 0.91)
})
test_that("Are the comparisons between DE tools sufficiently similar? (edger/basic)", {
    expect_gt(eb, 0.92)
})
test_that("Are the comparisons between DE tools sufficiently similar? (deseq/basic)", {
    expect_gt(db, 0.92)
})
combined_table <- combine_de_tables(hpgl_all, excel = FALSE)

expected_annotations <- c(
  "ensembltranscriptid", "ensemblgeneid",
  "description", "genebiotype",
  "cdslength", "chromosomename", "strand",
  "startposition", "endposition", "deseq_logfc",
  "deseq_adjp", "edger_logfc", "edger_adjp",
  "limma_logfc", "limma_adjp", "basic_num",
  "basic_den", "basic_numvar", "basic_denvar",
  "basic_logfc", "basic_t", "basic_p",
  "basic_adjp", "deseq_basemean", "deseq_lfcse",
  "deseq_stat", "deseq_p", "deseq_num", "deseq_den",
  "edger_logcpm", "edger_lr", "edger_p",
  "limma_ave", "limma_t", "limma_b",
  "limma_p", "limma_adjp_ihw", "deseq_adjp_ihw",
  "edger_adjp_ihw", "basic_adjp_ihw",
  "lfc_meta", "lfc_var", "lfc_varbymed", "p_meta",
  "p_var")
num_cols <- length(expected_annotations) - 1
expected <- c(7531, num_cols)
actual <- dim(combined_table[["data"]][[1]])
test_that("Has the untreated/treated combined table been filled in?", {
    expect_equal(expected, actual)
})

table <- "untreated_vs_treated"
sig_tables <- extract_significant_genes(combined_table,
                                        according_to = "all",
                                        excel = FALSE)
expected <- 90
actual <- nrow(sig_tables[["limma"]][["ups"]][[table]])
test_that("Are the limma significant ups expected?", {
    expect_gt(actual, expected)
})

expected <- 90
actual <- nrow(sig_tables[["limma"]][["downs"]][[table]])
test_that("Are the limma significant downs expected?", {
    expect_gt(actual, expected)
})

expected <- 90
actual <- nrow(sig_tables[["edger"]][["ups"]][[table]])
test_that("Are the edger significant ups expected?", {
    expect_gt(actual, expected)
})

expected <- 90
actual <- nrow(sig_tables[["edger"]][["downs"]][[table]])
test_that("Are the limma significant ups expected?", {
    expect_gt(actual, expected)
})

expected <- 90
actual <- nrow(sig_tables[["deseq"]][["ups"]][[table]])
test_that("Are the deseq significant ups expected?", {
    expect_gt(actual, expected)
})

expected <- 90
actual <- nrow(sig_tables[["deseq"]][["downs"]][[table]])
test_that("Are the deseq significant downs expected?", {
    expect_gt(actual, expected)
})

expected <- 40
actual <- nrow(sig_tables[["basic"]][["ups"]][[table]])
test_that("Are the basic significant ups expected?", {
    expect_gt(actual, expected)
})

expected <- 20
actual <- nrow(sig_tables[["basic"]][["downs"]][[table]])
test_that("Are the basic significant downs expected?", {
    expect_gt(actual, expected)
})

## I significantly changed the format of this function's output.
funkytown <- plot_num_siggenes(combined_table[["data"]][[1]])
## expected <- c(11.02373, 10.91238, 10.80103, 10.68968, 10.57833, 10.46698)
expected <- c(5.52222, 5.46644, 5.41066, 5.35488, 5.29910, 5.24332)
actual <- as.numeric(head(funkytown[["up_data"]][[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected, actual, tolerance = 0.02)
})

expected <- c("gg", "ggplot")
actual <- class(combined_excel[["plots"]][["treatment"]][["limma_scatter_plots"]][["scatter"]])
test_that("Do we get a pretty limma scatter plot?", {
    expect_equal(expected, actual)
})
actual <- class(combined_excel[["plots"]][["treatment"]][["deseq_scatter_plots"]][["scatter"]])
test_that("Do we get a pretty deseq scatter plot?", {
    expect_equal(expected, actual)
})
actual <- class(combined_excel[["plots"]][["treatment"]][["edger_scatter_plots"]][["scatter"]])
test_that("Do we get a pretty edger scatter plot?", {
    expect_equal(expected, actual)
})

table <- "treatment"
actual <- colnames(combined_excel[["data"]][[table]])
test_that("Do we get expected columns from the excel sheet?", {
    expect_equal(expected_annotations, actual)
})

## Test that we can extract the significant genes and get pretty graphs
significant_excel <- extract_significant_genes(combined_excel,
                                               excel = "excel_test_sig.xlsx")
test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("excel_test_sig.xlsx"))
})

## How many significant up genes did limma find?
actual <- dim(significant_excel[["limma"]][["ups"]][[table]])
expected_rows <- 90
expected_cols <- num_cols
test_that("Is the number of significant up genes as expected? (limma)", {
  expect_gt(actual[1], expected_rows)
  expect_equal(actual[2], expected_cols)
})

actual <- nrow(significant_excel[["deseq"]][["ups"]][[table]])
expected <- 90
test_that("Is the number of significant up genes as expected? (deseq)", {
    expect_gt(actual, expected)
})

actual <- nrow(significant_excel[["edger"]][["ups"]][[table]])
expected <- 90
test_that("Is the number of significant up genes as expected? (edger)", {
    expect_gt(actual, expected)
})

actual <- nrow(significant_excel[["limma"]][["downs"]][[table]])
expected <- 90
test_that("Is the number of significant down genes as expected? (limma)", {
    expect_gt(actual, expected)
})

actual <- nrow(significant_excel[["deseq"]][["downs"]][[table]])
expected <- 90
test_that("Is the number of significant down genes as expected? (deseq)", {
    expect_gt(actual, expected)
})

actual <- nrow(significant_excel[["edger"]][["downs"]][[table]])
expected <- 90
test_that("Is the number of significant down genes as expected? (edger)", {
  expect_gt(actual, expected)
})

actual <- class(significant_excel[["sig_bar_plots"]][["limma"]])[[1]]
expected <- "gg"
test_that("Are the significance bar plots generated? (limma)",  {
    expect_equal(expected, actual)
})

## Check to make sure that if we specify a direction for the comparison, that it is maintained.
forward_keepers <- list("treatment" = c("treated", "untreated"))
reverse_keepers <- list("treatment" = c("untreated", "treated"))
reverse_combined_excel <- combine_de_tables(hpgl_all, keepers = reverse_keepers, excel = FALSE)
forward_combined_excel <- combine_de_tables(hpgl_all, keepers = forward_keepers, excel = FALSE)
forward_fold_changes <- forward_combined_excel[["data"]][[table]][["limma_logfc"]]
expected <- sort(forward_fold_changes)
actual <- sort(reverse_combined_excel[["data"]][[table]][["limma_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["edger_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["edger_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["deseq_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["deseq_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (deseq)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["basic_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["basic_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (basic)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["limma_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["limma_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["edger_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["edger_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][[table]][["deseq_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][[table]][["deseq_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (deseq)", {
    expect_equal(expected, actual)
})

## Make sure that MA plots from combined tables are putting the logFCs in the right direction

## See that we can compare different analysis types
combined_sva <- combine_de_tables(hpgl_sva_result, excel = NULL, keepers = test_keepers)
sva_batch_test <- compare_de_results(combined_excel, combined_sva)
expected <- 0.71
actual <- sva_batch_test[["result"]][["limma"]][[table]][["logfc"]]
test_that("Do limma with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

expected <- 0.97
actual <- sva_batch_test[["result"]][["deseq"]][[table]][["logfc"]]
test_that("Do deseq with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

expected <- 0.97
actual <- sva_batch_test[["result"]][["edger"]][[table]][["logfc"]]
test_that("Do edger with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

## See if the intersection between limma, deseq, and edger is decent.
test_intersect <- intersect_significant(combined_sva, excel = NULL)
expected <- 85
actual <- nrow(test_intersect[["ups"]][[table]][["data"]][["all"]])
test_that("Do we get the expected number of agreed upon significant genes between edger/deseq/limma?", {
    expect_gt(actual, expected)
})

actual <- nrow(test_intersect[["downs"]][[table]][["data"]][["all"]])
expected <- 90
test_that("Ibid, but in the down direction?", {
    expect_gt(actual, expected)
})

actual <- sum(length(test_intersect[["ups"]][[table]][["limma"]]) +
              length(test_intersect[["ups"]][[table]][["edger"]]) +
              length(test_intersect[["ups"]][[table]][["deseq"]]))
expected <- 77
test_that("Are there very few genes observed without the others?", {
    expect_gt(actual, expected)
})

actual <- sum(length(test_intersect[["downs"]][[table]][["limma"]]) +
              length(test_intersect[["downs"]][[table]][["edger"]]) +
              length(test_intersect[["downs"]][[table]][["deseq"]]))
expected <- 300
test_that("Ibid, but down?", {
    expect_gt(actual, expected)
})

actual <- nrow(test_intersect[["ups"]][[table]][["data"]][["all"]])
expected <- 88
test_that("Do all methods have some genes in common? (up)", {
  expect_gt(actual, expected)
})

actual <- nrow(test_intersect[["downs"]][[table]][["data"]][["all"]])
expected <- 90
test_that("Do all methods have some genes in common? (down)", {
  expect_gt(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 29de_shared.R in ", elapsed,  " seconds.")
