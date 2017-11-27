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
normalized_expt <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant",
                                     convert="cbcbcpm", filter="cbcb", thresh=1))

## Interestingly, doParallel does not work when run from packrat.
test_keepers <- list("treatment" = c("treated", "untreated"))
hpgl_all <- sm(all_pairwise(pasilla_expt, parallel=FALSE,
                            keepers=test_keepers,
                            combined_excel="excel_test.xlsx"))
combined_excel <- hpgl_all[["combined"]]

test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("excel_test.xlsx"))
})

hpgl_sva_result <- sm(all_pairwise(normalized_expt, model_batch="sva", which_voom="limma",
                                   limma_method="robust", edger_method="long",
                                   edger_test="qlr", parallel=FALSE))

expected <- deseq[["hpgl_deseq"]][["all_tables"]][["untreated_vs_treated"]]
table_order <- rownames(expected)
actual <- hpgl_all[["deseq"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (DESeq2)?", {
    expect_equal(expected, actual)
})

expected <- edger[["hpgl_edger"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["edger"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (edgeR)?", {
    expect_equal(expected, actual)
})

expected <- limma[["hpgl_limma"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["limma"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (limma)?", {
    expect_equal(expected, actual)
})

expected <- basic[["hpgl_basic"]][["all_tables"]][["untreated_vs_treated"]]
expected <- expected[table_order, ]
actual <- hpgl_all[["basic"]][["all_tables"]][["untreated_vs_treated"]]
actual <- actual[table_order, ]
test_that("Do we get similar results to previous DE runs: (basic)?", {
    expect_equal(expected, actual)
})

le <- hpgl_all[["comparison"]][["comp"]][[1]]
ld <- hpgl_all[["comparison"]][["comp"]][[2]]
ed <- hpgl_all[["comparison"]][["comp"]][[3]]
lb <- hpgl_all[["comparison"]][["comp"]][[4]]
eb <- hpgl_all[["comparison"]][["comp"]][[5]]
db <- hpgl_all[["comparison"]][["comp"]][[6]]
test_that("Are the comparisons between DE tools sufficiently similar? (limma/edger)", {
    expect_gt(le, 0.96)
})
test_that("Are the comparisons between DE tools sufficiently similar? (limma/deseq)", {
    expect_gt(ld, 0.95)
})
test_that("Are the comparisons between DE tools sufficiently similar? (edger/deseq)", {
    expect_gt(ed, 0.97)
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

combined_table <- sm(combine_de_tables(hpgl_all, excel=FALSE))
expected <- c(10153, 45)
actual <- dim(combined_table[["data"]][[1]])
test_that("Has the untreated/treated combined table been filled in?", {
    expect_equal(expected, actual)
})

sig_tables <- sm(extract_significant_genes(combined_table,
                                           according_to="all",
                                           excel=FALSE))
expected <- 123
actual <- nrow(sig_tables[["limma"]][["ups"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 114
actual <- nrow(sig_tables[["limma"]][["downs"]][[1]])
test_that("Are the limma significant downs expected?", {
    expect_equal(expected, actual)
})

expected <- 141
actual <- nrow(sig_tables[["edger"]][["ups"]][[1]])
test_that("Are the edger significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 190
actual <- nrow(sig_tables[["edger"]][["downs"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 113
actual <- nrow(sig_tables[["deseq"]][["ups"]][[1]])
test_that("Are the deseq significant ups expected?", {
    expect_equal(expected, actual)
})

expected <- 109
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

funkytown <- plot_num_siggenes(combined_table[["data"]][[1]])
expected <- c(4.79556, 4.74712, 4.69868, 4.65024, 4.60180, 4.55336)
actual <- as.numeric(head(funkytown[["up_data"]][[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected, actual, tolerance=0.02)
})

expected <- c(-3.22938, -3.19676, -3.16414, -3.13152, -3.09890, -3.06628)
actual <- as.numeric(head(funkytown[["down_data"]][[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected, actual, tolerance=0.02)
})

## We previously checked that we can successfully combine tables, let us now ensure that plots get created etc.
## Check that there are some venn plots in the excel workbook:
## expected <- "recordedplot"
## actual <- class(combined_excel[["venns"]][["treatment"]][["up_noweight"]])
## test_that("Are venn plots getting generated for the excel sheets?", {
##     expect_equal(expected, actual)
## })

expected <- "gg"
actual <- class(combined_excel[["limma_plots"]][["treatment"]][["scatter"]])[[1]]
test_that("Do we get a pretty limma scatter plot?", {
    expect_equal(expected, actual)
})
actual <- class(combined_excel[["deseq_plots"]][["treatment"]][["scatter"]])[[1]]
test_that("Do we get a pretty deseq scatter plot?", {
    expect_equal(expected, actual)
})
actual <- class(combined_excel[["edger_plots"]][["treatment"]][["scatter"]])[[1]]
test_that("Do we get a pretty edger scatter plot?", {
    expect_equal(expected, actual)
})

expected <- c(
  "ensembltranscriptid", "ensemblgeneid", "version", "transcriptversion",
  "description", "genebiotype", "cdslength", "chromosomename", "strand",
  "startposition", "endposition", "limma_logfc", "limma_adjp", "deseq_logfc",
  "deseq_adjp", "edger_logfc", "edger_adjp", "limma_ave", "limma_t", "limma_b",
  "limma_p", "deseq_basemean", "deseq_lfcse", "deseq_stat", "deseq_p",
  "edger_logcpm", "edger_lr", "edger_p",
  "basic_nummed", "basic_denmed", "basic_numvar", "basic_denvar",
  "basic_logfc", "basic_t", "basic_p", "basic_adjp",
  "limma_adjp_fdr", "deseq_adjp_fdr", "edger_adjp_fdr", "basic_adjp_fdr",
  "lfc_meta", "lfc_var", "lfc_varbymed",
  "p_meta", "p_var")
actual <- colnames(combined_excel[["data"]][["treated_vs_untreated"]])
test_that("Do we get expected columns from the excel sheet?", {
    expect_equal(expected, actual)
})

## Test that we can extract the significant genes and get pretty graphs
significant_excel <- sm(extract_significant_genes(combined_excel,
                                                  excel="excel_test_sig.xlsx"))
test_that("Does combine_de_tables create an excel file?", {
    expect_true(file.exists("excel_test_sig.xlsx"))
})

## How many significant up genes did limma find?
actual <- dim(significant_excel[["limma"]][["ups"]][["treated_vs_untreated"]])
expected <- c(114, 45)
test_that("Is the number of significant up genes as expected? (limma)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel[["deseq"]][["ups"]][["treated_vs_untreated"]])
expected <- c(109, 45)
test_that("Is the number of significant up genes as expected? (deseq)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel[["edger"]][["ups"]][["treated_vs_untreated"]])
expected <- c(190, 45)
test_that("Is the number of significant up genes as expected? (edger)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel[["limma"]][["downs"]][["treated_vs_untreated"]])
expected <- c(123, 45)
test_that("Is the number of significant down genes as expected? (limma)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel[["deseq"]][["downs"]][["treated_vs_untreated"]])
expected <- c(113, 45)
test_that("Is the number of significant down genes as expected? (deseq)", {
    expect_equal(expected, actual)
})

actual <- dim(significant_excel[["edger"]][["downs"]][["treated_vs_untreated"]])
expected <- c(141, 45)
test_that("Is the number of significant down genes as expected? (edger)", {
    expect_equal(expected, actual)
})

actual <- class(significant_excel[["sig_bar_plots"]][["limma"]])[[1]]
expected <- "gg"
test_that("Are the significance bar plots generated? (limma)",  {
    expect_equal(expected, actual)
})

## Check to make sure that if we specify a direction for the comparison, that it is maintained.
forward_keepers <- list("treatment" = c("treated", "untreated"))
reverse_keepers <- list("treatment" = c("untreated", "treated"))
reverse_combined_excel <- sm(combine_de_tables(hpgl_all, keepers=reverse_keepers, excel=FALSE))
forward_combined_excel <- sm(combine_de_tables(hpgl_all, keepers=forward_keepers, excel=FALSE))
forward_fold_changes <- forward_combined_excel[["data"]][["treated_vs_untreated"]][["limma_logfc"]]
expected <- sort(forward_fold_changes)
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["limma_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["edger_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["edger_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["deseq_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["deseq_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (deseq)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["basic_logfc"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["basic_logfc"]] * -1)
test_that("When we reverse a combined_de_tables(), we get reversed results? (basic)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["limma_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["limma_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (limma)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["edger_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["edger_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (edger)", {
    expect_equal(expected, actual)
})

expected <- sort(forward_combined_excel[["data"]][["treated_vs_untreated"]][["deseq_adjp"]])
actual <- sort(reverse_combined_excel[["data"]][["untreated_vs_treated"]][["deseq_adjp"]])
test_that("When we reverse a combined_de_tables(), we get appropriate p-values? (deseq)", {
    expect_equal(expected, actual)
})

## Make sure that MA plots from combined tables are putting the logFCs in the right direction
forward_plot <- extract_de_plots(forward_combined_excel, type="limma")[["ma"]]
reverse_plot <- extract_de_plots(reverse_combined_excel, type="limma")[["ma"]]
expected <- sort(forward_plot[["df"]][["logfc"]])
actual <- sort(reverse_plot[["df"]][["logfc"]] * -1)
test_that("Plotting an MA plot from a combined DE table provides logFCs in the correct orientation?", {
    expect_equal(expected, actual)
})

## See that we can compare different analysis types
combined_sva <- sm(combine_de_tables(hpgl_sva_result,
                                     excel=NULL,
                                     keepers=test_keepers))
sva_batch_test <- sm(compare_de_results(combined_excel, combined_sva))
expected <- 0.96
actual <- sva_batch_test[["limma"]][["treated_vs_untreated"]][["logfc"]]
test_that("Do limma with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

expected <- 0.99
actual <- sva_batch_test[["deseq"]][["treated_vs_untreated"]][["logfc"]]
test_that("Do deseq with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

expected <- 0.99
actual <- sva_batch_test[["edger"]][["treated_vs_untreated"]][["logfc"]]
test_that("Do edger with combat and sva agree vis a vis logfc?", {
    expect_gt(actual, expected)
})

## See if the intersection between limma, deseq, and edger is decent.
test_intersect <- sm(intersect_significant(combined_sva, excel=NULL))
expected <- 82
actual <- nrow(test_intersect[["up_treated_vs_untreated"]][["led"]])
test_that("Do we get the expected number of agreed upon significant genes between edger/deseq/limma?", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["down_treated_vs_untreated"]][["led"]])
expected <- 88
test_that("Ibid, but in the down direction?", {
    expect_equal(actual, expected)
})
actual <- sum(nrow(test_intersect[["up_treated_vs_untreated"]][["l"]]) +
              nrow(test_intersect[["up_treated_vs_untreated"]][["e"]]) +
              nrow(test_intersect[["up_treated_vs_untreated"]][["d"]]))
expected <- 5
test_that("Are there very few genes observed without the others?", {
    expect_equal(actual, expected)
})
actual <- sum(nrow(test_intersect[["down_treated_vs_untreated"]][["l"]]) +
              nrow(test_intersect[["down_treated_vs_untreated"]][["e"]]) +
              nrow(test_intersect[["down_treated_vs_untreated"]][["d"]]))
expected <- 6
test_that("Ibid, but down?", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["up_treated_vs_untreated"]][["le"]])
expected <- 0
test_that("Do limma and edger have some genes in common? (up)", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["down_treated_vs_untreated"]][["le"]])
expected <- 0
test_that("Do limma and edger have some genes in common? (down)", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["up_treated_vs_untreated"]][["ld"]])
expected <- 0
test_that("Do limma and deseq have some genes in common? (up)", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["down_treated_vs_untreated"]][["ld"]])
expected <- 4
test_that("Do limma and deseq have some genes in common? (down)", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["up_treated_vs_untreated"]][["de"]])
expected <- 0
test_that("Do edger and deseq have some genes in common? (up)", {
    expect_equal(actual, expected)
})
actual <- nrow(test_intersect[["down_treated_vs_untreated"]][["de"]])
expected <- 0
test_that("Do edger and deseq have some genes in common? (down)", {
    expect_equal(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 29de_shared.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
