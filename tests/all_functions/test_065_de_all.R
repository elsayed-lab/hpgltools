start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("065de_all.R:
  123456789012345678901234567890123456789012345678901234567890123\n")
## 2017-12, exported functions in de_basic:
## basic_pairwise(), write_basic()
## 2017-12, exported functions in de_deseq:
## deseq2_pairwise(), write_deseq()
## 2017-12, exported functions in de_edger:
## edger_pairwise(), write_edger()
## 2017-12, exported functions in de_limma:
## hpgl_voomweighted(), hpgl_voom(), limma_pairwise(), write_limma()
## 2017-12, exported functions in de_shared:
## all_pairwise(), choose_model(), choose_dataset(), compare_de_results(),
## compare_led_tables(), compare_logfc_plots(), compare_significant_contrasts(),
## do_pairwise(), get_abundant_genes(), get_pairwise_gene_abundances(),
## get_sig_genes(), make_pairwise_contrasts(), semantic_copynumber_filter(),
## semantic_copynumber_extract(),
## 2017-12, exported functions in de_plots:
## significant_barplots(), extract_de_plots(),
## extract_coefficient_scatter(), de_venn(), plot_num_siggenes()
## 2017-12, exported functions in de_xlsx:
## combine_de_tables(), extract_abundant_genes()
## extract_significant_genes(), intersect_significant(),
## write_de_table()

if (! "Vennerable" %in% installed.packages()) {
    devtools::install_github("js229/Vennerable")
}

## All of these functions will depend on an expt to play with:
pombe_expt <- make_pombe_expt(annotation=FALSE)
pombe_subset <- subset_expt(
  pombe_expt,
  subset="minute == 0 | minute == 15 | minute == 30")

## Well, in the previous test, we created pombe_expt, so let us use it.
testing <- sm(basic_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 15
## 01
test_that("Basic performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 1
## 02
test_that("Basic got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Add a little fudge-factor to some of these tests.
actual <- sum(as.numeric(test[["p"]]) < 0.1)
expected <- 358
## 03
test_that("Basic got some expected results (p)?", {
  expect_equal(expected, actual, tolerance=3)
})

test <- write_basic(testing, excel="test_basic_pairwise.xlsx")
## 04
test_that("write_basic() did something?", {
  expect_true(file.exists("test_basic_pairwise.xlsx"))
})

testing <- sm(deseq_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 15
## 05
test_that("DESeq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 50
## 06
test_that("DESeq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 319
## 07
test_that("DESeq got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

written_test <- write_deseq(testing, excel="test_deseq_pairwise.xlsx")
## 08
test_that("write_deseq() did something?", {
  expect_true(file.exists("test_deseq_pairwise.xlsx"))
})

## edger_pairwise()
testing <- sm(edger_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 15
## 09
test_that("edgeR performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 61
## 10
test_that("edgeR got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["PValue"]]) < 0.1)
expected <- 328
## 11
test_that("edgeR got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

test <- write_edger(testing, excel="test_edger_pairwise.xlsx")
## 12
test_that("write_edger() did something?", {
  expect_true(file.exists("test_edger_pairwise.xlsx"))
})

## hpgl_voomweighted()
## hpgl_voom()
testing <- limma_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
## 13
test_that("limma performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 9
## 14
test_that("limma got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 421
## 15
test_that("limma got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

test <- write_limma(testing, excel="test_limma_pairwise.xlsx")
## 16
test_that("write_limma() did something?", {
  expect_true(file.exists("test_limma_pairwise.xlsx"))
})

test_condbatch <- all_pairwise(pombe_subset, parallel=FALSE)
actual <- min(test_condbatch[["comparison"]][["comp"]])
expected <- 0.71
## 17
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

test_cond <- all_pairwise(pombe_subset, model_batch=FALSE)
actual <- min(test_cond[["comparison"]][["comp"]])
expected <- 0.76
## 18
test_that("all_pairwise() provided results reasonably similar (no batch in model)?", {
  expect_gt(actual, expected)
})

tmp <- normalize_expt(pombe_subset, filter=TRUE)
test_sva <- all_pairwise(tmp, model_batch="svaseq")
actual <- min(test_sva[["comparison"]][["comp"]])
expected <- 0.63
## 19
test_that("all_pairwise() provided results reasonably similar? (svaseq in model)", {
  expect_gt(actual, expected)
})

cond_model <- choose_model(pombe_subset, model_batch=FALSE)
expected <- "~ 0 + condition"
actual <- cond_model[["chosen_string"]]
## 20
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## I think forcing the expt to keep conditions as specifically ordered factors
## has caused some oddities in how the downstream model is getting made.
## I think I should therefore ensure fully that the conditions of the model
## match the conditions in the original design.
model_df <- as.data.frame(cond_model[["chosen_model"]])
test_df <- data.frame(row.names=names(pombe_subset[["conditions"]]))
for (cond in pombe_subset[["conditions"]]) {
  test_df[[cond]] <- 0
}
for (c in 1:length(pombe_subset[["conditions"]])) {
  name <- names(pombe_subset[["conditions"]])[c]
  value <- as.character(pombe_subset[["conditions"]][c])
  test_df[name, value] <- 1
}
## 21
test_that("choose_model provides a model which matches the design?", {
  expect_equal(model_df, test_df)
})

condbatch_model <- choose_model(pombe_subset)
expected <- "~ 0 + condition + batch"
actual <- condbatch_model[["chosen_string"]]
## 22
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## choose_dataset()
testing <- choose_dataset(pombe_subset)
expected <- c("libsize", "conditions", "batches", "data")
actual <- names(testing)
## 23
test_that("choose_dataset provides some expected output?", {
  expect_equal(expected, actual)
})

## we did test_condbatch, test_cond, test_sva
test_condbatch_combined <- combine_de_tables(test_condbatch)
## 24
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_condbatch_combined[["data"]]), 15)
})

## Same query, condition in model
test_cond_combined <- combine_de_tables(test_cond)
## 25
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_cond_combined[["data"]]), 15)
})

testing <- compare_de_results(test_condbatch_combined, test_cond_combined)
expected <- 135
actual <- length(unlist(testing[["result"]]))
## 26
test_that("compare_de_results provides some expected output?", {
  expect_equal(expected, actual)
})

expected <- 0.96
actual <- min(testing[["logfc"]])
## 27
test_that("compare_de_results provides some expected logfc comparisons?", {
  expect_gt(actual, expected)
})

testing <- correlate_de_tables(test_sva)
actual <- min(testing$comp)
expected <- 0.63
## 28
test_that("compare_led_tables provides some expected comparisons?", {
  expect_gt(actual, expected)
})

message("\nHuff huff, half way done.\n")
## Strange, I got a failure here when running make test
## but running manually everything seems to be working fine...
## 15 compare_significant_contrasts()
cb_sig <- extract_significant_genes(combined=test_condbatch_combined, excel=NULL)
actual <- dim(cb_sig$limma$ups[[1]])
expected <- c(182, 34)
## 2930
test_that("Did extract_significant_genes provide some sensible result?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

testing <- compare_significant_contrasts(
  cb_sig,
  contrasts=c("wt15_vs_mut0", "wt30_vs_mut0", "wt30_vs_mut15"))
shared_all <- testing[["up_intersections"]][["all"]]
actual <- length(shared_all)
expected <- 106
## 31
test_that("Did compare_significant_contrasts provide some sensible result?", {
  expect_equal(actual, expected)
})

## 3233
test_that("Did compare_significant_contrasts provide some plots?", {
  expect_equal(class(testing[["up_plot"]]), "recordedplot")
  expect_equal(class(testing[["down_plot"]]), "recordedplot")
})

## Saving this so we can use it for ontology searches later.
save(list="test_condbatch_combined", file="test_065_combined.rda", compress=TRUE)
save(list="cb_sig", file="test_065_significant.rda", compress=TRUE)


## do_pairwise()
## This is done by a bunch of other functions, I am not testing it.

testing <- get_abundant_genes(test_sva)
actual <- length(testing)
expected <- 6
## 34
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

actual <- names(head(testing[[1]]))
expected <- c("SPAC212.09c", "SPAC212.04c", "SPAC977.11",
              "SPAC977.13c", "SPAC977.15", "SPAC977.16c")
## 35
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

testing <- get_pairwise_gene_abundances(test_sva)
expected <- c(5720, 6)
actual <- dim(testing[["expression_values"]])
## 3637
test_that("Did get_pairwise_gene_abundances() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

testing <- get_sig_genes(table=test_sva$deseq$all_tables[[1]])
expected <- c(199, 6)
actual <- dim(testing[["up_genes"]])
## 3839
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})
expected <- c(183, 6)
actual <- dim(testing[["down_genes"]])
## 4041
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

pombe_model <- choose_model(pombe_subset)
testing <- make_pairwise_contrasts(model=pombe_model[["chosen_model"]],
                                   conditions=pombe_subset$conditions)
actual <- length(names(testing[["all_pairwise"]]))
expected <- 15
## 42
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})
## If we add back some experimental factors, we should get bigger
## models/contrast lists.
pombe_model <- choose_model(pombe_expt)
testing <- make_pairwise_contrasts(model=pombe_model[["chosen_model"]],
                                   conditions=pombe_expt$conditions)
actual <- length(names(testing[["all_pairwise"]]))
expected <- 66
## 43
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})

testing <- semantic_copynumber_filter(cb_sig[["limma"]],
                                      semantic="RNA",
                                      semantic_column="rownames")
table <- "wt30_vs_wt0"
pre <- nrow(cb_sig[["limma"]][["ups"]][[table]])
post1 <- nrow(testing[["ups"]][[table]])
expect_lt(post1, pre)
testing <- semantic_copynumber_filter(cb_sig[["limma"]],
                                      invert=FALSE,
                                      semantic="RNA",
                                      semantic_column="rownames")
post2 <- nrow(testing[["ups"]][[table]])
## 4445
test_that("Do we get expected results from semantic_copynumber_filter?", {
  expect_lt(post2, pre)
  expect_equal(pre, (post1 + post2))
})

testing <- significant_barplots(combined=test_condbatch_combined)
## 464748
test_that("significant_barplots() gave some plots?", {
  expect_equal(class(testing[["deseq"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["limma"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["edger"]]), c("gg", "ggplot"))
})

testing <- extract_de_plots(pairwise=test_sva)
## 4950
test_that("extract_de_plots() gave some plots?", {
  expect_equal(class(testing[["ma"]][["plot"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["volcano"]][["plot"]]), c("gg", "ggplot"))
})

testing <- extract_coefficient_scatter(output=test_sva)
## 51
test_that("extract_de_plots() gave some plots?", {
  expect_equal(class(testing[["scatter"]]), c("gg", "ggplot"))
})

testing <- de_venn(test_condbatch_combined[["data"]][[1]])
## 5253
test_that("de_venn() gave some plots?", {
  expect_equal(class(testing[["up_noweight"]]), "recordedplot")
  expect_equal(class(testing[["down_noweight"]]), "recordedplot")
})

testing <- plot_num_siggenes(test_condbatch_combined[["data"]][[1]])
expected <- c("gg", "ggplot")
## 57
test_that("plot_num_siggenes() gave some plots?", {
  expect_equal(class(testing[["up"]]), expected)
  expect_equal(class(testing[["down"]]), expected)
  expect_equal(class(testing[["pup"]]), expected)
  expect_equal(class(testing[["pdown"]]), expected)
})

testing <- extract_abundant_genes(test_sva, excel=NULL)
## 58
test_that("extract_abundant_genes() gave some stuff?", {
  expect_equal(200, length(testing[["abundances"]][["limma"]][[1]]))
})

testing <- extract_significant_genes(combined=test_condbatch_combined, excel=NULL)
actual <- dim(testing$limma$ups[[1]])
expected <- c(182, 34)
## 60
test_that("Did extract_significant_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

testing <- intersect_significant(combined=test_condbatch_combined, excel=NULL)
expected <- 384
## 61
test_that("Did intersect_significant() get some stuff?", {
  expect_equal(testing[["summary"]]["up", "all"], expected)
})

testing <- write_de_table(data=test_sva, type="deseq")
## 62
test_that("Did write_de_table() write something?", {
  expect_equal(testing, 1)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 065de_all.R in ", elapsed,  " seconds."))
