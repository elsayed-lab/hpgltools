start <- as.POSIXlt(Sys.time())
context("065de_all.R:\n")
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

## All of these functions will depend on an expt to play with:
pombe_expt <- make_pombe_expt(annotation=FALSE)
pombe_subset <- subset_expt(pombe_expt, subset="condition!='wt.180'&condition!='mut.180'")

## Well, in the previous test, we created pombe_expt, so let us use it.
## 01 basic_pairwise()
testing <- sm(basic_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 45
test_that("Basic performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt.0_vs_mut.0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 4
test_that("Basic got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["p"]]) < 0.1)
expected <- 363
test_that("Basic got some expected results (p)?", {
  expect_equal(expected, actual)
})

## 02 write_basic()
test <- sm(write_basic(testing, excel="test_basic_pairwise.xlsx"))
test_that("write_basic() did something?", {
  expect_true(file.exists("test_basic_pairwise.xlsx"))
})

## 03 deseq2_pairwise()
testing <- sm(deseq_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 45
test_that("DESeq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt.0_vs_mut.0"]]
if (is.null(test)) {
  test <- testing[["all_tables"]][["mut.0_vs_wt.0"]]
}
actual <- sum(test[["logFC"]] > 2)
expected <- 57
test_that("DESeq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 403
test_that("DESeq got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

## 04 write_deseq()
test <- sm(write_deseq(testing, excel="test_deseq_pairwise.xlsx"))
test_that("write_deseq() did something?", {
  expect_true(file.exists("test_deseq_pairwise.xlsx"))
})

## edger_pairwise()
testing <- sm(edger_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 45
test_that("edgeR performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt.0_vs_mut.0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 65
test_that("edgeR got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test$PValue) < 0.1)
expected <- 433
test_that("edgeR got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

## 05 write_edger()
test <- sm(write_edger(testing, excel="test_edger_pairwise.xlsx"))
test_that("write_edger() did something?", {
  expect_true(file.exists("test_edger_pairwise.xlsx"))
})

## hpgl_voomweighted()

## hpgl_voom()

## 06 limma_pairwise()
testing <- sm(limma_pairwise(pombe_subset))
actual <- length(testing[["contrasts_performed"]])
expected <- 45
test_that("limma performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt.0_vs_mut.0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 10
test_that("limma got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 462
test_that("limma got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

## 07 write_limma()
test <- sm(write_limma(testing, excel="test_limma_pairwise.xlsx"))
test_that("write_limma() did something?", {
  expect_true(file.exists("test_limma_pairwise.xlsx"))
})

## 08 all_pairwise()
## Setting parallel to FALSE so that I can make changes and test them
## immediately without having to reinstall.
test_condbatch <- sm(all_pairwise(pombe_subset))
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(min(test_condbatch[["comparison"]][["comp"]]), 0.75)
})

test_cond <- sm(all_pairwise(pombe_subset, model_batch=FALSE))
test_that("all_pairwise() provided results reasonably similar (no batch in model)?", {
  expect_gt(min(test_cond[["comparison"]][["comp"]]), 0.75)
})

tmp <- sm(normalize_expt(pombe_subset, filter=TRUE))
test_sva <- sm(all_pairwise(tmp, model_batch="svaseq"))
test_that("all_pairwise() provided results reasonably similar? (svaseq in model)", {
  expect_gt(min(test_sva[["comparison"]][["comp"]]), 0.65)
})

## 09 choose_model()
cond_model <- sm(choose_model(pombe_subset, model_batch=FALSE))
expected <- "~ 0 + condition"
actual <- cond_model[["chosen_string"]]
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## I think forcing the expt to keep conditions as specifically ordered factors
## has caused some oddities in how the downstream model is getting made.
## I think I should therefore ensure fully that the conditions of the model
## match the conditions in the original design.
model_df <- as.data.frame(cond_model[["chosen_model"]])
test_df <- data.frame(row.names=names(input[["conditions"]]))
for (cond in pombe_subset[["conditions"]]) {
  test_df[[cond]] <- 0
}
for (c in 1:length(pombe_subset[["conditions"]])) {
  name <- names(pombe_subset[["conditions"]])[c]
  value <- as.character(pombe_subset[["conditions"]][c])
  test_df[name, value] <- 1
}
test_that("choose_model provides a model which matches the design?", {
  expect_equal(model_df, test_df)
})

condbatch_model <- sm(choose_model(pombe_subset))
expected <- "~ 0 + condition + batch"
actual <- condbatch_model[["chosen_string"]]
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## 10 choose_dataset()
testing <- sm(choose_dataset(pombe_subset))
expected <- c("libsize", "conditions", "batches", "data")
actual <- names(testing)
test_that("choose_dataset provides some expected output?", {
  expect_equal(expected, actual)
})

## 11 combine_de_tables()
## 12 compare_de_results()
## we did test_condbatch, test_cond, test_sva
test_condbatch <- sm(combine_de_tables(test_condbatch))
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_condbatch[["data"]]), 45)
})
test_cond <- sm(combine_de_tables(test_cond))
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_cond[["data"]]), 45)
})

testing <- sm(compare_de_results(test_condbatch, test_cond))
expected <- 405
actual <- length(unlist(testing[["result"]]))
test_that("compare_de_results provides some expected output?", {
  expect_equal(expected, actual)
})

expected <- 0.98
actual <- min(testing[["logfc"]])
test_that("compare_de_results provides some expected logfc comparisons?", {
  expect_gt(actual, expected)
})

## 13 compare_led_tables()
testing <- sm(compare_led_tables(limma=test_sva[["limma"]],
                                 deseq=test_sva[["deseq"]],
                                 edger=test_sva[["edger"]],
                                 basic=test_sva[["basic"]]))
actual <- min(testing$comp)
expected <- 0.65
test_that("compare_led_tables provides some expected comparisons?", {
  expect_gt(actual, expected)
})

## 14 compare_logfc_plots()
testing <- sm(compare_logfc_plots(test_condbatch))
actual <- length(testing)
expected <- 45
test_that("Did compare_logfc_plots return something useful?", {
  expect_equal(actual, expected)
})

message("\nHuff huff, half way done.\n")

## 15 compare_significant_contrasts()
cb_sig <- sm(extract_significant_genes(combined=test_condbatch, excel=NULL))
testing <- sm(compare_significant_contrasts(cb_sig))
expected <- c(70, 34)
actual <- dim(testing[["shared_up"]])
test_that("Did compare_significant_contrasts provide some sensible result?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

test_that("Did compare_significant_contrasts provide some plots?", {
  expect_equal(class(testing[["up_venn_plot"]]), "recordedplot")
  expect_equal(class(testing[["down_venn_plot"]]), "recordedplot")
})

## do_pairwise()
## This is done by a bunch of other functions, I am not testing it.

## 16 get_abundant_genes()
testing <- get_abundant_genes(test_sva)
actual <- length(testing)
expected <- length(testing)
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

actual <- names(head(testing[[1]]))
expected <- c("SPAC212.09c", "SPAC212.04c", "SPAC977.11",
              "SPAC977.13c", "SPAC977.15", "SPAC977.16c")
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

## 17 get_pairwise_gene_abundances()
testing <- get_pairwise_gene_abundances(test_sva)
expected <- c(5801, 10)
actual <- dim(testing[["expression_values"]])
test_that("Did get_pairwise_gene_abundances() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## 18 get_sig_genes()
testing <- sm(get_sig_genes(table=test_sva$deseq$all_tables[[1]]))
expected <- c(209, 6)
actual <- dim(testing[["up_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})
expected <- c(334, 6)
actual <- dim(testing[["down_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## 19 make_pairwise_contrasts()
pombe_model <- sm(choose_model(pombe_subset))
testing <- sm(make_pairwise_contrasts(model=pombe_model[["chosen_model"]],
                                      conditions=pombe_subset$conditions))
actual <- length(names(testing[["all_pairwise"]]))
expected <- 45
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})
## If we add back some experimental factors, we should get bigger
## models/contrast lists.
pombe_model <- sm(choose_model(pombe_expt))
testing <- sm(make_pairwise_contrasts(model=pombe_model[["chosen_model"]],
                                      conditions=pombe_expt$conditions))
actual <- length(names(testing[["all_pairwise"]]))
expected <- 66
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})

## 20 semantic_copynumber_filter()
testing <- sm(semantic_copynumber_filter(de_list=cb_sig$limma,
                                         semantic="RNA",
                                         semantic_column="rownames"))
table <- "wt.120_vs_wt.0"
pre <- nrow(cb_sig[["limma"]][["ups"]][[table]])
post1 <- nrow(testing[["ups"]][[table]])
expect_lt(post1, pre)
testing <- sm(semantic_copynumber_filter(de_list=cb_sig$limma,
                                         invert=FALSE,
                                         semantic="RNA",
                                         semantic_column="rownames"))
post2 <- nrow(testing[["ups"]][[table]])
expect_lt(post2, pre)
expect_equal(pre, (post1 + post2))

## 21 significant_barplots()
testing <- sm(significant_barplots(combined=test_condbatch))
test_that("significant_barplots() gave some plots?", {
  expect_equal(class(testing[["deseq"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["limma"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["edger"]]), c("gg", "ggplot"))
})

## 22 extract_de_plots()
testing <- extract_de_plots(pairwise=test_sva)
test_that("extract_de_plots() gave some plots?", {
  expect_equal(class(testing[["ma"]][["plot"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["volcano"]][["plot"]]), c("gg", "ggplot"))
})

## 23 extract_coefficient_scatter()
testing <- sm(extract_coefficient_scatter(output=test_sva))
test_that("extract_de_plots() gave some plots?", {
  expect_equal(class(testing[["scatter"]]), c("gg", "ggplot"))
})

## 24 de_venn()
testing <- de_venn(test_condbatch[["data"]][[1]])
test_that("de_venn() gave some plots?", {
  expect_equal(class(testing[["up_noweight"]]), "recordedplot")
  expect_equal(class(testing[["down_noweight"]]), "recordedplot")
})

## 25 plot_num_siggenes()
testing <- plot_num_siggenes(test_condbatch[["data"]][[1]])
test_that("plot_num_siggenes() gave some plots?", {
  expect_equal(class(testing[["up"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["down"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["p"]]), c("gg", "ggplot"))
})

## 26 extract_abundant_genes()
testing <- sm(extract_abundant_genes(test_sva, excel=NULL))
test_that("extract_abundant_genes() gave some stuff?", {
  expect_equal(200, length(testing[["abundances"]][["limma"]][[1]]))
})

## 27 extract_significant_genes()
testing <- sm(extract_significant_genes(combined=test_condbatch, excel=NULL))
actual <- dim(testing$limma$ups[[1]])
expected <- c(80, 34)
test_that("Did extract_significant_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## 28 intersect_significant(),
testing <- sm(intersect_significant(combined=test_condbatch, excel=NULL))
test_that("Did intersect_significant() get some stuff?", {
  expect_equal(length(testing), 90)
})

## 29 write_de_table()
testing <- sm(write_de_table(data=test_sva, type="deseq"))
test_that("Did write_de_table() write something?", {
  expect_equal(testing, 1)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 065de_all.R in ", elapsed,  " seconds."))
