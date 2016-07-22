library(testthat)
library(hpgltools)

context("Do the combined differential expression searches work?")

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

normalized_expt <- s_p(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm"))$result
hpgl_result <- s_p(all_pairwise(normalized_expt, model_batch=TRUE))$result

previous_deseq <- deseq$hpgl_deseq$all_tables[["untreated_vs_treated"]]
previous_edger <- edger$hpgl_edger$all_tables[["untreated_vs_treated"]]
previous_limma <- limma$hpgl_limma$all_tables[["untreated_vs_treated"]]
previous_basic <- basic$hpgl_basic$all_tables[["untreated_vs_treated"]]

this_deseq <- hpgl_result$deseq$all_tables$untreated_vs_treated
this_edger <- hpgl_result$edger$all_tables$untreated_vs_treated
this_limma <- hpgl_result$limma$all_tables$untreated_vs_treated
this_basic <- hpgl_result$basic$all_tables$untreated_vs_treated

d <- e <- b <- l <- NULL
d <- expect_equal(previous_deseq, this_deseq)
e <- expect_equal(previous_edger, this_edger)
l <- expect_equal(previous_limma, this_limma)
b <- expect_equal(previous_basic, this_basic)

test_that("Do we get similar results to previous DE runs: (DESeq2)?", {
    expect_equal(previous_deseq, this_deseq)
})

test_that("Do we get similar results to previous DE runs: (edgeR)?", {
    expect_equal(previous_edger, this_edger)
})

test_that("Do we get similar results to previous DE runs: (limma)?", {
    expect_equal(previous_limma, this_limma)
})

test_that("Do we get similar results to previous DE runs: (basic)?", {
    expect_equal(previous_basic, this_basic)
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

combined_table <- s_p(combine_de_tables(hpgl_result, excel=NULL))$result
expected_size <- c(10153,42)
actual_size <- dim(combined_table$data[[1]])
test_that("Has the untreated/treated combined table been filled in?", {
    expect_equal(expected_size, actual_size)
})

sig_tables <- s_p(extract_significant_genes(combined_table, according_to="all", excel=NULL))$result
expected_limma_ups <- 134
actual_limma_ups <- nrow(sig_tables[["limma"]][["ups"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected_limma_ups, actual_limma_ups)
})

expected_limma_downs <- 126
actual_limma_downs <- nrow(sig_tables[["limma"]][["downs"]][[1]])
test_that("Are the limma significant downs expected?", {
    expect_equal(expected_limma_downs, actual_limma_downs)
})

expected_edger_ups <- 115
actual_edger_ups <- nrow(sig_tables[["edger"]][["ups"]][[1]])
test_that("Are the edger significant ups expected?", {
    expect_equal(expected_edger_ups, actual_edger_ups)
})

expected_edger_downs <- 152
actual_edger_downs <- nrow(sig_tables[["edger"]][["downs"]][[1]])
test_that("Are the limma significant ups expected?", {
    expect_equal(expected_edger_downs, actual_edger_downs)
})

expected_deseq_ups <- 66
actual_deseq_ups <- nrow(sig_tables[["deseq"]][["ups"]][[1]])
test_that("Are the deseq significant ups expected?", {
    expect_equal(expected_deseq_ups, actual_deseq_ups)
})

expected_deseq_downs <- 51
actual_deseq_downs <- nrow(sig_tables[["deseq"]][["downs"]][[1]])
test_that("Are the deseq significant downs expected?", {
    expect_equal(expected_deseq_downs, actual_deseq_downs)
})

expected_basic_ups <- 37
actual_basic_ups <- nrow(sig_tables[["basic"]][["ups"]][[1]])
test_that("Are the basic significant ups expected?", {
    expect_equal(expected_basic_ups, actual_basic_ups)
})

expected_basic_downs <- 22
actual_basic_downs <- nrow(sig_tables[["basic"]][["downs"]][[1]])
test_that("Are the basic significant downs expected?", {
    expect_equal(expected_basic_downs, actual_basic_downs)
})

funkytown <- plot_num_siggenes(combined_table$data[[1]])
expected_head_up_fc <- c(4.8609, 4.8118, 4.7627, 4.7136, 4.6645, 4.6154)
actual_head_up_fc <- as.numeric(head(funkytown$up_data[[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected_head_up_fc, actual_head_up_fc)
})

expected_head_down_fc <- c(-3.26997, -3.23694, -3.20391, -3.17088, -3.13785, -3.10482)
actual_head_down_fc <- as.numeric(head(funkytown$down_data[[1]]))
test_that("Can we monitor changing significance (up_fc)?", {
    expect_equal(expected_head_down_fc, actual_head_down_fc)
})
