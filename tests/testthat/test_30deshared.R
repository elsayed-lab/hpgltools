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

normalized_expt <- normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm")
hpgl_result <- all_pairwise(normalized_expt, model_batch=TRUE)

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
test_that("Do we get similar results to previous DE runs?", {
    expect_equal(previous_deseq, this_deseq)
    expect_equal(previous_edger, this_edger)
    expect_equal(previous_limma, this_limma)
    expect_equal(previous_basic, this_basic)
})

le <- hpgl_result$comparison$comp[[1]]
ld <- hpgl_result$comparison$comp[[2]]
ed <- hpgl_result$comparison$comp[[3]]
lb <- hpgl_result$comparison$comp[[4]]
eb <- hpgl_result$comparison$comp[[5]]
db <- hpgl_result$comparison$comp[[6]]
test_that("Are the comparisons between DE tools sufficiently similar?", {
    expect_gt(le, 0.95)
    expect_gt(ld, 0.80)
    expect_gt(ed, 0.75)
    expect_gt(lb, 0.88)
    expect_gt(eb, 0.83)
    expect_gt(db, 0.68)
})

##combined_table <- combine_de_tables(hpgl_result, excel=NULL)
##sig_tables <- extract_significant_genes(combined_table, excel=NULL)
##sig_up_genes <- sig_tables$limma$ups
