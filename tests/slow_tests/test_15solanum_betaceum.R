start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("15solanum_betaceum.R: Test some Pseudomonas data.\n")

sb_annot <- load_trinotate_annotations(
  trinotate=system.file("sb/trinotate_head.csv.xz", package="hpgltools"))

expected <- c(136152, 33)
actual <- dim(sb_annot)
test_that("Do we get the expected trinotate annotation data?", {
  expect_equal(expected, actual)
})

## Try loading some data generated with Sandra and make sure it still works.
untarred <- utils::untar(tarfile=system.file("sb/preprocessing.tar.xz", package="hpgltools"),
                         compressed="xz")

## A little context: Sandra is working with an organism with no genome and no love.
## We therefore passed all the reads to trinity and trinotate and used the resulting,
## entirely bad, putative transcriptome as an input for kallisto and rsem.
## Given that as input, we think we can hunt down some meaningful genes.
## Let us see if that is true.

sb_expt <- sm(create_expt(metadata="preprocessing/kept_samples.xlsx", gene_info=sb_annot))
sb_count <- sm(normalize_expt(sb_expt, filter=TRUE, thresh=4))
sb_plots <- sm(graph_metrics(sb_count))
sb_norm <- sm(normalize_expt(sb_count, norm="quant", convert="cpm", transform="log2"))
sb_nplots <- sm(graph_metrics(sb_norm))
