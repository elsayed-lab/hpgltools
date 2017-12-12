start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("20annotation_gff.R: Test functions in annotation_genbank.r")
## 2017-12, exported functions in annotation_gff:
##   gff2irange(), load_gff_annotations(), pattern_count_genome()
##   sequence_attributes(), sum_exons()
## I deleted make_tooltips(), that was stupid.
## I moved get_gff_gene_lengths() to get_genelengths() and made it less stupid.

pa_gff <- system.file("paeruginosa_pa14.gff", package="hpgltools")
pa_gff <- system.file("paeruginosa_pa14.fasta", package="hpgltools")

## gff2irange()
pa_irange <- sm(gff2irange(gff=pa_gff))
test_that("Do we get suitable irange data?", {
  expect_equal("GRanges object with 11946 ranges and 11 metadata columns",
               summary(pa_irange))
})

## load_gff_annotations()
pa_annot <- load_gff_annotations(gff=pa_gff)
test_that("Do we get some gff data for Pseudomonas?", {
  expect_equal(11946, nrow(pa_annot))
  expect_equal(16, ncol(pa_annot))
})

## pattern_count_genome()
pa_tas <- pattern_count_genome(fasta=pa_fasta, gff=pa_gff)
expected <- c(26, 16, 20, 39, 14, 14)
actual <- head(pa_tas[["number"]])
test_that("Do we get sensible numbers of TAs in the pseudomonas genome?", {
  expect_equal(expected, actual)
})

## sequence_attributes()
pa_attribs_genes <- sequence_attributes(fasta=pa_fasta, gff=pa_gff)
expected <- c(0.62589, 0.37411, 0.4757282, 0.5242718)
actual <- as.numeric(pa_attribs_genes["gene1650835", ])
test_that("Do we get sensible gene attributes by gene?", {
  expect_equal(expected, actual, tolerance=0.001)
})

pa_attribs_genome <- sequence_attributes(fasta=pa_fasta)
expected <- c(0.6629220, 0.3370763, 0.4998674, 0.5001309)
actual <- as.numeric(pa_attribs_genome)
test_that("Do we get sensible gene attributes by genome?", {
  expect_equal(expected, actual, tolerance=0.001)
})

## sum_exons()
## I need a gff to test this with

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 20annotation_gff.R in ", elapsed,  " seconds."))
