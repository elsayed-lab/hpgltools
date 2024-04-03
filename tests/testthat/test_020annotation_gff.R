start <- as.POSIXlt(Sys.time())
context("020annotation_gff.R")
## 2017-12, exported functions in annotation_gff:
##   gff2irange(), load_gff_annotations(), pattern_count_genome()
##   sequence_attributes(), sum_exons()
## I deleted make_tooltips(), that was stupid.
## I moved get_gff_gene_lengths() to get_genelengths() and made it less stupid.

pa_gff <- system.file("share/paeruginosa_pa14.gff", package = "hpgldata")
pa_fasta <- system.file("share/paeruginosa_pa14.fasta", package = "hpgldata")

## gff2irange()
pa_irange <- gff2irange(pa_gff)
## 01
test_that("Do we get suitable irange data?", {
  expect_equal("GRanges object with 11946 ranges and 11 metadata columns",
               summary(pa_irange))
})

## load_gff_annotations()
pa_annot <- load_gff_annotations(pa_gff)
## 0203
test_that("Do we get some gff data for Pseudomonas?", {
  expect_equal(11946, nrow(pa_annot))
  expect_equal(16, ncol(pa_annot))
})

## pattern_count_genome()
pa_tas <- pattern_count_genome(pa_fasta, gff = pa_gff)
expected <- c(26, 16, 20, 39, 14, 14)
actual <- head(pa_tas[["number"]])
## 04
test_that("Do we get sensible numbers of TAs in the pseudomonas genome?", {
  expect_equal(expected, actual)
})

## sequence_attributes()
pa_attribs_genes <- sequence_attributes(pa_fasta, gff = pa_gff)
expected <- c(0.62589, 0.37411, 0.4757282, 0.5242718)
actual <- as.numeric(pa_attribs_genes["gene1650835", ])
## 05
test_that("Do we get sensible gene attributes by gene?", {
  expect_equal(expected, actual, tolerance = 0.001)
})

pa_attribs_genome <- sequence_attributes(pa_fasta)
expected <- c(0.6629220, 0.3370763, 0.4998674, 0.5001309)
actual <- as.numeric(pa_attribs_genome)
## 06
test_that("Do we get sensible gene attributes by genome?", {
  expect_equal(expected, actual, tolerance = 0.001)
})

## sum_exons()
## I need a gff to test this with

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 020annotation_gff.R in ", elapsed,  " seconds.")
