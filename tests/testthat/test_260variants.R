start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("260variants.R:
")

## The functions in variants.R deal with the outputs from samtools mpileup,
## bcftools, and vcfutils.  I have a small perl function which reads the bcf
## data and parses it into a table of high-confidence variants where each row is
## identified by position, reference, and new identity; columns are some of the
## summary information about the given variant.  This may then be distilled into
## an expressionset-like matrix of positions/identity and coverage.
## Note: the cyoa function defines high-confidence by default as greater than
## 80% agreement for 5 counts or more per position/variant.  Thus if there are 5
## hits and 4 agree, that is not sigificant, but 5 of 6 or 5 of 5 are.

## Therefore, it requires an expressionset to hold the metadata and provide the structure.
gff_file <- system.file("share/clbr/clbrener_8.1_complete_genes.gff.gz", package = "hpgltools")
gff_annotations <- load_gff_annotations(gff_file, type = "gene")
rownames(gff_annotations) <- gff_annotations[["ID"]]
meta <- system.file("share/clbr/clbr_samples_combined.xlsx", package = "hpgltools")
untarred <- utils::untar(tarfile = system.file("share/clbr/clbr_counts.tar.xz", package = "hpgltools"))
all_expt <- create_expt(metadata = meta, file_column = "clbrenerfile", gene_info = gff_annotations)

expected <- 25100
test_that("We have a functional expressionset?", {
  expect_equal(expected, nrow(exprs(all_expt)))
})

## Given that expressionset, note that there is a metadata column 'bcffile'
## This of course refers to the bcf data files, much like 'clbrenerfile' refers
## to the count tables for the clbrener haplotype.
untarred <- utils::untar(tarfile = system.file("share/clbr/vcfutils_output.tar.xz",
                                               package = "hpgltools"))
## Type in this context may be either percent or counts, this just defines the
## column to extract from the bcf file.
snp_expt <- count_expt_snps(all_expt, type = "percent", annot_column = "bcffile")
expected <- 8295
test_that("Do we have a decent number of variants?", {
  expect_equal(expected, nrow(exprs(snp_expt)))
})

## get_snp_sets uses Vennerable to cross reference variants against some column
## in the pData.  In this case, we are looking for the group of sets of
## unions/intersections with respect to condition.
snp_conditions <- pData(all_expt)[["condition"]]
expt_conditions <- pData(snp_expt)[["condition"]]
actual <- levels(snp_conditions)
expected <- c("CLBr.A96", "CLBr.A60", "CLBr.Tryp", "CLBr.Epi")
## Thus, get_snp_sets on condition will look for variants unique to the A96
## timepoint, A60 timepoint..., A96+A60, A96+Tryp..., the threes, and common to
## all; just like a Venn diagram!
test_that("Do we have sensible pData?", {
  expect_equal(snp_conditions, expt_conditions)
})

## For complext experimental designs, this can take quite a long time.
snp_sets <- get_snp_sets(snp_expt, factor = "condition")
actual <- snp_sets[["set_names"]][["0011"]]
expected <- "CLBr.Tryp, CLBr.Epi"
test_that("Do we get sensible set names from get_snp_sets?", {
  expect_equal(actual, expected)
})

## We can extract from this data structure lots of random information, like
## How many variants are shared among all conditions on chromosome 10P
expected <- 25
actual <- snp_sets$chr_data[["TcChr10-P"]]$venn
actual <- length(actual@IntersectionSets[["1111"]])
test_that("Do we get the expected number of variants on chromosome 10 shared among all conditions?", {
  expect_equal(actual, expected)
})

snp_gene_summary <- sm(snps_vs_genes(all_expt, snp_sets, expt_name_col = "seqnames"))
expected <- 18
actual <- snp_gene_summary[["summary_by_gene"]][["TcCLB.507505.10"]]
test_that("Do we observe the expected variants in a specific gene?", {
  expect_equal(actual, expected)
})

## Here we can ask for variants specific to samples with given condition(s)
## specific to a gene
snp_genes <- sm(snps_intersections(all_expt, snp_sets,
                                   chr_column = "seqnames"))
actual <- 11
## Thus, we expect 11 variant positions found only in the 3 Tryp samples
## in gene TcCLB.510483.360
expected <- snp_genes[["gene_summaries"]][["CLBr.Tryp"]][["TcCLB.510483.360"]]
test_that("Do we observe the expected variants in a specific gene under a specific condition?",
{
  expect_equal(actual, expected)
})

## The Epis are the pink samples, Tryps are purple, A60 are orange, A96 are green
##snp_norm <- normalize_expt(snp_expt, filter = TRUE, transform = "log2")
##test <- plot_corheat(snp_norm)
##test2 <- plot_sample_heatmap(snp_norm, Rowv = FALSE)

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 260variants.R in ", elapsed,  " seconds.")
