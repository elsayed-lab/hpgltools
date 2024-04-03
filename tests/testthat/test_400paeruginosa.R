start <- as.POSIXlt(Sys.time())
context("400paeruginosa.R: Test some Pseudomonas data.")

pa_gff <- load_gff_annotations(gff = system.file("share/paeruginosa_pa14.gff", package = "hpgldata"))
## This does provide some of what we need:
## seqnames, start, end, width, strand, source, type, score, phase, ID(thisone!),
## Name, Dbxref, Alias, name, Parent(or this one!), locus
good_parents <- pa_gff[["Parent"]] != ""
pa_gff <- pa_gff[good_parents, ]
pa_gff <- data.table::as.data.table(pa_gff)
pa_gff[["rownames"]] <- pa_gff[["Parent"]]
##pa_annotations <- as.data.frame(merge(x = pa_annotations, y = pa_gff, by.x = "sysName", by.y = "locus"))
##rownames(pa_annotations) <- pa_annotations[["Parent"]]

pa_expt <- create_expt(
  metadata = system.file("share/pa_samples.xlsx", package = "hpgldata"),
  countdir = system.file("share/counts", package = "hpgldata"),
  gene_info = pa_gff,
  title = "Pseudomonas aeruginosa RNAseq data of two strains and two time points.")

expected <- c(5979, 12)
actual <- dim(exprs(pa_expt))
## 01
test_that("Created pseudomonas expt of the correct count table size?", {
  expect_equal(expected, actual)
})

expected <- c(5979, 16)
actual <- dim(fData(pa_expt))
## 02
test_that("Created pseudomonas expt of the correct annotation table size?", {
  expect_equal(expected, actual)
})

expected <- c(12, 27)
actual <- dim(pData(pa_expt))
## 03
test_that("Created pseudomonas expt of the correct metadata size?", {
  expect_equal(expected, actual)
})

expected <- "Pseudomonas aeruginosa RNAseq data of two strains and two time points."
actual <- pa_expt[["title"]]
## 04
test_that("The expt has a title.", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 01paeruginosa.R in ", elapsed,  " seconds.")
