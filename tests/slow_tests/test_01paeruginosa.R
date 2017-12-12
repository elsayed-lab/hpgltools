start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("01paeruginosa.R: Test some Pseudomonas data.\n")

pa_ids <- sm(get_microbesonline_ids("PA14"))
pa_id <- pa_ids[1, 1]
pa_annotations <- sm(load_microbesonline_annotations(ids=pa_id)[[1]])
pa_annotations <- data.table::as.data.table(pa_annotations)
## This provides the following columns:
## locusId, accession, GI, scaffoldId, start, stop, strand, sysName
## name, desc, COG, COGFun, COGDesc, TIGRFam, TIGRRoles, GO, EC, ECDesc
## Sadly, none of them match the gene16... rownames in the count tables.

pa_gff <- sm(load_gff_annotations(gff=system.file("paeruginosa_pa14.gff", package="hpgltools")))
## This does provide some of what we need:
## seqnames, start, end, width, strand, source, type, score, phase, ID(thisone!),
## Name, Dbxref, Alias, name, Parent(or this one!), locus
good_parents <- pa_gff[["Parent"]] != ""
pa_gff <- pa_gff[good_parents, ]
pa_gff <- data.table::as.data.table(pa_gff)
pa_annotations <- as.data.frame(merge(x=pa_annotations, y=pa_gff, by.x="sysName", by.y="locus"))
rownames(pa_annotations) <- pa_annotations[["Parent"]]

pa_expt <- sm(create_expt(metadata=system.file("pa_samples.xlsx", package="hpgltools"),
                          countdir=system.file("counts", package="hpgltools"),
                          gene_info=pa_annotations,
                          title="Pseudomonas aeruginosa RNAseq data of two strains and two time points."))

expected <- c(5979, 12)
actual <- dim(exprs(pa_expt))
test_that("Created pseudomonas expt of the correct count table size?", {
  expect_equal(expected, actual)
})

expected <- c(5979, 33)
actual <- dim(fData(pa_expt))
test_that("Created pseudomonas expt of the correct annotation table size?", {
  expect_equal(expected, actual)
})

expected <- c(12, 27)
actual <- dim(pData(pa_expt))
test_that("Created pseudomonas expt of the correct metadata size?", {
  expect_equal(expected, actual)
})

expected <- "Pseudomonas aeruginosa RNAseq data of two strains and two time points."
actual <- pa_expt[["title"]]
test_that("The expt has a title.", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 01paeruginosa.R in ", elapsed,  " seconds."))
