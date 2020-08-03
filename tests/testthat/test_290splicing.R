start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("290alt_splicing.R:\n")

mtb_rmats_untarred <- utils::untar(tarfile=system.file("mtb_rmats.tar.xz", package="hpgltools"))

basedir <- "mtb_rmats/rmats_hg38_91"
## Define the prefixes
skipped_exon <- "SE"
alt_5p_splice_site <- "A5SS"
alt_3p_splice_site <- "A3SS"
mutually_exclusive <- "MXE"
retained_introns <- "RI"

## Define the suffixes
## Described in the documentation as 'with only reads that span splicing
## junctions'.
spanning_suffix <- ".MATS.JC.txt"
## Described in the documentation as 'with reads that span splicing junctions
## and reads on target (striped regions in the home page figure.)."
both_suffix <- ".MATS.JCEC.txt"

spanning_a3ss <- file.path(basedir, paste0(alt_3p_splice_site, spanning_suffix))
both_a3ss <- file.path(basedir, paste0(alt_3p_splice_site, both_suffix))
spanning_a5ss <- file.path(basedir, paste0(alt_5p_splice_site, spanning_suffix))
both_a5ss <- file.path(basedir, paste0(alt_5p_splice_site, both_suffix))
spanning_mxe <- file.path(basedir, paste0(mutually_exclusive, spanning_suffix))
both_mxe <- file.path(basedir, paste0(mutually_exclusive, both_suffix))
spanning_se <- file.path(basedir, paste0(skipped_exon, spanning_suffix))
both_se <- file.path(basedir, paste0(skipped_exon, both_suffix))
spanning_ri <- file.path(basedir, paste0(retained_introns, spanning_suffix))
both_ri <- file.path(basedir, paste0(retained_introns, both_suffix))

rmats_plots <- plot_rmats(se=both_se,
                          a5ss=both_a5ss,
                          a3ss=both_a3ss,
                          mxe=both_mxe,
                          ri=both_ri)

expected <- "gg"
actual <- class(rmats_plots[["ma"]])[1]
test_that("We get a plot of some rMATS data?", {
  expect_equal(expected, actual)
})

suppa_rmats_untarred <- utils::untar(tarfile=system.file("mtb_suppa.tar.xz", package="hpgltools"))
basedir <- "mtb_suppa/preprocessing/outputs/suppa_hg38_91"
dpsi_file <- file.path(basedir, "uninf_inf_diffsplice.dpsi")
tpm_file <- file.path(basedir, "uninf_inf_diffsplice_avglogtpm.tab")
events_file <- file.path(basedir, "ensembl_hg38.events.ioe")
psi_file <- file.path(basedir, "uninf_inf_diffsplice.psivec")

suppa_plots <- plot_suppa(dpsi_file, tpm_file,
                          events=events_file,
                          label_type="Skipping exon",
                          psi=psi_file)
actual <- class(suppa_plots[["ma"]])[1]
test_that("We get a plot of some suppa data?", {
  expect_equal(expected, actual)
})

annotation <- load_biomart_annotations()[["annotation"]]
stuff <- write_suppa_table(suppa_plots[["data"]], annotations=annotation)
actual <- nrow(stuff)
expected <- 72000
test_that("Can we write out suppa data with annotations?", {
  expect_gt(actual, expected)
})


end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 290alt_splicing.R in ", elapsed,  " seconds."))
