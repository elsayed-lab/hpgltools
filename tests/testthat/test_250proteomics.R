start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("250proteomics.R:
")

## Available functions:  add_conditional_nas(), extract_mayu_pps_fdr(),
## extract_scan_data(), extract_mzXML_scans(), extract_mzML_scans(),
## extract_msraw_data(), extract_peprophet_data(), extract_pyprophet_data(),
## impute_expt(), mean_by_bioreplicate(), read_thermo_xlsx(), s2s_all_filters(),
## subset_pyprophet_data(), and gather_masses().

## Plotting functions:
## plot_intensity_mz(), plot_mzxml_boxplot(), plot_pyprophet_counts(),
## plot_pyprophet_xy(), plot_pyprophet_distribution(), plot_pyprophet_protein(),
## plot_pyprophet_points(), plot_peprophet_data(), plot_cleaved(),
## cleavage_histogram()


## I would like to test out my various proteomics functions.  Unfortunately,
## many of those functions require rather large input files.  So, I made a
## tarfile of the smallest of them in order to play with some of the
## functionality.

## As a reminder to myself, the general process followed is: (check out my
## 02_preprocessing*.Rmd for more details))
##  1.  Convert the thermo-fisher proprietary output files to mzXML/mzML using
##      MSreader on a windows virtual machine.
##  2.  Create a spectral library in the pqp format.  This may either be created
##      from an existing DDA experiment, or from some downloaded data.  Let us
##      assume the latter for the moment, then the command is
##      'TargetedFileConverter'
##  3.  Check over the mzXML/mzML files to ensure that there are no problems so
##      significant that they will cause later steps to fail. (That is what the
##      first few functions help with)
##  4.  Write appropriate window files for OpenSwathWorkflow, the function
##      'extract_msraw_data()' will do that for you, but it requires the rather
##      large mzXML files which I don't want to include here.  So unless I find
##      a data package containing some, I guess you will have to trust me that
##      those parsers actually work.
##  5.  Pass the transition library and raw data files to OpenSwathWorkFlow.
##  6.  Run pyprophet on the openswath outputs to normalize the error rates.
##  7.  Examine the pyprophet outputs to see if there are problems.
##      extract_pyprophet_data() and related functions help with that.
##  8.  Invoke SWATH2stats to filter the data and create matrices in the formats
##      expected by MSstats etc.
##  9.  Conversely, take those matrices and pass them to hpgltools.
##  10. Conversely, pass the mzXML data directly to encylopeDIA and extract the
##      output matrices.

## For most of the above choices, I have functions in the hpgltools to help.

meta <- system.file("share/mtb_prot/dia_samples.ods", package = "hpgltools")
untarred <- utils::untar(tarfile = system.file("share/mtb_prot/sb_prot.tar.xz",
                                               package = "hpgltools"))
##mtb_expt <- create_expt(meta = meta)

## As the name implies, this function uses the diascored column in the metadata
## and reads the pyprophet output files indicated therein; it parses them into a
## list of data tables which are later used for plotting/examining.
## In that list is 'failed, colors, metadata, and sample_data'.  Failed is
## comprised of the files which failed to read properly, colors are for
## plotting, metadata is a copy of the metadata, and sample_data is the fun.
pyprophet_fun <- extract_pyprophet_data(metadata = meta,
                                        pyprophet_column = "diascored")
test_that("Did extract_pyprotphet_data have failures?", {
  expect_equal(NULL, pyprophet_fun[["failed"]])
})
expected <- 35
test_that("Did extract_pyprophet_data provide the expected number of elements?", {
  expect_equal(expected, length(pyprophet_fun[["sample_data"]]))
})
## These are the columns provided by pyprophet.  They are columns provided by
## openswath with some extra scores appended by pyprophet.
expected <- c(
    "transition_group_id", "decoy", "run_id", "filename", "rt",  "assay_rt",
    "delta_rt", "irt", "assay_irt", "delta_irt", "id", "sequence",
    "fullpeptidename", "charge", "mz", "intensity", "aggr_prec_peak_area",
    "aggr_prec_peak_apex", "leftwidth", "rightwidth", "peak_group_rank",
    "d_score", "m_score", "aggr_peak_area", "aggr_peak_apex",
    "aggr_fragment_annotation", "proteinname", "mass", "seqlength")
test_that("Did extract_pyprotphet_data provide the expected columns?", {
  expect_equal(expected, colnames(pyprophet_fun[["sample_data"]][[1]]))
})
## The various plot_pyprophet functions may be used on pretty much any of the
## above columns, some of them are more useful than others for diagnosing
## problems.

mass_plot <- plot_pyprophet_distribution(pyprophet_fun, column = "mass")
test_that("Does plot_pyprotphet_distribution return some plots?", {
  ## Yeah, so I couldn't decide which type of plot was best for representing this...
  expect_equal(class(mass_plot[["violin"]])[1], "gg")
  expect_equal(class(mass_plot[["boxplot"]])[1], "gg")
  expect_equal(class(mass_plot[["dotboxplot"]])[1], "gg")
  expect_equal(class(mass_plot[["density"]])[1], "gg")
})

## plot_pyprophet_counts has some parameters to ignore/include decoys
## as well as different types: 'protein_count' which counts how many proteins
## were found, 'intensity' which sums the intensities observed, 'count' which
## counts how many identifications were observed, otherwise it sums whatever
## columnname was provided, assuming it was numeric.
peptide_identifications <- plot_pyprophet_counts(pyprophet_fun, keep_decoys = FALSE,
                                                 type = "count")
test_that("Does plot_pyprophet_counts return some plots?", {
  expect_equal(class(peptide_identifications[["plot"]])[1], "gg")
})

## widths with respect to counts are a surprisingly reliable way to find
## problematic samples.
pyprophet_lwidths <- plot_pyprophet_xy(pyprophet_fun, x_type = "count", y_type = "leftwidth")
test_that("Does plot_pyprophet_xy return a plot?", {
  expect_equal(class(pyprophet_lwidths)[1], "gg")
})

## The above are global metrics, we can plot some metrics of individual
## proteins.  Any column from above may be used.
intensities_esxG <- plot_pyprophet_protein(pyprophet_fun, scale = "log",
                                           title = "esxG Intensities",
                                           column = "intensity", protein = "Rv0287")
test_that("Does plot_pyprophet_protein return a plot?", {
  expect_equal(class(intensities_esxG)[1], "gg")
})

## The final step of pyprophet is to export the individual matrices into a
## single table, that is the input for SWATH2stats.
## I have one function which helps with that, s2s_all_filters() which just
## invokes all of the SWATH2stats filters with some hopefully sane defaults.

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 250proteomics.R in ", elapsed,  " seconds.")
