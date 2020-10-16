start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("280tnseq.R:
")

## I am going to just copy some of my tasks from here:
## https://github.com/abelew/sagalacticae_2019
## I should therefore just be able to repeat them in my tests
## to ensure that my tnseq functions still work.

a909_microbes <- load_microbesonline_annotations(species="A909")
expected <- 2200
actual <- nrow(a909_microbes)
test_that("We downloaded annotations for strain a909?", {
  expect_gt(actual, expected)
})

## I copied some of the reference data into inst/gbs_tnseq
gff_file <- system.file("gbs_tnseq/sagalactiae_a909.gff", package="hpgltools")
a909_gff <- load_gff_annotations(gff_file)
expected <- 4400
actual <- nrow(a909_gff)
test_that("We acquired annotations from the gff file used to map the data?", {
  expect_gt(actual, expected)
})

## Combine the two annotation sources
a909_microbes <- as.data.frame(a909_microbes)
rownames(a909_gff) <- make.names(a909_gff[["locus_tag"]], unique=TRUE)
## I am going to only pay attention to the first annotation for each locus tag from microbesonline.
a909_microbes[["sysName"]] <- make.names(a909_microbes[["sysName"]], unique=TRUE)
a909_annot <- merge(a909_gff, a909_microbes, by.x="old_locus_tag", by.y="sysName")
rownames(a909_annot) <- make.names(a909_annot[["locus_tag"]], unique=TRUE)
## Rename the merged start/strand columns
colnames(a909_annot)[3] <- "start"
colnames(a909_annot)[6] <- "strand"
## And drop the duplicate columns
a909_annot[["start.y"]] <- NULL
a909_annot[["strand.y"]] <- NULL

expected <- 4000
actual <- nrow(a909_annot)
test_that("We merged the annotation data?", {
  expect_gt(actual, expected)
})

a909_counts <- utils::untar(tarfile=system.file("gbs_tnseq/gbs_essentiality_counts.tar.xz",
                                                package="hpgltools"))
metadata <- system.file("gbs_tnseq/sagalactiae_samples.xlsx", package="hpgltools")
a909_expt <- create_expt(metadata=metadata, batch=FALSE, gene_info=a909_annot,
                         file_column="a909_filename")
expected <- 2000
actual <- nrow(exprs(a909_expt))
test_that("We created an expressionset?", {
  expect_gt(actual, expected)
})

## Grab copies of some of the essentiality results
a909_wig <- utils::untar(tarfile=system.file("gbs_tnseq/gbs_essentiality_wig.tar.xz",
                                             package="hpgltools"))
a909_csv <- utils::untar(tarfile=system.file("gbs_tnseq/gbs_essentiality.tar.xz",
                                             package="hpgltools"))

saturation <- tnseq_saturation(
  "preprocessing/01/outputs/essentiality_sagalactiae_a909/trimmed_ca-v0M1.wig",
  adjust=2)
test_that("tnseq_saturation returns expected outputs?", {
  expect_equal("gg", class(saturation[["plot"]])[1])
})
expected <- 21000
test_that("We expect more than 21000 TAs with more than 16 hits:", {
  expect_gt(saturation[["gt_16"]], expected)
})
expected <- 67
test_that("We expect an average of 67ish hits per TA:", {
  expect_gt(saturation[["hits_summary"]][["Mean"]], expected)
})

ess_plts <- plot_essentiality(
  "preprocessing/01/outputs/essentiality_sagalactiae_a909/mh_ess-trimmed_ca-v0M1_gene_tas_m1.csv")
test_that("plot_essentiality returns expected outputs?", {
  expect_equal("gg", class(ess_plts[["zbar"]])[1])
  expect_equal("gg", class(ess_plts[["span_plot"]])[1])
})


plt <- sm(tnseq_multi_saturation(meta=pData(a909_expt), meta_column="a909esswig"))
test_that("tnseq_multi_saturation returns some fun?", {
  expect_equal("gg", class(plt[["plot"]])[1])
  expect_equal("gg", class(plt[["ggstats"]])[1])
})

## Perform my 'fitness' analysis; which is just a normal differential expression analysis
a909_de <- all_pairwise(a909_expt, model_batch=FALSE)
test_that("all_pairwise returned?", {
  expect_equal("all_pairwise", class(a909_de)[1])
})

## Make a couple tables out of that:
a909_contrasts <- list(
  "low_vs_control" = c("cal_low", "control"),
  "high_vs_control" = c("cal_high", "control"))
a909_tables <- combine_de_tables(
  a909_de, keepers=a909_contrasts,
  excel="a909_tables.xlsx")
test_that("all_pairwise returned?", {
  expect_equal("combined_de", class(a909_tables)[1])
})

a909_sig <- extract_significant_genes(
    a909_tables, excel="a909_sig.xlsx")
expected <- 29
actual <- a909_sig[["summary_df"]]["low_vs_control", "edger_change_counts_up"]
test_that("Did we get the expected number of up genes between low Ca+ and control according to EdgeR?", {
  expect_equal(expected, actual)
})

colors <- c("990000", "008800", "000000", "0000AA")
names(colors) <- c("E", "NE", "S", "U")
low_df <- a909_tables[["data"]][["low_vs_control"]]
high_df <- a909_tables[["data"]][["high_vs_control"]]

circos_cfg <- circos_prefix(annotation=a909_annot, name="a909")
a909_fasta <- system.file("gbs_tnseq/sagalactiae_a909.fasta", package="hpgltools")
a909_kary <- circos_karyotype(circos_cfg, fasta=a909_fasta)
a909_plus_minus <- circos_plus_minus(circos_cfg, width=0.06, thickness=40)
a909_low <- circos_hist(circos_cfg, low_df, colname="deseq_logfc", basename="low",
                        outer=a909_plus_minus, fill_color="vvdpgreen", width=0.06, thickness=0.1)
a909_suffix <- circos_suffix(circos_cfg)
made <- sm(circos_make(circos_cfg, target="a909"))

test_that("circos provided an imagemap output?", {
  expect_true(file.exists("circos/a909.html"))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 280tnseq.R in ", elapsed,  " seconds."))
