## ----options, include = FALSE-------------------------------------------------
## These are the options I tend to favor
library("hpgltools")
## tt <- devtools::load_all("~/hpgltools")
knitr::opts_knit$set(progress = TRUE,
                     verbose = TRUE,
                     width = 90,
                     echo = TRUE)
knitr::opts_chunk$set(error = TRUE,
                      fig.width = 8,
                      fig.height = 8,
                      dpi = 96)
old_options <- options(digits = 4,
                       stringsAsFactors = FALSE,
                       knitr.duplicate.label = "allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))
set.seed(1)
rmd_file <- "c-03_fission_differential_expression.Rmd"

## ----setup--------------------------------------------------------------------
library(hpgltools)
tt <- sm(library(fission))
tt <- data(fission)

## ----spombe_annotations-------------------------------------------------------
pombe_annotations <- load_biomart_annotations(
    host = "fungi.ensembl.org",
    trymart = "fungal_mart",
    trydataset = "spombe_eg_gene",
    gene_requests = c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                      "hgnc_symbol", "description", "gene_biotype"),
    species = "spombe", overwrite = TRUE)
pombe_mart <- pombe_annotations[["mart"]]
annotations <- pombe_annotations[["annotation"]]
rownames(annotations) <- make.names(gsub(pattern = "\\.\\d+$",
                                         replacement = "",
                                         x = rownames(annotations)), unique = TRUE)

## ----data_import--------------------------------------------------------------
## Extract the meta data from the fission dataset
meta <- as.data.frame(fission@colData)
## Make conditions and batches
meta[["condition"]] <- paste(meta$strain, meta$minute, sep = ".")
meta[["batch"]] <- meta[["replicate"]]
meta[["sample.id"]] <- rownames(meta)
## Grab the count data
fission_data <- fission@assays[["data"]][["counts"]]
## This will make an experiment superclass called 'expt' and it contains
## an ExpressionSet along with any arbitrary additional information one might want to include.
## Along the way it writes a Rdata file which is by default called 'expt.Rdata'
fission_expt <- create_expt(metadata = meta,
                            count_dataframe = fission_data,
                            gene_info = annotations)

## ----simple_subset------------------------------------------------------------
fun_data <- subset_expt(fission_expt,
                        subset = "condition=='wt.120'|condition=='wt.30'")
fun_filt <- normalize_expt(fun_data, filter = "simple")
fun_norm <- sm(normalize_expt(fun_filt, batch = "limma", norm = "quant",
                              transform = "log2", convert = "cpm"))

## ----simple_limma-------------------------------------------------------------
limma_comparison <- sm(limma_pairwise(fun_data))
names(limma_comparison[["all_tables"]])
summary(limma_comparison[["all_tables"]][["wt30_vs_wt120"]])
scatter_wt_mut <- extract_coefficient_scatter(limma_comparison, type = "limma",
                                              x = "wt30", y = "wt120")
scatter_wt_mut[["scatter"]]
scatter_wt_mut[["both_histogram"]][["plot"]] +
  ggplot2::scale_y_continuous(limits = c(0, 0.20))
ma_wt_mut <- extract_de_plots(limma_comparison, type = "limma")
ma_wt_mut[["ma"]][["plot"]]
ma_wt_mut[["volcano"]][["plot"]]

## ----simple_deseq2------------------------------------------------------------
deseq_comparison <- sm(deseq2_pairwise(fun_data))
summary(deseq_comparison[["all_tables"]][["wt30_vs_wt120"]])
scatter_wt_mut <- extract_coefficient_scatter(deseq_comparison, type = "deseq",
                                              x = "wt30", y = "wt120", gvis_filename = NULL)
scatter_wt_mut[["scatter"]]
plots_wt_mut <- extract_de_plots(deseq_comparison, type = "deseq")
plots_wt_mut[["ma"]][["plot"]]
plots_wt_mut[["volcano"]][["plot"]]

## ----simple_edger1------------------------------------------------------------
edger_comparison <- sm(edger_pairwise(fun_data, model_batch = TRUE))
plots_wt_mut <- extract_de_plots(edger_comparison, type = "edger")
scatter_wt_mut <- extract_coefficient_scatter(edger_comparison, type = "edger",
                                              x = "wt30", y = "wt120", gvis_filename = NULL)
scatter_wt_mut[["scatter"]]
plots_wt_mut[["ma"]][["plot"]]
plots_wt_mut[["volcano"]][["plot"]]

## ----simple_basic-------------------------------------------------------------
basic_comparison <- sm(basic_pairwise(fun_data))
summary(basic_comparison$all_tables$wt30_vs_wt120)
scatter_wt_mut <- extract_coefficient_scatter(basic_comparison, type = "basic",
                                              x = "wt30", y = "wt120")
scatter_wt_mut[["scatter"]]
plots_wt_mut <- extract_de_plots(basic_comparison, type = "basic")
plots_wt_mut[["ma"]][["plot"]]
plots_wt_mut[["volcano"]][["plot"]]

