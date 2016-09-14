## ----options, include=FALSE----------------------------------------------
## These are the options I tend to favor
library("hpgltools")
knitr::opts_knit$set(
    progress = TRUE,
    verbose = TRUE,
    width = 90,
    echo = TRUE)
knitr::opts_chunk$set(
    error = TRUE,
    fig.width = 8,
    fig.height = 8,
    dpi = 96)
options(
    digits = 4,
    stringsAsFactors = FALSE,
    knitr.duplicate.label = "allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
set.seed(1)
rmd_file <- "c-03_fission_differential_expression.Rmd"

## ----setup, include=TRUE-------------------------------------------------
## These first 4 lines are not needed once hpgltools is installed.
## source("http://bioconductor.org/biocLite.R")
## biocLite("devtools")
## library(devtools)
## install_github("elsayed-lab/hpgltools")
library(hpgltools)
require.auto("fission")
tt <- sm(library(fission))
tt <- data(fission)
knitr::opts_knit$set(progress=TRUE, verbose=TRUE, error=TRUE,  fig.width=7, fig.height=7)

## ----data_import---------------------------------------------------------
## Extract the meta data from the fission dataset
meta <- as.data.frame(fission@colData)
## Make conditions and batches
meta$condition <- paste(meta$strain, meta$minute, sep=".")
meta$batch <- meta$replicate
meta$sample.id <- rownames(meta)
## Grab the count data
fission_data <- fission@assays$data$counts
## This will make an experiment superclass called 'expt' and it contains
## an ExpressionSet along with any arbitrary additional information one might want to include.
## Along the way it writes a Rdata file which is by default called 'expt.Rdata'
fission_expt <- create_expt(metadata=meta, count_dataframe=fission_data)

## ----simple_subset-------------------------------------------------------
fun_data <- expt_subset(fission_expt, subset="condition=='wt.120'|condition=='wt.30'")
fun_norm <- sm(normalize_expt(fun_data, batch="limma", norm="quant", transform="log2", convert="cpm"))

## ----simple_limma--------------------------------------------------------
limma_comparison <- sm(limma_pairwise(fun_norm))
names(limma_comparison$all_tables)
summary(limma_comparison$all_tables$wt.30_vs_wt.120)
scatter_wt_mut <- limma_coefficient_scatter(limma_comparison, x="wt.30", y="wt.120", gvis_filename=NULL)
scatter_wt_mut$scatter
scatter_wt_mut$both_histogram
ma_wt_mut <- limma_ma(limma_comparison)
ma_wt_mut
## So sad
ma_wt_mut <- limma_ma(limma_comparison, p_col="P.Value", fc=0.5)
ma_wt_mut

## ----simple_deseq2-------------------------------------------------------
deseq_comparison <- sm(deseq2_pairwise(fun_data))
summary(deseq_comparison$all_tables$wt.30_vs_wt.120)
scatter_wt_mut <- deseq_coefficient_scatter(deseq_comparison, x="wt.30", y="wt.120", gvis_filename=NULL)
scatter_wt_mut$scatter
ma_wt_mut <- deseq_ma(deseq_comparison)
ma_wt_mut
ma_wt_mut <- deseq_ma(deseq_comparison, fc=0.5)
ma_wt_mut

## ----simple_edger--------------------------------------------------------
edger_comparison <- sm(edger_pairwise(fun_data, model_batch=TRUE))
scatter_wt_mut <- edger_coefficient_scatter(edger_comparison, x="wt.30", y="wt.120", gvis_filename=NULL)
scatter_wt_mut$scatter
ma_wt_mut <- edger_ma(edger_comparison)
ma_wt_mut

## ----simple_basic--------------------------------------------------------
basic_comparison <- sm(basic_pairwise(fun_data))
summary(basic_comparison$all_tables$wt.30_vs_wt.120)
scatter_wt_mut <- basic_coefficient_scatter(basic_comparison)
ma_wt_mut <- basic_ma(basic_comparison)

## ----simple_all----------------------------------------------------------
all_comparisons <- sm(all_pairwise(fun_data, model_batch=TRUE))
all_combined <- sm(combine_de_tables(all_comparisons))
head(all_combined)
sig_genes <- sm(extract_significant_genes(all_combined, excel=NULL))
head(sig_genes)

## ----ontology_setup------------------------------------------------------
limma_results <- limma_comparison$all_tables
## The set of comparisons performed
names(limma_results)
table <- limma_results$wt.30_vs_wt.120
dim(table)
gene_names <- rownames(table)

updown_genes <- get_sig_genes(table, p=0.1, fc=0.8, p_column="P.Value")
tt <- require.auto("GenomicFeatures")
tt <- require.auto("biomaRt")
ensembl_pombe <- biomaRt::useMart("fungal_mart", dataset="spombe_eg_gene", host="fungi.ensembl.org")
pombe_filters <- biomaRt::listFilters(ensembl_pombe)
head(pombe_filters, n=20) ## 11 looks to be my guy

possible_pombe_attributes <- biomaRt::listAttributes(ensembl_pombe)
##pombe_goids <- biomaRt::getBM(attributes=c('pombase_gene_name', 'go_accession'), filters="biotype", values=gene_names, mart=ensembl_pombe)
pombe_goids <- biomaRt::getBM(attributes=c('pombase_transcript', 'go_accession'), values=gene_names, mart=ensembl_pombe)
colnames(pombe_goids) <- c("ID","GO")

pombe <- sm(GenomicFeatures::makeTxDbFromBiomart(biomart="fungal_mart",
                                                 dataset="spombe_eg_gene",
                                                 host="fungi.ensembl.org"))
pombe_transcripts <- as.data.frame(GenomicFeatures::transcriptsBy(pombe))
lengths <- pombe_transcripts[, c("group_name","width")]
colnames(lengths) <- c("ID","width")
## Something useful I didn't notice before:
## makeTranscriptDbFromGFF()  ## From GenomicFeatures, much like my own gff2df()
gff_from_txdb <- GenomicFeatures::asGFF(pombe)
## why is GeneID: getting prefixed to the IDs!?
gff_from_txdb$ID <- gsub(x=gff_from_txdb$ID, pattern="GeneID:", replacement="")
written_gff <- rtracklayer::export.gff3(gff_from_txdb, con="pombe.gff")

## ----test_goseq----------------------------------------------------------
summary(updown_genes)
funkytown <- simple_goseq(de_genes=updown_genes$down_genes, go_db=pombe_goids, length_db=lengths)

