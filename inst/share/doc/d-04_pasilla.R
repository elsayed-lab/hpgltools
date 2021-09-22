## ----options, include=FALSE---------------------------------------------------
library("hpgltools")
knitr::opts_knit$set(progress=TRUE,
                     verbose=TRUE,
                     width=90,
                     echo=TRUE)
knitr::opts_chunk$set(error=TRUE,
                      fig.width=8,
                      fig.height=8,
                      dpi=96)
old_options <- options(digits=4,
                       stringsAsFactors=FALSE,
                       knitr.duplicate.label="allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
set.seed(1)
ver <- "20170820"
rmd_file <- "d-04_pasilla.Rmd"

## ----load_data----------------------------------------------------------------
## I use sm to keep functions from printing too much (well, anything really)
tt <- sm(library(hpgltools))
tt <- sm(library(pasilla))
tt <- sm(data(pasillaGenes))

## ----biomart------------------------------------------------------------------
## Try loading some annotation information for this species.
gene_info_lst <- sm(load_biomart_annotations(species="dmelanogaster",
                                             host="useast.ensembl.org"))
gene_info <- gene_info_lst[["annotation"]]
info_idx <- gene_info[["gene_biotype"]] == "protein_coding"
gene_info <- gene_info[info_idx, ]
rownames(gene_info) <- make.names(gene_info[["ensembl_gene_id"]], unique=TRUE)
head(gene_info)

## ----load_counts--------------------------------------------------------------
## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata[["sampleid"]] <- rownames(metadata)

## Make sure it is still possible to create an expt
pasilla_expt <- sm(create_expt(count_dataframe=counts, metadata=metadata,
                               savefile="pasilla", gene_info=gene_info))

## ----graph_metrics, fig.show="hide"-------------------------------------------
pasilla_metrics <- sm(graph_metrics(pasilla_expt, ma=TRUE, qq=TRUE))
summary(pasilla_metrics)

## ----print_graphs-------------------------------------------------------------
pasilla_metrics$libsize
## The library sizes range from 8-21 million reads, this might be a problem for
## some analyses.
pasilla_metrics$nonzero
## Ergo, the lower abundance libraries have more genes of counts == 0 (bottom
## left).
pasilla_metrics$boxplot
## And a boxplot downshifts them (but not that much because it decided to put
## the data on the log scale).
pasilla_metrics$density
## Similarly, one can see those samples are a bit lower with respect to density

## Unless the data is very well behaved, the rest of the plots are not likely to
## look good until the data is normalized, nonetheless, lets see
pasilla_metrics$corheat
pasilla_metrics$disheat
pasilla_metrics$pc_plot
## So the above 3 plots are pretty much the worst case scenario for this data.

## ----normalize, fig.show="hide"-----------------------------------------------
norm <- default_norm(pasilla_expt, transform="log2")
norm_metrics <- graph_metrics(norm)

## ----show_norm----------------------------------------------------------------
norm_metrics$corheat
norm_metrics$smc
norm_metrics$disheat
norm_metrics$smd
## some samples look a little troublesome here.
norm_metrics$pc_plot

## ----perform_pairwise, fig.show="hide"----------------------------------------
pasilla_pairwise <- sm(all_pairwise(pasilla_expt))
pasilla_tables <- sm(combine_de_tables(
  pasilla_pairwise,
  excel="pasilla_tables.xlsx"))
pasilla_sig <- sm(extract_significant_genes(
  pasilla_tables,
  excel="pasilla_sig.xlsx"))
pasilla_ab <- sm(extract_abundant_genes(
  pasilla_pairwise,
  excel="pasilla_abundant.xlsx"))

## ----de_pictures--------------------------------------------------------------
pasilla_tables$plots[[1]][["deseq_ma_plots"]]$plot
pasilla_tables$plots[[1]][["edger_ma_plots"]]$plot
pasilla_tables$plots[[1]][["limma_ma_plots"]]$plot

## ----goseq_test---------------------------------------------------------------
up_genes <- pasilla_sig[["deseq"]][["ups"]][[1]]
down_genes <- pasilla_sig[["deseq"]][["downs"]][[1]]
pasilla_go <- load_biomart_go(species="dmelanogaster")$go
pasilla_length <- fData(pasilla_expt)[, c("ensembl_gene_id", "cds_length")]
colnames(pasilla_length) <- c("ID", "length")
pasilla_goseq <- simple_goseq(sig_genes=up_genes, go_db=pasilla_go, length_db=pasilla_length)
pasilla_goseq$pvalue_plots$bpp_plot_over

pasilla_goseq <- simple_goseq(sig_genes=down_genes, go_db=pasilla_go, length_db=pasilla_length)
pasilla_goseq$pvalue_plots$bpp_plot_over

test <- simple_goseq(sig_genes=names(pasilla_ab$abundances$deseq$treated), go_db=pasilla_go, length_db=pasilla_length)
test$pvalue_plots$bpp_plot_over

## ----saveme-------------------------------------------------------------------
pander::pander(sessionInfo())
message(paste0("This is hpgltools commit: ", get_git_commit()))

