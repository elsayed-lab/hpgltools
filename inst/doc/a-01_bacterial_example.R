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
rmd_file <- "01_hpgltools_bacterial_examples.Rmd"

## ----rendering, include=FALSE, eval=FALSE--------------------------------
#  ## This block is used to render a document from within it.
#  rmarkdown::render(rmd_file)
#  
#  ## An extra renderer for pdf output
#  rmarkdown::render(rmd_file, output_format="pdf_document", output_options=c("skip_html"))
#  
#  ## Or to save/load large Rdata files.
#  hpgltools:::saveme()
#  hpgltools:::loadme()
#  rm(list=ls())

## ----loading_data--------------------------------------------------------
library(hpgltools)
data_file <- system.file("cdm_expt.rda", package="hpgltools")
cdm <- new.env()
load(data_file, envir=cdm)
rm(data_file)

ls()
## 2 variables should exist now: rmd_file in case I want to knitr this file, cdm which is a list including the data required to make an expressionset.

expt <- create_expt(count_dataframe=cdm$cdm_counts, metadata=cdm$cdm_metadata, gene_info=cdm$gene_info)
## The gff information is in 'annotations'
## The experiment is in most_v0M1
## Here is the meta-data! (well, the first 6 lines anyway).
knitr::kable(head(expt$design))
summary(expt)

## ----graph_original------------------------------------------------------
raw_metrics <- sm(graph_metrics(expt, qq=TRUE))
## View a raw library size plot
raw_metrics$libsize
## Or boxplot to see the data distribution
raw_metrics$boxplot
## The warning is because it automatically uses a log scale and there are some 0 count genes.
## Perhaps you prefer density plots
raw_metrics$density
## quantile/quantile plots compared to the median of all samples
raw_metrics$qqrat
## Here we can see some samples are differently 'shaped' compared to the median than others
## There are other plots one may view, but this data set is a bit too crowded as is.
## The following summary shows the other available plots:
summary(raw_metrics)

## ----subset_data---------------------------------------------------------
head(expt$design)
## elt stands for: "early/late in thy"
batch_a <- expt_subset(expt, subset="batch=='a'")
batch_b <- expt_subset(expt, subset="batch=='b'")

a_metrics <- graph_metrics(batch_a)
a_metrics$pcaplot
b_metrics <- graph_metrics(batch_b)
b_metrics$pcaplot

## ----normalize_subset----------------------------------------------------
## doing nothing to the data except log2 transforming it has a surprisingly large effect
norm_test <- normalize_expt(expt, transform="log2")
plot_pca(norm_test)$plot
## Looks like 2 samples got switched
## a quantile normalization alone affect some, but not all of the data
norm_test <- sm(normalize_expt(expt, norm="quant"))
plot_pca(norm_test)$plot
## cpm alone brings out some samples, too
norm_test <- sm(normalize_expt(expt, convert="cpm"))
plot_pca(norm_test)$plot
## low count filtering has some effect, too
norm_test <- sm(normalize_expt(expt, filter="pofa"))
plot_pca(norm_test)$plot
## how about if we mix and match methods?
norm_test <- sm(normalize_expt(expt, transform="log2", convert="cpm", norm="quant", batch="combat_scale", filter=TRUE, batch_step=4, low_to_zero=TRUE))
plot_pca(norm_test)$plot
## The different batch effect testing methods have a pretty widely ranging effect on the clustering
## play with them by changing the batch= parameter to:
## "limma", "sva", "svaseq", "limmaresid", "ruvg", "combat", combatmod"
pca_test <- plot_pca(norm_test)
head(pca_test$res)
## Thus we see a dramatic decrease in variance accounted for
## by batch after applying limma's 'removebatcheffect'
## (see batch.R2 here vs. above)

## Some metrics are not very useful on (especially quantile) normalized data
norm_graphs <- sm(graph_metrics(norm_test))
norm_graphs$smc
norm_graphs$disheat  ## svaseq's batch correction seems to draw out the signal quite nicely.
## It is worth noting that the wt, early log, thy, replicate c samples are still a bit weird.
norm_graphs$pcaplot

## ----de_test-------------------------------------------------------------
spyogenes_de <- sm(all_pairwise(expt))
spyogenes_tables <- sm(combine_de_tables(spyogenes_de))
spyogenes_sig <- sm(extract_significant_genes(spyogenes_tables))
knitr::kable(head(spyogenes_sig$limma$ups[[1]]))

## ----circos--------------------------------------------------------------
microbe_ids <- as.character(get_microbesonline_ids("pyogenes MGAS5005"))
mgas_df <- sm(get_microbesonline_annotation(microbe_ids[[1]])[[1]])
mgas_df$sysName <- gsub(pattern="Spy_", replacement="Spy", x=mgas_df$sysName)
rownames(mgas_df) <- make.names(mgas_df$sysName, unique=TRUE)


## First make a template configuration
circos_test <- circos_prefix()
## Fill it in with the data for s.pyogenes
circos_kary <- circos_karyotype("mgas", length=1895017)
## Fill in the gene category annotations by gene-strand
circos_plus <- sm(circos_plus_minus(mgas_df, circos_test))

circos_limma_hist <- sm(circos_hist(spyogenes_de$limma$all_tables[[1]], mgas_df, circos_test, outer=circos_plus))
circos_deseq_hist <- sm(circos_hist(spyogenes_de$deseq$all_tables[[1]], mgas_df, circos_test, outer=circos_limma_hist))
circos_edger_hist <- sm(circos_hist(spyogenes_de$edger$all_tables[[1]], mgas_df, circos_test, outer=circos_deseq_hist))
circos_suffix(cfgout=circos_test)
## circos_made <- sm(circos_make(target="mgas"))
## For some reason this fails weirdly when not run interactively.

## ----sysinfo, results='asis'---------------------------------------------
pander::pander(sessionInfo())

