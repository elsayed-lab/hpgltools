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

## ----setup, include=TRUE-------------------------------------------------
## These first 4 lines are not needed once hpgltools is installed.
## source("http://bioconductor.org/biocLite.R")
## biocLite("devtools")
## library(devtools)
## install_github("elsayed-lab/hpgltools")
library(hpgltools)
require.auto("fission")
library(fission)
data(fission)

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
fission_expt <- create_expt(meta, count_dataframe=fission_data)

## ----norm_explore--------------------------------------------------------
## First make a bar plot of the library sizes in the experiment.
## Notice that the colors were auto-chosen by create_expt() and they should
## be maintained throughout this process
fis_libsize <- plot_libsize(fission_expt)
fis_libsize
## Here we see that the wild type replicate 3 sample for 15 minutes has fewer non-zero genes than all its friends.
fis_nonzero <- plot_nonzero(fission_expt, labels="boring", title="nonzero vs. cpm")
fis_nonzero

## ----pca-----------------------------------------------------------------
## Something in this is causing a build loop on travis...

## I am no longer certain that code is maintained, what remains from it I might pull in to my own.
## require.auto("kokrah/cbcbSEQ")  ## Install Kwame's cbcbSEQ
## Unsurprisingly, the raw data doesn't cluster well at all...
fis_rawpca <- plot_pca(fission_expt, expt_labels=fission_expt$condition)
fis_rawpca$plot
## So, normalize the data
norm_expt <- normalize_expt(fission_expt, transform="log2", norm="quant", convert="cpm")
## And try the pca again
fis_normpca <- plot_pca(norm_expt, plot_labels="normal", title="normalized pca")
fis_normpca$plot

normbatch_expt <- normalize_expt(fission_expt, transform="log2", norm="quant", convert="cpm", batch="sva")
fis_normbatchpca <- plot_pca(normbatch_expt, title="Normalized PCA with batch effect correction.")
fis_normbatchpca$plot
## ok, that caused the 0, 60, 15, and 30 minute samples to cluster nicely
## the 120 and 180 minute samples are still a bit tight

## pca_information provides some more information about the call to
## fast.svd that went into making the pca plot
fis_info <- pca_information(norm_expt, expt_factors=c("condition","batch","strain","minute"), num_components=6)
## The r^2 table shows that quite a lot of the variance in the data is explained by condition
head(fis_info$rsquared_table)
## We can look at the correlation between the principle components and the factors in the experiment
## in this case looking at condition/batch vs the first 4 components.
fis_info$pca_cor
## And p-values to lend some credence(or not to those assertions)
fis_info$anova_p

## Try again with batch removed data
batchnorm_expt <- normalize_expt(fission_expt, batch="limma", norm="quant", transform="log2", convert="cpm")
fis_batchnormpca <- plot_pca(batchnorm_expt, plot_title="limma corrected pca")
fis_batchnormpca$plot
test_pca <- pca_information(batchnorm_expt, expt_factors=c("condition","batch","strain","minute"), num_components=6)

## ----distributions-------------------------------------------------------
plot_boxplot(fission_expt)
sf_expt <- normalize_expt(fission_expt, norm="sf")
plot_boxplot(sf_expt)
tm_expt <- normalize_expt(fission_expt, norm="tmm")
plot_boxplot(tm_expt)
rle_expt <- normalize_expt(fission_expt, norm="rle")
plot_boxplot(rle_expt)
up_expt <- normalize_expt(fission_expt, norm="upperquartile")
plot_boxplot(up_expt)

plot_density(norm_expt)
plot_density(sf_expt)
plot_density(tm_expt)

compare_12 <- plot_qq_plot(fission_expt, x=1, y=2)
compare_12$log

## ----clustering----------------------------------------------------------
plot_corheat(norm_expt)
plot_corheat(batchnorm_expt)
plot_disheat(norm_expt)
plot_disheat(batchnorm_expt)

