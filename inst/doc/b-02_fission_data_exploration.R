## ----options, include=FALSE----------------------------------------------
## These are the options I tend to favor
library("hpgltools")
devtools::load_all("~/hpgltools")
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
rmd_file <- "b-02_fission_data_exploration.Rmd"

## ----setup---------------------------------------------------------------
if (! "BiocInstaller" %in% installed.packages()) {
  source("http://bioconductor.org/biocLite.R")
}
if (! "devtools" %in% installed.packages()) {
  biocLite("devtools")
}
if (! "hpgltools" %in% installed.packages()) {
  devtools::install_github("elsayed-lab/hpgltools")
}
if (! "fission" %in% installed.packages()) {
  biocLite("fission")
}
## I use the function 'sm()' to quiet loud functions.
tt <- sm(library(fission))
tt <- sm(data(fission))

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

## ----norm_explore--------------------------------------------------------
## First make a bar plot of the library sizes in the experiment.
## Notice that the colors were auto-chosen by create_expt() and they should
## be maintained throughout this process
fis_libsize <- plot_libsize(fission_expt)
fis_libsize$plot
## Here we see that the wild type replicate 3 sample for 15 minutes has fewer non-zero genes than all its friends.
fis_nonzero <- plot_nonzero(fission_expt, labels="boring", title="nonzero vs. cpm")
fis_nonzero$plot

fis_density <- plot_density(fission_expt)
fis_density$plot
fis_density$condition_summary
fis_density$batch_summary
fis_density$sample_summary

## ----pca-----------------------------------------------------------------
## Something in this is causing a build loop on travis...
## Unsurprisingly, the raw data doesn't cluster well at all...
##fis_rawpca <- plot_pca(fission_expt, expt_names=fission_expt$condition, cis=NULL)
fis_rawpca <- plot_pca(fission_expt)
fis_rawpca$plot
## So, normalize the data
norm_expt <- sm(normalize_expt(fission_expt, transform="log2", norm="quant", convert="cpm"))
## And try the pca again
fis_normpca <- plot_pca(norm_expt, plot_labels="normal", title="normalized pca", cis=NULL)
fis_normpca$plot

summary(fis_normpca)

## ----test_3d-------------------------------------------------------------
testing <- plot_pca(fission_expt, num_pc=3)
silly <- make_3d_pca(testing)
silly$plot

## ----normalized_pca------------------------------------------------------
normbatch_expt <- sm(normalize_expt(fission_expt, transform="log2", norm="quant",
                                    convert="cpm", batch="sva"))
fis_normbatchpca <- plot_pca(normbatch_expt,
                             title="Normalized PCA with batch effect correction.", cis=NULL)
fis_normbatchpca$plot
## ok, that caused the 0, 60, 15, and 30 minute samples to cluster nicely
## the 120 and 180 minute samples are still a bit tight

## pca_information provides some more information about the call to
## fast.svd that went into making the pca plot
fis_info <- pca_information(norm_expt,
                            expt_factors=c("condition","batch","strain","minute"),
                            num_components=6)
## The r^2 table shows that quite a lot of the variance in the data is explained by condition
head(fis_info$rsquared_table)
## We can look at the correlation between the principle components and the factors in the experiment
## in this case looking at condition/batch vs the first 4 components.
fis_info$pca_cor
## And p-values to lend some credence(or not to those assertions)
fis_info$anova_p

## Try again with batch removed data
batchnorm_expt <- sm(normalize_expt(fission_expt, batch="limma", norm="quant",
                                    transform="log2", convert="cpm"))
fis_batchnormpca <- plot_pca(batchnorm_expt, plot_title="limma corrected pca")
fis_batchnormpca$plot
test_pca <- pca_information(batchnorm_expt,
                            expt_factors=c("condition","batch","strain","minute"),
                            num_components=6)

## ----distributions-------------------------------------------------------
fission_boxplot <- sm(plot_boxplot(fission_expt))
fission_boxplot
sf_expt <- sm(normalize_expt(fission_expt, norm="sf"))
fission_boxplot <- sm(plot_boxplot(sf_expt))
fission_boxplot
tm_expt <- sm(normalize_expt(fission_expt, norm="tmm"))
fission_boxplot <- sm(plot_boxplot(tm_expt))
fission_boxplot
rle_expt <- sm(normalize_expt(fission_expt, norm="rle"))
fission_boxplot <- sm(plot_boxplot(rle_expt))
fission_boxplot
up_expt <- sm(normalize_expt(fission_expt, norm="upperquartile"))
fission_boxplot <- sm(plot_boxplot(up_expt))
fission_boxplot

fission_density <- plot_density(norm_expt)
fission_density$plot
fission_density <- plot_density(sf_expt)
fission_density$plot
fission_density <- plot_density(tm_expt)
fission_density$plot

compare_12 <- plot_single_qq(fission_expt, x=1, y=2)
compare_12$log

## ----clustering----------------------------------------------------------
fission_cor <- plot_corheat(norm_expt)
fission_cor$plot
fission_cor <- plot_corheat(batchnorm_expt)
fission_cor$plot
fission_dis <- plot_disheat(norm_expt)
fission_dis$plot
fission_dis <- plot_disheat(batchnorm_expt)
fission_dis$plot

## ----variancePartition---------------------------------------------------
test_varpart <- varpart(fission_expt, predictor=NULL, factors=c("condition","batch"))
test_varpart$percent_plot
test_varpart$partition_plot

## Here, let us test the variance contributed by strain, time, and replicate.
test_varpart <- varpart(fission_expt, predictor=NULL, factors=c("condition", "strain", "minute", "replicate"))
test_varpart$percent_plot
test_varpart$partition_plot

## ----sysinfo, results='asis'---------------------------------------------
pander::pander(sessionInfo())

