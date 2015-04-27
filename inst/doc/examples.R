## ----setup---------------------------------------------------------------

## These first 4 lines are not needed once hpgltools is installed.
##source("http://bioconductor.org/biocLite.R")
##biocLite("devtools")
##library(devtools)
##install_github("elsayed-lab/hpgltools")
library(hpgltools)
autoloads_all()
require.auto("fission")
library(fission)
data(fission)
opts_knit$set(progress=TRUE, verbose=TRUE, error=TRUE,  fig.width=7, fig.height=7)



## ----data_import---------------------------------------------------------

## Extract the meta data from the fission dataset
meta = as.data.frame(fission@colData)
## Make conditions and batches
meta$condition = paste(meta$strain, meta$minute, sep=".")
meta$batch = meta$replicate
meta$sample.id = rownames(meta)
## Write it out in the format expected by my toys
write.csv(meta, file="fission.csv")
## Grab the count data
fission_data = fission@assays$data$counts
## This will make an experiment superclass called 'expt' and it contains
## an ExpressionSet along with any arbitrary additional information one might want to include.
## Along the way it writes a Rdata file which is by default called 'expt.Rdata'
fission_expt = create_expt("fission.csv", count_dataframe=fission_data)


## ----norm_explore--------------------------------------------------------

## First make a bar plot of the library sizes in the experiment.
## Notice that the colors were auto-chosen by create_expt() and they should
## be maintained throughout this process
fis_libsize = hpgl_libsize(expt=fission_expt)
fis_libsize

## Here we see that the wild type replicate 3 sample for 15 minutes has fewer non-zero genes than all its friends.
fis_nonzero = hpgl_nonzero(expt=fission_expt, labels="boring", title="nonzero vs. cpm")
fis_nonzero


## ----pca-----------------------------------------------------------------

## Unsurprisingly, the raw data doesn't cluster well at all...
fis_rawpca = hpgl_pca(expt=fission_expt, labels=fission_expt$condition)
fis_rawpca$plot

## So, normalize the data
norm_expt = normalize_expt(fission_expt)
## And try the pca again
fis_normpca = hpgl_pca(expt=norm_expt, labels=norm_expt$condition)
fis_normpca$plot
## ok, that caused the 0, 60, 15, and 30 minute samples to cluster nicely
## the 120 and 180 minute samples are still a bit tight

## pca_information provides some more information about the call to
## fast.svd that went into making the pca plot
fis_info = pca_information(df=exprs(norm_expt$expressionset), design=norm_expt$design, factors=c("condition","batch"), num_components=4)
## The r^2 table shows that quite a lot of the variance in the data is explained by condition
head(fis_info$rsquared_table)
## We can look at the correlation between the principle components and the factors in the experiment
## in this case looking at condition/batch vs the first 4 components.
fis_info$pca_cor
## And p-values to lend some credence(or not to those assertions)
fis_info$anova_p


