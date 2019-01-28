start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("130model_tsne.R:
  12\n")
## 2018-03, exported functions in model_tsne:
## plot_tsne_genes(), plot_tsne(), tsne_res()

## This file is somewhat deprecated, as I now just treat tsne as
## another dimension reduction technique.

pombe_expt <- make_pombe_expt(annotation=FALSE)

## 01 plot_tsne_genes()
## I have no idea what to do with the result of plot_tsne_genes()
## It seems like I should take the coordinates provided by Rtsne and do a kmeans
## clustering or something on them in order to get the gene IDs associated with
## a cluster and then see if they are a part of some specific group or fall into
## a group of co-expressed genes or something.  But until I figure out something
## like that, these are not very useful.
##stuff <- plot_tsne_genes(pombe_expt)

## 02 plot_tsne()
## ok, so this is not actually doing tsne plots because I wanted to test
## something with pca.  It turns out that under the right conditions, we can get
## a nice pca split by time.  But the legend for time is messed up.
pombe_strain <- set_expt_conditions(expt=pombe_expt, fact="strain")
pombe_norm <- normalize_expt(pombe_strain, filter=TRUE, batch="fsva",
                              transform="log2", num_surrogates=1)

tsne_stuff <- plot_tsne(pombe_norm, size_column="minute")
##tsne_stuff$plot
expected <- c(36, 12)
actual <- dim(tsne_stuff[["table"]])
test_that("Did we get an expected tsne table?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## 03 tsne_res()

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 115model_pca.R in ", elapsed,  " seconds."))
