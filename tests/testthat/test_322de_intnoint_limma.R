start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

## I want to understand why I am getting different p-values between
## a model with an intercept and one without.  This should not be I think.

pasilla <- new.env()
load("pasilla.rda", envir=pasilla)
pasilla <- pasilla[["expt"]]

## First, let us invoke limma with batch in the model with no intercept.
## The following is mostly copy/pasted from the limma manual.

## Section 15.4 of the manual:
## Ok, so the first thing I am seeing is that one might not always want to use voom().
## But since all my work so far uses voom(), I will do these test with that.

counts <- exprs(pasilla)
design <- pData(pasilla)

noint_model <- stats::model.matrix(~ condition + batch, data=design)
int_model <- stats::model.matrix(~ 0 + condition + batch, data=design)

noint_voom <- limma::voom(counts, noint_model, plot=TRUE, normalize.method="quantile")
int_voom <- limma::voom(counts, int_model, plot=TRUE, normalize.method="quantile")

noint_fit <- limma::lmFit(noint_voom, noint_model)
int_fit <- limma::lmFit(int_voom, int_model)
int_contrast_matrix <- limma::makeContrasts(untreat_vs_treated=conditionuntreated-conditiontreated, levels=int_model)
int_contrasts <- limma::contrasts.fit(int_fit, int_contrast_matrix)

noint_eb <- suppressWarnings(limma::eBayes(noint_fit))
int_eb <- limma::eBayes(int_contrasts)

noint_table <- limma::topTable(noint_eb, coef="conditiontreated", n=Inf, adjust="BH")
int_table <- limma::topTable(int_eb, coef="untreat_vs_treated", n=Inf, adjust="BH")

## https://hopstat.wordpress.com/2014/06/26/be-careful-with-using-model-design-in-r/
## Toward the end of this document, I found that the big difference is in the setting
## of the intercept attribute to the model.  It looks like there is a test in summary.lm
## which checks for the 'intercept' attribute and changes the sum of squares accordingly.

## Ok, so I have a big misconception:
## https://support.bioconductor.org/p/13601/
## "If you specify the model with ~ 0, then you are saying that you don't want an intercept."
## I thought the 0 (or -1 in the example) _is_ the intercept, thus:
## y = b + mx <--> y ~ 0 + m     In this case b is the y intercept.

interaction_model <- stats::model.matrix(~ condition + batch + condition:batch, data=design)
interaction_voom <- limma::voom(counts, interaction_model, plot=TRUE, normalize.method="quantile")
interaction_fit <- limma::lmFit(interaction_voom, interaction_model)
interaction_eb <- suppressWarnings(limma::eBayes(interaction_fit))
interaction_table <- limma::topTable(interaction_eb, coef="conditiontreated", n=Inf, adjust="BH")
head(interaction_table)
