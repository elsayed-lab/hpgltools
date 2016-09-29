## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE

## #' @import scales grid gtable
## #' @importFrom plyr defaults
## #' @importFrom stats setNames
## NULL

## Some notes!
##
## Something I should have realized:
## From: hcorrada@gmail.com [mailto:hcorrada@gmail.com] On Behalf Of Hector Corrada Bravo
## To: Najib M. El-Sayed
## Subject: Re: fold change
## Sorry for the delay in answering, I missed this. Using 'topTable' you can get confidence
## intervals for the logFC estimate. If you want the standard error itself you can get by dividing
## the logFC by the t statistic.

#' Pipe operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' dopar
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %dopar%
#' @rdname dopar
#' @keywords internal
#' @export
#' @importFrom foreach %dopar%
NULL

#' hpgltools: a suite of tools to make our analyses easier
#'
#' This provides a series of helpers for working with sequencing data
#'
#' It falls under a few main topics
#'
#' \itemize{
#' \item Data exploration, look for trends in sequencing data and identify batch effects or skewed distributions.
#' \item Differential expression analyses, use DESeq2/limma/EdgeR in a hopefully robust and flexible fashion.
#' \item Ontology analyses, use goseq/clusterProfiler/topGO/GOStats/gProfiler in hopefully robust ways.
#' \item Perform some simple TnSeq analyses.
#' }
#'
#' To see examples of this inaction, check out the vignettes:
#' \code{browseVignettes(package = 'hpgltools')}
#'
#' @docType package
#' @name hpgltools
#' @importFrom ggplot2 aes aes_string ggplot
#' @import data.table
#' @import foreach
#' @import grDevices
#' @import graphics
#' @import knitr
#' @import magrittr
#' @import methods
#' @import stats
#' @import utils
NULL
