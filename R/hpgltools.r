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

#' data.table's funky column assignment operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name :=
#' @rdname column_assignment
#' @keywords internal
#' @export
#' @importFrom data.table :=
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
#' \item Data exploration, look for trends in sequencing data and identify batch
#'       effects or skewed distributions.
#' \item Differential expression analyses, use DESeq2/limma/EdgeR in a hopefully
#'       robust and flexible fashion.
#' \item Ontology analyses, use goseq/clusterProfiler/topGO/GOStats/gProfiler in
#'       hopefully robust ways.
#' \item Perform some simple TnSeq analyses.
#' }
#'
#' To see examples of this in action, check out the vignettes:
#' \code{browseVignettes(package = 'hpgltools')}
#'
#' @docType package
#' @name hpgltools
#' @importFrom Biobase exprs pData fData notes sampleNames
#' @importFrom data.table data.table
#' @importFrom dplyr filter group_by n summarise
#' @importFrom foreach foreach
#' @importFrom ggplot2 aes aes_string ggplot
#' @importFrom glue glue glue_data
#' @importFrom grDevices recordPlot
#' @importFrom rlang abort sym
#' @importFrom stats
#'  aggregate as.dendrogram as.formula ave biplot coef coefficients complete.cases
#'  cor cor.test density dist dnorm formula glm hclust lm lowess median
#'  model.matrix na.omit order.dendrogram p.adjust p.adjust.methods pnorm
#'  princomp quantile relevel reorder resid residuals rnbinom sd setNames
#'  t.test var
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
NULL

#' The following sets the ggplot2 default text size.
base_size <- 16

#' Set the xlsx table style
table_style <- "TableStyleMedium9"

#' R CMD check is super annoying about :::.
#'
#' In a fit of pique, I did a google search to see if anyone else has been
#' annoyed in the same way as was I.  Yihui Xie was, and in his email to r-devel
#' in 2013 he proposed a game of hide-and-seek; which I am repeating here.
#'
#' This just implements ::: as an infix operator that will not trip check.
#'
#' @param pkg on the left hand side
#' @param fun on the right hand side
`%:::%` <- function(pkg, fun) {
  get(fun, envir = asNamespace(pkg), inherits = FALSE)
}

getMaintainer <- "GenomicFeatures" %:::% ".getMaintainer"
getMetaDataValue <- "GenomicFeatures" %:::% ".getMetaDataValue"
getTxDbVersion <- "GenomicFeatures" %:::% ".getTxDbVersion"
normAuthor <- "GenomicFeatures" %:::% ".normAuthor"
sortCols <- "variancePartition" %:::% ".sortCols"
aprior <- "sva" %:::% "aprior"
bprior <- "sva" %:::% "bprior"
it.sol <- "sva" %:::% "it.sol"
int.eprior <- "sva" %:::% "int.eprior"
getGOLevel <- "clusterProfiler" %:::% "getGOLevel"
