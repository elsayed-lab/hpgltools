#' Subset an expt
#'
#' Working on improving my understanding of how R sets up functions by type.
#' @param expt expt to subset
#' @param i Set of genes to keep
#' @param j Set of samples to keep.
#' @param ... Parameters to pass to subset_genes/subset_expt.
#' @export
`[.expt` <- function(expt, i, j, ...) {
  if (!missing(i)) {
    message("Subsetting on features.")
    expt <- subset_genes(expt, ids = i, method = "keep", ...)
  }
  if (!missing(j)) {
    message("Subsetting on samples.")
    expt <- subset_expt(expt, ids = j, ...)
  }
  if (missing(i) && missing(j)) {
    message("No subset was provided, returning the original.")
  }
  return(expt)
}

#' Simplifying subset on metadata.
#'
#' @param object an expt
#' @param i Column to extract
#' @export
setMethod(
  "[[", signature = c(x = "expt", i = "character"),
  definition = function(x, i, j, ...) {
    message("Got here?")
    pData(object)[[i]]
  })
setMethod(
  "[[", signature = c(x = "expt", i = "ANY"),
  definition = function(x, i) {
    message("Got here?")
    pData(x)[[i]]
  })

setMethod(
  "[[", signature = c("expt", "ANY", "missing"),
  definition = function(x, i, j, ...) {
    message("Got here?")
    pData(x)[[i]]
  })
setReplaceMethod(
  "[[", c("expt", "ANY", "missing"),
  function(x, i, j, ..., value) {
    pData(x)[[i, ...]] <- value
    x
  })
