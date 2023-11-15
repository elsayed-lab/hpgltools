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
    message("Subsetting on genes.")
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
