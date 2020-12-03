#' Call various count filters.
#'
#' This calls the various filtering functions in genefilter along with
#' suggestions made in our lab meetings; defaulting to the threshold based
#' filter suggested by Hector.
#'
#' @param count_table Some counts to filter.
#' @param method Filtering method to apply (cbcb, pofa, kofa, cv right now).
#' @param p Used by genefilter's pofa().
#' @param A Also for pofa().
#' @param k Used by genefilter's kofa().
#' @param cv_min Used by genefilter's cv().
#' @param cv_max Also used by cv().
#' @param thresh Minimum threshold across samples for cbcb.
#' @param min_samples Minimum number of samples for cbcb.
#' @param ... More options might be needed, especially if I fold cv/p/etc into ...
#' @return Data frame of filtered counts.
#' @seealso \pkg{genefilter}
#' @examples
#' \dontrun{
#'  new <- filter_counts(old)
#' }
#' @export
filter_counts <- function(count_table, method="cbcb", p=0.01, A=1, k=1,
                          cv_min=0.01, cv_max=1000, thresh=1, min_samples=2, ...) {
  arglist <- list(...)
  if (tolower(method) == "povera") {
    type <- "pofa"
  } else if (tolower(method) == "kovera") {
    type <- "kofa"
  }
  if (isTRUE(method)) {
    filter <- "cbcb"
  }
  filtered_counts <- NULL
  switchret <- switch(
    method,
    "cbcb" = {
      filtered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                            min_samples=min_samples)
    },
    "hpgl" = {
      filtered_counts <- hpgl_filter_counts(count_table, threshold=thresh,
                                            min_samples=min_samples)
    },
    "pofa" = {
      filtered_counts <- genefilter_pofa_counts(count_table, p=p, A=A)
    },
    "kofa" = {
      filtered_counts <- genefilter_kofa_counts(count_table, k=k, A=A)
    },
    "cv" = {
      filtered_counts <- genefilter_cv_counts(count_table, cv_min=cv_min,
                                              cv_max=cv_max)
    },
    "simple" = {
      filtered_counts <- simple_filter_counts(count_table, threshold=thresh)
    },
    {
      message("The requested filter did not match anything, defaulting to 'cbcb'.")
      filtered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                            min_samples=min_samples, ...)
    }
  ) ## Ending the switch
  return(filtered_counts)
}

#' Filter low-count genes from a data set using cpm data and a threshold.
#'
#' This was a function written by Kwame Okrah and perhaps also Laura Dillon to remove low-count
#' genes.  It drops genes based on a cpm threshold and number of samples.
#'
#' @param count_table Data frame of (pseudo)counts by sample.
#' @param threshold Lower threshold of counts for each gene.
#' @param min_samples Minimum number of samples.
#' @param libsize  Table of library sizes.
#' @return Dataframe of counts without the low-count genes.
#' @seealso \pkg{edgeR}
#' @examples
#' \dontrun{
#'  filtered_table <- cbcb_filter_counts(count_table)
#' }
#' @export
cbcb_filter_counts <- function(count_table, threshold=1, min_samples=2, libsize=NULL) {
  log2CPM <- function(qcounts, libsize=NULL) {
    if (is.null(libsize)) {
      libsize <- colSums(qcounts, na.rm=TRUE)
    }
    count_table <- t(log2(t(qcounts + 0.5) / (libsize + 1) * 1e+06))
    retlist <- list(
      "count_table" = count_table,
      "libsize" = libsize)
    return(retlist)
  }
  ##cpms <- edgeR::cpm(count_table)
  l2cpm <- log2CPM(count_table, libsize=libsize)
  cpms <- 2 ^ l2cpm[["count_table"]]
  keep <- rowSums(cpms > threshold, na.rm=TRUE) >= min_samples
  num_before <- nrow(count_table)
  count_table <- count_table[keep, ]
  message(sprintf("Removing %d low-count genes (%d remaining).",
                  num_before - nrow(count_table), nrow(count_table)))

  libsize <- l2cpm[["libsize"]]
  counts <- list(
    "count_table" = count_table,
    "libsize" = libsize)
  return(counts)
}

#' Filter low-count genes from a data set using cpm data and a threshold.
#'
#' This is identical to cbcb_filter_counts except it does not do the somewhat tortured
#' log2CPM() but instead just uses a 4 cpm non-log threshold.  It should therefore give
#' basically the same result, but without the shenanigans.
#'
#' @param count_table Data frame of (pseudo)counts by sample.
#' @param threshold Lower threshold of counts for each gene.
#' @param min_samples Minimum number of samples.
#' @param libsize  Table of library sizes.
#' @param ...  Arguments passed to cpm and friends.
#' @return Dataframe of counts without the low-count genes.
#' @seealso \pkg{edgeR}
#' @examples
#' \dontrun{
#'  filtered_table <- cbcb_filter_counts(count_table)
#' }
#' @export
hpgl_filter_counts <- function(count_table, threshold=2, min_samples=2, libsize=NULL, ...) {
  neg_idx <- count_table < 0
  neg_sum <- sum(neg_idx)
  if (sum(neg_sum) > 0) {
    warning("Found ", neg_sum, " negative entries, setting them to 0.")
    count_table[neg_idx] <- 0
  }
  cpms <- edgeR::cpm(count_table)
  keep <- rowSums(cpms > threshold) >= min_samples
  num_before <- nrow(count_table)
  count_table <- count_table[keep, ]
  num_after <- nrow(count_table)
  removed_rows <- num_before - num_after
  message("Removing ", removed_rows, " low-count genes (",
          num_after, " remaining).")
  libsize <- colSums(count_table)
  counts <- list(
    "count_table" = count_table,
    "libsize" = libsize)
  return(counts)
}

#' Filter low-count genes from a data set only using a simple threshold and number of samples.
#'
#' This was a function written by Kwame Okrah and perhaps also Laura Dillon to remove low-count
#' genes.  It drops genes based on a threshold and number of samples.
#'
#' @param count_table Data frame of (pseudo)counts by sample.
#' @param threshold Lower threshold of counts for each gene.
#' @return Dataframe of counts without the low-count genes.
#' @seealso \pkg{edgeR}
#' @examples
#' \dontrun{
#'  filtered_table <- simple_filter_counts(count_table)
#' }
#' @export
simple_filter_counts <- function(count_table, threshold=2) {
  num_before <- nrow(count_table)
  sums <- rowSums(count_table)
  keepers <- (sums >= threshold)
  count_table <- count_table[keepers, ]

  message(sprintf("Removing %d low-count genes (%d remaining).",
                  num_before - nrow(count_table), nrow(count_table)))

  libsize <- colSums(count_table)
  counts <- list("count_table" = count_table,
                 "libsize" = libsize)
  return(counts)
}

#' Filter low-count genes from a data set using genefilter's pOverA().
#'
#' I keep thinking this function is pofa... oh well.  Of the various tools in
#' genefilter, this one to me is the most intuitive.  Take the ratio of
#' counts/samples and make sure it is >= a score.
#'
#' @param count_table Input data frame of counts by sample.
#' @param p Minimum proportion of each gene's counts/sample to be greater than a
#'   minimum(A).
#' @param A Minimum number of counts in the above proportion.
#' @return Dataframe of counts without the low-count genes.
#' @seealso \pkg{genefilter}
#'  \code{\link[genefilter]{pOverA}}
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_pofa_counts(count_table)
#' }
#' @export
genefilter_pofa_counts <- function(count_table, p=0.01, A=100) {
  ## genefilter has functions to work with expressionsets directly, but I think
  ## I will work merely with tables in this.
  num_before <- nrow(count_table)

  if ("ExpressionSet" %in% class(count_table)) {
    counts <- exprs(count_table)
  }
  test <- genefilter::pOverA(p=p, A=A)
  filter_list <- genefilter::filterfun(test)
  answer <- genefilter::genefilter(count_table, filter_list)
  count_table <- count_table[answer, ]

  removed <- num_before - nrow(count_table)
  message("Removing ", removed, " low-count genes (", nrow(count_table), " remaining).")

  libsize <- colSums(count_table)
  counts <- list(count_table=count_table, libsize=libsize)
  return(counts)
}

#' Filter genes from a dataset outside a range of variance.
#'
#' This function from genefilter removes genes surpassing a variance cutoff.  It
#' is not therefore a low-count filter per se.
#'
#' @param count_table Input data frame of counts by sample.
#' @param cv_min Minimum coefficient of variance.
#' @param cv_max Maximum coefficient of variance.
#' @return Dataframe of counts without the high/low variance genes.
#' @seealso \pkg{genefilter}
#'  \code{\link[genefilter]{kOverA}}
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_cv_counts <- function(count_table, cv_min=0.01, cv_max=1000) {
  ## genefilter has functions to work with expressionsets directly, but I think
  ## I will work merely with tables in this.
  num_before <- nrow(count_table)

  if ("ExpressionSet" %in% class(count_table)) {
    counts <- exprs(count_table)
  }
  test <- genefilter::cv(cv_min, cv_max)
  filter_list <- genefilter::filterfun(test)
  answer <- genefilter::genefilter(count_table, filter_list)
  count_table <- count_table[answer, ]

  message(sprintf("Removing %d low-count genes (%d remaining).",
                  num_before - nrow(count_table), nrow(count_table)))
  libsize <- colSums(count_table)
  counts <- list(count_table=count_table, libsize=libsize)
  return(counts)
}

#' Filter low-count genes from a data set using genefilter's kOverA().
#'
#' This is the most similar to the function suggested by Hector I think.
#'
#' @param count_table Input data frame of counts by sample.
#' @param k Minimum number of samples to have >A counts.
#' @param A Minimum number of counts for each gene's sample in kOverA().
#' @return Dataframe of counts without the low-count genes.
#' @seealso \pkg{genefilter}
#'  \code{\link[genefilter]{kOverA}}
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_kofa_counts <- function(count_table, k=1, A=1) {
  ## genefilter has functions to work with expressionsets directly, but I think
  ## I will work merely with tables in this.
  num_before <- nrow(count_table)

  if ("ExpressionSet" %in% class(count_table)) {
    counts <- exprs(count_table)
  }
  test <- genefilter::kOverA(k=k, A=A)
  filter_list <- genefilter::filterfun(test)
  answer <- genefilter::genefilter(count_table, filter_list)
  count_table <- count_table[answer, ]

  message(sprintf("Removing %d low-count genes (%d remaining).",
                  num_before - nrow(count_table), nrow(count_table)))
  libsize <- colSums(count_table)
  counts <- list(count_table=count_table, libsize=libsize)
  return(counts)
}

## EOF
