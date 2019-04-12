#' Perform a simple transformation of a count table (log2)
#'
#' the add argument is only important if the data was previously cpm'd because
#' that does a +1, thus this will avoid a double+1 on the data.
#'
#' @param count_table  A matrix of count data
#' @param design  Sometimes the experimental design is also required.
#' @param transform   A type of transformation to perform: log2/log10/log.
#' @param base   Other log scales?
#' @param ...  Options I might pass from other functions are dropped into
#'   arglist.
#' @return dataframe of transformed counts.
#' @seealso \pkg{limma}
#' @examples
#' \dontrun{
#'  filtered_table = transform_counts(count_table, transform='log2', converted='cpm')
#' }
#' @export
transform_counts <- function(count_table, design=NULL, transform="raw",
                             base=NULL, ...) {
  arglist <- list(...)
  ## Short circuit this if we are going with raw data.


  switchret <- switch(
    transform,
    "raw" = {
      libsize <- colSums(count_table)
      counts <- list(count_table=count_table, libsize=libsize)
      return(counts)
    },
    "round" = {
      ## Also short-circuit it if we are asked to round the data.
      count_table <- round(count_table)
      less_than <- count_table < 0
      count_table[less_than] <- 0
      libsize <- colSums(count_table)
      counts <- list(count_table=count_table, libsize=libsize)
      return(counts)
    },
    "voom" = {
      libsize <- colSums(count_table)
      model_choice <- choose_model(count_table,
                                   conditions=design[["condition"]],
                                   batches=design[["batch"]], ...)
      fun_model <- model_choice[["chosen_model"]]
      voomed <- limma::voom(counts=count_table, design=fun_model)
      counts <- list(count_table=voomed[["E"]], libsize=libsize)
      return(counts)
    },
    "voomweight" = {
      libsize <- colSums(count_table)
      model_choice <- choose_model(count_table, design[["condition"]],
                                   design[["batch"]])
      fun_model <- model_choice[["chosen_model"]]
      voomed <- try(limma::voomWithQualityWeights(counts=count_table, design=fun_model))
      if (class(voomed) == "try-error") {
        warning("voomwithqualityweights failed.  Falling back to voom.")
        voomed <- limma::voom(counts=count_table, design=fun_model)
      }
      counts <- list(count_table=voomed[["E"]], libsize=libsize)
      return(counts)
    },
    {
    }
  ) ## Ending the switch statement

  ## If we are performing a transformation, then the minimum value I want is 1
  ## before performing the logn
  less_zero_counts <- count_table < 0
  less_zero <- sum(less_zero_counts, na.rm=TRUE)
  if (less_zero > 0) {
    message("transform_counts: Found ", less_zero, " values less than 0.")
  }

  num_zero_counts <- count_table == 0
  num_zero <- sum(num_zero_counts, na.rm=TRUE)
  if (num_zero > 0) {
    message("transform_counts: Found ", num_zero,
            " values equal to 0, adding 1 to the matrix.")
    count_table <- count_table + 1
  }

  if (!is.null(base)) {
    count_table <- (log(count_table) / log(base))
  } else if (transform == "log2") {
    count_table <- log2(count_table)
  } else if (transform == "log10") {
    count_table <- log10(count_table)
  } else if (transform == "log") {
    ## Natural log
    count_table <- log(count_table)  ## Apparently log1p does this.
  } else {
    message("Did not recognize the transformation, leaving the table.
 Recognized transformations include: 'log2', 'log10', 'log'
")
  }
  count_table <- as.matrix(count_table)
  ## As a final check, remove any NaNs produced due to some shenanigans.
  ## This logic was removed, which is causing some unintended consequences.
  ## We should consider further how to deal with this.
  num_before <- nrow(count_table)
  nans <- is.nan(count_table)
  ## print(head(nans))
  nans_sum <- sum(nans)
  ##if (nans_sum > 0) {
    ##message("Setting nan containing elements to 0.")
    ##count_table[nans] <- 0
    ##nans <- nans == 0
    ##count_table <- count_table[!nans, ]
    ##message(sprintf("Removing %d NaN containing rows (%d remaining).",
    ##                num_before - nrow(count_table), nrow(count_table)))
  ##}

  libsize <- colSums(count_table)
  counts <- list(
    "count_table" = count_table,
    "libsize" = libsize)
  return(counts)
}

## EOF
