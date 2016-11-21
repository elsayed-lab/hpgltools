#' Perform a simple transformation of a count table (log2)
#'
#' the add argument is only important if the data was previously cpm'd because that does a +1, thus
#' this will avoid a double+1 on the data.
#'
#' @param count_table  A matrix of count data
#' @param transform   A type of transformation to perform: log2/log10/log
#' @param base   for other log scales
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return dataframe of logx(counts)
#' @examples
#' \dontrun{
#' filtered_table = transform_counts(count_table, transform='log2', converted='cpm')
#' }
#' @export
transform_counts <- function(count_table, design=NULL, transform="raw",
                             base=NULL, ...) {
    arglist <- list(...)
    ## Short circuit this if we are going with raw data.
    if (transform == "raw") {
        libsize <- colSums(count_table)
        counts <- list(count_table=count_table, libsize=libsize)
        return(counts)
    }
    ## Also short-circuit it if we are asked to round the data.
    if (transform == "round") {
        count_table <- round(count_table)
        less_than <- count_table < 0
        count_table[less_than] <- 0
        libsize <- colSums(count_table)
        counts <- list(count_table=count_table, libsize=libsize)
        return(counts)
    }
    if (transform == "voom") {
        libsize <- colSums(count_table)
        model_choice <- choose_model(count_table, design[["condition"]],
                                     design[["batch"]])
        fun_model <- model_choice[["chosen_model"]]
        voomed <- limma::voom(counts=count_table, design=fun_model)
        counts <- list(count_table=voomed[["E"]], libsize=libsize)
        return(counts)
    }
    if (transform == "voomweight") {
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
    }

    ## If we are performing a transformation, then the minimum value I want is 1 before performing the logn
    less_zero <- sum(count_table < 0)
    if (less_zero > 0) {
        message(paste0("transform_counts: Found ", less_zero, " values less than 0."))
    }

    num_zero <- sum(count_table == 0)
    if (num_zero > 0) {
        message(paste0("transform_counts: Found ", num_zero, " values equal to 0, adding 1 to the matrix."))
        count_table <- count_table + 1
    }

    if (!is.null(base)) {
        count_table <- (log(count_table) / log(base))
    } else if (transform == "log2") {
        count_table <- log2(count_table)
    } else if (transform == "log10") {
        count_table <- log10(count_table)
    } else if (transform == "log") {  ## Natural log
        count_table <- log(count_table)  ## Apparently log1p does this.
    } else {
        message("Did not recognize the transformation, leaving the table.
 Recognized transformations include: 'log2', 'log10', 'log'
")
    }

    ## As a final check, remove any NaNs produced due to some shenanigans.
    num_before <- nrow(count_table)
    nans <- rowSums(is.nan(x=count_table))
    if (sum(nans) > 0) {
        nans <- nans == 0
        count_table <- count_table[nans, ]
        message(sprintf("Removing %d NaN containing rows (%d remaining).",
                        num_before - nrow(count_table), nrow(count_table)))
    }

    libsize <- colSums(count_table)
    counts <- list(
        "count_table" = count_table,
        "libsize" = libsize)
    return(counts)
}

## EOF
