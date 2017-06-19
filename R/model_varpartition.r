#' A shortcut for replotting the percent plots from variancePartition.
#'
#' In case I wish to look at different numbers of genes from variancePartition and/or
#' different columns to sort from.
#'
#' @param varpart_output  List returned by varpart()
#' @param n  How many genes to plot.
#' @param column  The df column to use for sorting.
#' @param decreasing  high->low or vice versa?
#' @return  The percent variance bar plots from variancePartition!
#' @seealso \pkg{variancePartition}
#'  \code{\link[variancePartition]{plotPercentBars}}
#' @export
replot_varpart_percent <- function(varpart_output, n=30, column=NULL, decreasing=TRUE) {
    sorted <- varpart_output[["sorted_df"]]
    if (!is.null(column)) {
        if (column %in% colnames(sorted)) {
            sorted <- sorted[ order(sorted[[column]], decreasing=decreasing), ]
        } else {
            message(paste0("The column ", column, "is not in the sorted data frame returned by varpart()."))
            message("Leaving the data frame alone.")
        }
    }
    new_plot <- variancePartition::plotPercentBars(sorted[1:n, ])
    return(new_plot)
}

#' Use variancePartition to try and understand where the variance lies in a data set.
#'
#' variancePartition is the newest toy introduced by Hector.
#'
#' @param expt  Some data
#' @param predictor  Non-categorical predictor factor with which to begin the model.
#' @param factors  Character list of columns in the experiment design to query
#' @param cpus  Number cpus to use
#' @param genes  Number of genes to count.
#' @param parallel  use doParallel?
#' @return partitions  List of plots and variance data frames
#' @seealso \pkg{doParallel} \pkg{variancePartition}
#' @export
varpart <- function(expt, predictor="condition", factors=c("batch"),
                    cpus=6, genes=40, parallel=TRUE) {
    cl <- NULL
    para <- NULL
    if (isTRUE(parallel)) {
        cl <- parallel::makeCluster(cpus)
        para <- doParallel::registerDoParallel(cl)
    }
    num_batches <- length(levels(as.factor(expt[["batches"]])))
    if (num_batches == 1) {
        message("varpart sees only 1 batch, adjusting the model accordingly.")
        factors <- factors[!grepl(pattern="batch", x=factors)]
    }
    model_string <- "~ "
    if (!is.null(predictor)) {
        model_string <- paste0(model_string, predictor, " +")
    }
    for (fact in factors) {
        model_string <- paste0(model_string, " (1|", fact, ") +")
    }
    model_string <- gsub(pattern="\\+$", replacement="", x=model_string)
    message(paste0("Attempting mixed linear model with: ", model_string))
    my_model <- as.formula(model_string)
    norm <- sm(normalize_expt(expt, filter=TRUE))
    data <- Biobase::exprs(norm[["expressionset"]])
    design <- expt[["design"]]
    message("Fitting the expressionset to the model, this is slow.")
    message("(Eg. Take the projected run time and mulitply by 3-6 and round up.)")
    ##my_fit <- try(variancePartition::fitVarPartModel(data, my_model, design))
    ##message("Extracting the variances.")
    ##my_extract <- try(variancePartition::extractVarPart(my_fit))
    my_extract <- try(variancePartition::fitExtractVarPartModel(data, my_model, design))
    if (class(my_extract) == "try-error") {
        stop("An error like 'vtv downdated' may be because there are too many 0s, try and filter the data and rerun.")
    }
    my_sorted <- variancePartition::sortCols(my_extract)
    my_sorted <- my_sorted[ order(my_sorted[["condition"]], decreasing=TRUE), ]
    percent_plot <- variancePartition::plotPercentBars(my_sorted[1:genes, ])
    partition_plot <- variancePartition::plotVarPart(my_sorted)
    if (isTRUE(parallel)) {
        para <- parallel::stopCluster(cl)
    }
    ret <- list(
        "model_used" = my_model,
        "percent_plot" = percent_plot,
        "partition_plot" = partition_plot,
        "sorted_df" = my_sorted,
        "fitted_df" = my_extract)
    return(ret)
}

#' Attempt to use variancePartition's fitVarPartModel() function.
#'
#' Note the word 'attempt'.  This function is so ungodly slow that it probably will never be used.
#'
#' @param expt  Input expressionset.
#' @param factors  Set of factors to query
#' @param cpus  Number of cpus to use in doParallel.
#' @return  Summaries of the new model,  in theory this would be a nicely batch-corrected data set.
#' @seealso \pkg{variancePartition}
varpart_summaries <- function(expt, factors=c("condition", "batch"), cpus=6) {
    cl <- parallel::makeCluster(cpus)
    doParallel::registerDoParallel(cl)
    model_string <- paste0("~ ")
    for (fact in factors) {
        model_string <- paste0(model_string, " (1|", fact, ") + ")
    }
    model_string <- gsub(pattern="\\+ $", replacement="", x=model_string)
    my_model <- as.formula(model_string)
    norm <- sm(normalize_expt(expt, filter=TRUE))
    data <- Biobase::exprs(norm[["expressionset"]])
    design <- expt[["design"]]
    summaries <- variancePartition::fitVarPartModel(data, my_model, design, fxn=summary)
    return(summaries)
}

## EOF
