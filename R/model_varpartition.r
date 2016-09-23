#' Use variancePartition to gather some plots about a troubling dataset
#'
#' variancePartition is the newest toy introduced by Hector
#'
#' @param expt  Some data
#' @param predictor  Non-categorical predictor factor with which to begin the model.
#' @param factors  Character list of columns in the experiment design to query
#' @param cpus  Number cpus to use
#' @param genes Number of genes to count
#' @return partitions  List of plots and variance data frames
#' @export
varpart <- function(expt, predictor="condition", factors=c("batch"), cpus=6, genes=20, parallel=TRUE) {
    cl <- NULL
    if (isTRUE(parallel)) {
        cl <- parallel::makeCluster(cpus)  ## I am keeping 2 processors to myself, piss off, R.
        doParallel::registerDoParallel(cl)
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
    ##my_fit <- try(variancePartition::fitVarPartModel(data, my_model, design))
    ##message("Extracting the variances.")
    ##my_extract <- try(variancePartition::extractVarPart(my_fit))
    my_extract <- try(variancePartition::fitExtractVarPartModel(data, my_model, design))
    if (class(my_extract) == "try-error") {
        stop("An error like 'vtv downdated' may be because there are too many 0s, try and filter the data and rerun.")
    }
    my_sorted <- variancePartition::sortCols(my_extract)
    percent_plot <- variancePartition::plotPercentBars(my_sorted[1:genes, ])
    partition_plot <- variancePartition::plotVarPart(my_sorted)
    if (isTRUE(parallel)) {
        parallel::stopCluster(cl)
    }
    ret <- list(
        "model_used" = my_model,
        "percent_plot" = percent_plot,
        "partition_plot" = partition_plot,
        "sorted_df" = my_sorted,
        "fitted_df" = my_extract)
    return(ret)
}


varpart_summaries <- function(expt, factors=c("condition","batch"), cpus=6) {
    cl <- parallel::makeCluster(cpus)  ## I am keeping 2 processors to myself, piss off, R.
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
