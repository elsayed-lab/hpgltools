#' Use variancePartition to gather some plots about a troubling dataset
#'
#' variancePartition is the newest toy introduced by Hector
#'
#' @param expt  Some data
#' @param factors  Character list of columns in the experiment design to query
#' @param cpus  Number cpus to use
#' @param genes Number of genes to count
#' @return partitions  List of plots and variance data frames
varpart <- function(expt, factors=c("condition","batch"), cpus=6, genes=20) {
    cl <- parallel::makeCluster(cpus)  ## I am keeping 2 processors to myself, piss off, R.
    doParallel::registerDoParallel(cl)
    model_string <- paste0("~ ")
    for (fact in factors) {
        model_string <- paste0(model_string, "(1|", fact, ") + ")
    }
    model_string <- gsub(pattern="+ ", replacement="", x=model_string)
    my_model <- as.formula(model_string)
    norm_par <- normalize_expt(expt, filter=TRUE)
    data <- Biobase::exprs(norm_par[["expressionset"]])
    design <- expt[["design"]]
    my_fit <- variancePartition::fitExtractVarPartModel(data, my_model, design)
    my_sorted <- variancePartition::sortCols(my_fit)
    percent_plot <- variancePartition::plotPercentBars(my_sorted[1:genes, ])
    partition_plot <- variancePartition::plotVarPart(my_sorted)
    ret <- list(
        "percent_plot" = percent_plot,
        "partition_plot" = partition_plot,
        "sorted_df" = my_sorted,
        "fitted_df" = my_fit)
    return(ret)
}
