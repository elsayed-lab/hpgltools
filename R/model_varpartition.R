## model_varpartition.r: Use VariancePartition to analyze the measureability of
## variance in a data set given different experimental models.  Variance
## Partition is fun, but weird.  This file seeks to simplify and standardize
## these methods.

#' A shortcut for replotting the percent plots from variancePartition.
#'
#' In case I wish to look at different numbers of genes from variancePartition
#' and/or different columns to sort from.
#'
#' @param varpart_output List returned by varpart()
#' @param n How many genes to plot.
#' @param column The df column to use for sorting.
#' @param decreasing high->low or vice versa?
#' @return The percent variance bar plots from variancePartition!
#' @seealso [variancePartition]
#' @export
replot_varpart_percent <- function(varpart_output, n = 30, column = NULL, decreasing = TRUE) {
  sorted <- varpart_output[["sorted_df"]]
  if (!is.null(column)) {
    if (column %in% colnames(sorted)) {
      sorted <- sorted[order(sorted[[column]], decreasing = decreasing), ]
    } else {
      message("The column ", column,
              "is not in the sorted data frame returned by varpart(). ",
              "Leaving the data frame alone.")
    }
  }
  new_plot <- variancePartition::plotPercentBars(sorted[1:n, ])
  retlist <- list(
      "resorted" = sorted,
      "plot" = new_plot)
  class(retlist) <- "reordered_varpart"
  return(retlist)
}

#' Use variancePartition to try and understand where the variance lies in a data set.
#'
#' The arguments and usage of variancePartition are a bit opaque.  This function
#' attempts to fill in reasonable values and simplify its invocation.
#'
#' @param expt Some data
#' @param predictor Non-categorical predictor factor with which to begin the
#'  model.
#' @param factors Character list of columns in the experiment design to query
#' @param chosen_factor When checking for sane 'batches', what column to
#'  extract from the design?
#' @param do_fit Perform a fitting using variancePartition?
#' @param cor_gene Provide a set of genes to look at the correlations, defaults
#'  to the first gene.
#' @param cpus Number cpus to use
#' @param genes Number of genes to count.
#' @param parallel Use doParallel?
#' @param strict_filter Perform a strict filtering of the results via median_by_factor and dropping
#'  any genes with a 0.
#' @param mixed Used a mixed model?
#' @param modify_expt Add annotation columns with the variance/factor?
#' @return List of plots and variance data frames
#' @seealso [variancePartition] DOI:10.1186/s12859-016-1323-z.
#' @export
simple_varpart <- function(expt, predictor = NULL, factors = c("condition", "batch"),
                           chosen_factor = "batch", do_fit = FALSE, cor_gene = 1,
                           cpus = NULL, genes = 40, parallel = TRUE, strict_filter = TRUE,
                           mixed = FALSE, modify_expt = TRUE) {
  cl <- NULL
  para <- NULL
  ## One is not supposed to use library() in packages, but it needs to do all
  ## sorts of foolish attaching.
  ## tt <- sm(library("variancePartition"))
  lib_result <- sm(requireNamespace("variancePartition"))
  att_result <- sm(try(attachNamespace("variancePartition"), silent = TRUE))
  lib_result <- sm(requireNamespace("BiocParallel"))
  att_result <- sm(try(attachNamespace("BiocParallel"), silent = TRUE))
  if (isTRUE(parallel)) {
    cl <- NULL
    if (is.null(cpus)) {
      cpus <- parallel::detectCores() - 4
      if (cpus < 1) {
        cpus <- 1
      }
      cl <- parallel::makeCluster(cpus)
    } else {
      cl <- parallel::makeCluster(cpus)
    }
    para <- doParallel::registerDoParallel(cl)
    ## multi <- BiocParallel::MulticoreParam()
  }
  design <- pData(expt)
  #for (f in factors) {
  #  design[[f]] <- as.factor(design[[f]])
  #}
  num_batches <- length(levels(as.factor(design[[chosen_factor]])))
  if (num_batches == 1) {
    message("varpart sees only 1 batch, adjusting the model accordingly.")
    factors <- factors[!grepl(pattern = chosen_factor, x = factors)]
  }

  model_string <- "~ "
  if (isTRUE(mixed)) {
    if (!is.null(predictor)) {
      model_string <- glue("{model_string}{predictor} + ")
    }
    for (fact in factors) {
      model_string <- glue("{model_string}(1|{fact}) + ")
    }
  } else {
    for (fact in factors) {
      model_string <- glue("{model_string}{fact} + ")
    }
  }
  model_string <- gsub(pattern = "\\+ $", replacement = "", x = model_string)
  if (isTRUE(mixed)) {
    mesg("Attempting mixed linear model with: ", model_string)
  } else {
    mesg("Attempting regular linear model with: ", model_string)
  }
  my_model <- as.formula(model_string)
  ## I think the simple filter is insufficient and I need there to be
  ## no genes with 0 counts in any one condition.
  norm <- sm(normalize_expt(expt, filter = "simple"))
  if (isTRUE(strict_filter)) {
    test <- sm(median_by_factor(norm, fact = "condition", fun = "mean"))
    all_condition_gt_zero_idx <- rowSums(test[["medians"]] == 0) == 0
    kept_gt <- rownames(exprs(norm))[all_condition_gt_zero_idx]
    norm <- norm[kept_gt, ]
  }
  data <- exprs(norm)

  design_sub <- design[, factors]
  mesg("Fitting the expressionset to the model, this is slow.")
  my_extract <- try(variancePartition::fitExtractVarPartModel(data, my_model, design_sub))
  ## my_extract <- try(variancePartition::fitVarPartModel(data, my_model, design))
  if ("try-error" %in% class(my_extract)) {
    mesg("A couple of common errors:
An error like 'vtv downdated' may be because there are too many 0s, filter the data and rerun.
An error like 'number of levels of each grouping factor must be < number of observations' means
that the factor used is not appropriate for the analysis - it really only works for factors
which are shared among multiple samples.")
    message("Retrying with only condition in the model.")
    my_model <- as.formula("~ condition")
    my_extract <- try(variancePartition::fitExtractVarPartModel(data, my_model, design))
    if ("try-error" %in% class(my_extract)) {
      message("Attempting again with only condition failed.")
      stop()
    }
  }
  chosen_column <- predictor
  if (is.null(predictor)) {
    chosen_column <- factors[[1]]
    mesg("Placing factor: ", chosen_column, " at the beginning of the model.")
  }

  ## A new dataset has some NAs!
  na_idx <- is.na(my_extract)
  if (sum(na_idx) > 0) {
    warning("There are ", sum(na_idx), " NAs in this data, something may be wrong.")
    message("There are ", sum(na_idx), " NAs in this data, something may be wrong.")
    message("Converting NAs to 0.")
    my_extract[na_idx] <- 0
  }

  my_sorted <- sortCols(my_extract)
  order_idx <- order(my_sorted[[chosen_column]], decreasing = TRUE)
  my_sorted <- my_sorted[order_idx, ]
  ## Recent error noticed when checking that variances sum to 1
  ## This is because sometimes we have smaller data sets
  if (genes > ncol(my_sorted)) {
    genes <- ncol(my_sorted)
  }

  percent_plot <- variancePartition::plotPercentBars(my_sorted[1:genes, ])
  partition_plot <- variancePartition::plotVarPart(my_sorted)

  fitting <- NULL
  stratify_batch_plot <- NULL
  stratify_condition_plot <- NULL
  if (isTRUE(do_fit)) {
    ## Try fitting with lmer4
    fitting <- variancePartition::fitVarPartModel(exprObj = data,
                                                  formula = my_model, data = design)
    last_fact <- factors[length(factors)]
    idx <- order(design[[chosen_column]], design[[last_fact]])
    ##first <- variancePartition::plotCorrStructure(fitting, reorder = idx)
    test_strat <- data.frame(Expression = data[3, ],
                             condition = design[[chosen_column]],
                             batch = design[[last_fact]])
    batch_expression <- as.formula("Expression ~ batch")
    cond_expression <- as.formula("Expression ~ condition")
    stratify_batch_plot <- variancePartition::plotStratify(batch_expression, test_strat)
    stratify_condition_plot <- variancePartition::plotStratify(cond_expression, test_strat)
  }

  if (isTRUE(parallel)) {
    para <- parallel::stopCluster(cl)
  }

  ret <- list(
    "model_string" = model_string,
    "model_used" = my_model,
    "percent_plot" = percent_plot,
    "partition_plot" = partition_plot,
    "sorted_df" = my_sorted,
    "fitted_df" = my_extract,
    "fitting" = fitting,
    "stratify_batch_plot" = stratify_batch_plot,
    "stratify_condition_plot" = stratify_condition_plot)
  if (isTRUE(modify_expt) && nrow(fData(expt)) > 0) {
    new_expt <- expt
    tmp_annot <- fData(new_expt)
    tmp_annot[["Row.names"]] <- NULL
    added_data <- my_sorted
    colnames(added_data) <- glue("variance_{colnames(added_data)}")
    ## Note that we are getting these variance numbers from data which was filtered
    ## thus we need all.x to get the IDs to match up.
    tmp_annot <- merge(tmp_annot, added_data, by = "row.names", all.x = TRUE)
    rownames(tmp_annot) <- tmp_annot[["Row.names"]]
    tmp_annot[["Row.names"]] <- NULL
    annot_order <- rownames(exprs(new_expt))
    tmp_annot <- tmp_annot[annot_order, ]
    ## Make it possible to use a generic expressionset, though maybe this is
    ## impossible for this function.
    fData(new_expt) <- tmp_annot
    ret[["modified_expt"]] <- new_expt
  }
  class(ret) <- "varpart"
  return(ret)
}

#' Attempt to use variancePartition's fitVarPartModel() function.
#'
#' Note the word 'attempt'.  This function is so ungodly slow that it probably
#' will never be used.
#'
#' @param expt Input expressionset.
#' @param factors Set of factors to query
#' @param cpus Number of cpus to use in doParallel.
#' @return Summaries of the new model,  in theory this would be a nicely
#'  batch-corrected data set.
#' @seealso [variancePartition]
varpart_summaries <- function(expt, factors = c("condition", "batch"), cpus = 6) {
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  model_string <- "~ "
  for (fact in factors) {
    model_string <- glue("{model_string} (1|{fact}) + ")
  }
  model_string <- gsub(pattern = "\\+ $", replacement = "", x = model_string)
  my_model <- as.formula(model_string)
  norm <- sm(normalize_expt(expt, filter = TRUE))
  data <- exprs(norm)
  design <- expt[["design"]]
  summaries <- variancePartition::fitVarPartModel(data, my_model, design, fxn = summary)
  return(summaries)
}

## EOF
