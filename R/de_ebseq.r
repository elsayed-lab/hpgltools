#' Set up model matrices contrasts and do pairwise comparisons of all conditions using EBSeq.
#'
#' Invoking EBSeq is confusing, this should help.
#'
#' @param input  Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions  Factor of conditions in the experiment.
#' @param patterns  Set of expression patterns to query.
#' @param ng_vector I think this is for isoform quantification, but am not yet
#'   certain.
#' @param rounds  Number of iterations for doing the multi-test
#' @param target_fdr  Definition of 'significant'
#' @param norm  Normalization method to use.
#' @param ... Extra arguments currently unused.
#' @export
ebseq_pairwise <- function(input=NULL, conditions=NULL, patterns=NULL,
                           ng_vector=NULL, rounds=10, target_fdr=0.05,
                           method="pairwise_subset", norm="median", ...) {
  arglist <- list(...)
  if (method == "pairwise_subset") {
    result <- ebseq_pairwise_subset(input, conditions=conditions,
                                    ng_vector=ng_vector, rounds=rounds,
                                    target_fdr=target_fdr, norm=norm,
                                    ...)
    return(result)
  }

  message("Starting EBSeq comparisons.")
  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force=force)
  ## Now that I understand pData a bit more, I should probably remove the conditions/batches slots
  ## from my expt classes.
  design <- pData(input)
  conditions <- design[["condition"]]
  batches <- design[["batches"]]
  data <- as.matrix(input_data[["data"]])
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  multi <- FALSE
  if (length(condition_levels) < 2) {
    stop("You have fewer than 2 conditions.")
  } else if (length(condition_levels) == 2) {
    message("Invoking ebseq with 2-condition parameters.")
    result <- ebseq_two(data, conditions,
                        ng_vector=ng_vector, rounds=rounds,
                        target_fdr=target_fdr, norm=norm)
  } else if (length(condition_levels) > 5) {
    if (is.null(patterns)) {
      stop("Beyond 5 conditions generates too many patterns, please provide a pattern matrix, or 'all_same'.")
      }
    message("Invoking ebseq with parameters for preset patterns.")
    result <- ebseq_many(data, conditions, patterns=patterns,
                         ng_vector=ng_vector, rounds=rounds,
                         target_fdr=target_fdr, norm=norm)
  } else {
    message("Invoking ebseq with parameters suitable for a few conditions.")
    result <- ebseq_few(data, conditions, patterns=patterns,
                        ng_vector=ng_vector, rounds=rounds,
                        target_fdr=target_fdr, norm=norm)
  }

  return(result)
}


ebseq_pairwise_subset <- function(input=NULL, conditions=NULL, patterns=NULL,
                                  ng_vector=NULL, rounds=10, target_fdr=0.05,
                                  method="pairwise_subset", norm="median", ...) {
  message("Starting EBSeq pairwise subset.")
  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force=force)
  ## Now that I understand pData a bit more, I should probably remove the conditions/batches slots
  ## from my expt classes.
  design <- pData(input)
  conditions <- design[["condition"]]
  batches <- design[["batches"]]
  data <- as.matrix(input_data[["data"]])
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  lenminus <- length(condition_levels) - 1
  retlst <- list()
  for (c in 1:lenminus) {
    c_name <- condition_levels[c]
    nextc <- c + 1
    for (d in nextc:length(condition_levels)) {
      d_name <- condition_levels[d]
      contrast_name <- paste0(c_name, "_vs_", d_name)
      message("Contrasting ", c_name, " and ", d_name, ".")
      pair <- subset_expt(
        expt=input,
        subset=paste0("condition=='",c_name,"' | condition=='", d_name, "'"))
      a_result <- ebseq_pairwise(pair, conditions=NULL,
                                 ng_vector=ng_vector,
                                 rounds=rounds, method="normal",
                                 target_fdr=target_fdr, norm=norm)
      retlst[[contrast_name]] <- a_result
    }
  }
  return(retlst)
}

ebseq_many <- function(data, conditions, patterns="all_same") {
  if (patterns == "all_same") {
    patterns <- data.frame(row.names="Pattern1")
    for (i in conditions) {
      patterns[1, i] <- 1
    }
  }
  normalized <- ebseq_size_factors(data, norm)
}

ebseq_size_factors <- function(data, norm=NULL) {
  ## Set up a null normalization vector
  normalized <- rep(x=1, times=ncol(data))
  names(normalized) <- colnames(data)
  ## If the parameter passed matches an ebseq function, use it, if it doesn't,
  ## we still have the null.
  if (norm == "median") {
    normalized <- EBSeq::MedianNorm(data)
  } else if (norm == "quantile") {
    normalized <- EBSeq::QuantileNorm(data)
  } else if (norm == "rank") {
    normalized <- EBSeq::RankNorm(data)
  }
  return(normalized)
}

ebseq_few <- function(data, conditions,
                      patterns=NULL, ng_vector=NULL, rounds=10,
                      target_fdr=0.05, norm="median") {

  ## Reminder about the meanings of 'patterns':
  ## Each row is a set of mean expression levels
  ## Each column is a condition
  ## Therefore, if a row has all three columns with a '1', then this pattern
  ## signifies a scenario in which all conditions have the same mean
  ## expression.
  ## If a row is '1', '1', '2', '2' then there are two sets of the same
  ## expression, 1 and 2.
  ## etc etc.
  if (is.null(patterns)) {
    patterns <- EBSeq::GetPatterns(conditions)
  } else if (patterns == "all_same") {
    patterns <- data.frame(row.names="Pattern1")
    for (i in conditions) {
      patterns[1, i] <- 1
    }
  }

  normalized <- ebseq_size_factors(data, norm)
  eb_output <- EBSeq::EBMultiTest(
                        Data=data, NgVector=ng_vector,
                        Conditions=conditions, AllParti=patterns,
                        sizeFactors=normalized, maxround=rounds)
  posteriors <- EBSeq::GetMultiPP(eb_output)
  fold_changes <- EBSeq::GetMultiFC(eb_output)

  pp_df <- as.data.frame(eb_output[["PPMat"]])
  ## Drop the pattern which looks for all the same
  interesting_patterns <- as.data.frame(patterns[-1, ])

  table_lst <- list()
  for (i in 1:ncol(fold_changes[["FCMat"]])) {
    column <- colnames(fold_changes[["FCMat"]])[i]
    contrast <- gsub(pattern="Over", replacement="_vs_", x=column)
    numerator <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\1", x=contrast)
    denominator <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\2", x=contrast)
    chosen_pattern_idx <- interesting_patterns[[numerator]] ==
      interesting_patterns[[denominator]]
    ## This should provide the name of the pattern which, when I query the
    ## PPMatrix will provide the likelihood that the two are the same
    ## I can therefore take the 1-that to get a likelihood different.
    chosen_pattern <- rownames(interesting_patterns)[chosen_pattern_idx]

    table <- data.frame(row.names=rownames(fold_changes[["FCMat"]]))
    table[["ebseq_FC"]] <- fold_changes[["FCMat"]][, i]
    table[["ebseq_logfc"]] <- fold_changes[["Log2FCMat"]][, i]
    table[["ebseq_postfc"]] <- fold_changes[["Log2PostFCMat"]][, i]
    table[["ebseq_mean"]] <- fold_changes[["CondMeans"]][, i]
    table[["pp_all_same"]] <- pp_df[["Pattern1"]]
    table[["pp_pair_diff"]] <- pp_df[[chosen_pattern]]
    table[["putative_p_value"]] <- 1 - table[["pp_pair_diff"]]
    table_lst[[contrast]] <- table
  }

  return(table_lst)
}

ebseq_two <- function(data, conditions,
                      ng_vector=NULL, rounds=10,
                      target_fdr=0.05, norm="median") {
  normalized <- ebseq_size_factors(data, norm=norm)
  eb_output <- EBSeq::EBTest(
                        Data=data, NgVector=NULL, Conditions=conditions,
                        sizeFactors=normalized, maxround=rounds)
  posteriors <- EBSeq::GetPP(eb_output)
  fold_changes <- EBSeq::PostFC(eb_output)
  eb_result <- EBSeq::GetDEResults(eb_output, FDR=target_fdr)
  table <- data.frame(row.names=names(fold_changes[["PostFC"]]))
  table[["ebseq_FC"]] <- fold_changes[["RealFC"]]
  table[["ebseq_logfc"]] <- log2(table[["ebseq_FC"]])
  table[["ebseq_postfc"]] <- fold_changes[["PostFC"]]
  table[["ebseq_mean"]] <- as.numeric(eb_output[["MeanList"]][[1]])
  table <- merge(table, as.data.frame(eb_result[["PPMat"]]),
                 by="row.names", all.x=TRUE)
  rownames(table) <- table[["Row.names"]]
  table <- table[, -1]
  return(table)
}

#' Writes out the results of a ebseq search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from ebseq and friends.
#'
#' @param data  Output from ebseq_pairwise()
#' @param ...  Options for writing the xlsx file.
#' @seealso \pkg{EBSeq} \link{write_xls}
#' @examples
#' \dontrun{
#'  finished_comparison = ebseq_pairwise(expressionset)
#'  data_list = write_ebseq(finished_comparison)
#' }
#' @export
write_ebseq <- function(data, ...) {
  result <- write_de_table(data, type="ebseq", ...)
  return(result)
}

## EOF
