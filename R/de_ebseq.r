#' Set up model matrices contrasts and do pairwise comparisons of all conditions using EBSeq.
#'
#' Invoking EBSeq is confusing, this should help.
#'
#' @param input  Dataframe/vector or expt class containing data, normalization state, etc.
#' @param patterns  Set of expression patterns to query.
#' @param ng_vector I think this is for isoform quantification, but am not yet
#'   certain.
#' @param rounds  Number of iterations for doing the multi-test
#' @param target_fdr  Definition of 'significant'
#' @param method  The default ebseq methodology is to create the set of all
#'   possible 'patterns' in the data; for data sets which are more than
#'   trivially complex, this is not tenable, so this defaults to subsetting the
#'   data into pairs of conditions.
#' @param norm  Normalization method to use.
#' @param ... Extra arguments currently unused.
#' @export
ebseq_pairwise <- function(input=NULL, patterns=NULL,
                           ng_vector=NULL, rounds=10, target_fdr=0.05,
                           method="pairwise_subset", norm="median", ...) {
  arglist <- list(...)

  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force=force)
  design <- pData(input)
  conditions <- design[["condition"]]
  batches <- design[["batches"]]
  data <- as.matrix(input_data[["data"]])
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  if (method == "pairwise_subset") {
    result <- ebseq_pairwise_subset(input,
                                    ng_vector=ng_vector, rounds=rounds,
                                    target_fdr=target_fdr, norm=norm,
                                    ...)
  } else {
    message("Starting single EBSeq invocation.")

    multi <- FALSE
    if (length(condition_levels) < 2) {
      stop("You have fewer than 2 conditions.")
    } else if (length(condition_levels) == 2) {
      message("Invoking ebseq with 2-condition parameters.")
      result <- ebseq_two(input,
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
  }

  retlist <- list(
    "all_tables" = result,
    "conditions" = conditions,
    "conditions_table" = conditions_table,
    "method" = "ebseq"
  )
  return(retlist)
}

ebseq_pairwise_subset <- function(input, ng_vector=NULL, rounds=10, target_fdr=0.05,
                                  model_batch=FALSE, method="pairwise_subset",
                                  model_cond=TRUE, model_intercept=FALSE, alt_model=NULL,
                                  conditions=NULL, norm="median", ...) {
  message("Starting EBSeq pairwise subset.")
  ## Now that I understand pData a bit more, I should probably remove the conditions/batches slots
  ## from my expt classes.
  design <- pData(input)
  conditions <- design[["condition"]]
  batches <- design[["batches"]]
  data <- exprs(input)
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  model_choice <- choose_model(input, conditions=conditions,
                               batch=batches,
                               model_batch=FALSE,
                               model_cond=TRUE,
                               model_intercept=FALSE,
                               alt_model=NULL, ...)
  model_data <- model_choice[["chosen_model"]]
  apc <- make_pairwise_contrasts(model_data, conditions, do_identities=FALSE, ...)
  contrasts_performed <- c()
  bar <- utils::txtProgressBar(style=3)
  retlst <- list()
  for (c in 1:length(apc[["names"]])) {
    pct_done <- c / length(apc[["names"]])
    utils::setTxtProgressBar(bar, pct_done)
    name  <- apc[["names"]][[c]]
    a_name <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\1", x=name)
    b_name <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\2", x=name)
    utils::setTxtProgressBar(bar, pct_done)
    pair <- sm(subset_expt(
      expt=input,
      subset=paste0("condition=='", b_name, "' | condition=='", a_name, "'")))
    pair_data <- exprs(pair)
    conditions <- pair[["conditions"]]
    a_result <- ebseq_two(pair_data, conditions,
                          numerator=b_name,
                          denominator=a_name,
                          ng_vector=ng_vector,
                          rounds=rounds,
                          target_fdr=target_fdr,
                          norm=norm)
    retlst[[name]] <- a_result
  }
  close(bar)
  return(retlst)
}

ebseq_many <- function(data, conditions, patterns="all_same",
                       ng_vector=ng_vector, rounds=rounds,
                       target_fdr=target_fdr, norm=norm) {

  if (patterns == "all_same") {
    patterns <- data.frame(row.names="Pattern1")
    for (i in conditions) {
      patterns[1, i] <- 1
    }
  }
  normalized <- ebseq_size_factors(data, norm)
  ## Not yet implemented.
  return(NULL)
}

ebseq_size_factors <- function(data_mtrx, norm=NULL) {
  ## Set up a null normalization vector
  normalized <- rep(x=1, times=ncol(data_mtrx))
  names(normalized) <- colnames(data_mtrx)
  ## If the parameter passed matches an ebseq function, use it, if it doesn't,
  ## we still have the null.
  if (norm == "median") {
    normalized <- EBSeq::MedianNorm(data_mtrx, alternative=TRUE)
  } else if (norm == "quantile") {
    normalized <- EBSeq::QuantileNorm(data_mtrx)
  } else if (norm == "rank") {
    normalized <- EBSeq::RankNorm(data_mtrx)
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
    table[["logFC"]] <- fold_changes[["Log2FCMat"]][, i]
    table[["ebseq_postfc"]] <- fold_changes[["Log2PostFCMat"]][, i]
    table[["ebseq_mean"]] <- fold_changes[["CondMeans"]][, i]
    table[["pp_all_same"]] <- pp_df[["Pattern1"]]
    table[["pp_pair_diff"]] <- pp_df[[chosen_pattern]]
    table[["putative_p_value"]] <- 1 - table[["pp_pair_diff"]]
    table_lst[[contrast]] <- table
  }

  retlst <- list(
    "all_tables" = table_lst,
    "conditions" = conditions,
    "method" = "ebseq")
  return(retlst)
}

ebseq_two <- function(pair_data, conditions,
                      numerator=2, denominator=1,
                      ng_vector=NULL, rounds=10,
                      target_fdr=0.05, norm="median") {
  normalized <- ebseq_size_factors(pair_data, norm=norm)
  message("Starting EBTest of ", numerator, " vs. ", denominator, ".")
  eb_output <- sm(EBSeq::EBTest(
                           Data=pair_data, NgVector=NULL, Conditions=conditions,
                           sizeFactors=normalized, maxround=rounds))
  posteriors <- EBSeq::GetPP(eb_output)
  fold_changes <- EBSeq::PostFC(eb_output)
  eb_result <- EBSeq::GetDEResults(eb_output, FDR=target_fdr)
  table <- data.frame(row.names=names(fold_changes[["PostFC"]]))
  table[["ebseq_FC"]] <- fold_changes[["RealFC"]]
  table[["logFC"]] <- log2(table[["ebseq_FC"]])
  table[["ebseq_postfc"]] <- fold_changes[["PostFC"]]
  table[["ebseq_mean"]] <- as.numeric(eb_output[["MeanList"]][[1]])
  table <- merge(table, as.data.frame(eb_result[["PPMat"]]),
                 by="row.names", all.x=TRUE)
  rownames(table) <- table[["Row.names"]]
  table <- table[, -1]
  ## This is incorrect I think, but being used as a placeholder until I figure out how to
  ## properly adjust a set prior probabilities.
  message("Copying the ppee values as an ajusted p-value until I figure out how to deal with them.")
  table[["ebseq_adjp"]] <- table[["PPEE"]]

  ## Finally, make sure the 'direction' matches my conception of numerator/denominator.
  eb_direction <- fold_changes[["Direction"]]
  eb_numerator <- gsub(pattern="^(.*) Over (.*)$", replacement="\\1", x=eb_direction)
  eb_denominator <- gsub(pattern="^(.*) Over (.*)$", replacement="\\2", x=eb_direction)
  ## message("EBD: ", eb_denominator, " EBN: ", eb_numerator, " D: ", denominator, " N: ", numerator)
  if (! eb_numerator == denominator) {
    ## message("Flipping the table FC and logFC.")
    table[["ebseq_FC"]] <- 1 / table[["ebseq_FC"]]
    table[["logFC"]] <- -1 * table[["logFC"]]
    table[["ebseq_postfc"]] <- 1 / table[["ebseq_postfc"]]
  }
  return(table)
}

## EOF
