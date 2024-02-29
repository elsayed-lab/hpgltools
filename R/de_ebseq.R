## de_ebseq.r: Simplify the invocation of EBSeq, the Bayesian statistical method
## of performing differential expression analysis.  EBSeq provides a third
## family of DE method: DESeq/EdgeR provides explicitly negative binomial
## distribution aware methods, while limma assumes microarrayish.  Both these
## families use various statistical models.  In contrast, EBSeq uses prior
## probabilities and sampling to gather its values for each gene.

#' Set up model matrices contrasts and do pairwise comparisons of all conditions
#' using EBSeq.
#'
#' Invoking EBSeq is confusing, this should help.
#'
#' @param input Dataframe/vector or expt class containing data, normalization
#'  state, etc.
#' @param patterns Set of expression patterns to query.
#' @param ng_vector I think this is for isoform quantification, but am not yet
#'  certain.
#' @param rounds Number of iterations for doing the multi-test
#' @param target_fdr Definition of 'significant'
#' @param method The default ebseq methodology is to create the set of all
#'  possible 'patterns' in the data; for data sets which are more than
#'  trivially complex, this is not tenable, so this defaults to subsetting the
#'  data into pairs of conditions.
#' @param norm Normalization method to use.
#' @param conditions Not currently used, but passed from all_pairwise()
#' @param batches Not currently used, but passed from all_pairwise()
#' @param model_cond Not currently used, but passed from all_pairwise()
#' @param model_intercept Not currently used, but passed from all_pairwise()
#' @param alt_model Not currently used, but passed from all_pairwise()
#' @param model_batch Not currently used, but passed from all_pairwise()
#' @param keepers Perform a specific set of contrasts instead of all?
#' @param force Force ebseq to accept bad data (notably NA containing stuff from proteomics.
#' @param ... Extra arguments currently unused.
#' @return List containing tables from ebseq, the conditions tested, and the
#'  ebseq table of conditions.
#' @seealso [limma_pairwise()] [deseq_pairwise()] [edger_pairwise()] [basic_pairwise()]
#' @examples
#'  \dontrun{
#'   expt <- create_expt(metadata = "sample_sheet.xlsx", gene_info = annotations)
#'   ebseq_de <- ebseq_pairwise(input = expt)
#' }
#' @export
ebseq_pairwise <- function(input = NULL, patterns = NULL, conditions = NULL,
                           batches = NULL, model_cond = NULL, model_intercept = NULL,
                           alt_model = NULL, model_batch = NULL, keepers = NULL,
                           ng_vector = NULL, rounds = 10, target_fdr = 0.05,
                           method = "pairwise_subset", norm = "median",
                           force = FALSE,
                           ...) {
  arglist <- list(...)

  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force = force)
  design <- pData(input)
  conditions <- pData(input)[["condition"]]
  batches <- pData(input)[["batches"]]
  data <- as.matrix(input_data[["data"]])
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  if (method == "pairwise_subset") {
    result <- ebseq_pairwise_subset(input,
                                    ng_vector = ng_vector, rounds = rounds,
                                    target_fdr = target_fdr, norm = norm, force = force,
                                    ...)
  } else {
    mesg("Starting single EBSeq invocation.")

    multi <- FALSE
    if (length(condition_levels) < 2) {
      stop("You have fewer than 2 conditions.")
    } else if (length(condition_levels) == 2) {
      mesg("Invoking ebseq with 2-condition parameters.")
      result <- ebseq_two(input,
                          ng_vector = ng_vector, rounds = rounds,
                          target_fdr = target_fdr, norm = norm)
    } else if (length(condition_levels) > 5) {
      stop("Beyond 5 conditions generates too many patterns, ",
           "please provide a pattern matrix, or 'all_same'.")
    } else {
      mesg("Invoking ebseq with parameters suitable for a few conditions.")
      result <- ebseq_few(data, conditions, patterns = patterns,
                          ng_vector = ng_vector, rounds = rounds,
                          target_fdr = target_fdr, norm = norm)
    }
  }

  retlist <- list(
      "all_tables" = result,
      "conditions" = conditions,
      "conditions_table" = conditions_table,
      "method" = "ebseq"
  )
  class(retlist) <- c("ebseq_pairwise", "list")
  return(retlist)
}

#' Perform pairwise comparisons with ebseq, one at a time.
#'
#' This uses the same logic as in the various *_pairwise functions to invoke
#' the 'normal' ebseq pairwise comparison for each pair of conditions in an
#' expressionset.  It therefore avoids the strange logic inherent in the ebseq
#' multitest function.
#'
#' @param input Expressionset/expt to perform de upon.
#' @param ng_vector Passed on to ebseq, I forget what this does.
#' @param rounds Passed on to ebseq, I think it defines how many iterations to
#'  perform before return the de estimates
#' @param target_fdr If we reach this fdr before iterating rounds times, return.
#' @param model_batch Provided by all_pairwise()  I do not think a Bayesian
#'  analysis really cares about models, but if one wished to try to add a batch
#'  factor, this would be the place to do it.  It is currently ignored.
#' @param model_cond Provided by all_pairwise(), ibid.
#' @param model_intercept Ibid.
#' @param alt_model Ibid.
#' @param keepers Specify a set of contrasts to perform here.
#' @param conditions Factor of conditions in the data, used to define the
#'  contrasts.
#' @param norm EBseq normalization method to apply to the data.
#' @param force Flag used to force inappropriate data into the various methods.
#' @param ... Extra arguments passed downstream, noably to choose_model()
#' @return A pairwise comparison of the various conditions in the data.
#' @seealso [ebseq_pairwise()]
ebseq_pairwise_subset <- function(input, ng_vector = NULL, rounds = 10, target_fdr = 0.05,
                                  model_batch = FALSE, model_cond = TRUE,
                                  model_intercept = FALSE, alt_model = NULL, keepers = NULL,
                                  conditions = NULL, norm = "median", force = FALSE, ...) {
  mesg("Starting EBSeq pairwise subset.")
  ## Now that I understand pData a bit more, I should probably remove the
  ## conditions/batches slots from my expt classes.
  design <- pData(input)
  conditions <- design[["condition"]]
  batches <- design[["batches"]]
  data <- exprs(input)
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  model_choice <- choose_model(
      input, conditions = conditions, batches = batches,
      model_batch = FALSE, model_cond = TRUE, model_intercept = FALSE, alt_model = NULL,
      ...)
  model_data <- model_choice[["chosen_model"]]
  apc <- make_pairwise_contrasts(model_data, conditions, do_identities = FALSE,
                                 do_extras = FALSE, keepers = keepers,
                                 ...)
  contrasts_performed <- c()
  retlst <- list()
  for (c in seq_along(apc[["names"]])) {
    name  <- apc[["names"]][[c]]
    a_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\1", x = name)
    b_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\2", x = name)
    if (! a_name %in% input[["conditions"]]) {
      message("The contrast ", a_name, " is not in the results.")
      message("If this is not an extra contrast, then this is an error.")
      next
    }
    pair <- sm(subset_expt(
        expt = input,
        subset = glue("condition=='{b_name}' | condition=='{a_name}'")))
    pair_data <- exprs(pair)
    conditions <- pair[["conditions"]]
    a_result <- ebseq_two(pair_data, conditions, numerator = b_name, denominator = a_name,
                          ng_vector = ng_vector, rounds = rounds, target_fdr = target_fdr,
                          norm = norm, force = force)
    retlst[[name]] <- a_result
  }
  return(retlst)
}

#' Choose the ebseq normalization method to apply to the data.
#'
#' EBSeq provides three normaliation methods.  Median, Quantile, and Rank.
#' Choose among them here.
#'
#' @param data_mtrx This is exprs(expressionset)
#' @param norm The method to pass along.
#' @return a new matrix using the ebseq specific method of choice.
#' @seealso [EBSeq]
ebseq_size_factors <- function(data_mtrx, norm = NULL) {
  ## Set up a null normalization vector
  normalized <- rep(x = 1, times = ncol(data_mtrx))
  names(normalized) <- colnames(data_mtrx)
  ## If the parameter passed matches an ebseq function, use it, if it doesn't,
  ## we still have the null.
  if (norm == "median") {
    normalized <- EBSeq::MedianNorm(data_mtrx, alternative = TRUE)
  } else if (norm == "quantile") {
    normalized <- EBSeq::QuantileNorm(data_mtrx)
  } else if (norm == "rank") {
    normalized <- EBSeq::RankNorm(data_mtrx)
  } else {
    normalized <- data_mtrx
  }
  return(normalized)
}

#' Invoke EBMultiTest() when we do not have too many conditions to deal with.
#'
#' Starting at approximately 5 conditions, ebseq becomes too unwieldy to use
#' effectively. But, its results until then are pretty neat.
#'
#' @param data Expressionset/matrix
#' @param conditions Factor of conditions in the data to compare.
#' @param patterns Set of patterns as described in the ebseq documentation to query.
#' @param ng_vector Passed along to ebmultitest().
#' @param rounds Passed to ebseq.
#' @param target_fdr Passed to ebseq.
#' @param norm Normalization method to apply to the data.
#' @seealso [ebseq_pairwise()]
ebseq_few <- function(data, conditions,
                      patterns = NULL, ng_vector = NULL, rounds = 10,
                      target_fdr = 0.05, norm = "median") {

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
    patterns <- data.frame(row.names = "Pattern1")
    for (i in conditions) {
      patterns[1, i] <- 1
    }
  }

  normalized <- ebseq_size_factors(data, norm)
  eb_output <- EBSeq::EBMultiTest(
                          Data = data, NgVector = ng_vector,
                          Conditions = conditions, AllParti = patterns,
                          sizeFactors = normalized, maxround = rounds)
  posteriors <- EBSeq::GetMultiPP(eb_output)
  fold_changes <- EBSeq::GetMultiFC(eb_output)

  pp_df <- as.data.frame(eb_output[["PPMat"]])
  ## Drop the pattern which looks for all the same
  interesting_patterns <- as.data.frame(patterns[-1, ])

  table_lst <- list()
  ## FIXME: Change this to use a portion of make_pairwise_contrasts()
  ## so that numerators and denominators do not get flipped.
  for (i in seq_len(ncol(fold_changes[["FCMat"]]))) {
    column <- colnames(fold_changes[["FCMat"]])[i]
    contrast <- gsub(pattern = "Over", replacement = "_vs_", x = column)
    numerator <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\1", x = contrast)
    denominator <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\2", x = contrast)
    chosen_pattern_idx <- interesting_patterns[[numerator]] ==
      interesting_patterns[[denominator]]
    ## This should provide the name of the pattern which, when I query the
    ## PPMatrix will provide the likelihood that the two are the same
    ## I can therefore take the 1-that to get a likelihood different.
    chosen_pattern <- rownames(interesting_patterns)[chosen_pattern_idx]

    table <- data.frame(row.names = rownames(fold_changes[["FCMat"]]))
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

#' The primary function used in my EBSeq implementation.
#'
#' Most of the time, my invocation of ebseq will fall into this function.
#'
#' @param pair_data Matrix containing the samples comprising two experimental
#'  factors of interest.
#' @param conditions Factor of conditions in the data.
#' @param numerator Which factor has the numerator in the data.
#' @param denominator Which factor has the denominator in the data.
#' @param ng_vector Passed to ebseq.
#' @param rounds Passed to ebseq.
#' @param target_fdr Passed to ebseq.
#' @param norm Normalization method of ebseq to apply.
#' @param force Force inappropriate data into ebseq?
#' @return EBSeq result table with some extra formatting.
#' @seealso [ebseq_pairwise()]
ebseq_two <- function(pair_data, conditions,
                      numerator = 2, denominator = 1,
                      ng_vector = NULL, rounds = 10,
                      target_fdr = 0.05, norm = "median",
                      force = FALSE) {
  normalized <- ebseq_size_factors(pair_data, norm = norm)
  mesg("Starting EBTest of ", numerator, " vs. ", denominator, ".")
  if (isTRUE(force)) {
    mesg("Forcing out NA values by putting in the mean of all data.")
    ## Put NA values (proteomics) to the mean of the existing values in the hopes
    ## they will not mess anything up too badly.
    na_idx <- is.na(pair_data)
    pair_data[na_idx] <- mean(pair_data, na.rm = TRUE)
  }
  eb_output <- sm(EBSeq::EBTest(
    Data = pair_data, NgVector = NULL, Conditions = conditions,
    sizeFactors = normalized, maxround = rounds))
  posteriors <- EBSeq::GetPPMat(eb_output)
  fold_changes <- EBSeq::PostFC(eb_output)
  eb_result <- EBSeq::GetDEResults(eb_output, FDR = target_fdr)
  table <- data.frame(row.names = names(fold_changes[["PostFC"]]))
  table[["ebseq_FC"]] <- fold_changes[["RealFC"]]
  table[["logFC"]] <- log2(table[["ebseq_FC"]])
  table[["ebseq_c1mean"]] <- as.numeric(eb_output[["C1Mean"]][[1]])
  table[["ebseq_c2mean"]] <- as.numeric(eb_output[["C2Mean"]][[1]])
  table[["ebseq_mean"]] <- as.numeric(eb_output[["MeanList"]][[1]])
  table[["ebseq_var"]] <- as.numeric(eb_output[["VarList"]][[1]])
  table[["ebseq_postfc"]] <- fold_changes[["PostFC"]]
  table <- merge(table, as.data.frame(eb_result[["PPMat"]]),
                 by = "row.names", all.x = TRUE)
  rownames(table) <- table[["Row.names"]]
  table <- table[, -1]
  ## This is incorrect I think, but being used as a placeholder until I figure out how to
  ## properly adjust a set prior probabilities.
  mesg("Copying ppee values as ajusted p-values until I figure out how to deal with them.")
  table[["ebseq_adjp"]] <- table[["PPEE"]]

  ## Finally, make sure the 'direction' matches my conception of numerator/denominator.
  eb_direction <- fold_changes[["Direction"]]
  eb_numerator <- gsub(pattern = "^(.*) Over (.*)$", replacement = "\\1", x = eb_direction)
  eb_denominator <- gsub(pattern = "^(.*) Over (.*)$", replacement = "\\2", x = eb_direction)
  if (! eb_numerator == denominator) {
    table[["ebseq_FC"]] <- 1 / table[["ebseq_FC"]]
    table[["logFC"]] <- -1 * table[["logFC"]]
    table[["ebseq_postfc"]] <- 1 / table[["ebseq_postfc"]]
  }
  return(table)
}

## EOF
