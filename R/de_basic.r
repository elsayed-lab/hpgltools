#' The simplest possible differential expression method.
#'
#' Perform a pairwise comparison among conditions which takes
#' nothing into account.  It _only_ takes the conditions, a mean value/variance among
#' them, divides by condition, and returns the result.  No fancy nomalizations, no
#' statistical models, no nothing.  It should be the very worst method possible.
#' But, it should also provide a baseline to compare the other tools against, they should
#' all do better than this, always.
#'
#' Tested in test_27de_basic.R
#' This function was written after the corresponding functions in de_deseq.R, de_edger.R,
#' and de_limma.R.  Like those, it performs the full set of pairwise comparisons and returns
#' a list of the results.  As mentioned above, unlike those, it is purposefully stupid.
#'
#' @param input Count table by sample.
#' @param design Data frame of samples and conditions.
#' @param conditions  Not currently used, but passed from all_pairwise()
#' @param batches  Not currently used, but passed from all_pairwise()
#' @param model_cond  Not currently used, but passed from all_pairwise()
#' @param model_intercept Not currently used, but passed from all_pairwise()
#' @param alt_model Not currently used, but passed from all_pairwise()
#' @param model_batch Not currently used, but passed from all_pairwise()
#' @param force Force as input non-normalized data?
#' @param ... Extra options passed to arglist.
#' @return Df of pseudo-logFC, p-values, numerators, and denominators.
#' @seealso \pkg{limma} \pkg{DESeq2} \pkg{edgeR}
#' @examples
#' \dontrun{
#' stupid_de <- basic_pairwise(expt)
#' }
#' @export
basic_pairwise <- function(input=NULL, design=NULL, conditions=NULL, batches=NULL, model_cond=TRUE,
                           model_intercept=FALSE, alt_model=NULL, model_batch=FALSE,
                           force=FALSE, ...) {
  arglist <- list(...)
  if (!is.null(arglist[["input"]])) {
    input <- arglist[["input"]]
  }
  if (!is.null(arglist[["design"]])) {
    conditions <- arglist[["design"]]
  }
  if (!is.null(arglist[["force"]])) {
    batches <- arglist[["force"]]
  }
  message("Starting basic pairwise comparison.")
  input <- sanitize_expt(input)
  input_data <- choose_basic_dataset(input, force=force)
  design <- pData(input)
  conditions <- input_data[["conditions"]]
  batches <- input_data[["batches"]]
  data <- input_data[["data"]]

  conditions <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=conditions)
  batches <- gsub(pattern="^(\\d+)$", replacement="b\\1", x=batches)
  types <- levels(as.factor(conditions))
  num_conds <- length(types)
  ## These will be filled with num_conds columns and numRows(input) rows.
  median_table <- data.frame()
  variance_table <- data.frame()
  ## First use conditions to rbind a table of medians by condition.
  message("Basic step 1/3: Creating median and variance tables.")
  median_colnames <- c()
  for (c in 1:num_conds) {
    condition_name <- types[c]
    median_colnames <- append(median_colnames, condition_name)
    columns <- which(conditions == condition_name)
    if (length(columns) == 1) {
      med <- data.frame(data[, columns], stringsAsFactors=FALSE)
      var <- as.data.frame(matrix(NA, ncol=1, nrow=nrow(med)), stringsAsFactors=FALSE)
    } else {
      med_input <- data[, columns]
      med <- data.frame(Biobase::rowMedians(as.matrix(med_input)))
      colnames(med) <- c(condition_name)
      var <- as.data.frame(genefilter::rowVars(as.matrix(med_input)))
      colnames(var) <- c(condition_name)
    }
    if (c == 1) {
      median_table <- med
      variance_table <- var
    } else {
      median_table <- cbind(median_table, med)
      variance_table <- cbind(variance_table, var)
    }
  } ## end creation of median and variance tables.
  colnames(median_table) <- median_colnames
  colnames(variance_table) <- median_colnames
  rownames(median_table) <- rownames(data)
  rownames(variance_table) <- rownames(data)
  ## We have tables of the median values by condition
  ## Now perform the pairwise comparisons
  comparisons <- data.frame()
  tvalues <- data.frame()
  pvalues <- data.frame()
  num_done <- 0
  column_list <- c()
  total_contrasts <- length(levels(as.factor(conditions)))
  total_contrasts <- (total_contrasts * (total_contrasts + 1)) / 2
  message("Basic step 2/3: Performing ", total_contrasts, " comparisons.")

  model_choice <- sm(choose_model(
    input, conditions=conditions, batches=batches, model_batch=FALSE,
    model_cond=TRUE, model_intercept=FALSE, alt_model=NULL,
    ...))
  model_data <- model_choice[["chosen_model"]]
  ## basic_pairwise() does not support extra contrasts, but they may be passed through via ...
  apc <- make_pairwise_contrasts(model_data, conditions, do_identities=FALSE,
                                 ...)
  contrasts_performed <- c()
  bar <- utils::txtProgressBar(style=3)
  for (c in 1:length(apc[["names"]])) {
    num_done <- num_done + 1
    pct_done <- c / length(apc[["names"]])
    name  <- apc[["names"]][[c]]
    c_name <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\1", x=name)
    d_name <- gsub(pattern="^(.*)_vs_(.*)$", replacement="\\2", x=name)
    utils::setTxtProgressBar(bar, pct_done)
    contrasts_performed <- append(name, contrasts_performed)
    if (! c_name %in% colnames(median_table)) {
      message("The contrast ", name, " is not in the results.")
      message("If this is not an extra contrast, then this is an error.")
      next
    }
    division <- data.frame(
      median_table[, c_name] - median_table[, d_name])
    column_list <- append(column_list, name)
    colnames(division) <- name
    ## Lets see if I can make a dirty p-value
    xcols <- which(conditions == c_name)
    ycols <- which(conditions == d_name)
    xdata <- as.data.frame(data[, xcols])
    ydata <- as.data.frame(data[, ycols])

    t_data <- vector("list", nrow(xdata))
    p_data <- vector("list", nrow(xdata))
    for (j in 1:nrow(xdata)) {
      test_result <- try(t.test(xdata[j, ], ydata[j, ]), silent=TRUE)
      if (class(test_result) == "htest") {
        t_data[[j]] <- test_result[[1]]
        p_data[[j]] <- test_result[[3]]
      } else {
        t_data[[j]] <- 0
        p_data[[j]] <- 1
      }
    } ## Done calculating cheapo p-values

    if (c == 1) {
      comparisons <- division
      tvalues <- t_data
      pvalues <- p_data
    } else {
      comparisons <- cbind(comparisons, division)
      tvalues <- cbind(tvalues, t_data)
      pvalues <- cbind(pvalues, p_data)
    }
  } ## End for each contrast
  close(bar)

  ## Because of the way I made tvalues/pvalues into a list
  ## If only 1 comparison was performed, the resulting data structure never gets coerced into a
  ## data frame.  Therefore I am performing this check which, if a single comparison was done, adds
  ## a second column, performs the coercion, then strips it away.  This is a stupid way
  ## of doing what I want.
  if (num_done == 1) {
    tvalues <- cbind(tvalues, t_data)
    pvalues <- cbind(pvalues, p_data)
    tvalues <- as.data.frame(tvalues)
    pvalues <- as.data.frame(pvalues)
    tvalues <- tvalues[-1]
    pvalues <- pvalues[-1]
  }
  comparisons[is.na(comparisons)] <- 0
  tvalues[is.na(tvalues)] <- 0
  pvalues[is.na(pvalues)] <- 1
  rownames(comparisons) <- rownames(data)
  rownames(tvalues) <- rownames(data)
  rownames(pvalues) <- rownames(data)
  all_tables <- list()

  message("Basic step 3/3: Creating faux DE Tables.")
  for (e in 1:length(colnames(comparisons))) {
    colname <- colnames(comparisons)[[e]]
    fc_column <- comparisons[, e]
    t_column <- as.numeric(tvalues[, e])
    p_column <- as.numeric(pvalues[, e])
    fc_column[mapply(is.infinite, fc_column)] <- 0
    numer_denom <- strsplit(x=colname, split="_vs_")[[1]]
    numerator <- numer_denom[1]
    denominator <- numer_denom[2]
    fc_table <- data.frame(
      "numerator_median" = median_table[[numerator]],
      "denominator_median" = median_table[[denominator]],
      "numerator_var" = variance_table[[numerator]],
      "denominator_var" = variance_table[[denominator]],
      "t" = t_column,
      "p" = p_column,
      "logFC" = fc_column)
    fc_table[["adjp"]] <- stats::p.adjust(as.numeric(fc_table[["p"]]), method="BH")

    fc_table[["numerator_median"]] <- signif(x=fc_table[["numerator_median"]], digits=4)
    fc_table[["denominator_median"]] <- signif(x=fc_table[["denominator_median"]], digits=4)
    fc_table[["numerator_var"]] <- format(x=fc_table[["numerator_var"]], digits=4, scientific=TRUE)
    fc_table[["denominator_var"]] <- format(x=fc_table[["denominator_var"]], digits=4, scientific=TRUE)
    fc_table[["t"]] <- signif(x=fc_table[["t"]], digits=4)
    fc_table[["p"]] <- format(x=fc_table[["p"]], digits=4, scientific=TRUE)
    fc_table[["adjp"]] <- format(x=fc_table[["adjp"]], digits=4, scientific=TRUE)
    fc_table[["logFC"]] <- signif(x=fc_table[["logFC"]], digits=4)
    rownames(fc_table) <- rownames(data)
    all_tables[[e]] <- fc_table
  }
  message("Basic: Returning tables.")
  names(all_tables) <- colnames(comparisons)

  retlist <- list(
    "all_pairwise" = comparisons,
    "all_tables" = all_tables,
    "conditions_table" = table(conditions),
    "conditions" = conditions,
    "contrasts_performed" = contrasts_performed,
    "input_data" = data,
    "medians" = median_table,
    "method" = "basic",
    "variances" = variance_table)
  if (!is.null(arglist[["basic_excel"]])) {
    retlist[["basic_excel"]] <- write_basic(retlist, excel=arglist[["basic_excel"]])
  }
  return(retlist)
}

#' Attempt to ensure that input data to basic_pairwise() is suitable.
#'
#' basic_pairwise() assumes log2 data as input, use this to ensure that is true.
#'
#' @param input  An expressionset containing expt to test and/or modify.
#' @param force  If we want to try out other distributed data sets, force it in using me.
#' @param ... future options, I think currently unused.
#' @return data ready for basic_pairwise()
#' @seealso \pkg{Biobase}
#' @examples
#' \dontrun{
#'  ready <- choose_basic_dataset(expt)
#' }
choose_basic_dataset <- function(input, force=FALSE, ...) {
  ## arglist <- list(...)
  warn_user <- 0
  conditions <- input[["conditions"]]
  batches <- input[["batches"]]
  data <- as.data.frame(exprs(input))
  tran_state <- input[["state"]][["transform"]]
  if (is.null(tran_state)) {
    tran_state <- "raw"
  }
  conv_state <- input[["state"]][["conversion"]]
  ## Note that voom takes care of this for us.
  if (is.null(conv_state)) {
    conv_state <- "raw"
  }
  norm_state <- input[["state"]][["normalization"]]
  if (is.null(norm_state)) {
    norm_state <- "raw"
  }
  filt_state <- input[["state"]][["filter"]]
  if (is.null(filt_state)) {
    filt_state <- "raw"
  }

  ready <- input
  if (isTRUE(force)) {
    message("Leaving the data alone, regardless of normalization state.")
  } else {
    if (filt_state == "raw") {
      message("Basic step 0/3: Filtering data.")
      ready <- sm(normalize_expt(ready, filter=TRUE))
    }
    if (norm_state == "raw") {
      message("Basic step 0/3: Normalizing data.")
      ready <- sm(normalize_expt(ready, norm="quant"))
    }
    if (conv_state == "raw") {
      message("Basic step 0/3: Converting data.")
      ready <- sm(normalize_expt(ready, convert="cbcbcpm"))
    }

  }
  ## No matter what we do, it must be logged.
  if (tran_state == "raw") {
    message("Basic step 0/3: Transforming data.")
    ready <- sm(normalize_expt(ready, transform="log2"))
  }
  data <- as.data.frame(exprs(ready))
  rm(ready)
  retlist <- list(
    "conditions" = conditions,
    "batches" = batches,
    "data" = data)
  return(retlist)
}

#' Writes out the results of a basic search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from basic and friends.
#'
#' Tested in test_26basic.R
#'
#' @param data  Output from basic_pairwise()
#' @param ...  Options for writing the xlsx file.
#' @seealso \code{\link{write_de_table}}
#' @examples
#' \dontrun{
#'  finished_comparison <- basic_pairwise(expressionset)
#'  data_list <- write_basic(finished_comparison)
#' }
#' @export
write_basic <- function(data, ...) {
  result <- write_de_table(data, type="basic", ...)
  return(result)
}

## EOF
