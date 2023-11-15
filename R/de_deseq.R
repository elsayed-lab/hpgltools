## de_deseq.r: Some methods to standardize the inputs/outputs from DESeq2.
## DESeq2 has lots of fun options, so many that it can be a bit overwhelming.
## This seeks to simplify these invocations and ensure that it will work under
## most likely use-case scenarios.

#' Bring together some of the likelihood ratio test analyses.
#'
#' This function hopes to wrap up some of the ideas/methods for LRT.
#'
#' @param expt Input expressionset
#' @param interactor_column Potentially interacting metadata
#' @param interest_column Essentially the condition in other analyses.
#' @param transform DESeq2 transformation applied (vst or rlog).
#' @param factors Other factors of interest
#' @param cutoff Significance cutoff
#' @param minc Minimum number of elements for a group
#' @param interaction Use an interaction model?
#' @seealso DOI:10.1186/s13059-014-0550-8
#' @export
deseq_lrt <- function(expt, interactor_column = "visitnumber",
                      interest_column = "clinicaloutcome", transform = "rlog",
                      factors = NULL, cutoff = 0.05, minc = 3, interaction = TRUE) {
  reduced_string <- glue("~ {interest_column}")
  full_string <- glue("~ {interactor_column} + {interest_column}")
  if (isTRUE(interaction)) {
    reduced_string <- glue("~ {interactor_column} + {interest_column}")
    full_string <- glue("{full_string} + \\
 {interactor_column}:{interest_column}")
  }
  full_model <- as.formula(full_string)
  reduced_model <- as.formula(reduced_string)
  col_data <- pData(expt)
  if (!is.null(factors)) {
    for (f in factors) {
      col_data[[f]] <- as.factor(col_data[[f]])
    }
  }
  if (is.null(col_data[[interactor_column]])) {
    stop("There is no ", interactor_column, " column in the experiment metadata.")
  }
  if (class(col_data[[interactor_column]]) != "factor") {
    warning("The ", interactor_column,
            " should probably be a factor, set it with the 'factors' arg.")
  }
  if (is.null(col_data[[interest_column]])) {
    stop("There is no ", interest_column, " column in the experiment metadata.")
  }
  if (class(col_data[[interest_column]]) != "factor") {
    warning("The ", interest_column,
            " should probably be a factor, set it with the 'factors' arg.")
  }
  mesg("The full model is: ", as.character(full_model), ".")
  mesg("The truncated model is: ", as.character(reduced_model), ".")

  ## Lets check the rank of this design before making the DESeq data structure.
  deseq_test <- pData(expt)
  data_full_model <- stats::model.matrix.default(full_model, data = pData(expt))
  data_reduced_model <-  stats::model.matrix.default(reduced_model, data = pData(expt))
  full_model_columns <- ncol(data_full_model)
  reduced_model_columns <- ncol(data_reduced_model)
  full_model_rank <- qr(data_full_model)[["rank"]]
  reduced_model_rank <- qr(data_reduced_model)[["rank"]]
  if (full_model_rank < full_model_columns) {
    message("Creating the data set will fail because the resulting model is too low rank.")
    message("Consider trying without the interaction.")
  }

  deseq_input <- DESeq2::DESeqDataSetFromMatrix(countData = exprs(expt),
                                                colData = col_data,
                                                design = full_model)
  deseq_lrt <- DESeq2::DESeq(deseq_input, test = "LRT", reduced = reduced_model)
  deseq_lrt_table <- DESeq2::results(deseq_lrt)

  ## Copy-pasting from:
  ## https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
  ## Subset the LRT results to return genes with padj < 0.05
  padj <- NULL ## R CMD check
  lrt_significant <- deseq_lrt_table %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    tibble::as_tibble() %>%
    filter(padj <= cutoff)
  if (nrow(lrt_significant) == 0) {
    warning("There are no significant differences given the ", cutoff, " adjusted p-value.")
    lrt_significant <- deseq_lrt_table %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "gene") %>%
      tibble::as_tibble()
    message("Returning the full LRT table just so that you have something to look at.")
    retlist <- list(
        "deseq_result" = deseq_lrt,
        "deseq_table" = lrt_significant)
    return(retlist)
  }

  rlog_matrix <- matrix()
  if (transform == "vst") {
    rlog_matrix <- DESeq2::vst(deseq_input)
  } else {
    rlog_matrix <- DESeq2::rlog(deseq_input)
  }
  clustering_amounts <- rlog_matrix[lrt_significant[["gene"]], ]
  cluster_data <- try(DEGreport::degPatterns(assay(clustering_amounts), metadata = col_data,
                                             time = interactor_column, col = interest_column,
                                             minc = minc))
  if ("try-error" %in% class(cluster_data)) {
    ## On my container image, I don't have DISPLAY causing this to fail.
    retlist <- list(
      "deseq_result" = deseq_lrt,
      "deseq_table" = deseq_lrt_table)
    class(retlist) <- "deseq_lrt_noclusters"
    return(retlist)
  }

  group_lst <- NULL
  if (! "try-error" %in% class(cluster_data)) {
    cluster_df <- cluster_data[["df"]]
    cluster_df[["cluster"]] <- as.factor(cluster_df[["cluster"]])
    group_lst <- list()
    for (c in levels(cluster_df[["cluster"]])) {
      group_idx <- cluster_df[["cluster"]] == c
      group_lst[[c]] <- cluster_df[group_idx, ]
    }
  }
  retlist <- list(
      "deseq_result" = deseq_lrt,
      "deseq_table" = deseq_lrt_table,
      "cluster_data" = cluster_data,
      "group_list" = group_lst,
      "favorite_genes" = cluster_data[["df"]])
  class(retlist) <- "deseq_lrt"
  return(retlist)
}

#' deseq_pairwise()  Because I can't be trusted to remember '2'.
#'
#' This calls deseq2_pairwise(...) because I am determined to forget typing deseq2.
#'
#' @param ... I like cats.
#' @return stuff deseq2_pairwise results.
#' @seealso [deseq2_pairwise()]
#' @export
deseq_pairwise <- function(...) {
  deseq2_pairwise(...)
}

#' Set up model matrices contrasts and do pairwise comparisons of all conditions using DESeq2.
#'
#' Invoking DESeq2 is confusing, this should help.
#'
#' Like the other _pairwise() functions, this attempts to perform all pairwise
#' contrasts in the provided data set.  The details are of course slightly
#' different when using DESeq2.  Thus, this uses the function
#' choose_binom_dataset() to try to ensure that the incoming data is appropriate
#' for DESeq2 (if one normalized the data, it will attempt to revert to raw
#' counts, for example). It continues on to extract the conditions and batches
#' in the data, choose an appropriate experimental model, and run the DESeq
#' analyses as described in the manual.  It defaults to using an experimental
#' batch factor, but will accept a string like 'sva' instead, in which case it
#' will use sva to estimate the surrogates, and append them to the experimental
#' design.  The deseq_method parameter may be used to apply different DESeq2
#' code paths as outlined in the manual.  If you want to play with non-standard
#' data, the force argument will round the data and shoe-horn it into DESeq2.
#'
#' @param input Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Is condition in the experimental model?
#' @param model_batch Is batch in the experimental model?
#' @param model_intercept Use an intercept model?
#' @param alt_model Provide an arbitrary model here.
#' @param extra_contrasts Provide extra contrasts here.
#' @param annot_df Include some annotation information in the results?
#' @param force Force deseq to accept data which likely violates its assumptions.
#' @param deseq_method The DESeq2 manual shows a few ways to invoke it, I make
#'  2 of them available here.
#' @param ... Triple dots!  Options are passed to arglist.
#' @return List including the following information:
#'  run = the return from calling DESeq()
#'  denominators = list of denominators in the contrasts
#'  numerators = list of the numerators in the contrasts
#'  conditions = the list of conditions in the experiment
#'  coefficients = list of coefficients making the contrasts
#'  all_tables = list of DE tables
#' @seealso [DESeq2] [basic_pairwise()] [limma_pairwise()] [edger_pairwise()] [ebseq_pairwise()]
#'  DOI:10.1186/s13059-014-0550-8.
#' @examples
#' \dontrun{
#'  pretend = deseq2_pairwise(data, conditions, batches)
#' }
#' @export
deseq2_pairwise <- function(input = NULL, conditions = NULL,
                            batches = NULL, model_cond = TRUE,
                            model_batch = TRUE, model_intercept = FALSE,
                            alt_model = NULL, extra_contrasts = NULL,
                            annot_df = NULL, force = FALSE, keepers = NULL,
                            deseq_method = "long", fittype = "parametric", ...) {
  arglist <- list(...)

  mesg("Starting DESeq2 pairwise comparisons.")
  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force = force)
  ## Now that I understand pData a bit more, I should probably remove the
  ## conditions/batches slots from my expt classes.
  design <- pData(input)
  conditions <- input_data[["conditions"]]
  batches <- input_data[["batches"]]
  data <- input_data[["data"]]
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(as.factor(conditions))

  ## Make a model matrix which will have one entry for
  ## each of the condition/batches
  summarized <- NULL
  ## Moving the size-factor estimation into this if(){} block in order to accomodate sva-ish
  ## batch estimation in the model
  deseq_sf <- NULL

  ## A caveat because this is a point of confusion
  ## choose_model() returns a few models, including intercept and non-intercept versions
  ## of the same things.  However, if model_batch is passed as something like 'sva', then
  ## it will gather surrogate estimates from sva and friends and return those estimates.
  model_choice <- choose_model(input, conditions, batches, model_batch = model_batch,
                               model_cond = model_cond, model_intercept = model_intercept,
                               alt_model = alt_model,
                               ...)
  model_data <- model_choice[["chosen_model"]]
  model_including <- model_choice[["including"]]
  model_string <- model_choice[["chosen_string"]]
  ## This is redundant with the definition of design above.
  column_data <- pData(input)
  if (class(model_choice[["model_batch"]])[1] == "matrix") {
    ## The SV matrix from sva/ruv/etc are put into the model batch slot of the return from choose_model.
    ## Use them here if appropriate
    model_batch <- model_choice[["model_batch"]]
    column_data <- cbind(column_data, model_batch)
  }
  ## choose_model should now take all of the following into account
  ## Therefore the following 8 or so lines should not be needed any longer.
  model_string <- NULL
  if (!is.null(alt_model)) {
    mesg("DESeq2 step 1/5: Using a user-supplied model.")
    model_string <- model_choice[["chosen_string"]]
    if (is.null(model_string)) {
      model_string <- model_choice[["int_string"]]
    }
    column_data[["condition"]] <- as.factor(column_data[["condition"]])
    column_data[["batch"]] <- as.factor(column_data[["batch"]])
    summarized <- import_deseq(data, column_data,
                               model_string, tximport = input[["tximport"]][["raw"]])
    dataset <- DESeq2::DESeqDataSet(se = summarized, design = as.formula(model_string))
  } else if (isTRUE(model_batch) && isTRUE(model_cond)) {
    mesg("DESeq2 step 1/5: Including batch and condition in the deseq model.")
    ## summarized <- DESeqDataSetFromMatrix(countData = data, colData = pData(input),
    ##                                     design=~ 0 + condition + batch)
    ## conditions and batch in this context is information taken from pData()
    model_string <- model_choice[["chosen_string"]]
    column_data[["condition"]] <- as.factor(column_data[["condition"]])
    column_data[["batch"]] <- as.factor(column_data[["batch"]])
    summarized <- import_deseq(data, column_data,
                               model_string, tximport = input[["tximport"]][["raw"]])
    dataset <- DESeq2::DESeqDataSet(se = summarized, design = as.formula(model_string))
  } else if (isTRUE(model_batch)) {
    mesg("DESeq2 step 1/5: Including only batch in the deseq model.")
    model_string <- model_choice[["chosen_string"]]
    column_data[["batch"]] <- as.factor(column_data[["batch"]])
    summarized <- import_deseq(data, column_data,
                               model_string, tximport = input[["tximport"]][["raw"]])
    dataset <- DESeq2::DESeqDataSet(se = summarized, design = as.formula(model_string))
  } else if (class(model_batch)[1] == "matrix") {
    mesg("DESeq2 step 1/5: Including a matrix of batch estimates in the deseq model.")
    sv_model_string <- model_choice[["chosen_string"]]
    column_data[["condition"]] <- as.factor(column_data[["condition"]])
    for (i in seq_along(ncol(data))) {
      data[[i]] <- as.integer(data[[i]])
    }
    summarized <- import_deseq(data, column_data,
                               sv_model_string, tximport = input[["tximport"]][["raw"]])

    dataset <- DESeq2::DESeqDataSet(se = summarized, design = as.formula(sv_model_string))
  } else {
    mesg("DESeq2 step 1/5: Including only condition in the deseq model.")
    model_string <- model_choice[["chosen_string"]]
    column_data[["condition"]] <- as.factor(column_data[["condition"]])
    summarized <- import_deseq(data, column_data,
                               model_string, tximport = input[["tximport"]][["raw"]])
    dataset <- DESeq2::DESeqDataSet(se = summarized, design = as.formula(model_string))
  }

  normalized_counts <- NULL ## When performing the 'long' analysis, I pull
  ## out the normalized counts so that we can compare against other analyses (e.g. with Julieth)
  deseq_run <- NULL
  chosen_beta <- model_intercept
  if (deseq_method == "short") {
    mesg("DESeq steps 2-4 in one shot.")
    deseq_run <- try(DESeq2::DESeq(dataset, fitType = fittype,
                                   betaPrior = chosen_beta), silent = TRUE)
    if (class(deseq_run)[1] == "try-error") {
      message("A fitType of 'parametric' failed for this data, trying 'mean'.")
      deseq_run <- try(DESeq2::DESeq(dataset, fitType = "mean"), silent = TRUE)
      if (class(deseq_run)[1] == "try-error") {
        message("Both 'parametric' and 'mean' failed.  Trying 'local'.")
        deseq_run <- try(DESeq2::DESeq(dataset, fitType = "local"), silent = TRUE)
        if (class(deseq_run)[1] == "try-error") {
          warning("All fitting types failed.  This will end badly.")
        } else {
          mesg("Using a local fit seems to have worked.")
        }
      } else {
        mesg("Using a mean fitting seems to have worked.")
      }
    }
  } else {
    ## Eg. Using the long method of invoking DESeq.
    ## If making a model ~0 + condition -- then must set betaPrior = FALSE
    mesg("DESeq2 step 2/5: Estimate size factors.")
    deseq_sf <- DESeq2::estimateSizeFactors(dataset)
    mesg("DESeq2 step 3/5: Estimate dispersions.")
    deseq_disp <- try(DESeq2::estimateDispersions(deseq_sf, fitType = fittype), silent = TRUE)
    if (class(deseq_disp)[1] == "try-error") {
      mesg("Trying a mean fitting.")
      deseq_disp <- try(DESeq2::estimateDispersions(deseq_sf, fitType = "mean"), silent = TRUE)
      if (class(deseq_disp)[1] == "try-error") {
        warning("Both 'parametric' and 'mean' failed.  Trying 'local'.")
        deseq_disp <- try(DESeq2::estimateDispersions(deseq_sf,
                                                      fitType = "local"), silent = TRUE)
        if (class(deseq_disp)[1] == "try-error") {
          warning("All fitting types failed.  This will end badly.")
        } else {
          mesg("Using a local fit seems to have worked.")
        }
      } else {
        mesg("Using a mean fitting seems to have worked.")
      }
    } else {
      mesg("Using a parametric fitting seems to have worked.")
    }
    mesg("DESeq2 step 4/5: nbinomWaldTest.")
    deseq_run <- DESeq2::nbinomWaldTest(deseq_disp, betaPrior = chosen_beta, quiet = TRUE)
  }
  normalized_counts <- DESeq2::counts(deseq_disp)
  dispersions <- sm(try(DESeq2::plotDispEsts(deseq_run), silent = TRUE))
  dispersion_plot <- NULL
  if (class(dispersions)[1] != "try-error") {
    dispersion_plot <- grDevices::recordPlot()
  }

  ## possible options:  betaPrior = TRUE, betaPriorVar, modelMatrix = NULL
  ## modelMatrixType, maxit = 100, useOptim = TRUE useT = FALSE df useQR = TRUE
  ## deseq_run = DESeq2::nbinomLRT(deseq_disp)
  ## Set contrast= for each pairwise comparison here!

  ## DESeq does not use contrasts in a way similar to limma/edgeR
  ## Therefore we will create all sets of c/d using these for loops.
  mesg("DESeq2 step 5/5: Collecting the results.")
  denominators <- list()
  numerators <- list()
  result_list <- list()
  coefficient_list <- list()
  ## The following is an attempted simplification of the contrast formulae
  end <- length(condition_levels) - 1
  number_comparisons <- sum(seq_len(end))
  ## Something peculiar has happened, since making the condition levels
  ## ordered in the expts, deseq no longer necessarily orders it contrasts
  ## the same as limma/edger.
  ## This change in ordering is quite annoying.
  ## Therefore, I will invoke make_pairwise_contrasts() here
  ## rather than make all the contrasts myself, then use that ordering
  ## to handle DESeq's contrast method.
  apc <- make_pairwise_contrasts(model_data, conditions,
                                 extra_contrasts = extra_contrasts,
                                 do_identities = FALSE, keepers = keepers,
                                 ...)
  contrast_order <- apc[["names"]]
  contrast_strings <- apc[["all_pairwise"]]
  ## These two character lists are relevant because of the possibility that I
  ## will ask for a series of extra contrasts; in addition, I use them to ask
  ## for specific tables.
  contrasts <- c()
  contrasts_full <- c()
  ##total_contrasts <- length(condition_levels)
  ##total_contrasts <- (total_contrasts * (total_contrasts + 1)) / 2
  total_contrasts <- length(contrast_order)
  if (isTRUE(verbose)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (i in seq_along(contrast_order)) {
    contrast_name <- contrast_order[[i]]
    contrast_string <- contrast_strings[[i]]
    if (isTRUE(verbose)) {
      pct_done <- i / length(contrast_order)
      utils::setTxtProgressBar(bar, pct_done)
    }
    num_den_string <- strsplit(x = contrast_name, split = "_vs_")[[1]]
    if (length(num_den_string) == 0) {
      warning("This contrast does not appear to have a _vs_ in it, it cannot be used.")
    }
    num_name <- num_den_string[1]
    den_name <- num_den_string[2]
    denominators[[contrast_string]] <- den_name
    numerators[[contrast_string]] <- num_name
    contrasts <- append(contrast_name, contrasts)
    ## I am pretty sure this is not needed, but I will hold on to it for the moment.
    contrasts_full <- append(contrast_string, contrasts_full)
    if (! glue("condition{num_name}") %in% DESeq2::resultsNames(deseq_run)) {
      message("The contrast ", num_name, " is not in the results.")
      message("If this is not an extra contrast, then this is an error.")
      next
    }
    result <- as.data.frame(DESeq2::results(object = deseq_run,
                                            contrast = c("condition", num_name, den_name),
                                            format = "DataFrame"))
    ##result <- DESeq2::results(object = deseq_run,
    ##                          contrast = c("condition", num_name, den_name))
    result <- result[order(result[["log2FoldChange"]]), ]
    colnames(result) <- c("baseMean", "logFC", "lfcSE", "stat", "P.Value", "adj.P.Val")
    ## From here on everything is the same.
    result[is.na(result[["P.Value"]]), "P.Value"] <- 1 ## Some p-values come out as NA
    result[is.na(result[["adj.P.Val"]]), "adj.P.Val"] <- 1 ## Some p-values come out as NA
    result[["baseMean"]] <- signif(x = as.numeric(result[["baseMean"]]), digits = 4)
    result[["logFC"]] <- signif(x = as.numeric(result[["logFC"]]), digits = 4)
    result[["lfcSE"]] <- signif(x = as.numeric(result[["lfcSE"]]), digits = 4)
    result[["stat"]] <- signif(x = as.numeric(result[["stat"]]), digits = 4)
    result[["P.Value"]] <- signif(x = as.numeric(result[["P.Value"]]), digits = 4)
    result[["adj.P.Val"]] <- signif(x = as.numeric(result[["adj.P.Val"]]), digits = 4)
    result_nas <- is.na(result)
    result[result_nas] <- 0
    if (!is.null(annot_df)) {
      result <- merge(result, annot_df, by.x = "row.names", by.y = "row.names")
    }
    result_list[[contrast_name]] <- result
  }
  if (isTRUE(verbose)) {
    close(bar)
  }
  ## The logic here is a little tortuous.
  ## Here are some sample column names from an arbitrary coef() call:
  ## "Intercept" "SV1" "SV2" "SV3" "condition_mtc_wtu_vs_mtc_mtu"
  ## First of all, we don't care about the 'condition' prefix.
  ## In addition, if we want the coefficient for mtc_wtu, then we need to subtract
  ## the mtc_wtu_vs_mtc_mtu from the Intercept, which is annoying.
  ## The following lines will attempt to do these things and
  ## appropriately rename the columns.
  coefficient_df <- as.data.frame(coef(deseq_run))
  ## Here I will just simplify the column names.
  colnames(coefficient_df) <- gsub(
      pattern = "^condition", replacement = "", x = colnames(coefficient_df))
  colnames(coefficient_df) <- gsub(
      pattern = "^batch", replacement = "", x = colnames(coefficient_df))
  colnames(coefficient_df) <- gsub(
      pattern = "^_", replacement = "", x = colnames(coefficient_df))
  remaining_list <- colnames(coefficient_df)

  ## Create a list of all the likely column names, depending on how deseq was
  ## called this might be numerator_vs_denominator or numerator denominator.
  ## So, I just make a list of them all.
  num_den <- unique(c(names(numerators), names(denominators)))
  ## AFAICT, the intercept is the second half of the contrasts listed.
  ## So grab that contrast name out
  if ("Intercept" %in% remaining_list) {
    ## When there is a bunch of x_vs_y, then the intercept will be set to the _y
    ## And all other columns will be subtracted from it to get their coefficients.
    vs_indexes <- grepl(pattern = "_vs_", x = colnames(coefficient_df))
    if (sum(vs_indexes) > 0) {
      intercept_pairing <- strsplit(x = colnames(coefficient_df)[vs_indexes], split = "_vs_")
      ## This gives a list like: [[1]][1]: 'numerator' [[1]][2]: 'denominator'
      ## So grab the second element of an arbitrary list element.
      intercept_name <- intercept_pairing[[1]][2]
      ## Now grab a list of every other column
      not_intercepts_idx <- ! grepl(pattern = intercept_name, x = unlist(intercept_pairing))
      not_intercepts <- unlist(intercept_pairing)[not_intercepts_idx]
      for (count in seq_len(ncol(coefficient_df))) {
        column_name <- colnames(coefficient_df)[count]
        if (count == 1) {
          colnames(coefficient_df)[1] <- intercept_name
          next
        } else if (! vs_indexes[count]) {
          ## Then this does not have _vs_ in it, so skip.
          next
        } else {
          numerator <- strsplit(x = column_name, split = "_vs_")[[1]][1]
          coefficient_df[[count]] <- abs(coefficient_df[[1]] - coefficient_df[[count]])
          colnames(coefficient_df)[count] <- numerator
        }
      }
      ## End if the columns have _vs_ in them.
    } else {
      ## In this case, we just want the name of the condition which is not in the set
      ## of columns of the coefficient df.
      ## This is a bit more verbose that strictly it needs to be, but I hope it
      ## is clearer therefore. 1st, if a numerator/denominator is missing, then
      ## it is the intercept name.
      columns <- colnames(coefficient_df)
      missing_name_idx <- ! num_den %in% columns
      missing_name <- num_den[missing_name_idx]
      ## Those indexes found in the numerator+denominator list will be subtracted
      ## Previously this line was 'containing_names_idx <- columns %in% num_den'
      ## I think that was wrong, it should be this:
      containing_names_idx <- colnames(coefficient_df) %in% num_den
      containing_names <- columns[containing_names_idx]
      ## If the are not in the numerator+denominator list, then they must be the
      ## SVs (except the first column of course, that is the intercept.
      extra_names_idx <- ! columns %in% num_den
      extra_names <- columns[extra_names_idx]
      for (count in seq_len(ncol(coefficient_df))) {
        column_name <- colnames(coefficient_df)[count]
        if (count == 1) {
          colnames(coefficient_df)[1] <- missing_name
          next
        } else if (column_name %in% extra_names) {
          ## The SV columns or batch or whatever
          next
        } else {
          coefficient_df[[count]] <- abs(coefficient_df[[1]] - coefficient_df[[count]])
        }
      }
    } ## End both likely types of intercept columns.
  }

  ## Let us add the coefficients of each contrast to the result tables from deseq.
  for (i in seq_along(contrast_order)) {
    contrast_name <- contrast_order[[i]]
    if (is.null(result_list[[contrast_name]])) {
      ## This contrast was not performed, skip it.
      next
    }
    num_den <- strsplit(x = contrast_name, split = "_vs_")
    numerator_name <- num_den[[1]][[1]]
    denominator_name <- num_den[[1]][[2]]
    ## Reminder to self: If you expect a dataframe, are using a matrix,
    ## and ask for a simplifying subset ([[]]), then you will be sad with
    ## the unhelpful error 'subscript out of bounds'.  Don't forget this.
    if (!is.null(coefficient_df[[numerator_name]]) &&
          !is.null(coefficient_df[[denominator_name]])) {
      this_coef <- coefficient_df[, c(numerator_name, denominator_name)]
      colnames(this_coef) <- c("deseq_num", "deseq_den")
      result_list[[contrast_name]] <- merge(result_list[[contrast_name]], this_coef,
                                            by.x = "row.names", by.y = "row.names")
      rownames(result_list[[contrast_name]]) <- result_list[[contrast_name]][["Row.names"]]
      result_list[[contrast_name]][["Row.names"]] <- NULL
    } else {
      result_list[[contrast_name]][["deseq_num"]] <- NA
      result_list[[contrast_name]][["deseq_den"]] <- NA
    }
  }

  retlist <- list(
      "all_tables" = result_list,
      "batches" = batches,
      "batches_table" = batches_table,
      "coefficients" = coefficient_df,
      "conditions" = conditions,
      "conditions_table" = conditions_table,
      "contrasts_full" = contrasts_full,
      "contrasts_performed" = contrasts,
      "denominators" = denominators,
      "dispersion_plot" = dispersion_plot,
      "input_data" = input,
      "method" = "deseq",
      "model" = model_data,
      "model_string" = model_string,
      "normalized_counts" = normalized_counts,
      "numerators" = numerators,
      "deseq_dataset" = dataset,
      "run" = deseq_run)
  class(retlist) <- c("deseq_result", "list")
  if (!is.null(arglist[["deseq_excel"]])) {
    retlist[["deseq_excel"]] <- write_deseq(retlist, excel = arglist[["deseq_excel"]])
  }
  class(retlist) <- c("deseq_pairwise", "list")
  return(retlist)
}

#' Given a set of surrogate variables from sva and friends, try adding them to a DESeqDataSet.
#'
#' Sometimes sva returns a set of surrogate variable estimates which lead to
#' models which are invalid according to DESeq2.  This function will try before
#' buying and tell the user if the sva model additions are valid according to
#' DESeq.
#'
#' @param data DESeqDataSet to test out.
#' @param summarized Existing DESeq metadata to append svs.
#' @param svs Surrogates from sva and friends to test out.
#' @param num_sv Optionally, provide the number of SVs, primarily used if
#'  recursing in the hunt for a valid number of surrogates.
#' @seealso [sva] [RUVSeq] [all_adjusters()] [normalize_batch()]
#' @return DESeqDataSet with at least some of the SVs appended to the model.
deseq_try_sv <- function(data, summarized, svs, num_sv = NULL) {
  counts <- DESeq2::counts(data)
  passed <- FALSE
  if (is.null(num_sv)) {
    num_sv <- ncol(svs)
  }
  formula_string <- "as.formula(~ "
  for (count in seq_len(num_sv)) {
    colname <- glue("SV{count}")
    summarized[[colname]] <- svs[, count]
    formula_string <- glue("{formula_string} {colname} + ")
  }
  formula_string <- glue("{formula_string}condition)")
  new_formula <- eval(parse(text = formula_string))
  new_summarized <- summarized
  DESeq2::design(new_summarized) <- new_formula
  data_model <- stats::model.matrix.default(DESeq2::design(summarized),
                                            data = as.data.frame(summarized@colData))
  model_columns <- ncol(data_model)
  model_rank <- qr(data_model)[["rank"]]
  if (model_rank < model_columns) {
    mesg("Including ", num_sv, " will fail because the resulting model is too low rank.")
    num_sv <- num_sv - 1
    message("Trying again with ", num_sv, " surrogates.")
    message("You should consider rerunning the pairwise comparison with the number of
surrogates explicitly stated with the option surrogates = number.")
    ret <- deseq_try_sv(data, summarized, svs, (num_sv - 1))
  } else {
    ## If we get here, then the number of surrogates should work with DESeq2.
    ## Perhaps I should re-calculate the variables with the specific number of variables.
    new_dataset <- DESeq2::DESeqDataSet(se = new_summarized, design = new_formula)
    return(new_dataset)
  }
  return(ret)
}

#' Try to add data to DESeq in a flexible fashion.  This currently only handles
#' matrices, htseq data, and tximport data.
#'
#' This will hopefully make adding counts to a DESeq data set easier, as it
#' tries to handle the various arguments with minimal fuss.
#'
#' @param data Counts from htseq/mtrx/tximport/etc
#' @param column_data I think this is the sample names, I forget.
#' @param model_string Model describing the data by sample names.
#' @param tximport Where is this data coming from?
#' @seealso [DESeq2::DESeqDataSetFromMatrix]
import_deseq <- function(data, column_data, model_string,
                         tximport = NULL) {
  summarized <- NULL
  ## column_data_na_idx <- is.na(column_data)
  ## column_data[column_data_na_idx] <- "undefined"
  ## The default.

  ## DESeq explicitly limits the input of the data to the range of 2^32 integers.
  ## If one is insane and wants to dump intensity data into deseq, this might get violated.
  integer_limit <- .Machine[["integer.max"]]
  na_idx <- is.na(data)
  data[na_idx] <- 0
  too_big_idx <- data > integer_limit
  if (sum(too_big_idx) > 0) {
    warning("Converted down ", sum(too_big_idx),
            " elements because they are larger than the maximum integer size.")
    data[too_big_idx] <- integer_limit
  }

  if (is.null(tximport)) {
    summarized <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                                 colData = column_data,
                                                 design = as.formula(model_string))
  } else if (tximport[1] == "htseq") {
    ## We are not likely to use this.
    summarized <- DESeq2::DESeqDataSetFromHTSeqCount(countData = data,
                                                     colData = column_data,
                                                     design = as.formula(model_string))
  } else {
    ## This may be insufficient, it may require the full tximport result, while
    ## this may just be that result$counts, so be aware!!

    ## First make sure that if we subsetted the data, that is maintained from
    ## the data to the tximportted data
    keepers <- rownames(data)
    ## These recasts were because some matrix subsets were failing.
    abundances <- as.data.frame(tximport[["abundance"]])
    abundances <- abundances[keepers, ]
    counts <- as.data.frame(tximport[["counts"]])
    counts <- counts[keepers, ]
    lengths <- as.data.frame(tximport[["length"]])
    lengths <- lengths[keepers, ]
    ## Something about this does not feel right, I need to look over this some more.
    na_rows <- grepl(pattern = "^NA", x = rownames(counts))
    counts <- counts[!na_rows, ]
    lengths <- lengths[!na_rows, ]
    abundances <- abundances[!na_rows, ]
    ## That should no longer be a problem now that I explicitly change the
    ## rownames of all the tximport data to match the rownames from
    ## as.data.table() but I am leaving it for now.
    tximport[["abundance"]] <- as.matrix(abundances)
    tximport[["counts"]] <- as.matrix(counts)
    tximport[["length"]] <- as.matrix(lengths)
    summarized <- DESeq2::DESeqDataSetFromTximport(txi = tximport,
                                                   colData = column_data,
                                                   design = as.formula(model_string))
  }
  return(summarized)
}

#' Writes out the results of a deseq search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from deseq and friends.
#'
#' Tested in test_24deseq.R
#'
#' @param data Output from deseq_pairwise()
#' @param ... Options for writing the xlsx file.
#' @seealso [write_de_table()]
#' @examples
#' \dontrun{
#'  finished_comparison <- deseq2_pairwise(expressionset)
#'  data_list <- write_deseq(finished_comparison)
#' }
#' @export
write_deseq <- function(data, ...) {
  result <- write_de_table(data, type = "deseq", ...)
  return(result)
}

## EOF
