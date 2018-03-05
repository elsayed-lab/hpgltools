#' Set up a model matrix and set of contrasts to do pairwise comparisons using EdgeR.
#'
#' This function performs the set of possible pairwise comparisons using EdgeR.
#'
#' Tested in test_26de_edger.R
#' Like the other _pairwise() functions, this attempts to perform all pairwise contrasts in the
#' provided data set.  The details are of course slightly different when using EdgeR.  Thus, this
#' uses the function choose_binom_dataset() to try to ensure that the incoming data is appropriate
#' for EdgeR (if one normalized the data, it will attempt to revert to raw counts, for example).
#' It continues on to extract the conditions and batches in the data, choose an appropriate
#' experimental model, and run the EdgeR analyses as described in the manual.  It defaults to using
#' an experimental batch factor, but will accept a string like 'sva' instead, in which case it will
#' use sva to estimate the surrogates, and append them to the experimental design.  The edger_method
#' parameter may be used to apply different EdgeR code paths as outlined in the manual.  If you
#' want to play with non-standard data, the force argument will round the data and shoe-horn it into
#' EdgeR.
#'
#' @param input Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the experimental model?
#' @param model_batch Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept Use an intercept containing model?
#' @param alt_model Alternate experimental model to use?
#' @param extra_contrasts Add some extra contrasts to add to the list of pairwise contrasts.
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df Annotation information to the data tables?
#' @param force Force edgeR to accept inputs which it should not have to deal with.
#' @param edger_method  I found a couple/few ways of doing edger in the manual, choose with this.
#' @param ... The elipsis parameter is fed to write_edger() at the end.
#' @return List including the following information:
#'  contrasts = The string representation of the contrasts performed.
#'  lrt = A list of the results from calling glmLRT(), one for each contrast.
#'  contrast_list = The list of each call to makeContrasts()
#'  I do this to avoid running into the limit on # of contrasts addressable by topTags()
#'  all_tables = a list of tables for the contrasts performed.
#' @seealso \pkg{edgeR}
#' @examples
#' \dontrun{
#'  pretend = edger_pairwise(data, conditions, batches)
#' }
#' @export
edger_pairwise <- function(input=NULL, conditions=NULL,
                           batches=NULL, model_cond=TRUE,
                           model_batch=TRUE, model_intercept=FALSE,
                           alt_model=NULL, extra_contrasts=NULL,
                           annot_df=NULL, force=FALSE,
                           edger_method="long", ...) {
  arglist <- list(...)

  edger_test <- "lrt"
  if (!is.null(arglist[["edger_test"]])) {
    edger_test <- arglist[["edger_test"]]
  }
  message("Starting edgeR pairwise comparisons.")
  input <- sanitize_expt(input)
  input_data <- choose_binom_dataset(input, force=force)
  design <- pData(input)
  conditions <- input_data[["conditions"]]
  conditions_table <- table(conditions)
  batches <- input_data[["batches"]]
  batches_table <- table(batches)
  data <- input_data[["data"]]
  conditions <- as.factor(conditions)
  batches <- as.factor(batches)

  model_choice <- choose_model(input, conditions=conditions,
                               batches=batches,
                               model_batch=model_batch,
                               model_cond=model_cond,
                               model_intercept=model_intercept,
                               alt_model=alt_model, ...)
  ##model_choice <- choose_model(input, conditions, batches,
  ##                             model_batch=model_batch,
  ##                             model_cond=model_cond,
  ##                             model_intercept=model_intercept,
  ##                             alt_model=alt_model)
  model_including <- model_choice[["including"]]
  if (class(model_choice[["model_batch"]]) == "matrix") {
    model_batch <- model_choice[["model_batch"]]
  }
  model_data <- model_choice[["chosen_model"]]
  model_string <- model_choice[["chosen_string"]]

  ## I have a strong sense that the most recent version of edgeR changed its dispersion estimate code
  ## Here is a note from the user's guide, which may have been there previously and I merely did not notice:
  ## To estimate common dispersion, trended dispersions and tagwise dispersions in one run
  ## y <- estimateDisp(y, design)
  ## raw <- edgeR::DGEList(counts=data, group=conditions)
  ## norm <- edgeR::calcNormFactors(raw)
  norm <- import_edger(data, conditions, tximport=input[["tximport"]][["raw"]])
  message("EdgeR step 1/9: importing and normalizing data.")
  final_norm <- NULL
  if (edger_method == "short") {
    message("EdgeR steps 2 through 6/9: All in one!")
    final_norm <- edgeR::estimateDisp(norm, design=model_data)
  } else {
    state <- TRUE
    message("EdgeR step 2/9: Estimating the common dispersion.")
    disp_norm <- try(edgeR::estimateCommonDisp(norm))
    if (class(disp_norm) == "try-error") {
      warning("estimateCommonDisp() failed.  Trying again with estimateDisp().")
      state <- FALSE
    }
    if (isTRUE(state)) {
      message("EdgeR step 3/9: Estimating dispersion across genes.")
      tagdisp_norm <- try(edgeR::estimateTagwiseDisp(disp_norm))
      if (class(tagdisp_norm) == "try-error") {
        warning("estimateTagwiseDisp() failed.  Trying again with estimateDisp().")
        state <- FALSE
      }
    }
    if (isTRUE(state)) {
      message("EdgeR step 4/9: Estimating GLM Common dispersion.")
      glm_norm <- try(edgeR::estimateGLMCommonDisp(tagdisp_norm, model_data))
      if (class(glm_norm) == "try-error") {
        warning("estimateGLMCommonDisp() failed.  Trying again with estimateDisp().")
        state <- FALSE
      }
    }
    if (isTRUE(state)) {
      message("EdgeR step 5/9: Estimating GLM Trended dispersion.")
      glm_trended <- try(edgeR::estimateGLMTrendedDisp(glm_norm, model_data))
      if (class(glm_trended) == "try-error") {
        warning("estimateGLMTrendedDisp() failed.  Trying again with estimateDisp().")
        state <- FALSE
      }
    }
    if (isTRUE(state)) {
      message("EdgeR step 6/9: Estimating GLM Tagged dispersion.")
      final_norm <- try(edgeR::estimateGLMTagwiseDisp(glm_trended, model_data))
      if (class(final_norm) == "try-error") {
        warning("estimateGLMTagwiseDisp() failed.  Trying again with estimateDisp().")
        state <- FALSE
      }
    }

    ## If we had a failure along the way, redo using estimateDisp()
    if (!isTRUE(state)) {
      warning("There was a failure when doing the estimations.")
      message("There was a failure when doing the estimations, using estimateDisp().")
      final_norm <- edgeR::estimateDisp(norm, design=model_data, robust=TRUE)
    }
  }
  cond_fit <- NULL
  if (edger_test == "lrt") {
    message("EdgeR step 7/9: Running glmFit, switch to glmQLFit by changing the argument 'edger_test'.")
    cond_fit <- edgeR::glmFit(final_norm, design=model_data, robust=TRUE)
  } else {
    message("EdgeR step 7/9: Running glmQLFit, switch to glmFit by changing the argument 'edger_test'.")
    cond_fit <- edgeR::glmQLFit(final_norm, design=model_data, robust=TRUE)
  }

  message("EdgeR step 8/9: Making pairwise contrasts.")
  apc <- make_pairwise_contrasts(model_data, conditions,
                                 extra_contrasts=extra_contrasts,
                                 do_identities=FALSE, ...)
  contrast_string <- apc[["contrast_string"]]

  ## This section is convoluted because glmLRT only seems to take up to 7 contrasts at a time.
  ## As a result, I iterate through the set of possible contrasts one at a time and ask for each
  ## separately.
  contrast_list <- list()
  result_list <- list()
  lrt_list <- list()
  sc <- vector("list", length(apc[["names"]]))
  end <- length(apc[["names"]])
  bar <- utils::txtProgressBar(style=3)
  for (con in 1:length(apc[["names"]])) {
    name <- apc[["names"]][[con]]
    pct_done <- con / length(apc[["names"]])
    utils::setTxtProgressBar(bar, pct_done)
    sc[[name]] <- gsub(pattern=",", replacement="", apc[["all_pairwise"]][[con]])
    tt <- parse(text=sc[[name]])
    ctr_string <- paste0("tt = mymakeContrasts(", tt, ", levels=model_data)")
    eval(parse(text=ctr_string))
    contrast_list[[name]] <- tt
    lrt_list[[name]] <- NULL
    tt <- sm(requireNamespace("edgeR"))
    tt <- sm(try(attachNamespace("edgeR"), silent=TRUE))
    if (edger_test == "lrt") {
      lrt_list[[name]] <- edgeR::glmLRT(cond_fit, contrast=contrast_list[[name]])
    } else {
      lrt_list[[name]] <- edgeR::glmQLFTest(cond_fit, contrast=contrast_list[[name]])
    }
    res <- edgeR::topTags(lrt_list[[name]], n=nrow(data), sort.by="logFC")
    res <- as.data.frame(res)
    res[["logFC"]] <- signif(x=as.numeric(res[["logFC"]]), digits=4)
    res[["logCPM"]] <- signif(x=as.numeric(res[["logCPM"]]), digits=4)
    if (!is.null(res[["LR"]])) {
      res[["LR"]] <- signif(x=as.numeric(res[["LR"]]), digits=4)
    } else if (!is.null(res[["F"]])) {
      res[["F"]] <- signif(x=as.numeric(res[["F"]]), digits=4)
    }
    res[["PValue"]] <- signif(x=as.numeric(res[["PValue"]]), digits=4)
    res[["FDR"]] <- signif(x=as.numeric(res[["FDR"]]), digits=4)
    result_list[[name]] <- res
  } ## End for loop
  close(bar)

  dispersions <- sm(try(edgeR::plotBCV(y=final_norm), silent=TRUE))
  dispersion_plot <- NULL
  if (class(dispersions)[[1]] != "try-error") {
    dispersion_plot <- grDevices::recordPlot()
  }

  retlist <- list(
    "all_tables" = result_list,
    "batches" = batches,
    "batches_table" = batches_table,
    "conditions" = conditions,
    "conditions_table" = conditions_table,
    "contrast_list" = contrast_list,
    "contrasts" = apc,
    "contrast_string" = contrast_string,
    "contrasts_performed" = apc[["names"]],
    "dispersion_plot" = dispersion_plot,
    "input_data" = input,
    "lrt" = lrt_list,
    "method" = "edger",
    "model" = model_data,
    "model_string" = model_string)
  if (!is.null(arglist[["edger_excel"]])) {
    retlist[["edger_excel"]] <- write_edger(retlist, excel=arglist[["edger_excel"]])
  }
  return(retlist)
}

## Taken from the tximport manual with minor modification.
import_edger <- function(data, conditions, tximport=NULL) {
  if (is.null(tximport)) {
    raw <- edgeR::DGEList(counts=data, group=conditions)
    norm <- edgeR::calcNormFactors(raw)
  } else {
    raw <- tximport[["counts"]]
    norm_mat <- tximport[["length"]]
    norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
    offset <- log(edgeR::calcNormFactors(raw / norm_mat)) + log(colSums(raw / norm_mat))
    norm <- edgeR::DGEList(raw)
    norm$offset <- t(t(log(norm_mat)) + offset)
  }
  return(norm)
}

#' Writes out the results of a edger search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from edger and friends.
#'
#' Tested in test_26edger.R
#'
#' @param data  Output from deseq_pairwise()
#' @param ...  Options for writing the xlsx file.
#' @seealso \pkg{limma}
#'  \code{\link[limma]{toptable}} \code{\link{write_xls}}
#' @examples
#' \dontrun{
#'  finished_comparison <- edger_pairwise(expressionset)
#'  data_list <- write_edger(finished_comparison)
#' }
#' @export
write_edger <- function(data, ...) {
  result <- write_de_table(data, type="edger", ...)
  return(result)
}

## EOF
