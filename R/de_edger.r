#' Set up a model matrix and set of contrasts to do pairwise comparisons using EdgeR.
#'
#' This function performs the set of possible pairwise comparisons using EdgeR.
#'
#' @param input Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the experimental model?
#' @param model_batch Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept Use cell means or intercept?
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
#'   contrasts = The string representation of the contrasts performed.
#'   lrt = A list of the results from calling glmLRT(), one for each contrast.
#'   contrast_list = The list of each call to makeContrasts()
#'   I do this to avoid running into the limit on # of contrasts addressable by topTags()
#'   all_tables = a list of tables for the contrasts performed.
#' @seealso \pkg{edgeR} \code{\link[edgeR]{topTags}} \code{\link[edgeR]{glmLRT}}
#'   \code{\link{make_pairwise_contrasts}} \code{\link[edgeR]{DGEList}}
#'   \code{\link[edgeR]{calcNormFactors}} \code{\link[edgeR]{estimateTagwiseDisp}}
#'   \code{\link[edgeR]{estimateCommonDisp}} \code{\link[edgeR]{estimateGLMCommonDisp}}
#'   \code{\link[edgeR]{estimateGLMTrendedDisp}} \code{\link[edgeR]{glmFit}}
#' @examples
#' \dontrun{
#'  pretend = edger_pairwise(data, conditions, batches)
#' }
#' @export
edger_pairwise <- function(input=NULL, conditions=NULL,
                           batches=NULL, model_cond=TRUE,
                           model_batch=TRUE, model_intercept=TRUE,
                           alt_model=NULL, extra_contrasts=NULL,
                           annot_df=NULL, force=FALSE,
                           edger_method="long", ...) {
    arglist <- list(...)
    if (!is.null(arglist[["input"]])) {
        input <- arglist[["input"]]
    }
    if (!is.null(arglist[["conditions"]])) {
        conditions <- arglist[["conditions"]]
    }
    if (!is.null(arglist[["batches"]])) {
        batches <- arglist[["batches"]]
    }
    if (!is.null(arglist[["model_cond"]])) {
        model_cond <- arglist[["model_cond"]]
    }
    if (!is.null(arglist[["model_batch"]])) {
        model_batch <- arglist[["model_batch"]]
    }
    if (!is.null(arglist[["model_intercept"]])) {
        model_intercept <- arglist[["model_intercept"]]
    }
    if (!is.null(arglist[["alt_model"]])) {
        alt_model <- arglist[["alt_model"]]
    }
    if (!is.null(arglist[["extra_contrasts"]])) {
        extra_contrasts <- arglist[["extra_contrasts"]]
    }
    if (!is.null(arglist[["annot_df"]])) {
        annot_df <- arglist[["annot_df"]]
    }
    if (!is.null(arglist[["force"]])) {
        force <- arglist[["force"]]
    }
    if (!is.null(arglist[["edger_method"]])) {
        edger_method <- arglist[["edger_method"]]
    }
    test_type <- "lrt"
    if (!is.null(arglist[["test_type"]])) {
        test_type <- arglist[["test_type"]]
    }
    message("Starting edgeR pairwise comparisons.")
    input_data <- choose_binom_dataset(input, force=force)
    design <- Biobase::pData(input[["expressionset"]])
    conditions <- input_data[["conditions"]]
    batches <- input_data[["batches"]]
    data <- input_data[["data"]]

    model_choice <- choose_model(input, conditions, batches,
                                 model_batch=model_batch,
                                 model_cond=model_cond,
                                 model_intercept=model_intercept,
                                 alt_model=alt_model, ...)
    model_including <- model_choice[["including"]]
    if (class(model_choice[["model_batch"]]) == "matrix") {
        model_batch <- model_choice[["model_batch"]]
    }
    model_data <- model_choice[["chosen_model"]]

    ## I have a strong sense that the most recent version of edgeR changed its dispersion estimate code
    ## Here is a note from the user's guide, which may have been there previously and I merely did not notice:
    ## To estimate common dispersion, trended dispersions and tagwise dispersions in one run
    ## y <- estimateDisp(y, design)
    raw <- edgeR::DGEList(counts=data, group=conditions)
    message("EdgeR step 1/9: normalizing data.")
    norm <- edgeR::calcNormFactors(raw)
    final_norm <- NULL
    if (edger_method == "short") {
        message("EdgeR steps 2 through 6/9: All in one!")
        final_norm <- edgeR::estimateDisp(norm, design=model_data, robust=TRUE)
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
                warning("estimateGLMTagwiseDisp() failed.  Trying again with estimateDisp.()")
                state <- FALSE
            }
        }
        if (!isTRUE(state)) {
            final_norm <- edgeR::estimateDisp(norm, design=model_data, robust=TRUE)
        }
    }
    cond_fit <- NULL
    if (test_type == "lrt") {
        message("EdgeR step 7/9: Running glmFit, switch to glmQLFit by changing the argument 'test_type'.")
        cond_fit <- edgeR::glmFit(final_norm, design=model_data, robust=TRUE)
    } else {
        message("EdgeR step 7/9: Running glmQLFit, switch to glmFit by changing the argument 'test_type'.")
        cond_fit <- edgeR::glmQLFit(final_norm, design=model_data, robust=TRUE)
    }

    message("EdgeR step 8/9: Making pairwise contrasts.")
    apc <- make_pairwise_contrasts(model_data, conditions,
                                   extra_contrasts=extra_contrasts,
                                   do_identities=FALSE)

    ## This section is convoluted because glmLRT only seems to take up to 7 contrasts at a time.
    ## As a result, I iterate through the set of possible contrasts one at a time and ask for each
    ## separately.
    contrast_list <- list()
    result_list <- list()
    lrt_list <- list()
    sc <- vector("list", length(apc[["names"]]))
    end <- length(apc[["names"]])
    for (con in 1:length(apc[["names"]])) {
        name <- apc[["names"]][[con]]
        message(paste0("EdgeR step 9/9: ", con, "/", end, ": Creating table: ", name, ".")) ## correct
        sc[[name]] <- gsub(pattern=",", replacement="", apc[["all_pairwise"]][[con]])
        tt <- parse(text=sc[[name]])
        ctr_string <- paste0("tt = limma::makeContrasts(", tt, ", levels=model_data)")
        eval(parse(text=ctr_string))
        contrast_list[[name]] <- tt
        lrt_list[[name]] <- NULL
        if (test_type == "lrt") {
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
        res[["qvalue"]] <- tryCatch({
            ttmp <- as.numeric(res[["PValue"]])
            ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
            format(x=ttmp, digits=4, scientific=TRUE)
        },
        error=function(cond) {
            message(paste0("The qvalue estimation failed for ", name, "."))
            return(1)
        },
        finally={
        })
        result_list[[name]] <- res
    } ## End for loop
    final <- list(
        "model" = model_data,
        "contrasts" = apc,
        "contrasts_performed" = apc[["names"]],
        "lrt" = lrt_list,
        "contrast_list" = contrast_list,
        "all_tables" = result_list)
    return(final)
}

#' Writes out the results of a edger search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from edger and friends.
#'
#' Tested in test_24deseq.R
#'
#' @param data  Output from deseq_pairwise()
#' @param ...  Options for writing the xlsx file.
#' @seealso \link[deseq]{toptable} \link{write_xls}
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
