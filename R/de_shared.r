#' Perform limma, DESeq2, EdgeR pairwise analyses.
#'
#' This takes an expt object, collects the set of all possible pairwise comparisons, sets up
#' experimental models appropriate for the differential expression analyses, and performs them.
#'
#' Tested in test_29de_shared.R
#' This runs limma_pairwise(), deseq_pairwise(), edger_pairwise(), basic_pairwise() each in turn.
#' It collects the results and does some simple comparisons among them.
#'
#' @param input  Dataframe/vector or expt class containing count tables, normalization state, etc.
#' @param conditions  Factor of conditions in the experiment.
#' @param batches  Factor of batches in the experiment.
#' @param model_cond  Include condition in the model?  This is likely always true.
#' @param modify_p  Depending on how it is used, sva may require a modification of the p-values.
#' @param model_batch  Include batch in the model?  This may be true/false/"sva" or other methods
#'        supported by get_model_adjust().
#' @param model_intercept  Use an intercept model instead of cell means?
#' @param extra_contrasts  Optional extra contrasts beyone the pairwise comparisons.  This can be
#'        pretty neat, lets say one has conditions A,B,C,D,E and wants to do (C/B)/A and (E/D)/A or
#'        (E/D)/(C/B) then use this with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla =
#'        (E-D)-A, de_vs_cb = (E-D)-(C-B)".
#' @param alt_model  Alternate model to use rather than just condition/batch.
#' @param libsize  Library size of the original data to help voom().
#' @param annot_df  Annotations to add to the result tables.
#' @param parallel  Use dopar to run limma, deseq, edger, and basic simultaneously.
#' @param ...  Picks up extra arguments into arglist, currently only passed to write_limma().
#' @return A list of limma, deseq, edger results.
#' @examples
#' \dontrun{
#'  finished_comparison = eBayes(limma_output)
#'  data_list = all_pairwise(expt)
#' }
#' @export
all_pairwise <- function(input=NULL, conditions=NULL,
                         batches=NULL, model_cond=TRUE,
                         modify_p=FALSE, model_batch=TRUE,
                         model_intercept=TRUE, extra_contrasts=NULL,
                         alt_model=NULL, libsize=NULL,
                         annot_df=NULL, parallel=TRUE, ...) {
    arglist <- list(...)
    surrogates <- "be"
    if (!is.null(arglist[["surrogates"]])) {
        surrogates <- arglist[["surrogates"]]
    }
    if (is.null(model_cond)) {
        model_cond <- TRUE
    }
    if (is.null(model_batch)) {
        model_batch <- FALSE
    }
    if (is.null(model_intercept)) {
        model_intercept <- FALSE
    }
    null_model <- NULL
    sv_model <- NULL
    if (class(model_batch) == "character") {
        model_params <- get_model_adjust(input, estimate_type=model_batch, surrogates=surrogates)
        model_batch <- model_params[["model_adjust"]]
        null_model <- model_params[["null_model"]]
        sv_model <- model_batch
    }

    results <- list(
        "limma" = NULL,
        "deseq" = NULL,
        "edger" = NULL,
        "basic" = NULL)
    res <- NULL
    if (isTRUE(parallel)) {
        cl <- parallel::makeCluster(4)
        doParallel::registerDoParallel(cl)
        requireNamespace("parallel")
        requireNamespace("doParallel")
        requireNamespace("iterators")
        requireNamespace("foreach")
        res <- foreach(c=1:length(names(results)), .packages=c("hpgltools")) %dopar% {
            type <- names(results)[c]
            results[[type]] <- do_pairwise(type,
                                           input=input,
                                           conditions=conditions,
                                           batches=batches,
                                           model_cond=model_cond,
                                           model_batch=model_batch,
                                           model_intercept=model_intercept,
                                           extra_contrasts=extra_contrasts,
                                           alt_model=alt_model,
                                           libsize=libsize,
                                           annot_df=annot_df, ...)
        } ## End foreach() %dopar% { }
        parallel::stopCluster(cl)
        message("Finished running DE analyses, collecting outputs.")
        ## foreach returns the results in no particular order
        ## Therefore, I will reorder the results now and ensure that they are happy.
        for (r in 1:length(res)) {
            a_result <- res[[r]]
            type <- a_result[["type"]]
            results[[type]] <- a_result
        }
        rm(res)
        ## End performing parallel comparisons
    } else {
        for (type in names(results)) {
            results[[type]] <- do_pairwise(type,
                                           input=input,
                                           conditions=conditions,
                                           batches=batches,
                                           model_cond=model_cond,
                                           model_batch=model_batch,
                                           model_intercept=model_intercept,
                                           extra_contrasts=extra_contrasts,
                                           alt_model=alt_model,
                                           libsize=libsize,
                                           annot_df=annot_df, ...)
        }
    } ## End performing a serial comparison

    original_pvalues <- NULL
    ## Add in a little work to re-adjust the p-values in the situation where sva was used
    ## For the moment, DO NOT DO THIS BECAUSE YOU ARE TOO STUPID
    ## Only perform this f adjustment if you modify the data without making
    ## limma/deseq/edger aware of the modified model.  Ergo, if we feed sv_model to this
    ## function, then by definition, we do not want to use this function.
    ## Instead, the opposite is true
    modified_data <- FALSE  ## Thus we will use modified_data to note the data was modified by sva.
    if (is.null(sv_model) & isTRUE(modified_data)) {
        original_pvalues <- data.table::data.table(rownames=rownames(results[["edger"]][["all_tables"]][[1]]))
        message("Using the f.pvalue() function to modify the returned p-values of deseq/limma/edger.")
        ## This is from section 5 of the sva manual:  "Adjusting for surrogate values using the
        ## f.pvalue function. The following chunk of code is longer and more complex than I would
        ## like. This is because f.pvalue() assumes a pairwise comparison of a data set containing
        ## only two experimental factors. As a way to provide an example of _how_ to calculate
        ## appropriately corrected p-values for surrogate factor adjusted models, this is great;
        ## but when dealing with actual data, it falls a bit short.
        for (it in 1:length(results[["edger"]][["all_tables"]])) {
            name <- names(results[["edger"]][["all_tables"]])[it]
            message(paste0("Readjusting the p-values for comparison: ", name))
            namelst <- strsplit(x=name, split="_vs_")
            ## something like 'mutant'
            first <- namelst[[1]][[1]]
            ## something like 'wildtype', ergo the contrast was "mutant_vs_wildtype"
            second <- namelst[[1]][[2]]
            ## The comments that follow will use mutant and wildtype as examples

            ## I am going to need to extract the set of data for the samples in 'first' and 'second'
            ## I will need to also extract the surrogates for those samples from sv_model
            ## Then I rewrite null_model as the subset(null_model, samples included) and rewrite
            ## sv_model as subset(sv_model, samples_included) in that rewrite, there will just be
            ## conditions a, b where a and b are the subsets for first and second. Then the
            ## sv_model will be a for the first samples and b for the second. With that
            ## information, I should e able to feed sva::f.pvalue the appropriate information
            ## for it to run properly. The resulting pvalues will then be appropriate for
            ## backfilling the various tables from edger/limma/deseq.

            ## Get the samples from the limma comparison which are condition 'mutant'
            samples_first_idx <- results[["limma"]][["conditions"]] == first
            num_first <- sum(samples_first_idx)
            ## Subset the expressionset and make a new matrix of only the 'mutant' samples
            samples_first <- Biobase::exprs(input[["expressionset"]])[, samples_first_idx]
            ## Repeat for the 'wildtype' samples, when finished the m columns for samples_first
            ## 'mutant' and the n samples of samples_second will be 'wildtype'
            samples_second_idx <- results[["limma"]][["conditions"]] == second
            num_second <- sum(samples_second_idx)
            samples_second <- Biobase::exprs(input[["expressionset"]])[, samples_second_idx]
            ## Concatenate the 'mutant' and 'wildtype' samples by column
            included_samples <- cbind(samples_first, samples_second)
            ## Arbitrarily call them 'first' and 'second'
            colnames(included_samples) <- c(rep("first", times=num_first),
                                            rep("second", times=num_second))
            ## Do the same thing, but using the rows of the sva model adjustment
            first_sva <- sv_model[samples_first_idx, ]
            second_sva <- sv_model[samples_second_idx, ]
            ## But instead of just appending them, I need a matrix of 0s and 1s identifying which
            ## sv rows correspond to 'wildtype' and 'mutant'
            first_model <- append(rep(1, num_first), rep(0, num_second))
            second_model <- append(rep(0, num_first), rep(1, num_second))
            ## Once I have them, make the subset model matrix with append and cbind
            new_sv_model <- append(first_sva, second_sva)
            new_model <- cbind(first_model, second_model, new_sv_model)
            colnames(new_model) <- c("first","second","sv")
            ## The sva f.pvalue requires a null model of the appropriate size, create that here.
            new_null <- cbind(rep(1, (num_first + num_second)), new_sv_model)
            ## And give its columns suitable names
            colnames(new_null) <- c("null", "sv")
            ## Now all the pieces are in place, call f.pvalue().
            new_pvalues <- try(sva::f.pvalue(included_samples, new_model, new_null), silent=TRUE)
            ## For some things, f.pvalue returns NA, this is unacceptable.
            na_pvalues_idx <- is.na(new_pvalues)
            ## For the NaN pvalues, just set it to 1 under the assumption that something is fubar.
            new_pvalues[na_pvalues_idx] <- 2.2E-16
            ## For strange non-pairwise contrasts, the f.pvalue() call is expected to fail.
            if (class(new_pvalues) == "try-error") {
                new_pvalues <- NULL
                warning(paste0("Unable to adjust pvalues for: ", name))
                warning("If this was not for an extra contrast, then this is a serious problem.")
            } else {
                ## Most of the time it should succeed, so do a BH adjustment of the new values.
                new_adjp <- p.adjust(new_pvalues, method="BH")
                ## Now I need to fill in the tables with these new values.
                ## This section is a little complex.  In brief, it pulls the appropriate
                ## columns from each of the limma, edger, and deseq tables
                ## copies them to a new data.table named 'original_pvalues', and
                ## then replaces them with the just-calculated (adj)p-values.

                ## Start with limma, make no assumptions about table return-order
                limma_table_order <- rownames(results[["limma"]][["all_tables"]][[name]])
                reordered_pvalues <- new_pvalues[limma_table_order]
                reordered_adjp <- new_adjp[limma_table_order]
                ## Create a temporary dt with the old p-values and merge it into original_pvalues.
                tmpdt <- data.table::data.table(results[["limma"]][["all_tables"]][[name]][["P.Value"]])
                tmpdt[["rownames"]] <- rownames(results[["limma"]][["all_tables"]][[name]])
                ## Change the column name of the new data to reflect that it is from limma.
                colnames(tmpdt) <- c(paste0("limma_", name), "rownames")
                original_pvalues <- merge(original_pvalues, tmpdt, by="rownames")
                ## Swap out the p-values and adjusted p-values.
                results[["limma"]][["all_tables"]][[name]][["P.Value"]] <- reordered_pvalues
                results[["limma"]][["all_tables"]][[name]][["adj.P.Val"]] <- reordered_adjp

                ## Repeat the above verbatim, but for edger
                edger_table_order <- rownames(results[["edger"]][["all_tables"]][[name]])
                reordered_pvalues <- new_pvalues[edger_table_order]
                reordered_adjp <- new_adjp[edger_table_order]
                tmpdt <- data.table::data.table(results[["edger"]][["all_tables"]][[name]][["PValue"]])
                tmpdt[["rownames"]] <- rownames(results[["edger"]][["all_tables"]][[name]])
                colnames(tmpdt) <- c(paste0("edger_", name), "rownames")
                original_pvalues <- merge(original_pvalues, tmpdt, by="rownames")
                results[["edger"]][["all_tables"]][[name]][["PValue"]] <- reordered_pvalues
                results[["edger"]][["all_tables"]][[name]][["FDR"]] <- reordered_adjp

                ## Ibid.
                deseq_table_order <- rownames(results[["deseq"]][["all_tables"]][[name]])
                tmpdt <- data.table::data.table(results[["deseq"]][["all_tables"]][[name]][["P.Value"]])
                tmpdt[["rownames"]] <- rownames(results[["deseq"]][["all_tables"]][[name]])
                colnames(tmpdt) <- c(paste0("deseq_", name), "rownames")
                original_pvalues <- merge(original_pvalues, tmpdt, by="rownames")
                results[["deseq"]][["all_tables"]][[name]][["P.Value"]] <- reordered_pvalues
                results[["deseq"]][["all_tables"]][[name]][["adj.P.Val"]] <- reordered_adjp
            } ## End checking that f.pvalue worked.
        }  ## End foreach table
        original_pvalues <- as.data.frame(original_pvalues)
    } ## End checking if we should f-test modify the p-values

    result_comparison <- list()
    if (class(results[["limma"]]) == "list" &
        class(results[["deseq"]]) == "list" &
        class(results[["edger"]]) == "list" &
        class(results[["basic"]]) == "list") {
        result_comparison <- compare_tables(limma=results[["limma"]],
                                            deseq=results[["deseq"]],
                                            edger=results[["edger"]],
                                            basic=results[["basic"]],
                                            annot_df=annot_df, ...)
    }
    ## The first few elements of this list are being passed through into the return
    ## So that if I use combine_tables() I can report in the resulting tables
    ## some information about what was performed.
    ret <- list(
        "model_cond" = model_cond,
        "model_batch" = model_batch,
        "extra_contrasts" = extra_contrasts,
        "original_pvalues" = original_pvalues,
        "input" = input,
        "limma" = results[["limma"]],
        "deseq" = results[["deseq"]],
        "edger" = results[["edger"]],
        "basic" = results[["basic"]],
        "comparison" = result_comparison)
    return(ret)
}

#' Try out a few experimental models and return a likely working option.
#'
#' The _pairwise family of functions all demand an experimental model.  This tries to choose a
#' consistent and useful model for all for them.  This does not try to do multi-factor, interacting,
#' nor dependent variable models, if you want those do them yourself and pass them off as alt_model.
#'
#' Invoked by the _pairwise() functions.
#'
#' @param input  Input data used to make the model.
#' @param conditions  Factor of conditions in the putative model.
#' @param batches  Factor of batches in the putative model.
#' @param model_batch  Try to include batch in the model?
#' @param model_cond  Try to include condition in the model? (Yes!)
#' @param model_intercept  Use an intercept model instead of cell-means?
#' @param alt_model  Use your own model.
#' @param alt_string  String describing an alternate model.
#' @param intercept  Choose an intercept for the model as opposed to 0.
#' @param reverse  Reverse condition/batch in the model?  This shouldn't/doesn't matter but I wanted
#'     to test.
#' @param surrogates  Number of or method used to choose the number of surrogate variables.
#' @param ...  Further options are passed to arglist.
#' @return List including a model matrix and strings describing cell-means and intercept models.
#' @export
choose_model <- function(input, conditions, batches, model_batch=TRUE,
                         model_cond=TRUE, model_intercept=TRUE,
                         alt_model=NULL, alt_string=NULL,
                         intercept=0, reverse=FALSE,
                         surrogates="be", ...) {
    arglist <- list(...)
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    ## It would be much smarter to generate the models in the following if() {} blocks
    ## But I have it in my head to eventually compare results using different models.
    cond_int_string <- "~ 0 + condition"
    cond_int_model <- stats::model.matrix(~ 0 + conditions,
                                          contrasts.arg=list(conditions="contr.treatment"))
    batch_int_string <- "~ 0 + batch"
    batch_int_model <- try(stats::model.matrix(~ 0 + batches,
                                               contrasts.arg=list(batches="contr.treatment")),
                           silent=TRUE)
    condbatch_int_string <- "~ 0 + condition + batch"
    condbatch_int_model <- try(stats::model.matrix(~ 0 + conditions + batches,
                                                   contrasts.arg=list(conditions="contr.treatment",
                                                                      batches="contr.treatment")),
                               silent=TRUE)
    batchcond_int_string <- "~ 0 + batch + condition"
    batchcond_int_model <- try(stats::model.matrix(~ 0 + batches + conditions,
                                                   contrasts.arg=list(conditions="contr.treatment",
                                                                      batches="contr.treatment")),
                               silent=TRUE)
    cond_noint_string <- "~ condition"
    cond_noint_model <- try(stats::model.matrix(~ conditions,
                                                contrasts.arg=list(conditions="contr.treatment")),
                            silent=TRUE)
    batch_noint_string <- "~ batch"
    batch_noint_model <- try(stats::model.matrix(~ batches,
                                                 contrasts.arg=list(batches="contr.treatment")),
                             silent=TRUE)
    condbatch_noint_string <- "~ condition + batch"
    condbatch_noint_model <- try(stats::model.matrix(~ conditions + batches,
                                                     contrasts.arg=list(conditions="contr.treatment",
                                                                        batches="contr.treatment")),
                                 silent=TRUE)
    batchcond_noint_string <- "~ batch + condition"
    batchcond_noint_model <- try(stats::model.matrix(~ batches + conditions,
                                                     contrasts.arg=list(conditions="contr.treatment",
                                                                        batches="contr.treatment")),
                                 silent=TRUE)
    noint_model <- NULL
    int_model <- NULL
    noint_string <- NULL
    int_string <- NULL
    including <- NULL
    if (is.null(model_batch)) {
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    } else if (isTRUE(model_cond) & isTRUE(model_batch)) {
        if (class(condbatch_int_model) == "try-error") {
            message("The condition+batch model failed.  Does your experimental design support both condition and batch?")
            message("Using only a conditional model.")
            int_model <- cond_int_model
            noint_model <- cond_noint_model
            int_string <- cond_int_string
            noint_string <- cond_noint_string
            including <- "condition"
        } else if (isTRUE(reverse)) {
            int_model <- batchcond_int_model
            noint_model <- batchcond_noint_model
            int_string <- batchcond_int_string
            noint_string <- batchcond_noint_string
            including <- "batch+condition"
        } else {
            int_model <- condbatch_int_model
            noint_model <- condbatch_noint_model
            int_string <- condbatch_int_string
            noint_string <- condbatch_noint_string
            including <- "condition+batch"
        }
    } else if (class(model_batch) == "character") {
        ## Then calculate the estimates using get_model_adjust
        message("Extracting surrogate estimate from sva/ruv/pca and adding them to the model.")
        model_batch_info <- get_model_adjust(input, estimate_type=model_batch, surrogates=surrogates)
        ## Changing model_batch from 'sva' to the resulting matrix.
        ## Hopefully this will simplify things later for me.
        model_batch <- model_batch_info[["model_adjust"]]
        int_model <- stats::model.matrix(~ 0 + conditions + model_batch)
                                         ## contrasts.arg=list(conditions="contr.sum"))
        noint_model <- stats::model.matrix(~ conditions + model_batch)
                                           ## contrasts.arg=list(conditions="contr.sum"))
        int_string <- condbatch_int_string
        noint_string <- condbatch_noint_string
        including <- "condition+batchestimate"
    } else if (class(model_batch) == "numeric" | class(model_batch) == "matrix") {
        message("Including batch estimates from sva/ruv/pca in the model.")
        int_model <- stats::model.matrix(~ 0 + conditions + model_batch)
                                           ## contrasts.arg=list(conditions="contr.sum"))
        noint_model <- stats::model.matrix(~ conditions + model_batch)
                                           ## contrasts.arg=list(conditions="contr.sum"))
        int_string <- condbatch_int_string
        noint_string <- condbatch_noint_string
        including <- "condition+batchestimate"
    } else if (isTRUE(model_cond)) {
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    } else if (isTRUE(model_batch)) {
        int_model <- batch_int_model
        noint_model <- batch_noint_model
        int_string <- batch_int_string
        noint_string <- batch_noint_string
        including <- "batch"
    } else {
        ## Default to the conditional model
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    }

    tmpnames <- colnames(int_model)
    tmpnames <- gsub("data[[:punct:]]", "", tmpnames)
    tmpnames <- gsub("-", "", tmpnames)
    tmpnames <- gsub("+", "", tmpnames)
    ## The next lines ensure that conditions/batches which are all numeric will not cause weird
    ## errors for contrasts. Ergo, if a condition is something like '111', now it will be 'c111'
    ## Similarly, a batch '01' will be 'b01'
    tmpnames <- gsub("^conditions(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub("^batches(\\d+)$", replacement="b\\1", x=tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    tmpnames <- gsub("batches", "", tmpnames)
    colnames(int_model) <- tmpnames

    tmpnames <- colnames(noint_model)
    tmpnames <- gsub("data[[:punct:]]", "", tmpnames)
    tmpnames <- gsub("-", "", tmpnames)
    tmpnames <- gsub("+", "", tmpnames)
    ## The next lines ensure that conditions/batches which are all numeric will not cause weird
    ## errors for contrasts. Ergo, if a condition is something like '111', now it will be 'c111'
    ## Similarly, a batch '01' will be 'b01'
    tmpnames <- gsub("conditions^(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub("batches^(\\d+)$", replacement="b\\1", x=tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    tmpnames <- gsub("batches", "", tmpnames)
    colnames(noint_model) <- tmpnames

    chosen_model <- NULL
    chosen_string <- NULL
    if (isTRUE(model_intercept)) {
        message("Choosing the intercept containing model.")
        chosen_model <- int_model
        chosen_string <- int_string
    } else {
        chosen_model <- noint_model
        chosen_string <- noint_string
    }

    if (!is.null(alt_model)) {
        fun_model <- alt_model
        fun_string <- alt_string
        including <- "alt"
    }

    retlist <- list(
        "int_model" = int_model,
        "noint_model" = noint_model,
        "int_string" = int_string,
        "noint_string" = noint_string,
        "chosen_model" = chosen_model,
        "chosen_string" = chosen_string,
        "model_batch" = model_batch,
        "including" = including)
    return(retlist)
}

#' Choose a suitable data set for Edger/DESeq
#'
#' The _pairwise family of functions all demand data in specific formats.
#' This tries to make that consistent.
#'
#' Invoked by _pairwise().
#'
#' @param input  Expt input.
#' @param force  Force non-standard data?
#' @param choose_for  One of limma, deseq, edger, or basic.  Defines the requested data state.
#' @param ...  More options for future expansion.
#' @return List the data, conditions, and batches in the data.
#' @export
choose_dataset <- function(input, choose_for="limma", force=FALSE, ...) {
    arglist <- list(...)
    result <- NULL
    if (choose_for == "limma") {
        result <- choose_limma_dataset(input, force=force, ...)
    } else if (choose_for == "basic") {
        result <- choose_basic_dataset(input, force=force, ...)
    } else if (choose_for == "edger") {
        result <- choose_binom_dataset(input, force=force, ...)
    } else if (choose_for == "deseq") {
        result <- choose_binom_dataset(input, force=force, ...)
    } else {
        message("Unknown tool for which to choose a data set.")
        result <- list(
            "conditions" = input[["design"]][["condition"]],
            "batches" = input[["design"]][["batch"]],
            "data" = as.data.frame(Biobase::exprs(input[["expressionset"]])))
    }
    return(result)
}

#' A sanity check that a given set of data is suitable for analysis by limma.
#'
#' Take an expt and poke at it to ensure that it will not result in troubled limma results.
#'
#' @param input Expressionset containing expt object.
#' @param force Ingore warnings and use the provided data asis.
#' @param which_voom  Choose between limma'svoom, voomWithQualityWeights, or the hpgl equivalents.
#' @param ... Extra arguments passed to arglist.
#' @return dataset suitable for limma analysis
choose_limma_dataset <- function(input, force=FALSE, which_voom="limma", ...) {
    arglist <- list(...)
    input_class <- class(input)[1]
    data <- NULL
    warn_user <- 0
    libsize <- NULL
    ## It turns out, a more careful examination of how normalization affects the results, the above
    ## seems only to be true if the following are true:
    ## 1.  There are >2-3k features(genes/transcripts) with a full range of count values.
    ## 2.  One does not attempt to use sva, or at least one uses sva before messing with the
    ##     normalization state.
    ## 2a. #2 primarily applies if one is using quantile normalization, it looks like tmm/rle
    ##     does not have so profound an effect, and this effect is tightly bound with the state of
    ##     #1 -- in other words, if one has nice dense data with low->high counts in an even
    ##     distribution, then quantile+sva might be ok. But if that is not true, then one should
    ##     expect a poo-show result.
    ## For these reasons I am telling this function to revert to non-normalized data unless force
    ## is on, just like I do for edger/deseq.  I think to help this, I will add a parameter which
    ## allows one to to turn on/off normalization at the voom() step.

    if (input_class == "expt") {
        conditions <- input[["conditions"]]
        batches <- input[["batches"]]
        data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))

        tran_state <- input[["state"]][["transform"]]
        ## Note that voom will take care of this for us.
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
        data <- Biobase::exprs(input[["expressionset"]])
        if (isTRUE(force)) {
            message("Leaving the data alone, regardless of normalization state.")
            retlist <- list(
                "libsize" = libsize,
                "conditions" = conditions,
                "batches" = batches,
                "data" = data)
            return(retlist)
        }

        ## If we are using limma::voom*, then make sure we do things the limma way.
        ## If we use the hpgltools::hpgl_voom*, let the freak flags fly.
        if (grepl(pattern="limma", x=which_voom)) {
            ## Limma's voom requires we return log2(cpm()) to base 10.
            ## Otherwise it should accept pretty much anything.
            if (tran_state == "log2") {
                message("Using limma's voom, returning to base 10.")
                data <- 2 ^ data
            }
        }
    } else {
        data <- as.data.frame(input)
    }
    head(data)
    retlist <- list(
        "libsize" = libsize,
        "conditions" = conditions,
        "batches" = batches,
        "data" = data)
    return(retlist)
}

#' A sanity check that a given set of data is suitable for analysis by edgeR or DESeq2.
#'
#' Take an expt and poke at it to ensure that it will not result in troubled results.
#'
#' @param input Expressionset containing expt object.
#' @param force Ignore every warning and just use this data.
#' @param ... Extra arguments passed to arglist.
#' @return dataset suitable for limma analysis
choose_binom_dataset <- function(input, force=FALSE, ...) {
    arglist <- list(...)
    input_class <- class(input)[1]
    ## I think I would like to make this function smarter so that it will remove the log2 from
    ## transformed data.
    data <- NULL
    warn_user <- 0
    if (input_class == "expt") {
        conditions <- input[["conditions"]]
        batches <- input[["batches"]]
        data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))
        ## As I understand it, EdgeR fits a binomial distribution
        ## and expects data as integer counts, not floating point nor a log2 transformation
        ## Thus, having the 'normalization' state set to something other than 'raw' is a likely
        ## violation of its stated preferred/demanded input.  There are of course ways around this
        ## but one should not take them lightly, or ever.
        tran_state <- input[["state"]][["transform"]]
        if (is.null(tran_state)) {
            tran_state <- "raw"
        }
        conv_state <- input[["state"]][["conversion"]]
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
        if (norm_state == "round") {
            norm_state <- "raw"
        }

        if (isTRUE(force)) {
            ## Setting force to TRUE allows one to round the data to fool edger/deseq into accepting it
            ## This is a pretty terrible thing to do
            message("About to round the data, this is a pretty terrible thing to do. But if you,
like me, want to see what happens when you put non-standard data into deseq, then here you go.")
            data <- round(data)
            less_than <- data < 0
            data[less_than] <- 0
            warn_user <- 1
        } else if (norm_state != "raw" & tran_state != "raw" & conv_state != "raw") {
            ## These if statements may be insufficient to check for the appropriate input for deseq.
            data <- Biobase::exprs(input[["original_expressionset"]])
        } else if (norm_state != "raw" | tran_state != "raw") {
            ## This makes use of the fact that the order of operations in the normalization
            ## function is static. filter->normalization->convert->batch->transform.
            ## Thus, if the normalized state is not raw, we can look back either to the filtered
            ## or original data. The same is true for the transformation state.
            message("EdgeR/DESeq expect raw data as input, reverting to the count filtered data.")
            data <- input[["normalized"]][["intermediate_counts"]][["filter"]][["count_table"]]
            if (is.null(data)) {
                data <- input[["normalized"]][["intermediate_counts"]][["original"]]
            }
        } else {
            message("The data should be suitable for EdgeR/DESeq.")
            message("If EdgeR/DESeq freaks out, check the state of the count table and ensure that it is in integer counts.")
        }
        ## End testing if normalization has been performed
    } else {
        data <- as.data.frame(input)
    }
    retlist <- list(
        "conditions" = conditions,
        "batches" = batches,
        "data" = data)
    if (warn_user == 1) {
        warning("This data was inappropriately forced into integers.")
    }
    return(retlist)
}

#' See how similar are results from limma/deseq/edger.
#'
#' limma, DEseq2, and EdgeR all make somewhat different assumptions.
#' and choices about what makes a meaningful set of differentially.
#' expressed genes.  This seeks to provide a quick and dirty metric
#' describing the degree to which they (dis)agree.
#'
#' @param limma Data from limma_pairwise().
#' @param deseq Data from deseq2_pairwise().
#' @param edger Data from edger_pairwise().
#' @param basic Data from basic_pairwise().
#' @param include_basic include the basic data?
#' @param annot_df Include annotation data?
#' @param ... More options!
#' @return Heatmap showing how similar they are along with some
#'     correlations betwee the three players.
#' @seealso \code{\link{limma_pairwise}} \code{\link{edger_pairwise}} \code{\link{deseq2_pairwise}}
#' @examples
#' \dontrun{
#'  l = limma_pairwise(expt)
#'  d = deseq_pairwise(expt)
#'  e = edger_pairwise(expt)
#'  fun = compare_tables(limma=l, deseq=d, edger=e)
#' }
#' @export
compare_tables <- function(limma=NULL, deseq=NULL, edger=NULL, basic=NULL,
                           include_basic=TRUE, annot_df=NULL, ...) {
    arglist <- list(...)
    ## Fill each column/row of these with the correlation between tools for one contrast performed
    if (class(limma) == "list") {
        ## Then this was fed the raw output from limma_pairwise,
        ## lets assume the same is true for deseq/edger too and pull out the result tables.
        limma <- limma[["all_tables"]]
        deseq <- deseq[["all_tables"]]
        edger <- edger[["all_tables"]]
        basic <- basic[["all_tables"]]
    }
    len <- length(names(deseq))
    limma_vs_edger <- list()
    limma_vs_deseq <- list()
    limma_vs_basic <- list()
    edger_vs_deseq <- list()
    edger_vs_basic <- list()
    deseq_vs_basic <- list()
    limma_vs_edger_scatter <- list()
    limma_vs_deseq_scatter <- list()
    limma_vs_basic_scatter <- list()
    edger_vs_deseq_scatter <- list()
    edger_vs_basic_scatter <- list()
    deseq_vs_basic_scatter <- list()
    cc <- 0
    for (comp in names(deseq)) {
        ## assume all three have the same names() -- note that limma has more than the other two though
        cc <- cc + 1
        message(paste0("Comparing analyses ", cc, "/", len, ": ", comp))
        l <- data.frame(limma[[comp]])
        e <- data.frame(edger[[comp]])
        d <- data.frame(deseq[[comp]])
        b <- data.frame(basic[[comp]])
        le <- merge(l, e, by.x="row.names", by.y="row.names")
        le <- le[,c("logFC.x","logFC.y")]
        colnames(le) <- c("limma logFC", "edgeR logFC")
        lec <- stats::cor.test(x=le[, 1], y=le[, 2])[["estimate"]]
        les <- plot_scatter(le) + ggplot2::labs(title=paste0(comp, ": limma vs. edgeR.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        ld <- merge(l, d, by.x="row.names", by.y="row.names")
        ld <- ld[, c("logFC.x","logFC.y")]
        colnames(ld) <- c("limma logFC", "DESeq2 logFC")
        ldc <- stats::cor.test(ld[,1], ld[,2])[["estimate"]]
        lds <- plot_scatter(ld) + ggplot2::labs(title=paste0(comp, ": limma vs. DESeq2.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        lb <- merge(l, b, by.x="row.names", by.y="row.names")
        lb <- lb[, c("logFC.x","logFC.y")]
        colnames(lb) <- c("limma logFC", "basic logFC")
        lbc <- stats::cor.test(lb[,1], lb[,2])[["estimate"]]
        lbs <- plot_scatter(lb) + ggplot2::labs(title=paste0(comp, ": limma vs. basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        ed <- merge(e, d, by.x="row.names", by.y="row.names")
        ed <- ed[, c("logFC.x","logFC.y")]
        colnames(ed) <- c("edgeR logFC", "DESeq2 logFC")
        edc <- stats::cor.test(ed[,1], ed[,2])[["estimate"]]
        eds <- plot_scatter(ed) + ggplot2::labs(title=paste0(comp, ": edgeR vs. DESeq2.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        eb <- merge(e, b, by.x="row.names", by.y="row.names")
        eb <- eb[, c("logFC.x","logFC.y")]
        colnames(eb) <- c("edgeR logFC", "basic logFC")
        ebc <- stats::cor.test(eb[,1], eb[,2])[["estimate"]]
        ebs <- plot_scatter(eb) + ggplot2::labs(title=paste0(comp, ": edgeR vs. basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        db <- merge(d, b, by.x="row.names", by.y="row.names")
        db <- db[, c("logFC.x","logFC.y")]
        colnames(db) <- c("DESeq2 logFC", "basic logFC")
        dbc <- stats::cor.test(db[,1], db[,2])[["estimate"]]
        dbs <- plot_scatter(db) +
            ggplot2::labs(title=paste0(comp, ": DESeq2 vs basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        limma_vs_edger[[comp]] <- lec
        limma_vs_deseq[[comp]] <- ldc
        edger_vs_deseq[[comp]] <- edc
        limma_vs_basic[[comp]] <- lbc
        edger_vs_basic[[comp]] <- ebc
        deseq_vs_basic[[comp]] <- dbc
        limma_vs_edger_scatter[[comp]] <- les
        limma_vs_deseq_scatter[[comp]] <- lds
        edger_vs_deseq_scatter[[comp]] <- eds
        limma_vs_basic_scatter[[comp]] <- lbs
        edger_vs_basic_scatter[[comp]] <- ebs
        deseq_vs_basic_scatter[[comp]] <- dbs
    } ## End loop
    names(limma_vs_edger) <- names(deseq)
    names(limma_vs_deseq) <- names(deseq)
    names(edger_vs_deseq) <- names(deseq)
    names(limma_vs_basic) <- names(deseq)
    names(edger_vs_basic) <- names(deseq)
    names(deseq_vs_basic) <- names(deseq)
    names(limma_vs_edger_scatter) <- names(deseq)
    names(limma_vs_deseq_scatter) <- names(deseq)
    names(edger_vs_deseq_scatter) <- names(deseq)
    names(limma_vs_basic_scatter) <- names(deseq)
    names(edger_vs_basic_scatter) <- names(deseq)
    names(deseq_vs_basic_scatter) <- names(deseq)

    comparison_df <- rbind(as.numeric(limma_vs_edger), as.numeric(limma_vs_deseq))
    comparison_df <- rbind(comparison_df, as.numeric(edger_vs_deseq))
    if (isTRUE(include_basic)) {
        comparison_df <- rbind(comparison_df, as.numeric(limma_vs_basic))
        comparison_df <- rbind(comparison_df, as.numeric(edger_vs_basic))
        comparison_df <- rbind(comparison_df, as.numeric(deseq_vs_basic))
        rownames(comparison_df) <- c("le", "ld", "ed", "lb", "eb", "db")
    } else {
        rownames(comparison_df) <- c("le", "ld", "ed")
    }
    comparison_df <- as.matrix(comparison_df)
    colnames(comparison_df) <- names(deseq)
    heat_colors <- grDevices::colorRampPalette(c("white","black"))
    comparison_heatmap <- try(heatmap.3(comparison_df, scale="none",
                                        trace="none", keysize=1.5,
                                        linewidth=0.5, margins=c(9,9),
                                        col=heat_colors, dendrogram="none",
                                        Rowv=FALSE, Colv=FALSE,
                                        main="Compare DE tools"), silent=TRUE)
    heat <- NULL
    if (class(comparison_heatmap) != 'try-error') {
        heat <- recordPlot()
    }
    ret <- list(
        "limma_vs_edger" = limma_vs_edger,
        "limma_vs_deseq" = limma_vs_deseq,
        "limma_vs_basic" = limma_vs_basic,
        "edger_vs_deseq" = edger_vs_deseq,
        "edger_vs_basic" = edger_vs_basic,
        "deseq_vs_basic" = deseq_vs_basic,
        "limma_vs_edger_scatter" = limma_vs_edger_scatter,
        "limma_vs_deseq_scatter" = limma_vs_deseq_scatter,
        "limma_vs_basic_scatter" = limma_vs_basic_scatter,
        "edger_vs_deseq_scatter" = edger_vs_deseq_scatter,
        "edger_vs_basic_scatter" = edger_vs_basic_scatter,
        "deseq_vs_basic_scatter" = deseq_vs_basic_scatter,
        "comp" = comparison_df,
        "heat" = heat)
    return(ret)
}

#' Compare logFC values from limma and friends
#'
#' There are some peculiar discrepencies among these tools, what is up with that?
#'
#' @param combined_tables The combined tables from limma et al.
#' @return Some plots
#' @export
compare_logfc_plots <- function(combined_tables) {
    plots <- list()
    data <- NULL
    if (!is.null(combined_tables[["data"]])) {
        data <- combined_tables[["data"]]
    } else {
        data <- combined_tables
    }
    for (count in 1:length(data)) {
        tab <- data[[count]]
        le_data <- tab[, c("limma_logfc", "edger_logfc", "limma_adjp", "edger_adjp")]
        le <- sm(plot_linear_scatter(le_data, pretty_colors=FALSE)[["scatter"]])
        ld_data <- tab[, c("limma_logfc", "deseq_logfc", "limma_adjp", "deseq_adjp")]
        ld <- sm(plot_linear_scatter(ld_data, pretty_colors=FALSE)[["scatter"]])
        de_data <- tab[, c("deseq_logfc", "edger_logfc", "deseq_adjp", "edger_adjp")]
        de <- sm(plot_linear_scatter(de_data, pretty_colors=FALSE)[["scatter"]])
        lb_data <- tab[, c("limma_logfc", "basic_logfc", "limma_adjp", "basic_p")]
        lb <- sm(plot_linear_scatter(lb_data, pretty_colors=FALSE)[["scatter"]])
        db_data <- tab[, c("deseq_logfc", "basic_logfc", "deseq_adjp", "basic_p")]
        db <- sm(plot_linear_scatter(db_data, pretty_colors=FALSE)[["scatter"]])
        eb_data <- tab[, c("edger_logfc", "basic_logfc", "edger_adjp", "basic_p")]
        eb <- sm(plot_linear_scatter(eb_data, pretty_colors=FALSE)[["scatter"]])
        name <- names(data)[[count]]
        compared <- list(
            "le" = le, "ld" = ld, "de" = de,
            "lb" = lb, "db" = db, "eb" = eb)
        plots[[name]] <- compared
    }
    return(plots)
}

#' Test for infected/control/beads -- a placebo effect?
#'
#' The goal is therefore to find responses different than beads
#' The null hypothesis is (H0): (infected == uninfected) || (infected == beads)
#' The alt hypothesis is (HA): (infected != uninfected) && (infected != beads)
#'
#' @param contrast_fit  The result of lmFit.
#' @param coef1  The first coefficient to query.
#' @param coef2  And the second.
#' @param ...  Extra arguments are passed to arglist, but basically ignored.
disjunct_pvalues <- function(contrast_fit, coef1, coef2, ...) {
    arglist <- list(...)
    stat <- BiocGenerics::pmin(abs(contrast_fit[, coef1]), abs(contrast_fit[, coef2]))
    pval <- BiocGenerics::pmax(contrast_fit$p.val[, coef1], contrast_fit$p.val[, coef2])
    return(pval)
}

#' Generalize pairwise comparisons
#'
#' I want to multithread my pairwise comparisons, this is the first step in doing so.
#'
#' @param type  Which type of pairwise comparison to perform
#' @param ...  The set of arguments intended for limma_pairwise(), edger_pairwise(), and friends.
#' @return The result from limma/deseq/edger/basic
#' @export
do_pairwise <- function(type, ...) {
    arglist <- list(...)
    res <- NULL
    if (type == "limma") {
        res <- try(limma_pairwise(...))
    } else if (type == "edger") {
        res <- try(edger_pairwise(...))
    } else if (type == "deseq") {
        res <- try(deseq2_pairwise(...))
    } else if (type == "basic") {
        res <- try(basic_pairwise(...))
    }
    res[["type"]] <- type
    return(res)
}

#' Get a set of up/down differentially expressed genes.
#'
#' Take one or more criteria (fold change, rank order, (adj)p-value,
#' z-score from median FC) and use them to extract the set of genes
#' which are defined as 'differentially expressed.'  If no criteria
#' are provided, it arbitrarily chooses all genes outside of 1-z.
#'
#' @param table Table from limma/edger/deseq.
#' @param n Rank-order top/bottom number of genes to take.
#' @param z Number of z-scores >/< the median to take.
#' @param fc Fold-change cutoff.
#' @param p P-value cutoff.
#' @param fold Identifier reminding how to get the bottom portion of a
#'     fold-change (plusminus says to get the negative of the
#'     positive, otherwise 1/positive is taken).  This effectively
#'     tells me if this is a log fold change or not.
#' @param column Table's column used to distinguish top vs. bottom.
#' @param p_column Table's column containing (adjusted or not)p-values.
#' @return Subset of the up/down genes given the provided criteria.
#' @export
get_sig_genes <- function(table, n=NULL, z=NULL, fc=NULL, p=NULL,
                          column="logFC", fold="plusminus", p_column="adj.P.Val") {
    if (is.null(z) & is.null(n) & is.null(fc) & is.null(p)) {
        message("No n, z, p, nor fc provided, setting p to 0.05 and fc to 1.0.")
        p <- 0.05
        fc <- 1.0
    }
    up_genes <- table
    down_genes <- table

    if (!is.null(p)) {
        up_idx <- as.numeric(up_genes[[p_column]]) <= p
        ## Remember we have these reformatted as scientific
        up_genes <- up_genes[up_idx, ]
        down_idx <- as.numeric(down_genes[[p_column]]) <= p
        down_genes <- down_genes[down_idx, ]
        ## Going to add logic in case one does not ask for fold change
        ## In that case, a p-value assertion should still know the difference between up and down
        ## But it should also still know the difference between ratio and log changes
        if (fold == 'plusminus' | fold == 'log') {
            message(paste0("Assuming the fold changes are on the log scale and so taking >< 0"))
            ##up_idx <- up_genes[, column] > 0.0
            up_idx <- as.numeric(up_genes[[column]]) > 0.0
            up_genes <- up_genes[up_idx, ]
            down_idx <- as.numeric(down_genes[[column]]) < 0.0
            down_genes <- down_genes[down_idx, ]
        } else {
            ## plusminus refers to a positive/negative number of logfold changes from a logFC(1) = 0
            up_idx <- as.numeric(up_genes[[column]]) >= 1.0
            up_genes <- up_genes[up_idx, ]
            down_idx <- as.numeric(down_genes[[column]]) <= -1.0
            down_genes <- down_genes[down_idx, ]
        }
        message(paste0("After (adj)p filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After (adj)p filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(fc)) {
        up_idx <- as.numeric(up_genes[[column]]) >= fc
        up_genes <- up_genes[up_idx, ]
        if (fold == 'plusminus' | fold == 'log') {
            message(paste0("Assuming the fold changes are on the log scale and so taking -1.0 * fc"))
            ## plusminus refers to a positive/negative number of logfold changes from a logFC(1) = 0
            down_idx <- as.numeric(down_genes[[column]]) <= (fc * -1.0)
            down_genes <- down_genes[down_idx, ]
        } else {
            message(paste0("Assuming the fold changes are on a ratio scale and so taking 1/fc"))
            ## If it isn't log fold change, then values go from 0..x where 1 is unchanged
            down_idx <- as.numeric(down_genes[[column]]) <= (1.0 / fc)
            down_genes <- down_genes[down_idx, ]
        }
        message(paste0("After fold change filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After fold change filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(z)) {
        ## Take an arbitrary number which are >= and <= a value which is z zscores from the median.
        message(paste0("Getting the genes >= ", z, " z scores away from the median of all."))
        ## Use the entire table for the summary
        out_summary <- summary(as.numeric(table[[column]]))
        out_mad <- stats::mad(as.numeric(table[[column]]), na.rm=TRUE)
        up_median_dist <- out_summary["Median"] + (out_mad * z)
        down_median_dist <- out_summary["Median"] - (out_mad * z)
        ## But use the (potentially already trimmed) up/down tables for indexing
        up_idx <- as.numeric(up_genes[[column]]) >= up_median_dist
        up_genes <- up_genes[up_idx, ]
        down_idx <- as.numeric(down_genes[[column]]) <= down_median_dist
        down_genes <- down_genes[down_idx, ]
        message(paste0("After z filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After z filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(n)) {
        ## Take a specific number of genes at the top/bottom of the rank ordered list.
        message(paste0("Getting the top and bottom ", n, " genes."))
        upranked <- up_genes[order(as.numeric(up_genes[[column]]), decreasing=TRUE), ]
        up_genes <- head(upranked, n=n)
        downranked <- down_genes[order(as.numeric(down_genes[[column]])), ]
        down_genes <- head(downranked, n=n)
        message(paste0("After top-n filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After bottom-n filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }
    up_genes <- up_genes[order(as.numeric(up_genes[[column]]), decreasing=TRUE), ]
    down_genes <- down_genes[order(as.numeric(down_genes[[column]]), decreasing=FALSE), ]
    ret = list(
        "up_genes" = up_genes,
        "down_genes" = down_genes)
    return(ret)
}

#' Small hack of limma's exampleData() to allow for arbitrary data set
#' sizes.
#'
#' exampleData has a set number of genes/samples it creates. This
#' relaxes that restriction.
#'
#' @param ngenes How many genes in the fictional data set?
#' @param columns How many samples in this data set?
#' @return Matrix of pretend counts.
#' @seealso \pkg{limma}
#' @examples
#' \dontrun{
#'  pretend = make_exampledata()
#' }
#' @export
make_exampledata <- function (ngenes=1000, columns=5) {
    q0 <- stats::rexp(ngenes, rate = 1/250)
    is_DE <- stats::runif(ngenes) < 0.3
    lfc <- stats::rnorm(ngenes, sd = 2)
    q0A <- ifelse(is_DE, q0 * 2^(lfc / 2), q0)
    q0B <- ifelse(is_DE, q0 * 2^(-lfc / 2), q0)
    ##    true_sf <- c(1, 1.3, 0.7, 0.9, 1.6)
    true_sf <- abs(stats::rnorm(columns, mean=1, sd=0.4))
    cond_types <- ceiling(sqrt(columns))
    ##    conds <- c("A", "A", "B", "B", "B")
    ##x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
    conds <- sample(LETTERS[1:cond_types], columns, replace=TRUE)
    m <- t(sapply(seq_len(ngenes),
                  function(i) sapply(1:columns,
                                     function(j) rnbinom(1,
                                                         mu = true_sf[j] * ifelse(conds[j] == "A",
                                                                                  q0A[i], q0B[i]),
                                                         size = 1/0.2))))
    rownames(m) <- paste("gene", seq_len(ngenes), ifelse(is_DE, "T", "F"), sep = "_")
    example <- DESeq::newCountDataSet(m, conds)
    return(example)
}

#' Run makeContrasts() with all pairwise comparisons.
#'
#' In order to have uniformly consistent pairwise contrasts, I decided
#' to avoid potential human erors(sic) by having a function generate
#' all contrasts.
#'
#' @param model  Describe the conditions/batches/etc in the experiment.
#' @param conditions  Factor of conditions in the experiment.
#' @param do_identities  Include all the identity strings? Limma can
#'  use this information while edgeR can not.
#' @param do_pairwise Include all pairwise strings? This shouldn't
#'  need to be set to FALSE, but just in case.
#' @param extra_contrasts Optional string of extra contrasts to include.
#' @return List including the following information:
#'  all_pairwise_contrasts = the result from makeContrasts(...)
#'  identities = the string identifying each condition alone
#'  all_pairwise = the string identifying each pairwise comparison alone
#'  contrast_string = the string passed to R to call makeContrasts(...)
#'  names = the names given to the identities/contrasts
#' @seealso \link[limma]{makeContrasts}
#' @examples
#' \dontrun{
#'  pretend = make_pairwise_contrasts(model, conditions)
#' }
#' @export
make_pairwise_contrasts <- function(model, conditions, do_identities=TRUE,
                                    do_pairwise=TRUE, extra_contrasts=NULL) {
    tmpnames <- colnames(model)
    tmpnames <- gsub(pattern="data[[:punct:]]", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="-", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="+", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="conditions^(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub(pattern="conditions", replacement="", x=tmpnames)
    colnames(model) <- tmpnames
    conditions <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=conditions)
    condition_table <- table(conditions)
    identities <- list()
    contrast_string <- ""
    eval_strings <- list()
    for (c in 1:length(condition_table)) {
        identity_name <- names(condition_table[c])
        identity_string <- paste(identity_name, " = ", identity_name, ",", sep="")
        identities[identity_name] <- identity_string
    }
    ## If I also create a sample condition 'alice', and also perform a subtraction
    ## of 'alice' from 'bob', then the full makeContrasts() will be:
    ## makeContrasts(bob=bob, alice=alice, bob_vs_alice=(bob)-(alice), levels=design)
    ## The parentheses in this case are superfluous, but they remind me that I am finally
    ## doing some match, and also they remind me that we can do much more complex things like:
    ## ((bob-alice)-(jane-alice)) or whatever....
    all_pairwise <- list()
    identity_names <- names(identities)
    lenminus <- length(identities) - 1
    for (c in 1:lenminus) {
        c_name <- names(identities[c])
        nextc <- c + 1
        for (d in nextc:length(identities)) {
            d_name <- names(identities[d])
            minus_string <- paste(d_name, "_vs_", c_name, sep="")
            exprs_string <- paste(minus_string, "=", d_name, "-", c_name, ",", sep="")
            all_pairwise[minus_string] <- exprs_string
        }
    }
    ## At this point, I have strings which represent the definition of every
    ## sample condition as well as strings which represent every possible
    ## B-A where B comes somewhere after A in the model matrix.
    ## The goal now is to create the variables in the R environment
    ## and add them to makeContrasts()
    if (isTRUE(do_identities)) {
        eval_strings <- append(eval_strings, identities)
    }
    if (isTRUE(do_pairwise)) {
        eval_strings <- append(eval_strings, all_pairwise)
    }
    eval_names <- names(eval_strings)

    if (!is.null(extra_contrasts)) {
        extra_eval_strings <- strsplit(extra_contrasts, ",")
        extra_eval_names <- extra_eval_strings
        extra_eval_names <- stringi::stri_replace_all_regex(extra_eval_strings[[1]],
                                                            "^(\\s*)(\\w+)=.*$", "$2")
        for (i in 1:length(extra_eval_strings)) {
            new_name <- extra_eval_names[[i]]
            extra_contrast <- paste0(extra_eval_strings[[i]], ", ")
            eval_strings <- append(eval_strings, extra_contrast)
            eval_names <- append(eval_names, new_name)
            all_pairwise[new_name] <- extra_contrast
        }
        names(eval_strings) <- eval_names
    }
    ## for (f in 1:length(eval_strings)) {
    ##     eval_name = names(eval_strings[f])
    ##     message(paste("Setting ", eval_name, " with expression:<<", eval_strings[f], ">>", sep=""))
    ##     eval(parse(text=as.character(eval_strings[f])))
    ## }
    ## Now we have bob=(somestuff) in memory in R's environment
    ## Add them to makeContrasts()
    contrast_string <- paste0("all_pairwise_contrasts = limma::makeContrasts(")
    for (f in 1:length(eval_strings)) {
        ## eval_name = names(eval_strings[f])
        eval_string <- paste0(eval_strings[f])
        contrast_string <- paste(contrast_string, eval_string, sep=" ")
    }
    ## The final element of makeContrasts() is the design from voom()
    contrast_string <- paste0(contrast_string, " levels=model)")
    eval(parse(text=contrast_string))
    ## I like to change the column names of the contrasts because by default
    ## they are kind of obnoxious and too long to type

    colnames(all_pairwise_contrasts) <- eval_names
    result <- list(
        "all_pairwise_contrasts" = all_pairwise_contrasts,
        "identities" = identities,
        "identity_names" = identity_names,
        "all_pairwise" = all_pairwise,
        "contrast_string" = contrast_string,
        "names" = eval_names)
    return(result)
}

#' Remove multicopy genes from up/down gene expression lists.
#'
#' In our parasite data, there are a few gene types which are
#' consistently obnoxious.  Multi-gene families primarily where the
#' coding sequences are divergent, but the UTRs nearly identical.  For
#' these genes, our sequence based removal methods fail and so this
#' just excludes them by name.
#'
#' @param de_list  List of sets of genes deemed significantly
#'  up/down with a column expressing approximate count numbers.
#' @param max_copies  Keep only those genes with <= n putative
#'  copies.
#' @param use_files  Use a set of sequence alignments to define the copy numbers?
#' @param semantic  Set of strings with gene names to exclude.
#' @param semantic_column  Column in the DE table used to find the
#'  semantic strings for removal.
#' @return Smaller list of up/down genes.
#' @export
semantic_copynumber_filter <- function(de_list, max_copies=2, use_files=FALSE,
                                       semantic=c('mucin','sialidase','RHS','MASP','DGF','GP63'),
                                       semantic_column='1.tooltip') {
    removed_up <- list()
    removed_down <- list()
    table_type <- "significance"
    if (!is.null(de_list[["data"]])) {
        table_type <- "combined"
    }

    table_list <- NULL
    up_to_down <- 0
    if (table_type == "combined") {
        table_list <- de_list[["data"]]
    } else {
        table_list <- c(de_list[["ups"]], de_list[["downs"]])
        up_to_down <- length(de_list[["ups"]])
    }

    removed <- list()
    for (count in 1:length(table_list)) {
        tab <- table_list[[count]]
        table_name <- names(table_list)[[count]]
        removed[[table_name]] <- list()
        message(paste0("Working on ", table_name))
        if (isTRUE(use_files)) {
            file <- ""
            if (table_type == "combined") {
                file <- paste0("singletons/gene_counts/", table_name, ".fasta.out.count")
            } else {
                file <- paste0("singletons/gene_counts/up_", table_name, ".fasta.out.count")
            }
            tmpdf <- try(read.table(file), silent=TRUE)
            if (class(tmpdf) == 'data.frame') {
                colnames(tmpdf) = c("ID", "members")
                tab <- merge(tab, tmpdf, by.x="row.names", by.y="ID")
                rownames(tab) <- tab[["Row.names"]]
                tab <- tab[, -1, drop=FALSE]
                tab <- tab[count <= max_copies, ]
            }
        }  ## End using empirically defined groups of multi-gene families.
        for (string in semantic) {
            idx <- grep(pattern=string, x=tab[, semantic_column])
            num_removed <- length(idx)
            removed[[table_name]][[string]] <- num_removed
            if (num_removed > 0) {
                tab <- tab[-idx, ]
                message(paste0("Removing entries with string ", string,
                               ", found ", num_removed, "; table has ", nrow(tab),  " rows left."))
            } else {
                message("Found no entries of type ", string, ".")
            }
        }
        if (table_type == "combined") {
            de_list[["data"]][[count]] <- tab
        } else {
            if (count <= up_to_down) {
                de_list[["ups"]][[count]] <- tab
            } else {
                de_list[["downs"]][[count]] <- tab
            }
        }
    }
    if (table_type == "significance") {
        new_removed <- list()
        for (count in 1:length(removed)) {
            old_name <- names(removed)[[count]]
            new_name <- NULL
            if (count <= up_to_down) {
                new_name <- paste0("up_", old_name)
            } else {
                new_name <- paste0("down_", old_name)
            }
            new_removed[[new_name]] <- removed[[old_name]]
        }
        removed <- new_removed
        rm(new_removed)
    }
    de_list[["removed"]] <- removed
    return(de_list)
}

## EOF
