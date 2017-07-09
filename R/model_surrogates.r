## Going to try and recapitulate the analyses found at:
## https://github.com/jtleek/svaseq/blob/master/recount.Rmd
## and use the results to attempt to more completely understand batch effects in our data

#' Extract some surrogate estimations from a raw data set using sva, ruv, and/or pca.
#'
#' This applies the methodologies very nicely explained by Jeff Leek at
#' https://github.com/jtleek/svaseq/blob/master/recount.Rmd
#' and attempts to use them to acquire estimates which may be applied to an experimental model
#' by either EdgeR, DESeq2, or limma.  In addition, it modifies the count tables using these
#' estimates so that one may play with the modified counts and view the changes (with PCA or heatmaps
#' or whatever).  Finally, it prints a couple of the plots shown by Leek in his document.
#' In other words, this is entirely derivative of someone much smarter than me.
#'
#' @param input  Expt or data frame to manipulate.
#' @param design  If the data is not an expt, provide experimental design here.
#' @param estimate_type  One of: sva_supervised, sva_unsupervised, ruv_empirical, ruv_supervised,
#'  ruv_residuals, or pca.
#' @param surrogates  Choose a method for getting the number of surrogates, be or leek, or a number.
#' @param expt_state  Current state of the expt object (to check for log2, cpm, etc)
#' @param ... Parameters fed to arglist.
#' @return List including the adjustments for a model matrix, a modified count table, and 3 plots of
#'  the known batch, surrogates, and batch/surrogate.
#' @seealso \pkg{Biobase} \pkg{sva} \pkg{EDASeq} \pkg{RUVseq} \pkg{edgeR}
#' @export
get_model_adjust <- function(input, design=NULL, estimate_type="sva",
                             surrogates="be", expt_state=NULL, ...) {
    arglist <- list(...)
    my_design <- NULL
    my_data <- NULL
    transform_state <- "raw"
    log_data <- NULL
    log2_mtrx <- NULL
    base10_data <- NULL
    base10_mtrx <- NULL
    ## Gather all the likely pieces we can use
    ## Without the following requireNamespace(ruv)
    ## we get an error 'unable to find an inherited method for function RUVr'
    ruv_loaded <- try(require(package="ruv", quietly=TRUE))
    ## In one test, this seems to have been enough, but in another, perhaps not.

    filter <- "raw"
    if (!is.null(arglist[["filter"]])) {
        filter <- arglist[["filter"]]
    }
    convert <- "cpm"
    if (!is.null(arglist[["convert"]])) {
        convert <- arglist[["convert"]]
    }
    if (class(input) == "expt") {
        ## Gather all the likely pieces we can use
        my_design <- input[["design"]]
        my_data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))
        transform_state <- input[["state"]][["transform"]]
        base10_mtrx <- as.matrix(my_data)
        log_mtrx <- as.matrix(my_data)
        if (transform_state == "raw") {
            ## I think this was the cause of some problems.  The order of operations performed here
            ## was imperfect and could potentially lead to multiple different matrix sizes.
            base10_data <- sm(normalize_expt(input, convert=convert, filter=filter, thresh=1))
            base10_mtrx <- Biobase::exprs(base10_data[["expressionset"]])
            log_data <- sm(normalize_expt(base10_data, transform="log2"))
            log2_mtrx <- Biobase::exprs(log_data[["expressionset"]])
            rm(log_data)
            rm(base10_data)
        } else {
            log2_mtrx <- as.matrix(my_data)
            base10_mtrx <- as.matrix(2 ^ my_data) - 1
        }
    } else {
        if (is.null(design)) {
            stop("If an expt is not passed, then design _must_ be.")
        }
        message("Not able to discern the state of the data.")
        message("Going to use a simplistic metric to guess if it is log scale.")
        my_design <- design
        if (max(input) > 100) {
            transform_state <- "raw"
        } else {
            transform_state <- "log2"
        }
        my_data <- input
        base10_mtrx <- as.matrix(my_data)
        log_mtrx <- as.matrix(my_data)
        if (transform_state == "raw") {
            log_data <- sm(hpgl_norm(input, convert="cpm", transform="log2", filter=filter, thresh=1))
            log2_mtrx <- as.matrix(log_data[["count_table"]])
            ## base10_data <- sm(hpgl_norm(data, convert="cpm", filter=filter, thresh=1))
            ## base10_mtrx <- as.matrix(base10_data[["count_table"]])
            base10_mtrx <- (2 ^ log2_mtrx) - 1
            rm(log_data)
            ## rm(base10_data)
        } else {
            log2_mtrx <- as.matrix(input)
            base10_mtrx <- as.matrix(2 ^ input) - 1
        }
    }

    conditions <- droplevels(as.factor(my_design[["condition"]]))
    batches <- droplevels(as.factor(my_design[["batch"]]))
    conditional_model <- model.matrix(~ conditions, data=my_design)
    null_model <- conditional_model[, 1]
    chosen_surrogates <- 1
    if (is.null(surrogates)) {
        message("No estimate nor method to find surrogates was provided. Assuming you want 1 surrogate variable.")
    } else {
        if (class(surrogates) == "character") {
            ## num.sv assumes the log scale.
            if (surrogates != "be" & surrogates != "leek") {
                message("A string was provided, but it was neither 'be' nor 'leek', assuming 'be'.")
                chosen_surrogates <- sm(sva::num.sv(dat=log2_mtrx, mod=conditional_model))
            } else {
                chosen_surrogates <- sm(sva::num.sv(dat=log2_mtrx,
                                                    mod=conditional_model, method=surrogates))
            }
            message(paste0("The ", surrogates, " method chose ", chosen_surrogates, " surrogate variable(s)."))
        } else if (class(surrogates) == "numeric") {
            message(paste0("A specific number of surrogate variables was chosen: ", surrogates, "."))
            chosen_surrogates <- surrogates
        }
    }

    ## empirical controls can take either log or base 10 scale depending on 'control_type'
    control_type <- "norm"
    control_likelihoods <- NULL
    if (control_type == "norm") {
        control_likelihoods <- try(sm(sva::empirical.controls(dat=log2_mtrx,
                                                              mod=conditional_model,
                                                              mod0=null_model,
                                                              n.sv=chosen_surrogates,
                                                              type=control_type)))
    } else {
        control_likelihoods <- try(sm(sva::empirical.controls(dat=base10_mtrx,
                                                              mod=conditional_model,
                                                              mod0=null_model,
                                                              n.sv=chosen_surrogates,
                                                              type=control_type)))
    }
    if (class(control_likelihoods) == "try-error") {
        message("The most likely error in sva::empirical.controls() is a call to density in irwsva.build.
Setting control_likelihoods to zero and using unsupervised sva.")
        warning("It is highly likely that the underlying reason for this error is too many 0's in
the dataset, please try doing a filtering of the data and retry.")
        control_likelihoods <- 0
    }
    if (sum(control_likelihoods) == 0) {
        if (estimate_type == "sva_supervised") {
            message("Unable to perform supervised estimations, changing to unsupervised_sva.")
            estimate_type <- "sva_unsupervised"
        } else if (estimate_type == "ruv_supervised") {
            message("Unable to perform supervised estimations, changing to empirical_ruv.")
            estimate_type <- "ruv_empirical"
        }
    }

    ## I use 'sva' as shorthand fairly often
    if (estimate_type == "sva") {
        estimate_type <- "sva_unsupervised"
        message("Estimate type 'sva' is shorthand for 'sva_unsupervised'.")
        message("Other sva options include: sva_supervised and svaseq.")
    }
    if (estimate_type == "ruv") {
        estimate_type <- "ruv_empirical"
        message("Estimate type 'ruv' is shorthand for 'ruv_empirical'.")
        message("Other ruv options include: ruv_residual and ruv_supervised.")
    }

    surrogate_result <- NULL
    model_adjust <- NULL
    adjusted_counts <- NULL
    type_color <- NULL
    returned_counts <- NULL
    if (estimate_type == "sva_supervised") {
        message("Attempting sva supervised surrogate estimation.")
        type_color <- "red"
        supervised_sva <- sm(sva::ssva(log2_mtrx,
                                       controls=control_likelihoods,
                                       n.sv=chosen_surrogates))
        model_adjust <- as.matrix(supervised_sva[["sv"]])
        surrogate_result <- supervised_sva
    } else if (estimate_type == "fsva") {
        message("Attempting fsva surrogate estimation.")
        type_color <- "darkred"
        sva_object <- sm(sva::sva(log2_mtrx, conditional_model,
                                  null_model, n.sv=chosen_surrogates))
        fsva_result <- sva::fsva(log2_mtrx, conditional_model,
                                 sva_object, newdat=as.matrix(log2_mtrx),
                                 method="exact")
        ##method="exact"))
        model_adjust <- as.matrix(fsva_result[["newsv"]])
        surrogate_result <- fsva_result
    } else if (estimate_type == "svaseq") {
        message("This ignores the surrogates parameter and uses the be method to estimate surrogates.")
        type_color <- "dodgerblue"
        svaseq_result <- sm(sva::svaseq(base10_mtrx,
                                        conditional_model,
                                        null_model))
        surrogate_result <- svaseq_result
        model_adjust <- as.matrix(svaseq_result[["sv"]])
    } else if (estimate_type == "sva_unsupervised") {
        message("Attempting sva unsupervised surrogate estimation.")
        type_color <- "blue"
        if (min(rowSums(base10_mtrx)) == 0) {
            warning("sva will likely fail because some rowSums are 0.")
        }
        unsupervised_sva_batch <- sm(sva::sva(log2_mtrx,
                                              conditional_model,
                                              null_model,
                                              n.sv=chosen_surrogates))
        surrogate_result <- unsupervised_sva_batch
        model_adjust <- as.matrix(unsupervised_sva_batch[["sv"]])
    } else if (estimate_type == "pca") {
        message("Attempting pca surrogate estimation.")
        type_color <- "green"
        data_vs_means <- as.matrix(log2_mtrx - rowMeans(log2_mtrx))
        svd_result <- corpcor::fast.svd(data_vs_means)
        surrogate_result <- svd_result
        model_adjust <- as.matrix(svd_result[["v"]][, 1:chosen_surrogates])
    } else if (estimate_type == "ruv_supervised") {
        message("Attempting ruvseq supervised surrogate estimation.")
        type_color <- "black"
        ## Re-calculating the numer of surrogates with this modified data.
        surrogate_estimate <- sm(sva::num.sv(dat=log2_mtrx, mod=conditional_model))
        if (min(rowSums(base10_mtrx)) == 0) {
            warning("empirical.controls will likely fail because some rows are all 0.")
        }
        control_likelihoods <- sm(sva::empirical.controls(dat=log2_mtrx,
                                                          mod=conditional_model,
                                                          mod0=null_model,
                                                          n.sv=surrogate_estimate))
        ##ruv_result <- RUVSeq::RUVg(round(base10_mtrx),
        ##                           cIdx=as.logical(control_likelihoods),
        ##                           k=surrogate_estimate)
        ruv_result <- RUVSeq::RUVg(round(base10_mtrx),
                                   k=surrogate_estimate,
                                   cIdx=as.logical(control_likelihoods))

        surrogate_result <- ruv_result
        returned_counts <- ruv_result[["normalizedCounts"]]
        model_adjust <- as.matrix(ruv_result[["W"]])
    } else if (estimate_type == "ruv_residuals") {
        message("Attempting ruvseq residual surrogate estimation.")
        type_color <- "purple"
        ## Use RUVSeq and residuals
        ruv_input <- edgeR::DGEList(counts=base10_mtrx, group=conditions)
        norm <- edgeR::calcNormFactors(ruv_input)
        ruv_input <- try(edgeR::estimateDisp(norm, design=conditional_model, robust=TRUE))
        ruv_fit <- edgeR::glmFit(ruv_input, conditional_model)
        ruv_res <- residuals(ruv_fit, type="deviance")
        ruv_normalized <- EDASeq::betweenLaneNormalization(base10_mtrx, which="upper")
        ## This also gets mad if you pass it a df and not matrix
        controls <- rep(TRUE, dim(base10_mtrx)[1])
        ruv_result <- RUVSeq::RUVr(ruv_normalized, controls, k=chosen_surrogates, ruv_res)
        model_adjust <- as.matrix(ruv_result[["W"]])
    } else if (estimate_type == "ruv_empirical") {
        message("Attempting ruvseq empirical surrogate estimation.")
        type_color <- "orange"
        ruv_input <- edgeR::DGEList(counts=base10_mtrx, group=conditions)
        ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method="upperquartile")
        ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, conditional_model)
        ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, conditional_model)
        ruv_fit <- edgeR::glmFit(ruv_input_tag, conditional_model)
        ## Use RUVSeq with empirical controls
        ## The previous instance of ruv_input should work here, and the ruv_input_norm
        ## Ditto for _glm and _tag, and indeed ruv_fit
        ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
        ruv_lrt <- edgeR::glmLRT(ruv_fit, coef=2)
        ruv_control_table <- ruv_lrt[["table"]]
        ranked <- as.numeric(rank(ruv_control_table[["LR"]]))
        bottom_third <- (summary(ranked)[[2]] + summary(ranked)[[3]]) / 2
        ruv_controls <- ranked <= bottom_third  ## what is going on here?!
        ## ruv_controls = rank(ruv_control_table$LR) <= 400  ## some data sets fail with 400 hard-set
        ruv_result <- RUVSeq::RUVg(round(base10_mtrx), ruv_controls, k=chosen_surrogates)
        surrogate_result <- ruv_result
        model_adjust <- as.matrix(ruv_result[["W"]])
    } else {
        type_color <- "grey"
        ## If given nothing to work with, use supervised sva
        message(paste0("Did not understand ", estimate_type, ", assuming supervised sva."))
        supervised_sva <- sva::svaseq(base10_mtrx,
                                      conditional_model,
                                      null_model,
                                      controls=control_likelihoods,
                                      n.sv=chosen_surrogates)
        model_adjust <- as.matrix(supervised_sva[["sv"]])
        surrogate_result <- supervised_sva
    }

    ## This is the old code, potentially a source of my recent error, but I think it probably is not.
    ## new_model <- cbind(conditional_model, model_adjust)
    ## data_modifier <- solve(t(new_model) %*% new_model) %*% t(new_model)
    ## transformation <- (data_modifier %*% t(mtrx))
    ## conds <- ncol(conditional_model)
    ## new_counts <- mtrx - t(as.matrix(new_model[, -c(1:conds)]) %*% transformation[-c(1:conds), ])
    ## counts_from_surrogates currently resides in normalize_batch.R
    new_counts <- counts_from_surrogates(base10_mtrx, model_adjust, design=my_design)
    plotbatch <- as.integer(batches)
    plotcond <- as.numeric(conditions)
    x_marks <- 1:length(colnames(data))

    surrogate_plots <- NULL
    if (class(input) == "expt") {
        surrogate_plots <- plot_batchsv(input, model_adjust)
    }

    ret <- list(
        "surrogate_result" = surrogate_result,
        "null_model" = null_model,
        "model_adjust" = model_adjust,
        "new_counts" = new_counts,
        "sample_factor" = surrogate_plots[["sample_factor"]],
        "factor_svs" = surrogate_plots[["factor_svs"]],
        "svs_sample" = surrogate_plots[["svs_sample"]])
    return(ret)
}

#' Perform a comparison of the surrogate estimators demonstrated by Jeff Leek.
#'
#' This is entirely derivative, but seeks to provide similar estimates for one's own actual data
#' and catch corner cases not taken into account in that document (for example if the estimators
#' don't converge on a surrogate variable). This will attempt each of the surrogate estimators
#' described by Leek: pca, sva supervised, sva unsupervised, ruv supervised, ruv residuals, ruv
#' empirical. Upon completion it will perform the same limma expression analysis and plot the ranked
#' t statistics as well as a correlation plot making use of the extracted estimators against
#' condition/batch/whatever else. Finally, it does the same ranking plot against a linear fitting
#' Leek performed and returns the whole pile of information as a list.
#'
#' @param expt Experiment containing a design and other information.
#' @param extra_factors Character list of extra factors which may be included in the final plot of
#'  the data.
#' @param do_catplots Include the catplots?  They don't make a lot of sense yet, so probably no.
#' @param surrogates  Use 'be' or 'leek' surrogate estimates, or choose a number.
#' @return List of the results.
#' @seealso \code{\link{get_model_adjust}}
#' @export
compare_surrogate_estimates <- function(expt, extra_factors=NULL,
                                        do_catplots=FALSE, surrogates="be") {
    design <- expt[["design"]]
    pca_plots <- list()
    pca_plots[["null"]] <- plot_pca(expt)[["plot"]]

    pca_adjust <- get_model_adjust(expt, estimate_type="pca", surrogates=surrogates)
    pca_plots[["pca"]] <- plot_pca(pca_adjust[["new_counts"]],
                                   design=design,
                                   plot_colors=expt[["colors"]])[["plot"]]

    sva_supervised <- get_model_adjust(expt, estimate_type="sva_supervised", surrogates=surrogates)
    pca_plots[["svasup"]] <- plot_pca(sva_supervised[["new_counts"]],
                                      design=design,
                                      plot_colors=expt[["colors"]])[["plot"]]

    sva_unsupervised <- get_model_adjust(expt, estimate_type="sva_unsupervised", surrogates=surrogates)
    pca_plots[["svaunsup"]] <- plot_pca(sva_unsupervised[["new_counts"]],
                                        design=design,
                                        plot_colors=expt[["colors"]])[["plot"]]

    ruv_supervised <- get_model_adjust(expt, estimate_type="ruv_supervised", surrogates=surrogates)
    pca_plots[["ruvsup"]] <- plot_pca(ruv_supervised[["new_counts"]],
                                      design=design,
                                      plot_colors=expt[["colors"]])[["plot"]]

    ruv_residuals <- get_model_adjust(expt, estimate_type="ruv_residuals", surrogates=surrogates)
    pca_plots[["ruvresid"]] <- plot_pca(ruv_residuals[["new_counts"]],
                                        design=design,
                                        plot_colors=expt[["colors"]])[["plot"]]

    ruv_empirical <- get_model_adjust(expt, estimate_type="ruv_empirical", surrogates=surrogates)
    pca_plots[["ruvemp"]] <- plot_pca(ruv_empirical[["new_counts"]],
                                      design=design,
                                      plot_colors=expt[["colors"]])[["plot"]]

    first_svs <- data.frame(
        "condition" = as.numeric(as.factor(expt[["conditions"]])),
        "batch" = as.numeric(as.factor(expt[["batches"]])),
        "pca_adjust" = pca_adjust[["model_adjust"]][, 1],
        "sva_supervised" = sva_supervised[["model_adjust"]][, 1],
        "sva_unsupervised" = sva_unsupervised[["model_adjust"]][, 1],
        "ruv_supervised" = ruv_supervised[["model_adjust"]][, 1],
        "ruv_residuals" = ruv_residuals[["model_adjust"]][, 1],
        "ruv_empirical" = ruv_empirical[["model_adjust"]][, 1])
    batch_adjustments <- list(
        "condition" = as.factor(expt[["conditions"]]),
        "batch" = as.factor(expt[["batches"]]),
        "pca_adjust" = pca_adjust[["model_adjust"]],
        "sva_supervised" = sva_supervised[["model_adjust"]],
        "sva_unsupervised" = sva_unsupervised[["model_adjust"]],
        "ruv_supervised" = ruv_supervised[["model_adjust"]],
        "ruv_residuals" = ruv_residuals[["model_adjust"]],
        "ruv_empirical" = ruv_empirical[["model_adjust"]])
    batch_names <- c("condition", "batch", "pca", "sva_sup", "sva_unsup",
                     "ruv_sup", "ruv_resid", "ruv_emp")
    silly <- testthat::compare(batch_names, batch_names)  ## I want to try something silly
    if (!is.null(extra_factors)) {
        for (fact in extra_factors) {
            if (!is.null(design[, fact])) {
                batch_names <- append(x=batch_names, values=fact)
                first_svs[[fact]] <- as.numeric(as.factor(design[, fact]))
                batch_adjustments[[fact]] <- as.numeric(as.factor(design[, fact]))
            }
        }
    }
    correlations <- cor(first_svs)
    corrplot::corrplot(correlations, method="ellipse", type="lower", tl.pos="d")
    ret_plot <- grDevices::recordPlot()

    adjustments <- c("+ batch_adjustments$batch", "+ batch_adjustments$pca",
                     "+ batch_adjustments$sva_sup", "+ batch_adjustments$sva_unsup",
                     "+ batch_adjustments$ruv_sup", "+ batch_adjustments$ruv_resid",
                     "+ batch_adjustments$ruv_emp")
    adjust_names <- c("null", "batch", "pca", "sva_sup", "sva_unsup",
                      "ruv_sup", "ruv_resid", "ruv_emp")
    starter <- edgeR::DGEList(counts=Biobase::exprs(expt[["expressionset"]]))
    norm_start <- edgeR::calcNormFactors(starter)
    catplots <- vector("list", length(adjustments) + 1)  ## add 1 for a null adjustment
    names(catplots) <- adjust_names
    tstats <- list()

    ## First do a null adjust
    adjust <- ""
    counter <- 1
    num_adjust <- length(adjustments)
    message(paste0(counter, "/", num_adjust + 1, ": Performing lmFit(data) etc. with null in the model."))
    modified_formula <- as.formula(paste0("~ condition ", adjust))
    limma_design <- model.matrix(modified_formula, data=design)
    voom_result <- limma::voom(norm_start, limma_design, plot=FALSE)
    limma_fit <- limma::lmFit(voom_result, limma_design)
    modified_fit <- limma::eBayes(limma_fit)
    tstats[["null"]] <- abs(modified_fit[["t"]][, 2])
    ##names(tstats[["null"]]) <- as.character(1:dim(data)[1])
    ## This needs to be redone to take into account how I organized the adjustments!!!
    num_adjust <- length(adjustments)
    oldpar <- par(mar=c(5, 5, 5, 5))
    for (adjust in adjustments) {
        counter <- counter + 1
        message(paste0(counter, "/", num_adjust + 1, ": Performing lmFit(data) etc. with ", adjust, " in the model."))
        modified_formula <- as.formula(paste0("~ condition ", adjust))
        limma_design <- model.matrix(modified_formula, data=design)
        voom_result <- limma::voom(norm_start, limma_design, plot=FALSE)
        limma_fit <- limma::lmFit(voom_result, limma_design)
        modified_fit <- limma::eBayes(limma_fit)
        tstats[[adjust]] <- abs(modified_fit[["t"]][, 2])
        ##names(tstats[[counter]]) <- as.character(1:dim(data)[1])
        catplot_together <- NULL
        if (isTRUE(do_catplots)) {
            if (!isTRUE("ffpe" %in% .packages(all.available=TRUE))) {
                ## ffpe has some requirements which do not install all the time.
                require.auto("ffpe")
            }
            if (isTRUE("ffpe" %in% .packages(all.available=TRUE))) {
                catplots[[counter]] <- ffpe::CATplot(-rank(tstats[[adjust]]),
                                                     -rank(tstats[["null"]]),
                                                     maxrank=1000,
                                                     make.plot=TRUE)
            } else {
                catplots[[counter]] <- NULL
            }
            plot(catplots[["pca"]], ylim=c(0, 1), col="black",
                 lwd=3, type="l", xlab="Rank",
                 ylab="Concordance between study and different methods.")
            lines(catplots[["sva_sup"]], col="red", lwd=3, lty=2)
            lines(catplots[["sva_unsup"]], col="blue", lwd=3)
            lines(catplots[["ruv_sup"]], col="green", lwd=3, lty=3)
            lines(catplots[["ruv_resid"]], col="orange", lwd=3)
            lines(catplots[["ruv_emp"]], col="purple", lwd=3)
            legend(200, 0.5, legend=c("some stuff about methods used."), lty=c(1, 2, 1, 3, 1), lwd=3)
            catplot_together <- grDevices::recordPlot()
            newpar <- par(oldpar)
        } ## End checking whether to do catplots
    }

    ret <- list(
        "pca_adjust" = pca_adjust,
        "sva_supervised_adjust" = sva_supervised,
        "sva_unsupervised_adjust" = sva_unsupervised,
        "ruv_supervised_adjust" = ruv_supervised,
        "ruv_residual_adjust" = ruv_residuals,
        "ruv_empirical_adjust" = ruv_empirical,
        "adjustments" = batch_adjustments,
        "correlations" = correlations,
        "plot" = ret_plot,
        "pca_plots" = pca_plots,
        "catplots" = catplot_together)
    return(ret)
}

## EOF
