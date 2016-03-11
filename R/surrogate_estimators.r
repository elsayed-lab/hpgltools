## Time-stamp: <Fri Mar 11 12:38:35 2016 Ashton Trey Belew (abelew@gmail.com)>

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
#' @param raw_expt a raw experiment object
#' @param estimate_type one of sva_supervised, sva_unsupervised, ruv_empirical, ruv_supervised, ruv_residuals, or pca
#' @return a list including the adjustments for a model matrix, a modified count table, and 3 plots of the known batch, surrogates, and batch/surrogate.
#' @export
get_model_adjust <- function(raw_expt, estimate_type="sva_supervised", ...) {
    arglist <- list(...)
    chosen_surrogates <- NULL
    if (!is.null(arglist$surrogates)) {
        chosen_surrogates <- arglist$surrogates
    }
    ## Gather all the likely pieces we can use
    start_low <- suppressMessages(normalize_expt(raw_expt, filter_low=TRUE))
    design <- start_low$design
    data <- as.data.frame(Biobase::exprs(start_low$expressionset))
    mtrx <- as.matrix(data)
    l2_data <- Biobase::exprs(suppressMessages(normalize_expt(start_low, transform="log2")$expressionset))
    conditions <- as.factor(design[, "condition"])
    batches <- as.factor(design[, "batch"])
    conditional_model <- model.matrix(~ conditions, data=data)
    null_model <- conditional_model[, 1]
    be_surrogate_estimate <- sva::num.sv(dat=mtrx, mod=conditional_model, method="be")
    leek_surrogate_estimate <- sva::num.sv(dat=mtrx, mod=conditional_model, method="leek")
    chosen_estimate <- 1
    if (is.null(chosen_surrogates)) {
        if (as.numeric(be_surrogate_estimate) > 0 & as.numeric(leek_surrogate_estimate) > 0) {
            chosen_estimate <- be_surrogate_estimate
        } else if (be_surrogate_estimate > 0) {
            chosen_estimate <- be_surrogate_estimate
        } else {
            chosen_estimate <- leek_surrogate_estimate
        }
    }
    if (chosen_estimate <= 4 | chosen_estimate >= 1) {
        chosen_surrogates = chosen_estimate
    }

    control_likelihoods <- try(sva::empirical.controls(dat=mtrx, mod=conditional_model, mod0=null_model, n.sv=chosen_surrogates), silent=TRUE)
    if (class(control_likelihoods) == 'try-error') {
        control_likelihoods = 0
    }
    if (sum(control_likelihoods) == 0) {
        if (estimate_type == "sva_supervised") {
            message("Unable to perform supervised estimations, changing to unsupervised_sva.")
            estimate_type <- "sva_supervised"
        } else if (type == "ruv_supervised") {
            message("Unable to perform supervised estimations, changing to empirical_ruv.")
            estimate_type <- "ruv_empirical"
        }
    }

    model_adjust <- NULL
    adjusted_counts <- NULL
    type_color <- NULL
    if (estimate_type == "sva_supervised") {
        type_color <- "red"
        supervised_sva <- sva::svaseq(mtrx, conditional_model, null_model, controls=control_likelihoods)
        model_adjust <- supervised_sva$sv
    } else if (estimate_type == "sva_unsupervised") {
        type_color <- "blue"
        unsupervised_sva_batch <- sva::svaseq(mtrx, conditional_model, null_model)
        model_adjust <- unsupervised_sva_batch$sv
    } else if (estimate_type == "pca") {
        type_color <- "green"
        model_adjust <- corpcor::fast.svd(l2_data - rowMeans(l2_data))$v[, 1]
    } else if (estimate_type == "ruv_supervised") {
        type_color <- "black"
        surrogate_estimate <- sva::num.sv(dat=mtrx, mod=conditional_model)
        control_likelihoods <- sva::empirical.controls(dat=mtrx, mod=conditional_model, mod0=null_model, n.sv=surrogate_estimate)
        model_adjust <- RUVSeq::RUVg(mtrx, cIdx=as.logical(control_likelihoods), k=1)$W
    } else if (estimate_type == "ruv_residuals") {
        type_color <- "purple"
        ## Use RUVSeq and residuals
        conditional_model <- model.matrix(~ conditions)
        ruv_input <- edgeR::DGEList(counts=data, group=conditions)
        ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method="upperquartile")
        ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, ruv_design)
        ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, ruv_design)
        ruv_fit <- edgeR::glmFit(ruv_input_tag, ruv_design)
        ruv_res <- residuals(ruv_fit, type="deviance")
        ruv_normalized <- betweenLaneNormalization(mtrx, which="upper")  ## This also gets mad if you pass it a df and not matrix
        controls <- rep(TRUE, dim(data)[1])
        model_adjust <- RUVSeq::RUVr(ruv_normalized, controls, k=1, ruv_res)$W
    } else if (estimate_type == "ruv_empirical") {
        type_color <- "orange"
        conditional_model <- model.matrix(~ conditions)
        ruv_input <- edgeR::DGEList(counts=data, group=conditions)
        ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method="upperquartile")
        ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, ruv_design)
        ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, ruv_design)
        ruv_fit <- edgeR::glmFit(ruv_input_tag, ruv_design)
        ## Use RUVSeq with empirical controls
        ## The previous instance of ruv_input should work here, and the ruv_input_norm
        ## Ditto for _glm and _tag, and indeed ruv_fit
        ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
        ruv_lrt <- glmLRT(ruv_fit, coef=2)
        ruv_controls = rank(ruv_lrt$table$LR) <= 400  ## what is going on here?!
        model_adjust <- RUVSeq::RUVg(mtrx, ruv_controls, k=1)$W
    } else {
        type_color <- "black"
        ## If given nothing to work with, use supervised sva
        message(paste0("Did not understand ", type, ", assuming supervised sva."))
        supervised_sva <- sva::svaseq(mtrx, conditional_model, null_model, controls=control_likelihoods)
        model_adjust <- supervised_sva$sv
    }
    plotbatch <- as.integer(batches)
    plotcond <- as.numeric(conditions)
    x_marks <- 1:length(colnames(data))
    plot(plotbatch, type="p", pch=19, col="black", main=paste0("Known batches by sample"), xaxt="n", yaxt="n", xlab="Sample", ylab="Known batch")
    axis(1, at=x_marks, cex.axis=0.75, las=2, labels=as.character(colnames(data)))
    axis(2, at=plotbatch, cex.axis=0.75, las=2, labels=as.character(batches))
    batch_plot <- grDevices::recordPlot()
    plot(as.numeric(model_adjust), type="p", pch=19, col=type_color,
         xaxt="n", xlab="Sample", ylab="Surrogate estimate", main=paste0("Surrogates estimated by ", estimate_type))
    axis(1, at=x_marks, cex.axis=0.75, las=2, labels=as.character(colnames(data)))
    surrogate_plot <- grDevices::recordPlot()
    plot(model_adjust ~ plotbatch, pch=19, col=type_color, main=paste0(estimate_type, " vs. known batches."))
    batch_vs_adjust_plot <- grDevices::recordPlot()
    new_model <- cbind(conditional_model, model_adjust)
    data_modifier <- solve(t(new_model) %*% new_model) %*% t(new_model)
    transformation <- (data_modifier %*% t(mtrx))
    conds <- ncol(conditional_model)
    new_counts <- mtrx - t(as.matrix(new_model[, -c(1:conds)]) %*% transformation[-c(1:conds), ])
    ret <- list("model_adjust" = model_adjust,
                "new_counts" = new_counts,
                "batch_plot" = batch_plot,
                "surrogate_plot" = surrogate_plot,
                "batch_vs_adjust" = batch_vs_adjust_plot)
    return(ret)
}

#' Perform a comparison of the surrogate estimators demonstrated by Jeff Leek.
#' Once again this is entirely derivative, but seeks to provide similar estimates for one's own actual data
#' and catch corner cases not taken into account in that document (for example if the estimators
#' don't converge on a surrogate variable)
#'
#' This will attempt each of the surrogate estimators described by Leek: pca, sva supervised,
#' sva unsupervised, ruv supervised, ruv residuals, ruv empirical. Upon completion it will perform
#' the same limma expression analysis and plot the ranked t statistics as well as a correlation plot
#' making use of the extracted estimators against condition/batch/whatever else.
#' Finally, it does the same ranking plot against a linear fitting Leek performed and returns the
#' whole pile of information as a list.
#'
#' @param expt an experiment containing a design and other information
#' @param a character list of extra factors which may be included in the final plot of the data
#' @return a list of toys
#' @export
compare_surrogate_estimates <- function(expt, extra_factors=NULL) {
    message("1/6: Attempting pca surrogate estimation.")
    pca_adjust <- get_model_adjust(expt, type="pca")
    message("2/6: Attempting sva supervised surrogate estimation.")
    sva_supervised <- get_model_adjust(expt, type="supervised_sva")
    message("3/6: Attempting sva unsupervised surrogate estimation.")
    sva_unsupervised <- get_model_adjust(expt, type="unsupervised_sva")
    message("4/6: Attempting ruv supervised surrogate estimation.")
    ruv_supervised <- get_model_adjust(expt, type="ruv_supervised")
    message("5/6: Attempting ruv residual surrogate estimation.")
    ruv_residuals <- get_model_adjust(expt, type="ruv_residuals")
    message("6/6: Attempting ruv empirical surrogate estimation.")
    ruv_empirical <- get_model_adjust(expt, type="ruv_empirical")

    batch_names <- c("condition","batch","pca","sva_sup","sva_unsup","ruv_sup","ruv_resid","ruv_emp")
    batch_adjustments <- cbind(as.factor(expt$conditions),
                               as.factor(expt$batches),
                               pca_adjust$model_adjust,
                               sva_supervised$model_adjust,
                               sva_unsupervised$model_adjust,
                               ruv_supervised$model_adjust,
                               ruv_residuals$model_adjust,
                               ruv_empirical$model_adjust)
    if (!is.null(extra_factors)) {
        for (fact in extra_factors) {
            if (!is.null(expt$design[, fact])) {
                batch_names <- append(x=batch_names, values=fact)
                batch_adjustments <- cbind(batch_adjustments, as.factor(expt$design[, fact]))
            }
        }
    }
    colnames(batch_adjustments) <- batch_names
    correlations <- cor(batch_adjustments)
    par(mar=c(5,5,5,5))
    corrplot::corrplot(correlations, method="ellipse", type="lower", tl.pos="d")
    ret_plot <- grDevices::recordPlot()


    adjustments <- c("", "+ batch", "+ pca_adjust", "+ sva_supervised", "+ sva_unsupervised",
                     "+ ruv_supervised", "+ ruv_residuals", "+ ruv_empirival")
    starter <- edgeR::DGEList(counts=data)
    norm_start <- edgeR::calcNormFactors(starter)
    catplots <- vector("list", length(adjustments))
    tstats <- vector("list", length(adjustments))
    counter <- 0
    for (adjust in adjustments) {
        counter <- counter + 1
        design <- model.matrix(as.formula(paste0("~ condition", adjust)))
        voom_result <- limma::voom(norm_start, design, plot=FALSE)
        fit <- limma::lmFit(voom_result, design)
        fit <- eBayes(fit)
        tstats[[counter]] <- abs(fit$t[, 2])
        names(tstats[[counter]]) <- as.character(1:dim(data)[1])
        catplots[[counter]] <- CATplot(-rank(tstats[[counter]]), -rank(tstats[[1]]), maxrank=1000, make.plot=FALSE)
    }

    plot(catplots[[1]], ylim=c(0, 1), col="black", lwd=3, type="l", ylab="Concordance between study and different methods.", xlab="Rank")
    lines(catplots[[2]], col="red", lwd=3, lty=2)
    lines(catplots[[3]], col="blue", lwd=3)
    lines(catplots[[4]], col="green", lwd=3, lty=3)
    lines(catplots[[5]], col="orange", lwd=3)
    legend(200, 0.5, legend=c("some stuff about methods used."), col="colors", lty=c(1,2,1,3,1), lwd=3)
    catplot_together <- grDevices::recordPlot()

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
        "catplots" = catplot_together)
    return(ret)
}
