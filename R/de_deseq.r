#' Plot out 2 coefficients with respect to one another from deseq2.
#'
#' It can be nice to see a plot of two coefficients from a deseq2 comparison with respect to one
#' another. This hopefully makes that easy.
#'
#' @param output Set of pairwise comparisons provided by deseq_pairwise().
#' @param x Name or number of the x-axis coefficient column to extract.
#' @param y Name or number of the y-axis coefficient column to extract.
#' @param gvis_filename Filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline Add a trendline to the gvis plot?
#' @param tooltip_data Dataframe of gene annotations to be used in the gvis plot.
#' @param base_url When plotting interactive plots, have link-outs to this base url.
#' @return Ggplot2 plot showing the relationship between the two coefficients.
#' @seealso \link{plot_linear_scatter} \link{deseq2_pairwise}
#' @examples
#' \dontrun{
#'  pretty = coefficient_scatter(deseq_data, x="wt", y="mut")
#' }
#' @export
deseq_coefficient_scatter <- function(output, toptable=NULL, x=1, y=2, ## gvis_filename="limma_scatter.html",
                                      gvis_filename=NULL, gvis_trendline=TRUE, z=1.5,
                                      tooltip_data=NULL, base_url=NULL,
                                      color_low="#DD0000", color_high="#7B9F35", ...) {
    arglist <- list(...)
    qlimit <- 0.1
    if (!is.null(arglist[["qlimit"]])) {
        qlimit <- arglist[["qlimit"]]
    }
    fc_column <- "limma_logfc"
    if (!is.null(arglist[["fc_column"]])) {
        fc_column <- arglist[["fc_column"]]
    }
    p_column <- "limma_adjp"
    if (!is.null(arglist[["p_column"]])) {
        p_column <- arglist[["p_column"]]
    }
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    thenames <- names(output[["coefficients"]])
    message("This can do comparisons among the following columns in the deseq2 result:")
    message(toString(thenames))
    xname <- ""
    yname <- ""
    if (is.numeric(x)) {
        xname <- thenames[[x]]
    } else {
        xname <- x
    }
    if (is.numeric(y)) {
        yname <- thenames[[y]]
    } else {
        yname <- y
    }
    message(paste0("Actually comparing ", xname, " and ", yname, "."))
    first_df <- output[["coefficients"]][[xname]]
    first_df[["delta"]] <- log2(first_df[["baseMean"]]) + first_df[["log2FoldChange"]]
    second_df <- output[["coefficients"]][[yname]]
    second_df[["delta"]] <- log2(second_df[["baseMean"]]) + second_df[["log2FoldChange"]]
    first_col <- first_df[, c("baseMean", "log2FoldChange", "delta")]
    colnames(first_col) <- c("mean.1", "fc.1", xname)
    second_col <- second_df[, c("baseMean", "log2FoldChange", "delta")]
    colnames(second_col) <- c("mean.2", "fc.2", yname)
    coefficient_df <- merge(first_col, second_col, by="row.names")
    rownames(coefficient_df) <- coefficient_df[["Row.names"]]
    coefficient_df <- coefficient_df[-1]
    coefficient_df <- coefficient_df[, c(xname, yname)]
    coefficient_df[is.na(coefficient_df)] <- 0
    maxvalue <- max(coefficient_df) + 1.0
    plot <- suppressMessages(plot_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename,
                                                 gvis_trendline=gvis_trendline, first=xname, second=yname,
                                                 tooltip_data=tooltip_data, base_url=base_url,
                                                 pretty_colors=FALSE, color_low=color_low, color_high=color_high))
    plot[["scatter"]] <- plot[["scatter"]] +
        ggplot2::scale_x_continuous(limits=c(0, maxvalue)) +
        ggplot2::scale_y_continuous(limits=c(0, maxvalue))
    if (!is.null(toptable)) {
        theplot <- plot[["scatter"]] + ggplot2::theme_bw()
        sig <- get_sig_genes(toptable, z=z, column=fc_column, p_column=p_column)
        sigup <- sig[["up_genes"]]
        sigdown <- sig[["down_genes"]]
        up_index <- rownames(coefficients) %in% rownames(sigup)
        down_index <- rownames(coefficients) %in% rownames(sigdown)
        up_df <- as.data.frame(coefficients[up_index, ])
        down_df <- as.data.frame(coefficients[down_index, ])
        colnames(up_df) <- c("first", "second")
        colnames(down_df) <- c("first", "second")
        theplot <- theplot +
            ggplot2::geom_point(data=up_df, colour=color_high) +
            ggplot2::geom_point(data=down_df, colour=color_low)
        plot[["scatter"]] <- theplot
    }
    plot[["df"]] <- coefficient_df
    return(plot)
}



#' deseq_pairwise()  Because I can't be trusted to remember '2'.
#'
#' This calls deseq2_pairwise(...) because I am determined to forget typing deseq2.
#'
#' @param ... I like cats.
#' @return stuff deseq2_pairwise results.
#' @seealso \link{deseq2_pairwise}
#' @export
deseq_pairwise <- function(...) {
    message("Hey you, use deseq2 pairwise.")
    deseq2_pairwise(...)
}

#' Set up model matrices contrasts and do pairwise comparisons of all conditions using DESeq2.
#'
#' Invoking DESeq2 is confusing, this should help.
#'
#' @param input Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param alt_model Provide an arbitrary model here.
#' @param extra_contrasts Provide extra contrasts here.
#' @param model_cond Is condition in the experimental model?
#' @param model_batch Is batch in the experimental model?
#' @param annot_df Include some annotation information in the results?
#' @param force Force deseq to accept data which likely violates its assumptions.
#' @param ... triple dots!  Options are passed to arglist.
#' @return List including the following information:
#'   run = the return from calling DESeq()
#'   denominators = list of denominators in the contrasts
#'   numerators = list of the numerators in the contrasts
#'   conditions = the list of conditions in the experiment
#'   coefficients = list of coefficients making the contrasts
#'   all_tables = list of DE tables
#' @seealso \pkg{DESeq2} \code{\link[DESeq2]{results}} \code{\link[DESeq2]{estimateSizeFactors}}
#'          \code{\link[DESeq2]{estimateDispersions}} \code{\link[DESeq2]{nbinomWaldTest}}
#' @examples
#' \dontrun{
#' pretend = deseq2_pairwise(data, conditions, batches)
#' }
#' @export
deseq2_pairwise <- function(input, conditions=NULL, batches=NULL, model_cond=TRUE,
                            alt_model=NULL, extra_contrasts=NULL, model_intercept=FALSE,
                            model_batch=TRUE, annot_df=NULL, force=FALSE, ...) {
    arglist <- list(...)
    message("Starting DESeq2 pairwise comparisons.")
    input_data <- choose_dataset(input, force=force)
    conditions <- input_data[["conditions"]]
    batches <- input_data[["batches"]]
    data <- input_data[["data"]]

    condition_table <- table(conditions)
    condition_levels <- levels(as.factor(conditions))
    batch_levels <- levels(as.factor(batches))
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    summarized <- NULL
    ## Moving the size-factor estimation into this if(){} block in order to accomodate sva-ish batch estimation in the model
    deseq_sf <- NULL

    model_choice <- choose_model(conditions, batches,
                                 model_batch=model_batch,
                                 model_cond=model_cond,
                                 alt_model=alt_model)
    model_including <- model_choice[["including"]]
    ## choose_model should now take all of the following into account
    ## Therefore the following 8 or so lines should not be needed any longer.
    model_string <- NULL
    ## I keep forgetting that deseq wants the thing you care about last in the list of the model string.
    if (model_including == "batch") {
        model_string <- "~ condition"
    } else if (model_including == "condition+batch") {
        model_string <- "~ batch + condition"
    } else {
        model_string <- "~ condition"
    }

    if (isTRUE(model_batch) & isTRUE(model_cond)) {
        message("DESeq2 step 1/5: Including batch and condition in the deseq model.")
        ## summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ 0 + condition + batch)
        ## conditions and batch in this context is information taken from pData()
        column_data <- Biobase::pData(input[["expressionset"]])
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=column_data,
                                                     ##design=~ batch_levels + condition_levels)
                                                     design=as.formula(model_string))
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=as.formula(model_string))
    } else if (isTRUE(model_batch)) {
        message("DESeq2 step 1/5: Including only batch in the deseq model.")
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=as.formula(model_string))
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=as.formula(model_string))
    } else if (class(model_batch) == 'numeric') {
        message("DESeq2 step 1/5: Including batch estimates from sva/ruv/pca in the deseq model.")
        model_string <- "~ condition"
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=as.formula(model_string))
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=as.formula(model_string))
        dataset$SV1 <- model_batch
        DESeq2::design(dataset) <- as.formula(~ SV1 + condition)
    } else if (class(model_batch) == 'matrix') {
        message("DESeq2 step 1/5: Including a matrix of batch estimates from sva/ruv/pca in the deseq model.")
        model_string <- "~ condition"
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=as.formula(model_string))
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=as.formula(model_string))
        formula_string <- "as.formula(~ "
        for (count in 1:ncol(model_batch)) {
            colname <- paste0("SV", count)
            dataset[[colname]] <- model_batch[, 1]
            formula_string <- paste0(formula_string, " ", colname, " + ")
        }
        formula_string <- paste0(formula_string, "condition)")
        new_formula <- eval(parse(text=formula_string))
        DESeq2::design(dataset) <- new_formula
    } else {
        message("DESeq2 step 1/5: Including only condition in the deseq model.")
        model_string <- "~ condition"
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=as.formula(model_string))
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=as.formula(model_string))
    }
    ## If making a model ~0 + condition -- then must set betaPrior=FALSE
    ## dataset = DESeqDataSet(se=summarized, design=~ 0 + condition)
    message("DESeq2 step 2/5: Estimate size factors.")
    deseq_sf <- DESeq2::estimateSizeFactors(dataset)
    message("DESeq2 step 3/5: Estimate dispersions.")
    deseq_disp <- DESeq2::estimateDispersions(deseq_sf, quiet=TRUE)
    ## deseq_run = nbinomWaldTest(deseq_disp, betaPrior=FALSE)
    message("DESeq2 step 4/5: nbinomWaldTest.")
    ## deseq_run <- DESeq2::DESeq(deseq_disp)
    deseq_run = DESeq2::nbinomWaldTest(deseq_disp, quiet=TRUE)
    ## possible options:  betaPrior=TRUE, betaPriorVar, modelMatrix=NULL
    ## modelMatrixType, maxit=100, useOptim=TRUE useT=FALSE df useQR=TRUE
    ## deseq_run = DESeq2::nbinomLRT(deseq_disp)
    ## Set contrast= for each pairwise comparison here!
    denominators <- list()
    numerators <- list()
    result_list <- list()
    coefficient_list <- list()
    ## The following is an attempted simplification of the contrast formulae
    number_comparisons <- sum(1:(length(condition_levels) - 1))
    inner_count <- 0
    for (c in 1:(length(condition_levels) - 1)) {
        denominator <- condition_levels[c]
        nextc <- c + 1
        for (d in nextc:length(condition_levels)) {
            inner_count <- inner_count + 1
            numerator <- condition_levels[d]
            comparison <- paste0(numerator, "_vs_", denominator)
            message(paste0("DESeq2 step 5/5: ", inner_count, "/",
                           number_comparisons, ": Printing table: ", comparison))
            result <- as.data.frame(DESeq2::results(deseq_run,
                                                    contrast=c("condition", numerator, denominator),
                                                    format="DataFrame"))
            result <- result[order(result$log2FoldChange),]
            colnames(result) <- c("baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
            ## From here on everything is the same.
            result[is.na(result[["P.Value"]]), "P.Value"] = 1 ## Some p-values come out as NA
            result[is.na(result[["adj.P.Val"]]), "adj.P.Val"] = 1 ## Some p-values come out as NA
            result[["baseMean"]] <- signif(x=as.numeric(result[["baseMean"]]), digits=4)
            result[["logFC"]] <- signif(x=as.numeric(result[["logFC"]]), digits=4)
            result[["lfcSE"]] <- signif(x=as.numeric(result[["lfcSE"]]), digits=4)
            result[["stat"]] <- signif(x=as.numeric(result[["stat"]]), digits=4)
            result[["P.Value"]] <- signif(x=as.numeric(result[["P.Value"]]), digits=4)
            result[["adj.P.Val"]] <- signif(x=as.numeric(result[["adj.P.Val"]]), digits=4)

            result$qvalue <- tryCatch(
            {
                ## Nested expressions are way too confusing for me
                ttmp <- as.numeric(result[["P.Value"]])
                ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
                signif(x=ttmp, digits=4)
                ## as.numeric(format(signif(qvalue::qvalue(as.numeric(result$P.Value), robust=TRUE)$qvalues, 4), scientific=TRUE))
            },
            error=function(cond) {
                message(paste0("The qvalue estimation failed for ", comparison, "."))
                return(1)
            },
            ##warning=function(cond) {
            ##    message("There was a warning?")
            ##    message(cond)
            ##    return(1)
            ##},
            finally={
            })
            result_name <- paste0(numerator, "_vs_", denominator)
            denominators[[result_name]] <- denominator
            numerators[[result_name]] <- numerator
            result_list[[result_name]] <- result
            if (!is.null(annot_df)) {
                result <- merge(result, annot_df, by.x="row.names", by.y="row.names")
            }
        } ## End for each d
        ## Fill in the last coefficient (since the for loop above goes from 1 to n-1
        denominator <- names(condition_table[length(conditions)])
        ## denominator_name = paste0("condition", denominator)  ## maybe needed in 6 lines
    }  ## End for each c
    ## Now that we finished the contrasts, fill in the coefficient list with each set of values
    for (c in 1:(length(condition_levels))) {
        coef <- condition_levels[c]
        coef_name <- paste0("condition", coef)
        coefficient_list[[coef]] <- as.data.frame(
            DESeq2::results(deseq_run, contrast=as.numeric(coef_name == DESeq2::resultsNames(deseq_run))))
        message(paste0("Collected coefficients for: ", coef))
        ## coefficient_list[[denominator]] = as.data.frame(results(deseq_run, contrast=as.numeric(denominator_name == resultsNames(deseq_run))))
    }
    ret_list <- list(
        "run" = deseq_run,
        "denominators" = denominators,
        "numerators" = numerators,
        "conditions" = conditions,
        "coefficients" = coefficient_list,
        "all_tables" = result_list
    )
    return(ret_list)
}

## EOF
