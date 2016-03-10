## Time-stamp: <Tue Mar  8 13:19:01 2016 Ashton Trey Belew (abelew@gmail.com)>

#'   Plot out 2 coefficients with respect to one another from deseq2
#'
#' It can be nice to see a plot of two coefficients from a deseq2 comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param output the set of pairwise comparisons provided by deseq_pairwise()
#' @param x   the name or number of the first coefficient column to extract, this will be the x-axis of the plot
#' @param y   the name or number of the second coefficient column to extract, this will be the y-axis of the plot
#' @param gvis_filename   A filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline   add a trendline to the gvis plot?
#' @param tooltip_data   a dataframe of gene annotations to be used in the gvis plot
#' @param base_url   for gvis plots
#' @return a ggplot2 plot showing the relationship between the two coefficients
#' @seealso \link{hpgl_linear_scatter} \link{deseq2_pairwise}
#' @examples
#' \dontrun{
#'  pretty = coefficient_scatter(deseq_data, x="wt", y="mut")
#' }
#' @export
deseq_coefficient_scatter <- function(output, x=1, y=2, ## gvis_filename="limma_scatter.html",
                                      gvis_filename=NULL,
                                      gvis_trendline=TRUE, tooltip_data=NULL,
                                      base_url=NULL) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    message("This can do comparisons among the following columns in the deseq2 result:")
    thenames <- names(output$coefficients)
    message(thenames)
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
    first_df <- output$coefficients[[xname]]
    first_df$delta <- log2(first_df$baseMean) + first_df$log2FoldChange
    second_df <- output$coefficients[[yname]]
    second_df$delta <- log2(second_df$baseMean) + second_df$log2FoldChange
    first_col <- first_df[,c("baseMean","log2FoldChange","delta")]
    colnames(first_col) <- c("mean.1", "fc.1", xname)
    second_col <- second_df[,c("baseMean","log2FoldChange","delta")]
    colnames(second_col) <- c("mean.2", "fc.2", yname)
    coefficient_df <- merge(first_col, second_col, by="row.names")
    rownames(coefficient_df) <- coefficient_df$Row.names
    coefficient_df <- coefficient_df[-1]
    coefficient_df <- coefficient_df[,c(xname, yname, "mean.1", "mean.2")]
    coefficient_df[is.na(coefficient_df)] <- 0
    maxvalue <- max(coefficient_df) + 1.0
    plot <- hpgl_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename,
                                gvis_trendline=gvis_trendline, first=xname, second=yname,
                                tooltip_data=tooltip_data, base_url=base_url)
    plot$scatter <- plot$scatter +
        ggplot2::scale_x_continuous(limits=c(0, maxvalue)) +
        ggplot2::scale_y_continuous(limits=c(0, maxvalue))
    plot$df <- coefficient_df
    return(plot)
}

#' deseq_pairwise()  Because I can't be trusted to remember '2'
#'
#' This calls deseq2_pairwise(...) because I am determined to forget typing deseq2
#'
#' @param ... I like cats
#' @return stuff from deseq2_pairwise
#' @seealso \link{deseq2_pairwise}
#' @export
deseq_pairwise <- function(...) {
    message("Hey you, use deseq2 pairwise.")
    deseq2_pairwise(...)
}

#' Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using DESeq2.
#'
#' @param input  A dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions   A factor of conditions in the experiment
#' @param batches  A factor of batches in the experiment
#' @param model_cond   Have condition in the experimental model?
#' @param model_batch   Have batch in the experimental model?
#' @param annot_df   Include some annotation information in the results?
#' @param ... triple dots!
#' @return A list including the following information:
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
                            model_batch=FALSE, annot_df=NULL, ...) {
    arglist <- list(...)
    message("Starting DESeq2 pairwise comparisons.")
    input_class <- class(input)[1]
    if (input_class == 'expt') {
        conditions <- input$conditions
        batches <- input$batches
        data <- as.data.frame(Biobase::exprs(input$expressionset))
        if (!is.null(input$normalized)) {
            ## As I understand it, DESeq2 (and edgeR) fits a binomial distribution
            ## and expects data as floating point counts,
            ## not a log2 transformation.
            if (input$normalized[[1]] != "raw") {
                message("DESeq2 demands raw data as input, reverting to the original expressionset.")
                data <- Biobase::exprs(input$original_expressionset)
            } else if (!is.null(input$transform)) {
                if (input$transform == "log2") {
                    ##data = (2^data) - 1
                    data <- input$normalized$normalized_counts$count_table
                }
            }
        } ## End testing if normalization has been performed
    } else {
        data <- as.data.frame(input)
    }
    condition_table <- table(conditions)
    conditions <- levels(as.factor(conditions))
    batches <- levels(as.factor(batches))
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    summarized <- NULL
    if (isTRUE(model_batch) & isTRUE(model_cond)) {
        message("DESeq2 step 1/5: Including batch and condition in the deseq model.")
        ## summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ 0 + condition + batch)
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=~ batch + condition)
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=~ batch + condition)
    } else if (isTRUE(model_batch)) {
        message("DESeq2 step 1/5: Including only batch in the deseq model.")
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=~ batch)
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=~ batch)
    } else {
        message("DESeq2 step 1/5: Including only condition in the deseq model.")
        summarized <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                     colData=Biobase::pData(input$expressionset),
                                                     design=~ condition)
        dataset <- DESeq2::DESeqDataSet(se=summarized, design=~ condition)
    }
    ## If making a model ~0 + condition -- then must set betaPrior=FALSE
    ## dataset = DESeqDataSet(se=summarized, design=~ 0 + condition)
    message("DESeq2 step 2/5: Estimate size factors.")
    deseq_sf <- DESeq2::estimateSizeFactors(dataset)
    message("DESeq2 step 3/5: estimate Dispersions.")
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
    number_comparisons <- sum(1:(length(conditions) - 1))
    inner_count <- 0
    for (c in 1:(length(conditions) - 1)) {
        denominator <- conditions[c]
        nextc <- c + 1
        for (d in nextc:length(conditions)) {
            inner_count <- inner_count + 1
            numerator <- conditions[d]
            comparison <- paste0(numerator, "_vs_", denominator)
            message(paste0("DESeq2 step 5/5: ", inner_count, "/",
                           number_comparisons, ": Printing table: ", comparison))
            result <- as.data.frame(DESeq2::results(deseq_run,
                                                    contrast=c("condition", numerator, denominator),
                                                    format="DataFrame"))
            result <- result[order(result$log2FoldChange),]
            colnames(result) <- c("baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
            ## From here on everything is the same.
            result[is.na(result$P.Value), "P.Value"] = 1 ## Some p-values come out as NA
            result[is.na(result$adj.P.Val), "adj.P.Val"] = 1 ## Some p-values come out as NA
            result$baseMean <- signif(x=as.numeric(result$baseMean), digits=4)
            result$logFC <- signif(x=as.numeric(result$logFC), digits=4)
            result$lfcSE <- signif(x=as.numeric(result$lfcSE), digits=4)
            result$stat <- signif(x=as.numeric(result$stat), digits=4)
            result$P.Value <- format(x=as.numeric(result$P.Value), digits=4, scientific=TRUE)
            result$adj.P.Val <- format(x=as.numeric(result$adj.P.Val), digits=4, scientific=TRUE)

            result$qvalue <- tryCatch(
            {
                ## Nested expressions are way too confusing for me
                ttmp <- as.numeric(result$P.Value)
                ttmp <- qvalue::qvalue(ttmp)$qvalues
                format(x=ttmp, digits=4, scientific=TRUE)
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
    for (c in 1:(length(conditions))) {
        coef <- conditions[c]
        coef_name <- paste0("condition", coef)
        coefficient_list[[coef]] <- as.data.frame(DESeq2::results(deseq_run, contrast=as.numeric(coef_name == DESeq2::resultsNames(deseq_run))))
        message(paste0("Collected coefficients for: ", coef))
        ## coefficient_list[[denominator]] = as.data.frame(results(deseq_run, contrast=as.numeric(denominator_name == resultsNames(deseq_run))))
    }
    ret_list <- list(
        run=deseq_run,
        denominators=denominators,
        numerators=numerators,
        conditions=conditions,
        coefficients=coefficient_list,
        all_tables=result_list
    )
    return(ret_list)
}

## EOF
