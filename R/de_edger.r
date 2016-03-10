## Time-stamp: <Fri Feb  5 12:31:16 2016 Ashton Trey Belew (abelew@gmail.com)>

#'   Plot out 2 coefficients with respect to one another from edger
#'
#' It can be nice to see a plot of two coefficients from a edger comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param output the set of pairwise comparisons provided by edger_pairwise()
#' @param x   the name or number of the first coefficient column to extract, this will be the x-axis of the plot
#' @param y   the name or number of the second coefficient column to extract, this will be the y-axis of the plot
#' @param gvis_filename   A filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline   add a trendline to the gvis plot?
#' @param tooltip_data   a dataframe of gene annotations to be used in the gvis plot
#' @param base_url   for gvis plots
#' @return a ggplot2 plot showing the relationship between the two coefficients
#' @seealso \link{hpgl_linear_scatter} \link{edger_pairwise}
#' @examples
#' \dontrun{
#'  pretty = coefficient_scatter(limma_data, x="wt", y="mut")
#' }
#' @export
edger_coefficient_scatter <- function(output, x=1, y=2,
                                      gvis_filename=NULL,
                                      gvis_trendline=TRUE, tooltip_data=NULL,
                                      base_url=NULL) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    message("This can do comparisons among the following columns in the edger result:")
    thenames <- names(output$contrasts$identities)
    cat(thenames)
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
    ## It looks like the lrt data structure is redundant, so I will test that by looking at the apparent
    ## coefficients from lrt[[1]] and then repeating with lrt[[2]]
    coefficient_df <- output$lrt[[1]]$coefficients
    coefficient_df <- coefficient_df[, c(xname, yname)]
    if (max(coefficient_df) < 0) {
        coefficient_df <- coefficient_df * -1.0
    }

    plot <- hpgl_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename,
                                gvis_trendline=gvis_trendline, first=xname, second=yname,
                                tooltip_data=tooltip_data, base_url=base_url)
    maxvalue <- as.numeric(max(coefficient_df) + 1)
    print(maxvalue)
    plot$scatter <- plot$scatter +
        ggplot2::scale_x_continuous(limits=c(0, maxvalue)) +
        ggplot2::scale_y_continuous(limits=c(0, maxvalue))
    plot$df <- coefficient_df
    return(plot)
}

#' Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using EdgeR.
#'
#' @param input  a dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions   a factor of conditions in the experiment
#' @param batches   a factor of batches in the experiment
#' @param model_cond   Include condition in the experimental model?  This is pretty much always true.
#' @param model_batch   Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept   Use cell means or intercept? (I default to the former,
#'   but they work out the same)
#' @param alt_model   An alternate experimental model to use
#' @param extra_contrasts   some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df   Add some annotation information to the data tables?
#' @param ... The elipsis parameter is fed to write_edger() at the end.
#' @return A list including the following information:
#'   contrasts = The string representation of the contrasts performed.
#'   lrt = A list of the results from calling glmLRT(), one for each contrast.
#'   contrast_list = The list of each call to makeContrasts()
#'     I do this to avoid running into the limit on # of contrasts addressable by topTags()
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
edger_pairwise <- function(input, conditions=NULL, batches=NULL, model_cond=TRUE,
                          model_batch=FALSE, model_intercept=FALSE, alt_model=NULL,
                          extra_contrasts=NULL, annot_df=NULL, ...) {
    message("Starting edgeR pairwise comparisons.")
    input_class <- class(input)[1]
    if (input_class == 'expt') {
        conditions <- input$conditions
        batches <- input$batches
        data <- as.data.frame(Biobase::exprs(input$expressionset))
        ## As I understand it, edgeR fits a binomial distribution
        ## and expects data as floating point counts,
        ## not a log2 transformation.
        if (!is.null(input$transform)) {
            if (input$transform == "log2") {
                ##data = (2^data) - 1
                data <- input$normalized$normalized_counts$count_table
            } ## End checking for log2 normalized data
        } ## End checking for transformed data
    } else { ## End checking if this is an expt
        data <- as.data.frame(input)
    }
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    ## It would be much smarter to generate the models in the following if() {} blocks
    ## But I have it in my head to eventually compare results using different models.
    cond_model <- stats::model.matrix(~ 0 + conditions)
    batch_model <- try(stats::model.matrix(~ 0 + batches), silent=TRUE)
    condbatch_model <- try(stats::model.matrix(~ 0 + conditions + batches), silent=TRUE)
    cond_int_model <- try(stats::model.matrix(~ conditions), silent=TRUE)
    batch_int_model <- try(stats::model.matrix(~ batches), silent=TRUE)
    condbatch_int_model <- try(stats::model.matrix(~ conditions + batches), silent=TRUE)
    fun_model <- NULL
    fun_int_model <- NULL
    if (isTRUE(model_cond) & isTRUE(model_batch)) {
        fun_model <- condbatch_model
        fun_int_model <- condbatch_int_model
    } else if (isTRUE(model_cond)) {
        fun_model <- cond_model
        fun_int_model <- cond_int_model
    } else if (isTRUE(model_batch)) {
        fun_model <- batch_model
        fun_int_model <- batch_int_model
    } else {
        ## Default to the conditional model
        fun_model <- cond_model
        fun_int_model <- cond_int_model
    }
    if (isTRUE(model_intercept)) {
        fun_model <- fun_int_model
    }
    if (!is.null(alt_model)) {
        fun_model <- alt_model
    }
    tmpnames <- colnames(fun_model)
    tmpnames <- gsub("data[[:punct:]]", "", tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    colnames(fun_model) <- tmpnames

    ##tmpnames = colnames(condbatch_model)
    ##tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    ##tmpnames = gsub("conditions", "", tmpnames)
    ##colnames(cond_model) = tmpnames
    raw <- edgeR::DGEList(counts=data, group=conditions)
    message("EdgeR step 1/9: normalizing data.")
    norm <- edgeR::calcNormFactors(raw)
    message("EdgeR step 2/9: Estimating the common dispersion.")
    disp_norm <- edgeR::estimateCommonDisp(norm)
    message("EdgeR step 3/9: Estimating dispersion across genes.")
    tagdisp_norm <- edgeR::estimateTagwiseDisp(disp_norm)
    message("EdgeR step 4/9: Estimating GLM Common dispersion.")
    glm_norm <- edgeR::estimateGLMCommonDisp(tagdisp_norm, fun_model)
    message("EdgeR step 5/9: Estimating GLM Trended dispersion.")
    glm_trended <- edgeR::estimateGLMTrendedDisp(glm_norm, fun_model)
    message("EdgeR step 6/9: Estimating GLM Tagged dispersion.")
    glm_tagged <- edgeR::estimateGLMTagwiseDisp(glm_trended, fun_model)
    message("EdgeR step 7/9: Running glmFit.")
    cond_fit <- edgeR::glmFit(glm_tagged, design=fun_model)
    message("EdgeR step 8/9: Making pairwise contrasts.")
    apc <- make_pairwise_contrasts(fun_model, conditions, do_identities=FALSE)
    ## This is pretty weird because glmLRT only seems to take up to 7 contrasts at a time...
    contrast_list <- list()
    result_list <- list()
    lrt_list <- list()
    sc <- vector("list", length(apc$names))
    end <- length(apc$names)
    for (con in 1:length(apc$names)) {
        name <- apc$names[[con]]
        message(paste0("EdgeR step 9/9: ", con, "/", end, ": Printing table: ", name, ".")) ## correct
        sc[[name]] <- gsub(pattern=",", replacement="", apc$all_pairwise[[con]])
        tt <- parse(text=sc[[name]])
        ctr_string <- paste0("tt = limma::makeContrasts(", tt, ", levels=fun_model)")
        eval(parse(text=ctr_string))
        contrast_list[[name]] <- tt
        lrt_list[[name]] <- edgeR::glmLRT(cond_fit, contrast=contrast_list[[name]])
        res <- edgeR::topTags(lrt_list[[name]], n=nrow(data), sort.by="logFC")
        res <- as.data.frame(res)
        res$logFC <- signif(x=as.numeric(res$logFC), digits=4)
        res$logCPM <- signif(x=as.numeric(res$logCPM), digits=4)
        res$LR <- signif(x=as.numeric(res$LR), digits=4)
        res$PValue <- format(x=as.numeric(res$PValue), digits=4, scientific=TRUE)
        res$FDR <- format(x=as.numeric(res$FDR), digits=4, scientific=TRUE)
        res$qvalue <- tryCatch(
        {
            ##as.numeric(format(signif(
            ##    suppressWarnings(qvalue::qvalue(
            ##        as.numeric(res$PValue), robust=TRUE))$qvalues, 4),
            ##scientific=TRUE))
            ## ok I admit it, I am not smart enough for nested expressions
            ttmp <- as.numeric(res$PValue)
            ttmp <- qvalue::qvalue(ttmp)$qvalues
            format(x=ttmp, digits=4, scientific=TRUE)
        },
        error=function(cond) {
            message(paste0("The qvalue estimation failed for ", name, "."))
            return(1)
        },
        ##warning=function(cond) {
        ##    message("There was a warning?")
        ##    message(cond)
        ##    return(1)
        ##},
        finally={
        })
        result_list[[name]] <- res
    } ## End for loop
    final <- list(
        contrasts=apc,
        lrt=lrt_list,
        contrast_list=contrast_list,
        all_tables=result_list)
    return(final)
}

## EOF
