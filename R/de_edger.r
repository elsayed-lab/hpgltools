#' Make a MA plot of some limma output with pretty colors and shapes
#'
#' Yay pretty colors and shapes!
#'
#' @param output  The result from all_pairwise(), which should be changed to handle other invocations too.
#' @param table  Result from edger to use, left alone it chooses the first.
#' @param p_col  Column to use for p-value data.
#' @return a plot!
#' @seealso \link{plot_ma_de}
#' @examples
#'  \dontrun{
#'   prettyplot <- edger_ma(all_aprwise) ## [sic, I'm witty! and can speel]
#' }
#' @export
edger_ma <- function(output, table=NULL) {
    counts <- NULL
    de_genes <- NULL
    pval <- NULL
    output <- output[["edger"]]  ## Currently this only will work with the output from all_pairwise()
    possible_tables <- output[["contrasts"]][["names"]]
    if (is.null(table)) {
        table <- possible_tables[1]
    } else if (is.numeric(table)) {
        table <- possible_tables[table]
    }

    de_genes <- output[["all_tables"]][[table]]
    pval <- "qvalue"
    plot <- plot_ma_de(table=de_genes, expr_col="logCPM", fc_col="logFC", p_col="qvalue")
    return(plot)
}

#' Plot two coefficients with respect to one another from edgeR.
#'
#' It can be nice to see a plot of two coefficients from a edger comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param output Set of pairwise comparisons provided by edger_pairwise().
#' @param toptable  The table to use for extracting the logfc values.
#' @param x Name or number of the x-axis coefficient column to extract.
#' @param y Name or number of the y-axis coefficient column to extract.
#' @param gvis_filename Filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline Add a trendline to the gvis plot?
#' @param z  Make pretty colors for genes this number of z-scores from the median.
#' @param tooltip_data Dataframe of gene annotations to be used in the gvis plot.
#' @param base_url Add a linkout to gvis plots to this base url.
#' @param color_low  Color to use for low-logfc values.
#' @param color_high  Color to use for high-logfc values.
#' @param ... A few options may be added outside this scope and are left in the arglist, notably
#'     qlimit, fc_column, p_column.  I need to make a consistent decision about how to handle these
#'     not-always needed parameters, either always define them in the function body, or always put
#'     them in arglist(...), doing a little of both is stupid.
#' @return Ggplot2 plot showing the relationship between the two coefficients.
#' @seealso \link{plot_linear_scatter} \link{edger_pairwise}
#' @examples
#' \dontrun{
#'  pretty = coefficient_scatter(limma_data, x="wt", y="mut")
#' }
#' @export
edger_coefficient_scatter <- function(output, toptable=NULL, x=1, y=2,
                                      gvis_filename=NULL, gvis_trendline=TRUE, z=1.5,
                                      tooltip_data=NULL, base_url=NULL,
                                      color_low="#DD0000", color_high="#7B9F35", ...) {
    arglist <- list(...)
    qlimit <- 0.1
    if (!is.null(arglist[["qlimit"]])) {
        qlimit <- arglist[["qlimit"]]
    }
    fc_column <- "edger_logfc"
    if (!is.null(arglist[["fc_column"]])) {
        fc_column <- arglist[["fc_column"]]
    }
    p_column <- "edger_adjp"
    if (!is.null(arglist[["p_column"]])) {
        p_column <- arglist[["p_column"]]
    }
    thenames <- names(output[["contrasts"]][["identities"]])
    message("This can do comparisons among the following columns in the edger result:")
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
    ## It looks like the lrt data structure is redundant, so I will test that by looking at the apparent
    ## coefficients from lrt[[1]] and then repeating with lrt[[2]]
    coefficient_df <- output[["lrt"]][[1]][["coefficients"]]
    coefficient_df <- coefficient_df[, c(xname, yname)]
    if (max(coefficient_df) < 0) {
        coefficient_df <- coefficient_df * -1.0
    }

    plot <- sm(plot_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename,
                                   gvis_trendline=gvis_trendline, first=xname, second=yname,
                                   tooltip_data=tooltip_data, base_url=base_url,
                                   pretty_colors=FALSE, color_low=color_low, color_high=color_high))
    maxvalue <- as.numeric(max(coefficient_df) + 1)
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
edger_pairwise <- function(input, conditions=NULL, batches=NULL, model_cond=TRUE,
                          model_batch=TRUE, model_intercept=TRUE, alt_model=NULL,
                          extra_contrasts=NULL, annot_df=NULL, force=FALSE, edger_method="default", ...) {
    message("Starting edgeR pairwise comparisons.")
    input_data <- choose_dataset(input, force=force)
    design <- Biobase::pData(input[["expressionset"]])
    conditions <- design[["condition"]]
    batches <- design[["batch"]]
    data <- input_data[["data"]]

    fun_model <- choose_model(conditions, batches,
                              model_batch=model_batch,
                              model_cond=model_cond,
                              model_intercept=model_intercept,
                              alt_model=alt_model)
    fun_model <- fun_model[["chosen_model"]]

    ## I have a strong sense that the most recent version of edgeR changed its dispersion estimate code
    ## Here is a note from the user's guide, which may have been there previously and I merely did not notice:
    ## To estimate common dispersion, trended dispersions and tagwise dispersions in one run
    ## y <- estimateDisp(y, design)
    raw <- edgeR::DGEList(counts=data, group=conditions)
    message("EdgeR step 1/9: normalizing data.")
    norm <- edgeR::calcNormFactors(raw)
    ##message("EdgeR step 2/9: Estimating the common dispersion.")
    ##disp_norm <- edgeR::estimateCommonDisp(norm)
    ##message("EdgeR step 3/9: Estimating dispersion across genes.")
    ##tagdisp_norm <- edgeR::estimateTagwiseDisp(disp_norm)
    ##message("EdgeR step 4/9: Estimating GLM Common dispersion.")
    ##glm_norm <- edgeR::estimateGLMCommonDisp(tagdisp_norm, fun_model)
    ##message("EdgeR step 5/9: Estimating GLM Trended dispersion.")
    ##glm_trended <- edgeR::estimateGLMTrendedDisp(glm_norm, fun_model)
    ##message("EdgeR step 6/9: Estimating GLM Tagged dispersion.")
    ##glm_tagged <- edgeR::estimateGLMTagwiseDisp(glm_trended, fun_model)
    ##message("EdgeR step 7/9: Running glmFit.")
    ##cond_fit <- edgeR::glmFit(glm_tagged, design=fun_model)
    ##message("EdgeR step 8/9: Making pairwise contrasts.")
    ##apc <- make_pairwise_contrasts(fun_model, conditions, do_identities=FALSE)

    ## Try this instead:
    norm <- edgeR::estimateDisp(norm, design=fun_model, robust=TRUE)
    ##cond_fit <- edgeR::glmFit(norm, design=fun_model)
    cond_fit <- edgeR::glmQLFit(norm, design=fun_model, robust=TRUE)
    apc <- make_pairwise_contrasts(fun_model, conditions,
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
        message(paste0("EdgeR step 9/9: ", con, "/", end, ": Printing table: ", name, ".")) ## correct
        sc[[name]] <- gsub(pattern=",", replacement="", apc[["all_pairwise"]][[con]])
        tt <- parse(text=sc[[name]])
        ctr_string <- paste0("tt = limma::makeContrasts(", tt, ", levels=fun_model)")
        eval(parse(text=ctr_string))
        contrast_list[[name]] <- tt
        ##lrt_list[[name]] <- edgeR::glmLRT(cond_fit, contrast=contrast_list[[name]])
        lrt_list[[name]] <- edgeR::glmQLFTest(cond_fit, contrast=contrast_list[[name]])
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
        res[["qvalue"]] <- tryCatch(
        {
            ##as.numeric(format(signif(
            ##    suppressWarnings(qvalue::qvalue(
            ##        as.numeric(res$PValue), robust=TRUE))$qvalues, 4),
            ##scientific=TRUE))
            ## ok I admit it, I am not smart enough for nested expressions
            ttmp <- as.numeric(res[["PValue"]])
            ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
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
        "model" = fun_model,
        "contrasts" = apc,
        "lrt" = lrt_list,
        "contrast_list" = contrast_list,
        "all_tables" = result_list)
    return(final)
}

## EOF
