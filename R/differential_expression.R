## Time-stamp: <Fri Nov 20 14:47:16 2015 Ashton Trey Belew (abelew@gmail.com)>

## Test for infected/control/beads -- a placebo effect?
## The goal is therefore to find responses different than beads
## The null hypothesis is (H0): (infected == uninfected) || (infected == beads)
## The alt hypothesis is (HA): (infected != uninfected) && (infected != beads)
disjunct_tab = function(contrast_fit, coef1, coef2, ...) {
    stat = pmin(abs(contrast_fit[,coef1]), abs(contrast_fit[,coef2]))
    pval = pmax(contrast_fit$p.val[,coef1], contrast_fit$p.val[,coef2])
}
## An F-test only does inf==uninf && inf==bead
## So the solution is to separately perform the two subtests and subset for the set of genes for which both are true.
## However, if you do that, the f-statistics are a little screwey, but there are a few ways to handle it:
## Perform the two separate tests and perform the following combination of the stat and p-value:
##    stat = min(|inf-uninf|, |inf-bead|)  (logFC)
##    ^^pval^^ = max(pval(inf-uninf), pval(inf-beads))
##    adj.pval = p.adjust(^^pval^^, method='BH')
## ReportingTools hwriter

#' all_pairwise()  Wrap up limma/DESeq2/EdgeR pairwise analyses in one call.
#'
#' @param input  a dataframe/vector or expt class containing count tables, normalization state, etc.
#' @param conditions default=NULL  a factor of conditions in the experiment
#' @param batches default=NULL  a factor of batches in the experiment
#' @param model_cond default=TRUE  include condition in the model?  This is likely always true.
#' @param model_batch default=FALSE  include batch in the model?
#' @param model_intercept default=FALSE  use an intercept model instead of cell means?
#' @param extra_contrasts default=NULL some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param alt_model default=NULL an optional alternate model to use rather than just condition/batch
#' @param libsize default=NULL the library size of the original data to help voom()
#' @param ... The elipsis parameter is fed to write_limma() at the end.
#'
#' @return A list of limma, deseq, edger results.
#' @export
#' @examples
#' ## finished_comparison = eBayes(limma_output)
#' ## data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
all_pairwise = function(input, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=TRUE, model_intercept=FALSE, extra_contrasts=NULL, alt_model=NULL, libsize=NULL) {
    conditions = get0('conditions')
    batches = get0('batches')
    model_cond = get0('model_cond')
    model_batch = get0('model_batch')
    model_intercept = get0('model_intercept')
    extra_contrasts = get0('model_contrasts')
    alt_model = get0('alt_model')
    libsize = get0('libsize')
    if (is.null(model_cond)) {
        model_cond = TRUE
    }
    if (is.null(model_batch)) {
        model_batch = FALSE
    }
    if (is.null(model_intercept)) {
        model_intercept = FALSE
    }

    limma_result = limma_pairwise(input, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize)
    deseq_result = deseq2_pairwise(input, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize)
    edger_result = edger_pairwise(input, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize)
    basic_result = basic_pairwise(input, conditions)

    result_comparison = compare_tables(limma=limma_result, deseq=deseq_result, edger=edger_result, basic=basic_result)
    ret = list(limma=limma_result, deseq=deseq_result, edger=edger_result, basic=basic_result, comparison=result_comparison)
    return(ret)
}

#' combine_de_tables()  Combine portions of deseq/limma/edger table output
#'
#' This hopefully makes it easy to compare the outputs from limma/DESeq2/EdgeR on a table-by-table basis.
#'
#' @param all_pairwise_result  the output from all_pairwise()
#' @param table default='wt_vs_mut'  the name of a table comparison performed by deseq/limma/edger.
#'
#' @return a table combinine limma/edger/deseq outputs.
#' @seealso \code{\link{all_pairwise}}
#' @export
#' @examples
#' ## pretty = combine_de_tables(big_result, table='t12_vs_t0')
combine_de_tables = function(all_pairwise_result, table='wt_vs_mut', annot_df=NULL) {
    limma = all_pairwise_result$limma
    deseq = all_pairwise_result$deseq
    edger = all_pairwise_result$edger
    limma = limma$all_tables[[table]]
    colnames(limma) = c("limma_logfc","limma_ave","limma_t","limma_p","limma_adjp","limma_b","limma_q")
    limma = limma[,c("limma_logfc","limma_ave","limma_t","limma_b","limma_p","limma_adjp","limma_q")]
    deseq = deseq$all_tables[[table]]
    colnames(deseq) = c("deseq_basemean","deseq_logfc","deseq_lfcse","deseq_stat","deseq_p","deseq_adjp","deseq_q")
    deseq = deseq[,c("deseq_logfc","deseq_basemean","deseq_lfcse","deseq_stat","deseq_p","deseq_adjp","deseq_q")]
    edger = edger$all_tables[[table]]
    colnames(edger) = c("edger_logfc","edger_logcpm","edger_lr","edger_p","edger_adjp","edger_q")
    combined = merge(limma, deseq, by="row.names")
    combined = merge(combined, edger, by.x="Row.names", by.y="row.names")
    rownames(combined) = combined$Row.names
    combined = combined[-1]
    combined[is.na(combined)] = 0
    temp_fc = cbind(as.numeric(combined$limma_logfc), as.numeric(combined$edger_logfc), as.numeric(combined$deseq_logfc))
    temp_fc = preprocessCore::normalize.quantiles(as.matrix(temp_fc))
    combined$fc_meta = rowMeans(temp_fc, na.rm=TRUE)
    combined$fc_var = rowVars(temp_fc, na.rm=TRUE)
    combined$fc_varbymed = combined$fc_var / combined$fc_meta
    temp_p = cbind(as.numeric(combined$limma_p), as.numeric(combined$edger_p), as.numeric(combined$deseq_p))
    combined$p_meta = rowMeans(temp_p, na.rm=TRUE)
    combined$p_var = rowVars(temp_p, na.rm=TRUE)

    combined$q_meta = tryCatch(
    {
        format(signif(qvalue(combined$p_meta, robust=TRUE)$qvalues, 4), scientific=TRUE)
    },
    error=function(cond) {
        message(paste0("The meta qvalue estimation failed."))
        return(1)
    },
    warning=function(cond) {
        message("There was a warning!")
        message(cond)
        return(1)
    },
    finally={
    })
    if (!is.null(annot_df)) {
        combined = merge(annot_df, combined, by="row.names")
        rownames(combined) = combined$Row.names
        combined = combined[-1]
    }
    return(combined)
}

#' limma_coefficient_scatter()  Plot out 2 coefficients with respect to one another from limma
#'
#' It can be nice to see a plot of two coefficients from a limma comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param limma_output the set of pairwise comparisons provided by limma_pairwise()
#' @param x default=1  the name or number of the first coefficient column to extract, this will be the x-axis of the plot
#' @param y default=2  the name or number of the second coefficient column to extract, this will be the y-axis of the plot
#' @param gvis_filename default='limma_scatter.html'  A filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline default=TRUE  add a trendline to the gvis plot?
#' @param tooltip_data default=NULL  a dataframe of gene annotations to be used in the gvis plot
#'
#' @return a ggplot2 plot showing the relationship between the two coefficients
#' @seealso \code{\link{hpgl_linear_scatter}} \code{\link{limma_pairwise}}
#' @export
#' @examples
#' ## pretty = coefficient_scatter(limma_data, x="wt", y="mut")
limma_coefficient_scatter = function(output, toptable=NULL, x=1, y=2, gvis_filename="limma_scatter.html", gvis_trendline=TRUE, tooltip_data=NULL, flip=FALSE, base_url=NULL) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    print("This can do comparisons among the following columns in the limma result:")
    thenames = colnames(output$pairwise_comparisons$coefficients)
    print(thenames)
    xname=""
    yname=""
    if (is.numeric(x)) {
        xname = thenames[[x]]
    } else {
        xname = x
    }
    if (is.numeric(y)) {
        yname = thenames[[y]]
    } else {
        yname = y
    }
    ## This is just a shortcut in case I want to flip axes without thinking.
    if (isTRUE(flip)) {
        tmp = x
        tmpname = xname
        x = y
        xname = yname
        y = tmp
        yname = tmpname
        rm(tmp)
        rm(tmpname)
    }
    print(paste0("Actually comparing ", xname, " and ", yname, "."))
    coefficients = output$pairwise_comparisons$coefficients
    coefficients = coefficients[,c(x,y)]
    plot = hpgl_linear_scatter(df=coefficients, loess=TRUE, gvis_filename=gvis_filename, gvis_trendline=gvis_trendline, first=xname, second=yname, tooltip_data=tooltip_data, base_url=base_url, pretty_colors=FALSE)

    if (!is.null(toptable)) {
        theplot = plot$scatter + theme_bw()
        sig = limma_subset(toptable, z=1.5)
        sigup = sig$up
        sigdown = sig$down
        sigup = subset(sigup, qvalue < 0.1)
        sigdown = subset(sigdown, qvalue < 0.1)
        if (isTRUE(flip)) {
            tmp = sigup
            sigup = sigdown
            sigdown = tmp
            rm(tmp)
        }
        up_index = rownames(coefficients) %in% rownames(sigup)
        down_index = rownames(coefficients) %in% rownames(sigdown)
        up_df = as.data.frame(coefficients[up_index, ])
        down_df = as.data.frame(coefficients[down_index, ])
        colnames(up_df) = c("first","second")
        colnames(down_df) = c("first","second")
        theplot = theplot + geom_point(data=up_df, colour="#7B9F35") + geom_point(data=down_df, colour="#DD0000")
        plot$scatter = theplot
    }
    plot$df = coefficients
    return(plot)
}

#' coefficient_scatter()  Plot out 2 coefficients with respect to one another from limma
#'
#' It can be nice to see a plot of two coefficients from a limma comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param limma_output the set of pairwise comparisons provided by limma_pairwise()
#' @param x default=1  the name or number of the first coefficient column to extract, this will be the x-axis of the plot
#' @param y default=2  the name or number of the second coefficient column to extract, this will be the y-axis of the plot
#' @param gvis_filename default='limma_scatter.html'  A filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline default=TRUE  add a trendline to the gvis plot?
#' @param tooltip_data default=NULL  a dataframe of gene annotations to be used in the gvis plot
#'
#' @return a ggplot2 plot showing the relationship between the two coefficients
#' @seealso \code{\link{hpgl_linear_scatter}} \code{\link{limma_pairwise}}
#' @export
#' @examples
#' ## pretty = coefficient_scatter(limma_data, x="wt", y="mut")
deseq_coefficient_scatter = function(output, x=1, y=2, gvis_filename="limma_scatter.html", gvis_trendline=TRUE, tooltip_data=NULL, flip=FALSE, base_url=NULL) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    print("This can do comparisons among the following columns in the deseq2 result:")
    thenames = names(output$coefficients)
    print(thenames)
    xname=""
    yname=""
    if (is.numeric(x)) {
        xname = thenames[[x]]
    } else {
        xname = x
    }
    if (is.numeric(y)) {
        yname = thenames[[y]]
    } else {
        yname = y
    }
    ## This is just a shortcut in case I want to flip axes without thinking.
    if (isTRUE(flip)) {
        tmp = x
        tmpname = xname
        x = y
        xname = yname
        y = tmp
        yname = tmpname
        rm(tmp)
        rm(tmpname)
    }
    print(paste0("Actually comparing ", xname, " and ", yname, "."))
    first_df = output$coefficients[[xname]]
    first_df$delta = log2(first_df$baseMean) + first_df$log2FoldChange
    second_df = output$coefficients[[yname]]
    second_df$delta = log2(second_df$baseMean) + second_df$log2FoldChange
    first_col = first_df[,c("baseMean","log2FoldChange","delta")]
    colnames(first_col) = c("mean.1", "fc.1", xname)
    second_col = second_df[,c("baseMean","log2FoldChange","delta")]
    colnames(second_col) = c("mean.2", "fc.2", yname)
    coefficient_df = merge(first_col, second_col, by="row.names")
    rownames(coefficient_df) = coefficient_df$Row.names
    coefficient_df = coefficient_df[-1]
    coefficient_df = coefficient_df[,c(xname, yname, "mean.1", "mean.2")]
    coefficient_df[is.na(coefficient_df)] = 0

    plot = hpgl_linear_scatter(df=coefficient_df, loess=TRUE, gvis_filename=gvis_filename, gvis_trendline=gvis_trendline, first=xname, second=yname, tooltip_data=tooltip_data, base_url=base_url)
    plot$df = coefficient_df
    return(plot)
}

#' compare_tables()  See how similar are results from limma/deseq/edger.
#'
#' limma, DEseq2, and EdgeR all make somewhat different assumptions
#' and choices about what makes a meaningful set of differentially
#' expressed genes.  This seeks to provide a quick and dirty metric
#' describing the degree to which they (dis)agree.
#'
#' @param limma default=NULL  limma data from limma_pairwise()
#' @param deseq default=NULL  deseq data from deseq2_pairwise()
#' @param edger default=NULL  edger data from edger_pairwise()
#'
#' @return a heatmap showing how similar they are along with some
#' correlations betwee the three players.
#' @seealso \code{\link{limma_pairwise}} \code{\link{edger_pairwise}} \code{\link{deseq2_pairwise}}
#' @export
#' @examples
#' ## l = limma_pairwise(expt)
#' ## d = deseq_pairwise(expt)
#' ## e = edger_pairwise(expt)
#' fun = compare_tables(limma=l, deseq=d, edger=e)
compare_tables = function(limma=NULL, deseq=NULL, edger=NULL, basic=NULL, include_basic=TRUE) {
    ## Fill each column/row of these with the correlation between tools for one contrast performed
    if (class(limma) == "list") { ## Then this was fed the raw output from limma_pairwise, lets assume the same is true for deseq/edger too and pull out the result tables.
        limma = limma$all_tables
        deseq = deseq$all_tables
        edger = edger$all_tables
        basic = basic$all_tables
    }
    len = length(names(deseq))
    limma_vs_edger = list()
    limma_vs_deseq = list()
    limma_vs_basic = list()
    edger_vs_deseq = list()
    edger_vs_basic = list()
    deseq_vs_basic = list()
    cc = 0
    last = length(names(deseq))
    for (comp in names(deseq)) {  ## assume all three have the same names() -- note that limma has more than the other two though
        cc = cc + 1
        message(paste0(cc, "/", len, ": Comparing analyses for: ", comp))
        l = data.frame(limma[[comp]])
        e = data.frame(edger[[comp]])
        d = data.frame(deseq[[comp]])
        b = data.frame(basic[[comp]])
        le = merge(l, e, by.x="row.names", by.y="row.names")
        le = le[,c("logFC.x","logFC.y")]
        lec = cor.test(le[,1], le[,2])$estimate
        ld = merge(l, d, by.x="row.names", by.y="row.names")
        ld = ld[,c("logFC.x","logFC.y")]
        ldc = cor.test(ld[,1], ld[,2])$estimate
        lb = merge(l, b, by.x="row.names", by.y="row.names")
        lb = lb[,c("logFC.x","logFC.y")]
        lbc = cor.test(lb[,1], lb[,2])$estimate
        ed = merge(e, d, by.x="row.names", by.y="row.names")
        ed = ed[,c("logFC.x","logFC.y")]
        edc = cor.test(ed[,1], ed[,2])$estimate
        eb = merge(e, b, by.x="row.names", by.y="row.names")
        eb = eb[,c("logFC.x","logFC.y")]
        ebc = cor.test(eb[,1], eb[,2])$estimate
        db = merge(d, b, by.x="row.names", by.y="row.names")
        db = db[,c("logFC.x","logFC.y")]
        dbc = cor.test(db[,1], db[,2])$estimate
        limma_vs_edger[[comp]] = lec
        limma_vs_deseq[[comp]] = ldc
        edger_vs_deseq[[comp]] = edc
        limma_vs_basic[[comp]] = lbc
        edger_vs_basic[[comp]] = ebc
        deseq_vs_basic[[comp]] = dbc
    } ## End loop
    names(limma_vs_edger) = names(deseq)
    names(limma_vs_deseq) = names(deseq)
    names(edger_vs_deseq) = names(deseq)
    names(limma_vs_basic) = names(deseq)
    names(edger_vs_basic) = names(deseq)
    names(deseq_vs_basic) = names(deseq)

    comparison_df = rbind(as.numeric(limma_vs_edger), as.numeric(limma_vs_deseq))
    comparison_df = rbind(comparison_df, as.numeric(edger_vs_deseq))
    if (isTRUE(include_basic)) {
        comparison_df = rbind(comparison_df, as.numeric(limma_vs_basic))
        comparison_df = rbind(comparison_df, as.numeric(edger_vs_basic))
        comparison_df = rbind(comparison_df, as.numeric(deseq_vs_basic))
        rownames(comparison_df) = c("le", "ld", "ed", "lb", "eb", "db")
    } else {
        rownames(comparison_df) = c("le", "ld", "ed")
    }
    comparison_df = as.matrix(comparison_df)
    colnames(comparison_df) = names(deseq)
    heat_colors = colorRampPalette(c("white","black"))
    comparison_heatmap = try(heatmap.3(comparison_df, scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=heat_colors, dendrogram="none", Rowv=FALSE, Colv=FALSE, main="Compare DE tools"), silent=TRUE)
    heat=NULL
    if (class(comparison_heatmap) != 'try-error') {
        heat = recordPlot()
    }
    ret = list(limma_vs_edger=limma_vs_edger, limma_vs_deseq=limma_vs_deseq, limma_vs_basic=limma_vs_basic,
               edger_vs_deseq=edger_vs_deseq, edger_vs_basic=edger_vs_basic,
               deseq_vs_basic=deseq_vs_basic,
               comp=comparison_df, heat=heat)
    return(ret)
}

#' deseq_pairwise()  Because I can't be trusted to remember '2'
#'
#' This calls deseq2_pairwise(...) because I am determined to forget typing deseq2
#' @param look at deseq2_pairwise
#'
#' @return stuff from deseq2_pairwise
#'
#' @export
deseq_pairwise = function(...) {
    print("Hey you, use deseq2 pairwise.")
    deseq2_pairwise(...)
}

#' deseq2_pairwise()  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using DESeq2.
#'
#' @param input  a dataframe/vector or  expt class containing data, normalization state, etc.
#' @param conditions default=NULL  a factor of conditions in the experiment
#' @param batches default=NULL a factor of batches in the experiment
#'
#' @return A list including the following information:
#'   run = the return from calling DESeq()
#'   result = the return from calling results()
#'   result_mle = the return from calling results() with the MLE parameter.
#'   results = the return when a deseq result class is cast as a dataframe.
#'
#' @seealso \code{\link{topTags}} \code{\link{glmLRT}} \code{\link{makeContrasts}}
#' @export
#' @examples
#' ## pretend = deseq2_pairwise(data, conditions, batches)
deseq2_pairwise = function(input, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, excel=FALSE, csv=TRUE, annot_df=NULL, workbook="excel/deseq.xls", ...) {
    print("Starting DESeq2 pairwise comparisons.")
    input_class = class(input)[1]
    if (input_class == 'expt') {
        conditions = input$conditions
        batches = input$batches
        data = as.data.frame(exprs(input$expressionset))
        if (!is.null(input$norm)) {
            ## As I understand it, DESeq2 (and edgeR) fits a binomial distribution
            ## and expects data as floating point counts,
            ## not a log2 transformation.
            if (input$norm != "raw") {
                print("DESeq2 demands raw data as input, reverting to the original expressionset.")
                data = exprs(input$original_expressionset)
            } else if (!is.null(input$transform)) {
                if (input$transform == "log2") {
                    ##data = (2^data) - 1
                    data = input$normalized$normalized_counts$count_table
                }
            }
        } ## End testing if normalization has been performed
    } else {
        data = as.data.frame(input)
    }
    condition_table = table(conditions)
    batch_table = table(batches)
    conditions = levels(as.factor(conditions))
    batches = levels(as.factor(batches))
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches

    ## An interesting note about the use of formulae in DESeq:
    ## "you should put the variable of interest at the end of the formula and make sure the control level is the first level."
    ## Thus, all these formulae should have condition(s) at the end.
    ##cond_model = model.matrix(~ 0 + conditions)
    ##batch_model = try(model.matrix(~ 0 + batches), silent=TRUE)
    ##condbatch_model = try(model.matrix(~ 0 + batches + conditions), silent=TRUE)
    ##cond_int_model = try(model.matrix(~ conditions), silent=TRUE)
    ##condbatch_int_model = try(model.matrix(~ batches + conditions), silent=TRUE)
    ##tmpnames = colnames(cond_model)
    ##tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    ##tmpnames = gsub("conditions", "", tmpnames)
    ##colnames(cond_model) = tmpnames
    ##fun_model = NULL
    ##fun_int_model = NULL
    summarized = NULL
    if (isTRUE(model_batch) & isTRUE(model_cond)) {
        message("Attempting to include batch and condition in the model for DESeq.")
        ##        summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ 0 + condition + batch)
        summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ batch + condition)
        dataset = DESeqDataSet(se=summarized, design=~ batch + condition)
    } else if (isTRUE(model_batch)) {
        message("Attempting to include only batch in the deseq model, this will likely fail.")
        summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ batch)
        dataset = DESeqDataSet(se=summarized, design=~ batch)
    } else {
        message("Including only condition in the deseq model.")
        summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(input$expressionset), design=~ condition)
        dataset = DESeqDataSet(se=summarized, design=~ condition)
    }
    ## If making a model ~0 + condition -- then must set betaPrior=FALSE
    ##dataset = DESeqDataSet(se=summarized, design=~ 0 + condition)
    deseq_sf = estimateSizeFactors(dataset)
    deseq_disp = estimateDispersions(deseq_sf)
    ##deseq_run = nbinomWaldTest(deseq_disp, betaPrior=FALSE)
    ##deseq_run = nbinomWaldTest(deseq_disp)
    deseq_run = DESeq(deseq_disp)
    ## Set contrast= for each pairwise comparison here!
    result_names = resultsNames(deseq_run)
    denominators = list()
    numerators = list()
    result_list = list()
    result_mle_list = list()
    binom_list = list()
    coefficient_list = list()
##    if (isTRUE(model_cond)) {
##        condition_list = grep("^condition", resultsNames(deseq_run), value=TRUE)
##    } else if (isTRUE(model_batch)) {
##        condition_list = grep("^batch", resultsNames(deseq_run), value=TRUE)
##    }
##    for (c in 1:(length(condition_list) - 1)) {
##        denominator = names(condition_table[c])
##        ## This is where it will fall down if you only want to look at batch.
##        denominator_name = paste0("condition", denominator)
##        my_contrast = as.numeric(denominator_name == resultsNames(deseq_run))
##        if (sum(my_contrast) == 0) {
##            next
##        }
##        coefficient_list[[denominator]] = as.data.frame(results(deseq_run, contrast=my_contrast))
##        nextc = c + 1
##        for (d in nextc:length(condition_list)) {
##            numerator = names(condition_table[d])
##            message(paste0("DESeq2:", c, "/", d, ": Printing table: ", numerator, "_vs_", denominator))
##            result = as.data.frame(results(deseq_run, contrast=c("condition", numerator, denominator), format="DataFrame"))
##            result = result[order(result$log2FoldChange),]
##            colnames(result) = c("baseMean","logFC", "lfcSE","stat","P.Value","adj.P.Val")

    ## The following is an attempted simplification of the contrast formulae
    for (c in 1:(length(conditions) - 1)) {
        denominator = conditions[c]
        nextc = c + 1
        for (d in nextc:length(conditions)) {
            numerator = conditions[d]
            message(paste0("DESeq2:", c, "/", d, ": Printing table: ", numerator, "_vs_", denominator))
            result = as.data.frame(results(deseq_run, contrast=c("condition", numerator, denominator), format="DataFrame"))
            result = result[order(result$log2FoldChange),]
            colnames(result) = c("baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
            ## From here on everything is the same.
            result[is.na(result$P.Value), "P.Value"] = 1 ## Some p-values come out as NA
            result[is.na(result$adj.P.Val), "adj.P.Val"] = 1 ## Some p-values come out as NA
            result$qvalue = tryCatch(
                {
                    format(signif(qvalue(result$P.Value, robust=TRUE)$qvalues, 4), scientific=TRUE)
                },
                error=function(cond) {
                    message(paste0("The qvalue estimation failed for ", comparison, "."))
                    return(1)
                },
                warning=function(cond) {
                    message("There was a warning?")
                    message(cond)
                    return(1)
                },
                finally={
                }
            )
            result$P.Value = format(signif(result$P.Value, 4), scientific=TRUE)
            result$adj.P.Val = format(signif(result$adj.P.Val, 4), scientific=TRUE)
            result_name = paste0(numerator, "_vs_", denominator)
            denominators[[result_name]] = denominator
            numerators[[result_name]] = numerator
            result_list[[result_name]] = result

            if (!is.null(annot_df)) {
                result = merge(result, annot_df, by.x="row.names", by.y="row.names")
            }
            testdir = dirname(workbook)
            if (isTRUE(excel) | isTRUE(csv)) {
                if (!file.exists(testdir)) {
                    dir.create(testdir)
                    message(paste0("Creating directory: ", testdir, "for writing excel/csv data."))
                }
            }
            if (isTRUE(excel)) {
                try(write_xls(data=result, sheet=result_name, file=workbook, overwrite=TRUE))
            }
            if (isTRUE(csv)) {
                csv_filename = gsub(".xls$", "", workbook)
                csv_filename = paste0(csv_filename, "_", result_name, ".csv")
                write.csv(result, file=csv_filename)
            }
        } ## End for each d
        ## Fill in the last coefficient (since the for loop above goes from 1 to n-1
        denominator = names(condition_table[length(conditions)])
        denominator_name = paste0("condition", denominator)
    }  ## End for each c
    ## Now that we finished the contrasts, fill in the coefficient list with each set of values
    for (c in 1:(length(conditions))) {
        coef = conditions[c]
        coef_name = paste0("condition", coef)
        coefficient_list[[coef]] = as.data.frame(results(deseq_run, contrast=as.numeric(coef_name == resultsNames(deseq_run))))
##        coefficient_list[[denominator]] = as.data.frame(results(deseq_run, contrast=as.numeric(denominator_name == resultsNames(deseq_run))))
    }

    ret_list = list(
        run=deseq_run,
        denominators=denominators,
        numerators=numerators,
        conditions=conditions,
        coefficients=coefficient_list,
        all_tables=result_list
    )
    return(ret_list)
    ##    deseq_result = results(deseq_run)
    ##    deseq_mle_result = results(deseq_run, addMLE=TRUE)
    ##    deseq_df = data.frame(deseq_result[order(deseq_result$log2FoldChange),])
    ##    plotMA(deseq_df)
    ## identify(deseq_result$baseMean, deseq_result$log2FoldChange)
    ##    ma = recordPlot()
    ##d = plotCounts(dataset, gene=which.min(deseq_result$padj), intgroup="condition", returnData=TRUE)
    ##ggplot(d, aes(x=condition, y=count)) +
    ##    geom_point(position=position_jitter(w=0.1,h=0)) +
    ##    scale_y_log10(breaks=c(25,100,400))

    ##  mcols(deseq_result)$description
    ## DESeq can do multi-factors, but the important one goes last
    ## design(summarized) = formula(~ batch + condition)
    ## deseq_run = DESeq(summarized)

    ## Different contrasts may be applied with
    ## contrast_res = results(deseq_run, contrast=c("condition", "numerator", "denominator"))

    ## data transformations:
    ## rld = rlog(deseq_run)
    ## vsd = try(varianceStabilizingTransformation(deseq_run), silent=TRUE)
    ## rlogMat = assay(rld)
    ## vstMat = try(assay(vsd), silent=TRUE)
    ## par(mfrow=c(1,3))
    ## notAllZero <- (rowSums(counts(deseq_run))>0)
    ##    meanSdPlot(log2(counts(deseq_run, normalized=TRUE)[notAllZero,] + 1))
    ##    meanSdPlot(assay(rld[notAllZero,]))
    ##    meanSdPlot(assay(vsd[notAllZero,]))
}

#' edger_pairwise()  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using EdgeR.
#'
#' @param input  a dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions default=NULL  a factor of conditions in the experiment
#' @param batches default=NULL  a factor of batches in the experiment
#' @param model_cond default=TRUE  Include condition in the experimental model?  This is pretty much always true.
#' @param model_batch default=FALSE  Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept default=FALSE Use cell means or intercept? (I default to the former, but they work out the same)
#' @param extra_contrasts default=NULL  some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param ... The elipsis parameter is fed to write_edger() at the end.
#'
#' @return A list including the following information:
#'   results = A list of tables returned by 'topTags', one for each contrast.
#'   contrasts = The string representation of the contrasts performed.
#'   lrt = A list of the results from calling glmLRT(), one for each contrast.
#'   contrast_list = The list of each call to makeContrasts()
#'     I do this to avoid running into the limit on # of contrasts addressable by topTags()
#'
#' @seealso \code{\link{topTags}} \code{\link{glmLRT}} \code{\link{makeContrasts}}
#' @export
#' @examples
#' ## pretend = edger_pairwise(data, conditions, batches)
edger_pairwise = function(input, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, model_intercept=FALSE, alt_model=NULL, extra_contrasts=NULL, excel=FALSE, csv=TRUE, annotation=NULL, workbook="excel/edger.xls", ...) {
    print("Starting EdgeR pairwise comparisons.")
    input_class = class(input)[1]
    if (input_class == 'expt') {
        conditions = input$conditions
        batches = input$batches
        data = as.data.frame(exprs(input$expressionset))
        ## As I understand it, edgeR fits a binomial distribution
        ## and expects data as floating point counts,
        ## not a log2 transformation.
        if (!is.null(input$transform)) {
            if (input$transform == "log2") {
                ##data = (2^data) - 1
                data = input$normalized$normalized_counts$count_table
            } ## End checking for log2 normalized data
        } ## End checking for transformed data
    } else { ## End checking if this is an expt
        data = as.data.frame(input)
    }
    message("At this time, this only does conditional models.")
    condition_table = table(conditions)
    batch_table = table(batches)
    conditions = as.factor(conditions)
    batches = as.factor(batches)
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    ## It would be much smarter to generate the models in the following if() {} blocks
    ## But I have it in my head to eventually compare results using different models.
    cond_model = model.matrix(~ 0 + conditions)
    batch_model = try(model.matrix(~ 0 + batches), silent=TRUE)
    condbatch_model = try(model.matrix(~ 0 + conditions + batches), silent=TRUE)
    cond_int_model = try(model.matrix(~ conditions), silent=TRUE)
    condbatch_int_model = try(model.matrix(~ conditions + batches), silent=TRUE)
    fun_model = NULL
    fun_int_model = NULL
    if (isTRUE(model_cond) & isTRUE(model_batch)) {
        fun_model = condbatch_model
        fun_int_model = condbatch_int_model
    } else if (isTRUE(model_cond)) {
        fun_model = cond_model
        fun_int_model = cond_int_model
    } else if (isTRUE(model_batch)) {
        fun_model = batch_model
        fun_int_model = batch_int_model
    } else {
        ## Default to the conditional model
        fun_model = cond_model
        fun_int_model = cond_int_model
    }
    if (isTRUE(model_intercept)) {
        fun_model = fun_int_model
    }
    if (!is.null(alt_model)) {
        fun_model = alt_model
    }
    tmpnames = colnames(fun_model)
    tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    tmpnames = gsub("conditions", "", tmpnames)
    colnames(fun_model) = tmpnames

    ##tmpnames = colnames(condbatch_model)
    ##tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    ##tmpnames = gsub("conditions", "", tmpnames)
    ##colnames(cond_model) = tmpnames

    raw = DGEList(counts=data, group=conditions)
    message("Using EdgeR to normalize the data.")
    norm = calcNormFactors(raw)
    message("EdgeR: Estimating the common dispersion.")
    disp_norm = estimateCommonDisp(norm)
    message("EdgeR: Estimating dispersion across genes.")
    tagdisp_norm = estimateTagwiseDisp(disp_norm)
    message("EdgeR: Estimating GLM Common dispersion.")
    glm_norm = estimateGLMCommonDisp(tagdisp_norm, fun_model)
    message("EdgeR: Estimating GLM Trended dispersion.")
    glm_trended = estimateGLMTrendedDisp(glm_norm, fun_model)
    message("EdgeR: Estimating GLM Tagged dispersion.")
    glm_tagged = estimateGLMTagwiseDisp(glm_trended, fun_model)
    cond_fit = edgeR::glmFit(glm_tagged, design=fun_model)

    apc = make_pairwise_contrasts(fun_model, conditions, do_identities=FALSE)
    ## This is pretty weird because glmLRT only seems to take up to 7 contrasts at a time...
    contrast_list = list()
    result_list = list()
    lrt_list = list()
    single_contrasts = list()
    sc = vector("list", length(apc$names))
    end = length(apc$names)
    for (con in 1:length(apc$names)) {
        name = apc$names[[con]]
        message(paste0("EdgeR:", con, "/", end, ": Printing table: ", name, ".")) ## correct
        sc[[name]] = gsub(pattern=",", replacement="", apc$all_pairwise[[con]])
        tt = parse(text=sc[[name]])
        ctr_string = paste0("tt = makeContrasts(", tt, ", levels=fun_model)")
        eval(parse(text=ctr_string))
        contrast_list[[name]] = tt
        lrt_list[[name]] = edgeR::glmLRT(cond_fit, contrast=contrast_list[[name]])
        res = topTags(lrt_list[[name]], n=nrow(data), sort.by="logFC")
        res = as.data.frame(res)
        res$qvalue = tryCatch(
            {
                as.numeric(format(signif(qvalue(res$PValue, robust=TRUE)$qvalues, 4), scientific=TRUE))
            },
            error=function(cond) {
                message(paste("The qvalue estimation failed for ", comparison, ".", sep=""))
                return(1)
            },
            warning=function(cond) {
                message("There was a warning?")
                message(cond)
                return(1)
            },
            finally={
            }
        )
        res$PValue = format(signif(res$PValue, 4), scientific=TRUE)
        res$FDR = format(signif(res$FDR, 4), scientific=TRUE)
        result_list[[name]] = res
    }
    testdir = dirname(workbook)
    if (isTRUE(excel) | isTRUE(csv)) {
        if (!file.exists(testdir)) {
            dir.create(testdir)
            message(paste0("Creating directory: ", testdir, "for writing excel/csv data."))
        }
    }
    if (isTRUE(excel)) {
        try(write_xls(data=res, sheet=name, file=workbook, overwrite=TRUE))
    }
    if (isTRUE(csv)) {
        csv_filename = gsub(".xls$", "", workbook)
        csv_filename = paste0(csv_filename, "_", name, ".csv")
        write.csv(res, file=csv_filename)
    }

    final = list(
        contrasts=apc,
        lrt=lrt_list,
        contrast_list=contrast_list,
        all_tables=result_list)

    return(final)
}

#' hpgl_voom()  A slight modification of limma's voom() function.
#' Estimate mean-variance relationship between samples and generate
#' 'observational-level weights' in preparation for linear modelling
#' RNAseq data.  This particular implementation was primarily scabbed
#' from cbcbSEQ, but changes the mean-variance plot slightly and
#' attempts to handle corner cases where the sample design is
#' confounded by setting the coefficient to 1 for those samples rather
#' than throwing an unhelpful error.  Also, the Elist output gets a
#' 'plot' slot which contains the plot rather than just printing it.
#'
#' @param dataframe a dataframe of sample counts which have been
#' normalized and log transformed
#' @param model default=NULL  an experimental model defining batches/conditions/etc
#' @param libsize default=NULL  the size of the libraries (usually provided by
#' edgeR).
#' @param stupid default=FALSE  whether or not to cheat when the resulting matrix is not solvable.
#' @param logged default=FALSE  whether the input data is known to be logged.
#' @param converted default=FALSE  whether the input data is known to be cpm converted.
#'
#' @return an EList containing the following information:
#'   E = The normalized data
#'   weights = The weights of said data
#'   design = The resulting design
#'   lib.size = The size in pseudocounts of the library
#'   plot = A ggplot of the mean/variance trend with a blue loess fit and red trend fit
#'
#' @seealso \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{lmFit}}
#'
#' @export
#' @examples
#' ## funkytown = hpgl_voom(samples, model)
hpgl_voom = function(dataframe, model=NULL, libsize=NULL, stupid=FALSE, logged=FALSE, converted=FALSE) {
    out = list()
    if (is.null(libsize)) {
        libsize = colSums(dataframe, na.rm=TRUE)
    }
    if (converted == 'cpm') {
        converted = TRUE
    }
    if (!isTRUE(converted)) {
        message("The voom input was not cpm, converting now.")
        posed = t(dataframe + 0.5)
        dataframe = t(posed/(libsize + 1) * 1e+06)
        ##y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1000000)) ## from voom()
    }
    if (logged == 'log2') {
        logged = TRUE
    }
    if (isTRUE(logged)) {
        if (max(dataframe) > 1000) {
            warning("This data appears to not be logged, the lmfit will do weird things.")
        }
    } else {
        if (max(dataframe) < 400) {
            warning("This data says it was not logged, but the maximum counts seem small.")
            warning("If it really was log2 transformed, then we are about to double-log it and that would be very bad.")
        }
        message("The voom input was not log2, transforming now.")
        dataframe = log2(dataframe)
    }
    dataframe = as.matrix(dataframe)

    if (is.null(design)) {
        design = matrix(1, ncol(dataframe), 1)
        rownames(design) = colnames(dataframe)
        colnames(design) = "GrandMean"
    }
    linear_fit = limma::lmFit(dataframe, model, method="ls")
    if (is.null(linear_fit$Amean)) {
        linear_fit$Amean = rowMeans(dataframe, na.rm=TRUE)
    }
    sx = linear_fit$Amean + mean(log2(libsize + 1)) - log2(1e+06)
    sy = sqrt(linear_fit$sigma)
    if (is.na(sum(sy))) { ## 1 replicate
        return(NULL)
    }
    allzero = rowSums(dataframe) == 0
    stupid_NAs = is.na(sx)
    sx = sx[!stupid_NAs]
    stupid_NAs = is.na(sy)
    sy = sy[!stupid_NAs]
    if (any(allzero == TRUE, na.rm=TRUE)) {
        sx = sx[!allzero]
        sy = sy[!allzero]
    }
    fitted = gplots::lowess(sx, sy, f=0.5)
    f = stats::approxfun(fitted, rule=2)
    mean_var_df = data.frame(mean=sx, var=sy)
    mean_var_plot = ggplot2::ggplot(mean_var_df, aes(x=mean, y=var)) +
        geom_point() +
        xlab("Log2(count size + 0.5)") +
        ylab("Square root of the standard deviation.") +
        stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE, show_guide=FALSE) +
        scale_fill_gradientn(colours=colorRampPalette(c("white","black"))(256)) +
        geom_smooth(method="loess") +
        stat_function(fun=f, colour="red") +
        theme(legend.position="none")
    if (is.null(linear_fit$rank)) {
        message("Some samples cannot be balanced across the experimental design.")
        if (stupid) {
            ## I think this is telling me I have confounded data, and so
            ## for those replicates I will have no usable coefficients, so
            ## I say set them to 1 and leave them alone.
            linear_fit$coefficients[is.na(linear_fit$coef)] = 1
            fitted.values = linear_fit$coef %*% t(linear_fit$design)
        }
    } else if (linear_fit$rank < ncol(linear_fit$design)) {
        j = linear_fit$pivot[1:linear_fit$rank]
        fitted.values = linear_fit$coef[, j, drop=FALSE] %*% t(linear_fit$design[, j, drop=FALSE])
    } else {
        fitted.values = linear_fit$coef %*% t(linear_fit$design)
    }
    fitted.cpm = 2^fitted.values
    fitted.count = 1e-06 * t(t(fitted.cpm) * (libsize + 1))
    fitted.logcount = log2(fitted.count)
    w = 1/f(fitted.logcount)^4
    dim(w) = dim(fitted.logcount)
    rownames(w) = rownames(dataframe)
    colnames(w) = colnames(dataframe)
    out$E = dataframe
    out$weights = w
    out$design = model
    out$lib.size = libsize
    out$plot = mean_var_plot
    new("EList", out)
}

#' limma_pairwise()  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using voom/limma.
#'
#' @param input  a dataframe/vector or expt class containing count tables, normalization state, etc.
#' @param conditions default=NULL  a factor of conditions in the experiment
#' @param batches default=NULL  a factor of batches in the experiment
#' @param extra_contrasts default=NULL  some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param model_cond default=TRUE  include condition in the model?
#' @param model_batch default=FALSE  include batch in the model? This is hopefully TRUE.
#' @param model_intercept default=FALSE  perform a cell-means or intercept model?  A little more difficult for me to understand.  I have tested and get the same answer either way.
#' @param alt_model default=NULL  a separate model matrix instead of the normal condition/batch.
#' @param libsize default=NULL  I've recently figured out that libsize is far more important than I previously realized.  Play with it here.
#' @param ... The elipsis parameter is fed to write_limma() at the end.
#'
#' @return A list including the following information:
#'   macb = the mashing together of condition/batch so you can look at it
#'   macb_model = The result of calling model.matrix(~0 + macb)
#'   macb_fit =  The result of calling lmFit(data, macb_model)
#'   voom_result = The result from voom()
#'   voom_design = The design from voom (redundant from voom_result, but convenient)
#'   macb_table = A table of the number of times each condition/batch pairing happens
#'   cond_table = A table of the number of times each condition appears (the denominator for the identities)
#'   batch_table = How many times each batch appears
#'   identities = The list of strings defining each condition by itself
#'   all_pairwise = The list of strings defining all the pairwise contrasts
#'   contrast_string = The string making up the makeContrasts() call
#'   pairwise_fits = The result from calling contrasts.fit()
#'   pairwise_comparisons = The result from eBayes()
#'   limma_result = The result from calling write_limma()
#'
#' @seealso \code{\link{write_limma}}
#' @export
#' @examples
#' ## pretend = balanced_pairwise(data, conditions, batches)
limma_pairwise = function(input, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, model_intercept=FALSE, extra_contrasts=NULL, alt_model=NULL, libsize=NULL) {
    print("Starting limma pairwise comparison.")
    input_class = class(input)[1]
    if (input_class == 'expt') {
        conditions = input$conditions
        batches = input$batches
        data = exprs(input$expressionset)
        if (is.null(libsize)) {
            message("libsize was not specified, this parameter has profound effects on limma's result.")
            if (!is.null(input$best_libsize)) {
                message("Using the libsize from expt$best_libsize.")
                ## libsize = expt$norm_libsize
                libsize = input$best_libsize
            } else {
                message("Using the libsize from expt$normalized$normalized_counts.")
                libsize = input$normalized$normalized_counts$libsize
            }
        } else {
            message("libsize was specified.  This parameter has profound effects on limma's result.")
        }
    } else {  ## Not an expt class, data frame or matrix
        data = as.data.frame(input)
    }
    if (is.null(libsize)) {
        libsize = colSums(data)
    }
    condition_table = table(conditions)
    batch_table = table(batches)
    conditions = as.factor(conditions)
    batches = as.factor(batches)
    ## Make a model matrix which will have one entry for each of these condition/batches
    cond_model = model.matrix(~ 0 + conditions)  ## I am not putting a try() on this, because if it fails, then we are effed.
    batch_model = try(model.matrix(~ 0 + batches), silent=TRUE)
    condbatch_model = try(model.matrix(~ 0 + conditions + batches), silent=TRUE)
    batch_int_model = try(model.matrix(~ batches), silent=TRUE)
    cond_int_model = try(model.matrix(~ conditions), silent=TRUE)
    condbatch_int_model = try(model.matrix(~ conditions + batches), silent=TRUE)
    fun_model = NULL
    fun_int_model = NULL
    if (isTRUE(model_cond) & isTRUE(model_batch)) {
        fun_model = condbatch_model
        fun_int_model = condbatch_int_model
    } else if (isTRUE(model_cond)) {
        fun_model = cond_model
        fun_int_model = cond_int_model
    } else if (isTRUE(model_batch)) {
        fun_model = batch_model
        fun_int_model = batch_int_model
    } else {
        ## Default to the conditional model
        fun_model = cond_model
        fun_int_model = cond_int_model
    }
    if (isTRUE(model_intercept)) {
        fun_model = fun_int_model
    }
    if (!is.null(alt_model)) {
        fun_model = alt_model
    }
    tmpnames = colnames(fun_model)
    tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    tmpnames = gsub("-", "", tmpnames)
    tmpnames = gsub("+", "", tmpnames)
    tmpnames = gsub("conditions", "", tmpnames)
    colnames(fun_model) = tmpnames
    fun_voom = NULL

    ## voom() it, taking into account whether the data has been log2 transformed.
    logged = input$transform
    if (is.null(logged)) {
        print("I don't know if this data is logged, testing if it is integer.")
        if (is.integer(all_data)) {
            logged = FALSE
        } else {
            logged = TRUE
        }
    } else {
        if (logged == "raw") {
            logged = FALSE
        } else {
            logged = TRUE
        }
    }
    converted = input$convert
    if (is.null(converted)) {
        print("I cannot determine if this data has been converted, assuming no.")
        converted = FALSE
    } else {
        if (converted == "raw") {
            converted = FALSE
        } else {
            converted = TRUE
        }
    }
    ##fun_voom = voom(data, fun_model)
    ##fun_voom = hpgl_voom(data, fun_model, libsize=libsize)
    ##fun_voom = voomMod(data, fun_model, lib.size=libsize)
    fun_voom = hpgl_voom(data, fun_model, libsize=libsize, logged=logged, converted=converted)
    one_replicate = FALSE
    if (is.null(fun_voom)) {
        one_replicate = TRUE
        fun_voom = data
        fun_design = design
    } else {
        fun_design = fun_voom$design
    }

    ## Extract the design created by voom()
    ## This is interesting because each column of the design will have a prefix string 'macb' before the
    ## condition/batch string, so for the case of clbr_tryp_batch_C it will look like: macbclbr_tryp_batch_C
    ## This will be important in 17 lines from now.
    ## Do the lmFit() using this model
    fun_fit = lmFit(fun_voom, fun_model)
    ##fun_fit = lmFit(fun_voom)
    ## The following three tables are used to quantify the relative contribution of each batch to the sample condition.
    if (isTRUE(model_intercept)) {
        contrasts = "intercept"
        identities = NULL
        contrast_string = NULL
        all_pairwise = NULL
        all_pairwise_fits = fun_fit
    } else {
        contrasts = make_pairwise_contrasts(fun_model, conditions, extra_contrasts=extra_contrasts)
        all_pairwise_contrasts = contrasts$all_pairwise_contrasts
        identities = contrasts$identities
        contrast_string = contrasts$contrast_string
        all_pairwise = contrasts$all_pairwise
        ## Once all that is done, perform the fit
        ## This will first provide the relative abundances of each condition
        ## followed by the set of all pairwise comparisons.
        all_pairwise_fits = contrasts.fit(fun_fit, all_pairwise_contrasts)
    }
    all_tables = NULL
    if (isTRUE(one_replicate)) {
        all_pairwise_comparisons = all_pairwise_fits$coefficients
    } else {
        all_pairwise_comparisons = eBayes(all_pairwise_fits)
        all_tables = try(topTable(all_pairwise_comparisons, number=nrow(all_pairwise_comparisons)))
    }
    if (isTRUE(model_intercept)) {
        limma_result = all_tables
    } else {
        limma_result = try(write_limma(all_pairwise_comparisons, excel=FALSE))
    }
    result = list(
        input_data=data,
        conditions_table=condition_table,
        batches_table=batch_table,
        conditions=conditions,
        batches=batches,
        model=fun_model,
        fit=fun_fit,
        voom_result=fun_voom,
        voom_design=fun_design,
        identities=identities,
        all_pairwise=all_pairwise,
        contrast_string=contrast_string,
        pairwise_fits=all_pairwise_fits,
        pairwise_comparisons=all_pairwise_comparisons,
        single_table=all_tables,
        all_tables=limma_result)
    return(result)
}

#' limma_scatter()  Plot arbitrary data from limma
#'
#' @param all_pairwise_result  the result from calling balanced_pairwise()
#' @param first_table default=1  the first table from all_pairwise_result$limma_result to look at (may be a name or number)
#' @param first_column default='logFC'  the name of the column to plot from the first table
#' @param second_table default=2  the second table inside all_pairwise_result$limma_result (name or number)
#' @param second_column  a column to compare against
#' @param type A type of scatter plot (linear model, distance, vanilla)
#' @param ... so that you may feed it the gvis/tooltip information to make clicky graphs if so desired.
#'
#' @return a hpgl_linear_scatter() set of plots comparing the chosen columns
#' If you forget to specify tables to compare, it will try the first vs the second.
#' @seealso \code{\link{hpgl_linear_scatter}}, \code{\link{topTable}},
#'
#' @export
#' @examples
#' ## compare_logFC = limma_scatter(all_pairwise, first_table="wild_type", second_column="mutant", first_table="AveExpr", second_column="AveExpr")
#' ## compare_B = limma_scatter(all_pairwise, first_column="B", second_column="B")
limma_scatter = function(all_pairwise_result, first_table=1, first_column="logFC", second_table=2, second_column="logFC", type="linear_scatter", ...) {
    tables = all_pairwise_result$all_tables
    if (is.numeric(first_table)) {
        x_name = paste(names(tables)[first_table], first_column, sep=":")
    }
    if (is.numeric(second_table)) {
        y_name = paste(names(tables)[second_table], second_column, sep=":")
    }

    ## This section is a little bit paranoid
    ## I want to make absolutely certain that I am adding only the
    ## two columns I care about and that nothing gets reordered
    ## As a result I am explicitly pulling a single column, setting
    ## the names, then pulling the second column, then cbind()ing them.
    x_name = paste(first_table, first_column, sep=":")
    y_name = paste(second_table, second_column, sep=":")
    df = data.frame(x=tables[[first_table]][[first_column]])
    rownames(df) = rownames(tables[[first_table]])
    second_column_list = tables[[second_table]][[second_column]]
    names(second_column_list) = rownames(tables[[second_table]])
    df = cbind(df, second_column_list)
    colnames(df) = c(x_name, y_name)

    plots = NULL
    if (type == "linear_scatter") {
        plots = hpgl_linear_scatter(df, loess=TRUE, ...)
    } else if (type == "dist_scatter") {
        plots = hpgl_dist_scatter(df, ...)
    } else {
        plots = hpgl_scatter(df, ...)
    }
    plots[['dataframe']] = df
    return(plots)
}

#' limma_subset()  A quick and dirty way to pull the top/bottom genes from toptable()
#'
#' @param table  the original data from limma
#' @param n default=NULL  a number of genes to keep
#' @param z default=NULL  a number of z-scores from the mean
#'
#' If neither n nor z is provided, it assumes you want 1.5 z-scores from the median.
#'
#' @return a dataframe subset from toptable
#'
#' @seealso \code{\link{limma}}
#'
#' @export
#' @examples
#' ## subset = limma_subset(df, n=400)
#' ## subset = limma_subset(df, z=1.5)
limma_subset = function(table, n=NULL, z=NULL) {
    if (is.null(n) & is.null(z)) {
        z = 1.5
    }
    if (is.null(n)) {
        out_summary = summary(table$logFC)
        out_mad = mad(table$logFC, na.rm=TRUE)
        up_median_dist = out_summary["Median"] + (out_mad * z)
        down_median_dist = out_summary["Median"] - (out_mad * z)
        up_genes = subset(table, logFC >= up_median_dist)
        down_genes = subset(table, logFC <= down_median_dist)
    } else if (is.null(z)) {
        upranked = table[order(table$logFC, decreasing=TRUE),]
        up_genes = head(upranked, n=n)
        down_genes = tail(upranked, n=n)
    }
    ret_list = list(up=up_genes, down=down_genes)
    return(ret_list)
}

#' make_exampledata()  A small hack of limma's exampleData()
#' function to allow for arbitrary data set sizes.
#'
#' @param ngenes default=1000  how many genes in the fictional data set.
#' @param columns default=5  how many samples in this data set.
#'
#' @return a matrix of pretend counts
#' @seealso \code{\link{makeExampleData}}
#' @export
#' @examples
#' ## pretend = make_exampledata()
make_exampledata = function (ngenes=1000, columns=5) {
    q0 <- rexp(ngenes, rate = 1/250)
    is_DE <- runif(ngenes) < 0.3
    lfc <- rnorm(ngenes, sd = 2)
    q0A <- ifelse(is_DE, q0 * 2^(lfc/2), q0)
    q0B <- ifelse(is_DE, q0 * 2^(-lfc/2), q0)
    ##    true_sf <- c(1, 1.3, 0.7, 0.9, 1.6)
    true_sf = abs(rnorm(columns, mean=1, sd=0.4))
    cond_types = ceiling(sqrt(columns))
    ##    conds <- c("A", "A", "B", "B", "B")
    ##x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
    conds = sample(LETTERS[1:cond_types], columns, replace=TRUE)
    m <- t(sapply(seq_len(ngenes), function(i) sapply(1:columns, function(j) rnbinom(1,
        mu = true_sf[j] * ifelse(conds[j] == "A", q0A[i], q0B[i]),
        size = 1/0.2))))
    rownames(m) <- paste("gene", seq_len(ngenes), ifelse(is_DE, "T", "F"), sep = "_")
    newCountDataSet(m, conds)
}

#' make_pairwise_contrasts()  Run makeContrasts() with all pairwise comparisons.
#'
#' @param model  a model describing the conditions/batches/etc in the experiment
#' @param conditions  a factor of conditions in the experiment
#' @param do_identities default=TRUE  whether or not to include all the identity strings.
#' Limma can handle this, edgeR cannot.
#' @param do_pairwise default=TRUE  whether or not to include all the pairwise strings.
#' This shouldn't need to be set to FALSE, but just in case.
#' @param extra_contrasts default=NULL  an optional string of extra contrasts to include.
#'
#' @return A list including the following information:
#'   all_pairwise_contrasts = the result from makeContrasts(...)
#'   identities = the string identifying each condition alone
#'   all_pairwise = the string identifying each pairwise comparison alone
#'   contrast_string = the string passed to R to call makeContrasts(...)
#'   names = the names given to the identities/contrasts
#'
#' @seealso \code{\link{makeContrasts}}
#' @export
#' @examples
#' ## pretend = make_pairwise_contrasts(model, conditions)
make_pairwise_contrasts = function(model, conditions, do_identities=TRUE, do_pairwise=TRUE, extra_contrasts=NULL) {
    tmpnames = colnames(model)
    tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    tmpnames = gsub("-", "", tmpnames)
    tmpnames = gsub("+", "", tmpnames)
    tmpnames = gsub("conditions", "", tmpnames)
    colnames(model) = tmpnames
    condition_table = table(conditions)
    identities = list()
    contrast_string = ""
    eval_strings = list()
    for (c in 1:length(condition_table)) {
        identity_name = names(condition_table[c])
        identity_string = paste(identity_name, " = ", identity_name, ",", sep="")
        identities[identity_name] = identity_string
        print(paste("As a reference, the identity is: ", identity_string, sep=""))
    }
    ## If I also create a sample condition 'alice', and also perform a subtraction
    ## of 'alice' from 'bob', then the full makeContrasts() will be:
    ## makeContrasts(bob=bob, alice=alice, bob_vs_alice=(bob)-(alice), levels=design)
    ## The parentheses in this case are superfluous, but they remind me that I am finally
    ## doing some match, and also they remind me that we can do much more complex things like:
    ## ((bob-alice)-(jane-alice)) or whatever....
    all_pairwise = list()
    identity_names = names(identities)
    lenminus = length(identities) - 1
    for (c in 1:lenminus) {
        c_name = names(identities[c])
        nextc = c+1
        for (d in nextc:length(identities)) {
            d_name = names(identities[d])
            minus_string = paste(d_name, "_vs_", c_name, sep="")
            exprs_string = paste(minus_string, "=", d_name, "-", c_name, ",", sep="")
            all_pairwise[minus_string] = exprs_string
        }
    }
    ## At this point, I have strings which represent the definition of every
    ## sample condition as well as strings which represent every possible
    ## B-A where B comes somewhere after A in the model matrix.
    ## The goal now is to create the variables in the R environment
    ## and add them to makeContrasts()
    if (isTRUE(do_identities)) {
        eval_strings = append(eval_strings, identities)
    }
    if (isTRUE(do_pairwise)) {
        eval_strings = append(eval_strings, all_pairwise)
    }
    eval_names = names(eval_strings)
    if (!is.null(extra_contrasts)) {
        extra_eval_strings = strsplit(extra_contrasts, "\\n")
        extra_eval_names = extra_eval_strings
        require.auto("stringi")
        extra_eval_names = stri_replace_all_regex(extra_eval_strings[[1]], "^(\\s*)(\\w+)=.*$", "$2")
        eval_strings = append(eval_strings, extra_contrasts)
    }
##    for (f in 1:length(eval_strings)) {
##        eval_name = names(eval_strings[f])
##        print(paste("Setting ", eval_name, " with expression:<<", eval_strings[f], ">>", sep=""))
##        eval(parse(text=as.character(eval_strings[f])))
##    }
    ## Now we have bob=(somestuff) in memory in R's environment
    ## Add them to makeContrasts()
    contrast_string = paste("all_pairwise_contrasts = makeContrasts(")
    for (f in 1:length(eval_strings)) {
        eval_name = names(eval_strings[f])
        eval_string = paste(eval_strings[f], sep="")
        contrast_string = paste(contrast_string, eval_string, sep="   ")
    }
    ## The final element of makeContrasts() is the design from voom()
    contrast_string = paste(contrast_string, "levels=model)")
    eval(parse(text=contrast_string))
    ## I like to change the column names of the contrasts because by default
    ## they are kind of obnoxious and too long to type

    if (!is.null(extra_contrasts)) {
        eval_names = append(eval_names, extra_eval_names)
    }
    colnames(all_pairwise_contrasts) = eval_names
    result = list(
        all_pairwise_contrasts = all_pairwise_contrasts,
        identities = identities,
        all_pairwise = all_pairwise,
        contrast_string = contrast_string,
        names = eval_names
        )
    return(result)
}

#' simple_comparison()  Perform a simple experimental/control comparison
#' This is a function written primarily to provide examples for how to use limma.
#' It does the following:  1.  Makes a model matrix using condition/batch
#' 2.  Optionally uses sva's combat (from cbcbSEQ)  3.  Runs voom/lmfit
#' 4.  Sets the first element of the design to "changed" and the second to "control".
#' 5.  Performs a makeContrasts() of changed - control.  6.  Fits them
#' 7.  Makes histograms of the two elements of the contrast, cor.tests() them,
#' makes a histogram of the p-values, ma-plot, volcano-plot, writes out the results in
#' an excel sheet, pulls the up/down significant and p-value significant (maybe this should be
#' replaced with write_limma()? 8.  And returns a list containining these data and plots.
#'
#' @param subset  an experimental subset with two conditions to compare.
#' @param workbook default='simple_comparison.xls'  an excel workbook to which to write.
#' @param worksheet default='simple_comparison'  an excel worksheet to which to write.
#' @param basename default=NA  a url to which to send click evens in clicky volcano/ma plots.
#' @param batch default=TRUE  whether or not to include batch in limma's model.
#' @param combat default=FALSE  whether or not to use combatMod().
#' @param combat_noscale default=TRUE  whether or not to include combat_noscale (makes combat a little less heavy-handed).
#' @param pvalue_cutoff default=0.05  p-value definition of 'significant.'
#' @param logfc_cutoff default=0.6  fold-change cutoff of significance. 0.6 on the low end and therefore 1.6 on the high.
#' @param tooltip_data default=NULL  text descriptions of genes if one wants google graphs.
#' @param verbose default=FALSE  be verbose?
#'
#' @return A list containing the following pieces:
#'   amean_histogram = a histogram of the mean values between the two conditions
#'   coef_amean_cor = a correlation test between the mean values and coefficients (this should be a p-value of 1)
#'   coefficient_scatter = a scatter plot of condition 2 on the y axis and condition 1 on x
#'   coefficient_x = a histogram of the x axis
#'   coefficient_y = a histogram of the y axis
#'   coefficient_both = a histogram of both
#'   coefficient_lm = a description of the line described by y=slope(y/x)+b where
#'   coefficient_lmsummary = the r-squared and such information for the linear model
#'   coefficient_weights = the weights against the linear model, higher weights mean closer to the line
#'   comparisons = the result from eBayes()
#'   contrasts = the result from contrasts.fit()
#'   contrast_histogram = a histogram of the coefficients
#'   downsignificant = a subset from toptable() of the 'down-regulated' genes (< 1 Z from the mean)
#'   fit = the result from lmFit(voom_result)
#'   ma_plot = an ma plot using the voom$E data and p-values
#'   psignificant = a subset from toptable() of all genes with p-values <= pvalue_cutoff
#'   pvalue_histogram = a histogram of all the p-values
#'   table = everything from toptable()
#'   upsignificant = a subset from toptable() of 'up-regulated' genes (> 1 Z from the mean)
#'   volcano_plot = a volcano plot of x/y
#'   voom_data = the result from calling voom()
#'   voom_plot = a plot from voom(), redunant with voom_data
#'
#' @seealso \code{\link{hpgl_gvis_ma_plot}}, \code{\link{toptable}},
#' \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{hpgl_voom}},
#' \code{\link{lmFit}}, \code{\link{makeContrasts}},
#' \code{\link{contrasts.fit}}
#'
#' @export
#' @examples
#' ## model = model.matrix(~ 0 + subset$conditions)
#' ## simple_comparison(subset, model)
#' ## Currently this assumes that a variant of toptable was used which
#' ## gives adjusted p-values.  This is not always the case and I should
#' ## check for that, but I have not yet.
simple_comparison = function(subset, workbook="simple_comparison.xls", sheet="simple_comparison", basename=NA, batch=TRUE, combat=FALSE, combat_noscale=TRUE, pvalue_cutoff=0.05, logfc_cutoff=0.6, tooltip_data=NULL, verbose=FALSE, ...) {
    condition_model = stats::model.matrix(~ 0 + subset$condition)
    if (length(levels(subset$batch)) == 1) {
        print("There is only one batch! I can only include condition in the model.")
        condbatch_model = stats::model.matrix(~ 0 + subset$condition)
    } else {
        condbatch_model = stats::model.matrix(~ 0 + subset$condition + subset$batch)
    }
    if (isTRUE(batch)) {
        model = condbatch_model
    } else {
        model = condition_model
    }
    expt_data = as.data.frame(exprs(subset$expressionset))
    if (combat) {
#        expt_data = ComBat(expt_data, subset$batches, condition_model)
        expt_data = cbcbSEQ::combatMod(expt_data, subset$batches, subset$conditions)
    }
    expt_voom = hpgltools::hpgl_voom(expt_data, model, libsize=subset$original_libsize, logged=subset$transform, converted=subset$convert)
    lf = limma::lmFit(expt_voom)
    colnames(lf$coefficients)
    coefficient_scatter = hpgltools::hpgl_linear_scatter(lf$coefficients)
    colnames(lf$design)[1] = "changed"
    colnames(lf$coefficients)[1] = "changed"
    colnames(lf$design)[2] = "control"
    colnames(lf$coefficients)[2] = "control"
    ## Now make sure there are no weird characters in the column names...
    if (length(colnames(lf$design)) >= 3) {
        for (counter in 3:length(colnames(lf$design))) {
            oldname = colnames(lf$design)[counter]
            newname = gsub("\\$","_", oldname, perl=TRUE)
            colnames(lf$design)[counter] = newname
            colnames(lf$coefficients)[counter] = newname
        }
    }
    contrast_matrix = limma::makeContrasts(changed_v_control = changed - control, levels=lf$design)
    cond_contrasts = contrasts.fit(lf, contrast_matrix)
    hist_df = data.frame(values=cond_contrasts$coefficients)
    contrast_histogram = hpgltools::hpgl_histogram(hist_df)
    hist_df = data.frame(values=cond_contrasts$Amean)
    amean_histogram = hpgltools::hpgl_histogram(hist_df, fillcolor="pink", color="red")
    coef_amean_cor = cor.test(cond_contrasts$coefficients, cond_contrasts$Amean, exact=FALSE)
    cond_comparison = limma::eBayes(cond_contrasts)
    hist_df = data.frame(values=cond_comparison$p.value)
    pvalue_histogram = hpgltools::hpgl_histogram(hist_df, fillcolor="lightblue", color="blue")
    cond_table = limma::topTable(cond_comparison, number=nrow(expt_voom$E), coef="changed_v_control", sort.by="logFC")
    if (!is.na(basename)) {
        vol_gvis_filename = paste(basename, "volplot.html", sep="_")
        a_volcano_plot = hpgltools::hpgl_volcano_plot(cond_table, gvis_filename=vol_gvis_filename, tooltip_data=tooltip_data)
    } else {
        a_volcano_plot = hpgltools::hpgl_volcano_plot(cond_table)
    }
    if (!is.na(basename)) {
        ma_gvis_filename=paste(basename, "maplot.html", sep="_")
        an_ma_plot = hpgltools::hpgl_ma_plot(expt_voom$E, cond_table, gvis_filename=ma_gvis_filename, tooltip_data=tooltip_data)
    } else {
        an_ma_plot = hpgltools::hpgl_ma_plot(expt_voom$E, cond_table)
    }
    hpgltools::write_xls(cond_table, sheet, file=workbook, rowname="row.names")
    upsignificant_table = subset(cond_table, logFC >=  logfc_cutoff)
    downsignificant_table = subset(cond_table, logFC <= (-1 * logfc_cutoff))
#    psignificant_table = subset(cond_table, adj.P.Val <= pvalue_cutoff)
    psignificant_table = subset(cond_table, P.Value <= pvalue_cutoff)

    if (verbose) {
        message("The model looks like:")
        message(model)
        message("The mean:variance trend follows")
        plot(expt_voom$plot)
        message("Drawing a scatterplot of the genes.")
        message("The following statistics describe the relationship between:")
        print(coefficient_scatter$scatter)
        message(paste("Setting the column:", colnames(lf$design)[2], "to control"))
        message(paste("Setting the column:", colnames(lf$design)[1], "to changed"))
        message("Performing contrasts of the experimental - control.")
        message("Taking a histogram of the subtraction values.")
        print(contrast_histogram)
        message("Taking a histogram of the mean values across samples.")
        message("The subtraction values should not be related to the mean values.")
        print(coef_amean_cor)
        message("Making a table of the data including p-values and F-statistics.")
        message("Taking a histogram of the p-values.")
        print(pvalue_histogram)
        message("Printing a volcano plot of this data.")
        message("Printing an maplot of this data.")
        message(paste("Writing excel sheet:", sheet))
    }
    return_info = list(
        amean_histogram=amean_histogram,
        coef_amean_cor=coef_amean_cor,
        coefficient_scatter=coefficient_scatter$scatter,
        coefficient_x=coefficient_scatter$x_histogram,
        coefficient_y=coefficient_scatter$y_histogram,
        coefficient_both=coefficient_scatter$both_histogram,
        coefficient_lm=coefficient_scatter$lm_model,
        coefficient_lmsummary=coefficient_scatter$lm_summary,
        coefficient_weights=coefficient_scatter$lm_weights,
        comparisons=cond_comparison,
        contrasts=cond_contrasts,
        contrast_histogram=contrast_histogram,
        downsignificant=downsignificant_table,
        fit=lf,
        ma_plot=an_ma_plot,
        psignificant=psignificant_table,
        pvalue_histogram=pvalue_histogram,
        table=cond_table,
        upsignificant=upsignificant_table,
        volcano_plot=a_volcano_plot,
        voom_data=expt_voom,
        voom_plot=expt_voom$plot)
    return(return_info)
}

#' basic_pairwise()  Perform a pairwise comparison among conditions which takes
#' nothing into account.  It _only_ takes the conditions, a mean value/variance among
#' them, divides by condition, and returns the result.  No fancy nomalizations, no
#' statistical models, no nothing.  It should be the very worst method possible.
#' But, it should also provide a baseline to compare the other tools against, they should
#' all do better than this, always.
#'
#' @param input a count table by sample
#' @param conditions a data frame of samples and conditions
#'
#' @return I am not sure yet
#'
#' @seealso \code{\link{limma}} \code{\link{deseq2}} \code{\link{edger}}
basic_pairwise = function(input, design=NULL) {
    print("Starting basic pairwise comparison.")
    input_class = class(input)[1]
    if (input_class == 'expt') {
        conditions = input$conditions
        batches = input$batches
        data = exprs(input$expressionset)
    } else {  ## Not an expt class, data frame or matrix
        data = as.data.frame(input)
        conditions = as.factor(design$condition)
    }
    ## design = data.frame(sample=c("hpgl010","hpgl011","hpgl012","hpgl013"), condition=c("a","b","c","a"))

    condition_table = table(conditions)
    types = levels(conditions)
    num_conds = length(types)
    median_table = data.frame()  ## This will be filled with num_conds columns and numRows(input) rows.
    variance_table = data.frame() ## This will be filled with num_conds columns and numRows(input) rows.
    ## First use conditions to rbind a table of medians by condition.
    for (c in 1:num_conds) {
        condition_name = types[c]
        columns = which(conditions == condition_name)
        if (length(columns) == 1) {
            med = data.frame(data[,columns])
            var = as.data.frame(matrix(NA, ncol=1, nrow=nrow(med)))
        } else {
            med_input = data[,columns]
            med = data.frame(Biobase::rowMedians(as.matrix(med_input)))
            var = matrixStats::rowVars(as.matrix(med_input))
        }

        if (c == 1) {
            median_table = med
            variance_table = var
        } else {
            median_table = cbind(median_table, med)
            variance_table = cbind(variance_table, var)
        }
    } ## end creation of median table / variance table
    colnames(median_table) = types
    colnames(variance_table) = types
    rownames(median_table) = rownames(input)
    rownames(variance_table) = rownames(input)
    ## We have tables of the median values by condition
    ## Now perform the pairwise comparisons

    comparisons = data.frame()
    lenminus = num_conds - 1
    num_done = 0
    column_list = c()
    for (c in 1:lenminus) {
        c_name = types[c]
        nextc = c+1
        for (d in nextc:length(types)) {
            num_done = num_done + 1
            d_name = types[d]
            message(paste0("Basic:", d, "/", c, ": Performing division: ", d_name, "_vs_", c_name))
            division = data.frame(median_table[,d] / median_table[,c])
            comparison_name = paste0(d_name, "_vs_", c_name)
            column_list = append(column_list, comparison_name)
            if (num_done == 1) {
                comparisons = division
            } else {
                comparisons = cbind(comparisons, division)
            }
        } ## End for each d
    }
    colnames(comparisons) = column_list
    comparisons[is.na(comparisons)] = 1
    rownames(comparisons) = rownames(input)

    all_tables = list()
    for (e in 1:length(colnames(comparisons))) {
        colname = colnames(comparisons)[[e]]
        column = comparisons[,e]
        column[mapply(is.infinite, column)] = 1
        column[column == 0] = 1
        tmpdf = cbind(column, log2(column))
        colnames(tmpdf) = c(colname, "logFC")
        rownames(tmpdf) = rownames(data)
        all_tables[[e]] = tmpdf
    }
    names(all_tables) = colnames(comparisons)

    retlist = list(
        input_data=data,
        conditions_table=table(conditions),
        conditions=conditions,
        all_pairwise=comparisons,
        all_tables=all_tables,
        medians=median_table,
        variances=variance_table)
    return(retlist)
}

#' write_limma()  Writes out the results of a limma search using toptable()
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data  the output from eBayes()
#' @param adjust default='fdr'  the pvalue adjustment chosen.
#' @param n default=0  the number of entries to report, 0 says do them all.
#' @param coef default=NULL  which coefficients/contrasts to report, NULL says do them all.
#' @param workbook default='excel/limma.xls'  an excel filename into which to write the data, used for csv files too.
#' @param excel default=FALSE  write an excel workbook?
#' @param csv default=TRUE  write out csv files of the tables?
#' @param annotation default=NULL  an optional data frame including annotation information to include with the tables.
#'
#' @return a list of data frames comprising the toptable output for each coefficient,
#'    I also added a qvalue entry to these toptable() outputs.
#'
#' @seealso \code{\link{toptable}}. \code{\link{write_xls}}
#'
#' @export
#' @examples
#' ## finished_comparison = eBayes(limma_output)
#' ## data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
write_limma = function(data, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls", excel=FALSE, csv=TRUE, annotation=NULL) {
    testdir = dirname(workbook)
    if (n == 0) {
        n = dim(data$coefficients)[1]
    }
    if (is.null(coef)) {
        coef = colnames(data$contrasts)
    } else {
        coef = as.character(coef)
    }
    return_data = list()
    end = length(coef)
    for (c in 1:end) {
        comparison = coef[c]
        message(paste0("limma:", c, "/", end, ": Printing table: ", comparison, "."))
        data_table = topTable(data, adjust=adjust, n=n, coef=comparison)

        data_table$qvalue = tryCatch(
            {
                as.numeric(format(signif(qvalue(data_table$P.Value, robust=TRUE)$qvalues, 4), scientific=TRUE))
            },
            error=function(cond) {
                message(paste("The qvalue estimation failed for ", comparison, ".", sep=""))
                return(1)
            },
            warning=function(cond) {
                message("There was a warning?")
                message(cond)
                return(1)
            },
            finally={
            }
        )
        data_table$P.Value = as.numeric(format(signif(data_table$P.Value, 4), scientific=TRUE))
        data_table$adj.P.Val = as.numeric(format(signif(data_table$adj.P.Val, 4), scientific=TRUE))

        if (!is.null(annotation)) {
            data_table = merge(data_table, annotation, by.x="row.names", by.y="row.names")
            ###data_table = data_table[-1]
        }
        ## This write_xls runs out of memory annoyingly often
        if (isTRUE(excel) | isTRUE(csv)) {
            if (!file.exists(testdir)) {
                dir.create(testdir)
                message(paste0("Creating directory: ", testdir, " for writing excel/csv data."))
            }
        }
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwrite=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename = gsub(".xls$", "", workbook)
            csv_filename = paste0(csv_filename, "_", comparison, ".csv")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] = data_table
    }
    return(return_data)
}

#' get_sig_genes()  Get a set of up/down genes using the top/bottom n or >/< z scores away from the median.
#'
#' @param table  a table from limma/edger/deseq.
#' @param n default=NULL  a rank-order top/bottom number of genes to take.
#' @param z default=NULL  a number of z-scores >/< the median to take.
#' @param fc default=NULL  a number of fold-changes to take
#' @param fold default='plusminus'  an identifier reminding how to get the bottom portion of a fold-change (plusminus says to get the negative of the positive, otherwise 1/positive is taken).
#' @param column default='logFC'  a column to use to distinguish top/bottom
#'
#' @return a list of up/down genes
#' @export
get_sig_genes = function(table, n=NULL, z=NULL, fc=NULL, column='logFC', fold='plusminus') {
    if (is.null(z) & is.null(n) & is.null(fc)) {
        print("No n, z, nor fc provided, setting z to 1.")
        z = 1
    }
    if (!is.null(n)) {
        ## Take a specific number of genes at the top/bottom of the rank ordered list.
        print("Getting the n genes up and down.")
        upranked = table[order(table[,column], decreasing=TRUE),]
        up_genes = head(upranked, n=n)
        down_genes = tail(upranked, n=n)
    } else if (!is.null(z)) {
        ## Take an arbitrary number which are >= and <= a value which is z zscores from the median.
        print("Getting the genes >= z scores away from the median.")
        out_summary = summary(table[,column])
        out_mad = mad(table[,column], na.rm=TRUE)
        up_median_dist = out_summary["Median"] + (out_mad * z)
        down_median_dist = out_summary["Median"] - (out_mad * z)
        up_idx = table[,column] >= up_median_dist
        up_genes = table[up_idx,]
        down_idx = table[,column] <= down_median_dist
        down_genes = table[down_idx,]
        print(paste0("The up genes table has ", dim(up_genes)[1], " genes."))
        print(paste0("The down genes table has ", dim(down_genes)[1], " genes."))
    } else {
        ## Take an arbitrary number which are >= and <= a given fold value
        up_idx = table[,column] >= fc
        up_genes = table[up_idx,]
        if (fold == 'plusminus') {
            ## plusminus refers to a positive/negative number of logfold changes from a logFC(1) = 0
            down_idx = table[,column] <= (fc * -1)
            down_genes = table[down_idx,]
        } else {
            ## If it isn't log fold change, then values go from 0..x where 1 is unchanged
            down_idx = table[,column] <= (1 / fc)
            down_genes = table[down_idx,]
        }
        print(paste0("The up genes table has ", dim(up_genes)[1], " genes."))
        print(paste0("The down genes table has ", dim(down_genes)[1], " genes."))
    }
    ret = list(up_genes=up_genes, down_genes=down_genes)
    return(ret)
}

## EOF
