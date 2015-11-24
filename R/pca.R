# Time-stamp: <Tue Nov 24 17:39:22 2015 Ashton Trey Belew (abelew@gmail.com)>

#' hpgl_pca()  Make a ggplot PCA plot describing the samples' clustering.
#'
#' @param data  an expt set of samples.
#' @param design default=NULL  a design matrix and.
#' @param colors default=NULL  a color scheme.
#' @param title default=NULL  a title for the plot.
#' @param size default=5  size for the glyphs on the plot.
#' @param labels default=NULL  add labels?  Also, what type?  FALSE, "default", or "fancy".
#'
#' @return a list containing the following:
#'   pca = the result of fast.svd()
#'   plot = ggplot2 pca_plot describing the principle component analysis of the samples.
#'   table = a table of the PCA plot data
#'   res = a table of the PCA res data
#'   variance = a table of the PCA plot variance
#' This makes use of cbcbSEQ and prints the table of variance by component.
#' @seealso \code{\link{makeSVD}}, \code{\link{pcRes}},
#' \code{\link{geom_dl}}
#'
#' @export
#' @examples
#' ## pca_plot = hpgl_pca(expt=expt)
#' ## pca_plot
hpgl_pca = function(data, colors=NULL, design=NULL, title=NULL, labels=NULL, size=5, ...) {
    ## hpgl_env = environment()
    data_class = class(data)[1]
    names = NULL
    if (data_class == 'expt') {
        design = data$design
        colors = data$colors
        names = data$names
        data = exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        data = exprs(data)
    } else if (data_class == 'list') {
        data = data$count_table
        if (is.null(data)) {
            stop("The list provided contains no count_table element.")
        }
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data = as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    if (is.null(labels)) {
        if (is.null(names)) {
            labels = colnames(data)
        } else {
            labels = names
        }
    } else if (labels[1] == 'boring') {
        if (is.null(names)) {
            labels = colnames(data)
        } else {
            labels = names
        }
    }

    if (is.null(design)) {
        print("No design was provided.  Making one with 1 condition, 1 batch.")
        design = cbind(labels, 1)
        design = as.data.frame(design)
        design$condition = as.numeric(design$labels)
        colnames(design) = c("name","batch","condition")
    }
    pca = hpgltools:::makeSVD(data)  ## This is a part of cbcbSEQ
    included_batches = as.factor(as.character(design$batch))
    included_conditions = as.factor(as.character(design$condition))
    if (length(levels(included_conditions)) == 1 & length(levels(included_batches)) == 1) {
        warning("There is only one condition and one batch, it is impossible to get meaningful pcRes information.")
    } else if (length(levels(included_conditions)) == 1) {
        warning("There is only one condition, but more than one batch.   Going to run pcRes with the batch information.")
        pca_res = cbcbSEQ::pcRes(pca$v, pca$d, design$batch)
    } else if (length(levels(included_batches)) == 1) {
        print("There is just one batch in this data.")
        pca_res = cbcbSEQ::pcRes(pca$v, pca$d, design$condition)
    } else {
        pca_res = cbcbSEQ::pcRes(pca$v, pca$d, design$condition, design$batch)
    }
    pca_variance = round((pca$d ^ 2) / sum(pca$d ^ 2) * 100, 2)
    xl = sprintf("PC1: %.2f%% variance", pca_variance[1])
    yl = sprintf("PC2: %.2f%% variance", pca_variance[2])
    if (is.null(colors)) {
        colors = as.numeric(levels(as.factor(design$condition)))
    }

    pca_data = data.frame(SampleID=labels,
        condition=as.character(design$condition),
        batch=as.character(design$batch),
        batch_int = as.integer(design$batch),
        PC1=pca$v[,1],
        PC2=pca$v[,2],
        colors=colors,
        labels=as.character(labels))

    num_batches = length(levels(included_batches))
    pca_plot = NULL
    if (num_batches <= 5) {
        pca_plot = pca_plot_smallbatch(pca_data, size=size, first='PC1', second='PC2')
    } else {
        pca_plot = pca_plot_largebatch(pca_data, size=size, first='PC1', second='PC2')
    }
    pca_plot = pca_plot + xlab(xl) + ylab(yl) + theme_bw() + theme(legend.key.size=unit(0.5, "cm"))
    if (!is.null(labels)) {
        if (labels[[1]] == "fancy") {
            pca_plot = pca_plot + directlabels::geom_dl(aes(label=SampleID), method="smart.grid")
        } else if (labels[[1]] == "normal") {
            pca_plot = pca_plot + geom_text(aes(x=PC1, y=PC2, label=paste(design$condition, design$batch, sep="_")), angle=45, size=4, vjust=2)
        } else {
            pca_plot = pca_plot + geom_text(aes(x=PC1, y=PC2, label=labels), angle=45, size=4, vjust=2)
        }
    }

    if (!is.null(title)) {
        pca_plot = pca_plot + ggtitle(title)
    }
    pca_return = list(
        pca=pca, plot=pca_plot, table=pca_data, res=pca_res, variance=pca_variance)
    return(pca_return)
}


#' pca_plot_largebatch()  ggplot2 plots of PCA data with >= 6 batches.
#'
#' @param df  A dataframe of PC1/PC2 and other arbitrary data.
#' @param size default=5  The size of glyphs in the plot.
#'
#' @return a ggplot2 plot of principle components 1 and 2.
pca_plot_largebatch = function(df, size=5, first="PC1", second="PC2") {
    hpgl_env = environment()
    ## The following 8 lines were written when I forgot to add environment= to ggplot() and it started throwing errors
    ## because it could not find the variable 'first', so I made it explicitly impossible for them to be nonexistent.
    ## get0 will set the value to null if the object does not exist, I then check for null and set the default value.
    ## This is entirely too verbose, but kind of nice if I will find myself debugging a bunch of code in the future.
    first = get0('first')
    second = get0('second')
    if (is.null(first)) {
        first = 'PC1'
    }
    if (is.null(second)) {
        second = 'PC2'
    }
    num_batches = length(levels(factor(df$batch)))
    plot = ggplot(df, aes(x=get(first), y=get(second)), environment=hpgl_env) +
        ## geom_point(size=3, aes(shape=factor(df$batch), fill=condition, colour=colors)) +
        geom_point(size=size, aes(shape=batch, fill=condition, colour=colors)) +
        scale_fill_manual(name="Condition", guide="legend", labels=levels(as.factor(df$condition)), values=levels(as.factor(df$colors))) +
        scale_color_manual(name="Condition", guide="legend", labels=levels(as.factor(df$condition)), values=levels(as.factor(df$colors))) +
        guides(fill=guide_legend(override.aes=list(colour=levels(factor(df$colors)))), colour=guide_legend(override.aes="black")) +
        scale_shape_manual(values=c(1:num_batches), name="Batch") + theme_bw()
    return(plot)
}

#' pca_plot_smallbatch()  ggplot2 plots of PCA data with <= 5 batches.
#'
#' This uses hard-coded scale_shape_manual values 21-25 to have solid shapes in the plot.
#' @param df  A dataframe of PC1/PC2 and other arbitrary data.
#' @param size default=5  The size of glyphs in the plot.
#'
#' @return a ggplot2 plot of principle components 1 and 2.
pca_plot_smallbatch = function(df, size=5, first='PC1', second='PC2') {
    hpgl_env = environment()
    plot = ggplot(df, aes(x=get(first), y=get(second)), environment=hpgl_env) +
        geom_point(size=size, aes(shape=factor(batch), fill=condition), colour='black') +
        scale_fill_manual(name="Condition", guide="legend", labels=levels(as.factor(df$condition)), values=levels(as.factor(df$colors))) +
        ##scale_fill_manual(name="Condition", guide="legend", labels=condition, values=colors) +
        scale_shape_manual(name="Batch", labels=levels(as.factor(df$batch)), values=21:25) +
        guides(fill=guide_legend(override.aes=list(colour=levels(factor(df$colors)))), colour=guide_legend(override.aes="black"))
    return(plot)
}

#' factor_rsquared()  Collect the r^2 values from a linear model fitting between a singular
#' value decomposition and factor.
#'
#' @param svd_v  the V' V = I portion of a fast.svd call.
#' @param factor  an experimental factor from the original data.
#'
#' @return The r^2 values of the linear model as a percentage.
#'
#' @seealso \code{\link{fast.svd}}
factor_rsquared = function(svd_v, factor) {
    svd_lm = try(lm(svd_v ~ factor), silent=TRUE)
    if (class(svd_lm) == 'try-error') {
        result = 0
    } else {
        lm_summary = summary.lm(svd_lm)
        r_squared = lm_summary$r.squared
        result = round(r_squared * 100, 3)
    }
    return(result)
}

#' plot_pcs()  A quick and dirty PCA plotter of arbitrary components against one another.
#'
#' @param data  a dataframe of principle components PC1 .. PCN with any other arbitrary information.
#' @param first default='PC1'  principle component PCx to put on the x axis.
#' @param second default='PC2'  principle component PCy to put on the y axis.
#' @param variances default=NULL  a list of the percent variance explained by each component.
#' @param design default=NULL  the experimental design with condition batch factors.
#' @param title default=NULL  a title for the plot.
#' @param labels default=NULL  a parameter for the labels on the plot.
#'
#' @return a ggplot2 PCA plot
#'
#' @seealso \code{\link{ggplot2}}, \code{\link{geom_dl}}
#'
#' @export
#' @examples
#' ## pca_plot = plot_pcs(pca_data, first="PC2", second="PC4", design=expt$design)
plot_pcs = function(data, first="PC1", second="PC2", variances=NULL, design=NULL, title=NULL, labels=NULL) {
    hpgl_env = environment()
    batches = design$batch
    point_labels = factor(design$condition)
    if (is.null(title)) {
        title = paste(first, " vs. ", second, sep="")
    }
    ## colors = levels(as.factor(unlist(design$color)))
    ## num_batches = length(levels(factor(design$batch)))

    pca_plot = ggplot(data=as.data.frame(data), environment=hpgl_env, fill=factor(design$condititon)) +
        geom_point(aes(x=get(first), y=get(second), shape=batches, colour=data$colors), stat="identity", size=3) +
        scale_color_manual(values=levels(factor(data$colors)), name="Condition", labels=levels(as.factor(design$condition))) +
        scale_shape_manual(values=batches, name="Batch", guide=guide_legend(override.aes=aes(size=1))) +
        ggtitle(title) +
        theme_bw() +
        theme(legend.key.size=unit(0.5, "cm"))

    if (!is.null(variances)) {
        x_var_num = as.numeric(gsub("PC", "", first))
        y_var_num = as.numeric(gsub("PC", "", second))
        x_label = paste("PC", x_var_num, ": ", variances[[x_var_num]], "%  variance", sep="")
        y_label = paste("PC", y_var_num, ": ", variances[[y_var_num]], "%  variance", sep="")
        pca_plot = pca_plot + xlab(x_label) + ylab(y_label)
    }

    if (!is.null(labels)) {
        if (labels[[1]] == "fancy") {
            pca_plot = pca_plot + geom_dl(aes(x=get(first), y=get(second), label=point_labels), list("top.bumpup", cex=0.5))
        } else {
            pca_plot = pca_plot + geom_text(aes(x=get(first), y=get(second), label=point_labels), angle=45, size=4, vjust=2)
        }
    }
    return(pca_plot)
}

## An alternate to plotting rank order of svd$u
## The plotted_u1s and such below
## y-axis is z(i), x-axis is i
## z(i) = cumulative sum of $u squared
## z = cumsum((svd$u ^ 2))

#' u_plot()  Plot the rank order svd$u elements to get a view of how much the first genes contribute to the total variance by PC.
#'
#' @param plotted_us  a list of svd$u elements
#'
#' @return a recordPlot() plot showing the first 3 PCs by rank-order svd$u.
u_plot = function(plotted_us) {
    plotted_us = abs(plotted_us[,c(1,2,3)])
    plotted_u1s = plotted_us[order(plotted_us[,1], decreasing=TRUE),]
    plotted_u2s = plotted_us[order(plotted_us[,2], decreasing=TRUE),]
    plotted_u3s = plotted_us[order(plotted_us[,3], decreasing=TRUE),]
    ## allS <- BiocGenerics::rank(allS, ties.method = "random")
    ## plotted_us$rank = rank(plotted_us[,1], ties.method="random")
    plotted_u1s = cbind(plotted_u1s, rev(rank(plotted_u1s[,1], ties.method="random")))
    plotted_u1s = plotted_u1s[,c(1,4)]
    colnames(plotted_u1s) = c("PC1","rank")
    plotted_u1s = data.frame(plotted_u1s)
    plotted_u1s$ID = as.character(rownames(plotted_u1s))
    plotted_u2s = cbind(plotted_u2s, rev(rank(plotted_u2s[,2], ties.method="random")))
    plotted_u2s = plotted_u2s[,c(2,4)]
    colnames(plotted_u2s) = c("PC2","rank")
    plotted_u2s = data.frame(plotted_u2s)
    plotted_u2s$ID = as.character(rownames(plotted_u2s))
    plotted_u3s = cbind(plotted_u3s, rev(rank(plotted_u3s[,3], ties.method="random")))
    plotted_u3s = plotted_u3s[,c(3,4)]
    colnames(plotted_u3s) = c("PC3","rank")
    plotted_u3s = data.frame(plotted_u3s)
    plotted_u3s$ID = as.character(rownames(plotted_u3s))
    plotted_us = merge(plotted_u1s, plotted_u2s, by.x="rank", by.y="rank")
    plotted_us = merge(plotted_us, plotted_u3s, by.x="rank", by.y="rank")
    colnames(plotted_us) = c("rank","PC1","ID1","PC2","ID2","PC3","ID3")
    rm(plotted_u1s)
    rm(plotted_u2s)
    rm(plotted_u3s)
    ## top_threePC = head(plotted_us, n=20)
    plotted_us = plotted_us[,c("PC1","PC2","PC3")]
    plotted_us$ID = rownames(plotted_us)
    message("The more shallow the curves in these plots, the more genes responsible for this principle component.")
    plot(plotted_us)
    u_plot = recordPlot()
    return(u_plot)
}

#' pca_information()  Gather information about principle components.
#'
#' Calculate some information useful for generating PCA plots.
#'
#' pca_information seeks to gather together interesting information
#' to make principle component analyses easier, including: the results
#' from (fast.)svd, a table of the r^2 values, a table of the
#' variances in the data, coordinates used to make a pca plot for an
#' arbitrarily large set of PCs, correlations and fstats between
#' experimental factors and the PCs, and heatmaps describing these
#' relationships.  Finally, it will provide a plot showing how much of
#' the variance is provided by the top-n genes and (optionally) the
#' set of all PCA plots with respect to one another. (PCx vs. PCy)
#'
#' @section Warning:
#'  This function has gotten too damn big and needs to be split up.
#'
#' @param data  the data to analyze (usually exprs(somedataset)).
#' @param design default=NULL  a dataframe describing the experimental design, containing columns with
#'   useful information like the conditions, batches, number of cells, whatever...
#' @param factors default=c("condition","batch")  a character list of experimental conditions to query
#'   for R^2 against the fast.svd of the data.
#' @param components default=NULL  a number of principle components to compare the design factors against.
#'   If left null, it will query the same number of components as factors asked for.
#' @param plot_pcas default=FALSE  plot the set of PCA plots for every pair of PCs queried.
#' @param labels default="fancy"  how to label the glyphs on the plot.
#'
#' @return a list of fun pca information:
#'   svd_u/d/v: The u/d/v parameters from fast.svd
#'   rsquared_table: A table of the rsquared values between each factor and principle component
#'   pca_variance: A table of the pca variances
#'   pca_data: Coordinates for a pca plot
#'   pca_cor: A table of the correlations between the factors and principle components
#'   anova_fstats: the sum of the residuals with the factor vs without (manually calculated)
#'   anova_f: The result from performing anova(withfactor, withoutfactor), the F slot
#'   anova_p: The p-value calculated from the anova() call
#'   anova_sums: The RSS value from the above anova() call
#'   cor_heatmap: A heatmap from recordPlot() describing pca_cor.
#' @seealso \code{\link{fast.svd}}, \code{\link{lm}}
#'
#' @export
#' @examples
#' ## pca_info = pca_information(exprs(some_expt$expressionset), some_design, "all")
#' ## pca_info
pca_information = function(data, design=NULL, factors=c("condition","batch"), num_components=NULL, plot_pcas=FALSE, labels="fancy") {
    ## hpgl_env = environment()
    data_class = class(data)[1]
    if (data_class == 'expt') {
        design = data$design
        data = exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        data = exprs(data)
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data = as.matrix(data)
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    data = as.matrix(data)
    means = rowMeans(data)
    decomposed = fast.svd(data - means)
    positives = decomposed$d
    u = decomposed$u
    v = decomposed$v
    ## A neat idea from Kwame, rank order plot the U's in the svd version of:
    ## [Covariates] = [U][diagonal][V] for a given PC (usually/always PC1)
    ## The idea being: the resulting decreasing line should be either a slow even
    ## decrease if many genes are contributing to the given component
    ## Conversely, that line should drop suddenly if dominated by one/few genes.

    rownames(u) = rownames(data)
    rownames(v) = colnames(data)

    u_plot = u_plot(u)

    component_variance = round((positives^2) / sum(positives^2) * 100, 3)
    cumulative_pc_variance = cumsum(component_variance)

    ## Include in this table the fstatistic and pvalue described in rnaseq_bma.rmd
    component_rsquared_table = data.frame(
        prop_var = component_variance,
        cumulative_prop_var = cumulative_pc_variance)

    if (is.null(factors)) {
        factors = colnames(design)
    } else if (factors[1] == "all") {
        factors = colnames(design)
    }

    for (component in factors) {
        comp = factor(as.character(design[,component]), exclude=FALSE)
        column = apply(v, 2, factor_rsquared, factor=comp)
        component_rsquared_table[component] = column
    }
    pca_variance = round((positives ^ 2) / sum(positives ^2) * 100, 2)
    xl = sprintf("PC1: %.2f%% variance", pca_variance[1])
    print(xl)
    yl = sprintf("PC2: %.2f%% variance", pca_variance[2])
    print(yl)
    labels = rownames(design)

    pca_data = data.frame(SampleID=labels,
        condition=design$condition,
        batch=design$batch,
        batch_int=as.integer(design$batch),
        colors=design$color)
    pc_df = data.frame(SampleID=labels)
    rownames(pc_df) = make.names(labels)

    if (is.null(num_components)) {
        num_components = length(factors)
    }
    for (pc in 1:num_components) {
        name = paste("PC", pc, sep="")
        pca_data[name] = v[,pc]
        pc_df[name] = v[,pc]

    }
    pc_df = pc_df[-1]

    pca_plots = list()
    if (isTRUE(plot_pcas)) {
        for (pc in 1:num_components) {
            next_pc = pc + 1
            name = paste("PC", pc, sep="")
            for (second_pc in next_pc:num_components) {
                if (pc < second_pc & second_pc <= num_components) {
                    second_name = paste("PC", second_pc, sep="")
                    list_name = paste(name, "_", second_name, sep="")
                    ## Sometimes these plots fail because too many grid operations are happening.
                    tmp_plot = try(print(plot_pcs(pca_data, design=design, variances=pca_variance, first=name, second=second_name, labels=labels)))
                    pca_plots[[list_name]] = tmp_plot
                }
            }
        }
    }

    factor_df = data.frame(SampleID=labels)
    rownames(factor_df) = make.names(labels)
    for (factor in factors) {
        factor_df[factor] = as.numeric(as.factor(as.character(design[,factor])))
    }
    factor_df = factor_df[-1]

    ## fit_one = data.frame()
    ## fit_two = data.frame()
    cor_df = data.frame()
    anova_rss = data.frame()
    anova_sums = data.frame()
    anova_f = data.frame()
    anova_p = data.frame()
    anova_rss = data.frame()
    anova_fstats = data.frame()
    for (factor in factors) {
        for (pc in 1:num_components) {
            factor_name = names(factor_df[factor])
            pc_name = names(pc_df[pc])
            tmp_df = merge(factor_df, pc_df, by="row.names")
            rownames(tmp_df) = tmp_df[,1]
            tmp_df = tmp_df[-1]

            lmwithfactor_test = try(lm(formula=get(pc_name) ~ 1 + get(factor_name), data=tmp_df))
            lmwithoutfactor_test = try(lm(formula=get(pc_name) ~ 1, data=tmp_df))
            ## This fstat provides a metric of how much variance is removed by including this specific factor
            ## in the model vs not.  Therefore higher numbers tell us that adding that factor
            ## removed more variance and are more important.
            fstat = sum(residuals(lmwithfactor_test)^2) / sum(residuals(lmwithoutfactor_test)^2)
            ##1.  Perform lm(pc ~ 1 + factor) which is fit1
            ##2.  Perform lm(pc ~ 1) which is fit2
            ##3.  The Fstat is then defined as (sum(residuals(fit1)^2) / sum(residuals(fit2)^2))
            ##4.  The resulting p-value is 1 - pf(Fstat, (n-(#levels in the factor)), (n-1))  ## n is the number of samples in the fit
            ##5.  Look at anova.test() to see if this provides similar/identical information
            another_fstat = try(anova(lmwithfactor_test, lmwithoutfactor_test), silent=TRUE)
            if (class(another_fstat)[1] == 'try-error') {
                anova_sums[factor,pc] = 0
                anova_f[factor,pc] = 0
                anova_p[factor,pc] = 0
                anova_rss[factor,pc] = 0
            } else {
                anova_sums[factor,pc] = another_fstat$S[2]
                anova_f[factor,pc] = another_fstat$F[2]
                anova_p[factor,pc] = another_fstat$P[2]
                anova_rss[factor,pc] = another_fstat$RSS[1]
            }
            anova_fstats[factor,pc] = fstat

            cor_test = NULL
            tryCatch(
                {
                    cor_test = cor.test(tmp_df[,factor_name], tmp_df[,pc_name], na.rm=TRUE)
                },
                error=function(cond) {
                    message(paste("The correlation failed for ", factor_name, " and ", pc_name, ".", sep=""))
                    cor_test = 0
                },
                warning=function(cond) {
                    message(paste("The standard deviation was 0 for ", factor_name, " and ", pc_name, ".", sep=""))
                },
                finally={
                }
            ) ## End of the tryCatch
            if (class(cor_test) == 'try-error') {
                cor_df[factor,pc] = 0
            } else {
                cor_df[factor,pc] = cor_test$estimate
            }
        }
    }
    rownames(cor_df) = colnames(factor_df)
    colnames(cor_df) = colnames(pc_df)
    colnames(anova_sums) = colnames(pc_df)
    colnames(anova_f) = colnames(pc_df)
    colnames(anova_p) = colnames(pc_df)
    colnames(anova_rss) = colnames(pc_df)
    colnames(anova_fstats) = colnames(pc_df)

    cor_df = as.matrix(cor_df)
    ##    silly_colors = grDevices::colorRampPalette(brewer.pal(9, "Purples"))(100)
    silly_colors = grDevices::colorRampPalette(c("purple","black","yellow"))(100)
    cor_df = cor_df[complete.cases(cor_df),]
    pc_factor_corheat = heatmap.3(cor_df, scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=silly_colors, dendrogram = "none", Rowv=FALSE, Colv=FALSE, main="cor(factor, PC)")
    pc_factor_corheat = recordPlot()

    anova_f_colors = grDevices::colorRampPalette(c("blue","black","red"))(100)
    anova_f_heat = heatmap.3(as.matrix(anova_f), scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=anova_f_colors, dendrogram = "none", Rowv=FALSE, Colv=FALSE, main="anova fstats for (factor, PC)")
    anova_f_heat = recordPlot()

    anova_fstat_colors = grDevices::colorRampPalette(c("blue","white","red"))(100)
    anova_fstat_heat = heatmap.3(as.matrix(anova_fstats), scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=anova_fstat_colors, dendrogram = "none", Rowv=FALSE, Colv=FALSE, main="anova fstats for (factor, PC)")
    anova_fstat_heat = recordPlot()

    neglog_p = -1 * log(as.matrix(anova_p) + 1)
    anova_neglogp_colors = grDevices::colorRampPalette(c("blue","white","red"))(100)
    anova_neglogp_heat = heatmap.3(as.matrix(neglog_p), scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=anova_f_colors, dendrogram = "none", Rowv=FALSE, Colv=FALSE, main="-log(anova_p values)")
    anova_neglogp_heat = recordPlot()
    ## Another option: -log10 p-value of the ftest for this heatmap.
    ## covariate vs PC score
    ## Analagously: boxplot(PCn ~ batch)
    pca_list = list(
        pc1_trend=u_plot,
        svd_d=positives, svd_u=u, svd_v=v, rsquared_table=component_rsquared_table,
        pca_variance=pca_variance, pca_data=pca_data, anova_fstats=anova_fstats,
        anova_sums=anova_sums, anova_f=anova_f, anova_p=anova_p,
        pca_cor=cor_df,
        cor_heatmap=pc_factor_corheat,
        anova_f_heatmap=anova_f_heat, anova_fstat_heatmap=anova_fstat_heat, anova_neglogp_heatmaph=anova_neglogp_heat,
        pca_plots=pca_plots
    )
    return(pca_list)
}

#' pca_highscores()  Get the highest/lowest scoring genes for every principle component.
#'
#' This function uses princomp to acquire a principle component biplot
#' for some data and extracts a dataframe of the top n genes for each
#' component by score.
#'
#' @param df default=NULL  a dataframe of (pseudo)counts
#' @param conditions default=NULL  a factor or character of conditions in the experiment.
#' @param batches default=NULL  a factor or character of batches in the experiment.
#' @param n default=20  the number of genes to extract.
#'
#' @return a list including the princomp biplot, histogram, and tables
#' of top/bottom n scored genes with their scores by component.
#' @seealso \code{\link{princomp}}
#' @export
#' @examples
#' ## information = pca_highscores(df=df, conditions=cond, batches=bat)
#' ## information$pca_bitplot  ## oo pretty
pca_highscores = function(df=NULL, conditions=NULL, batches=NULL, n=20) {
    ## Another method of using PCA
    ## cond = as.factor(as.numeric(conditions))
    ## batch = as.factor(as.numeric(batches))
    another_pca = try(princomp(x=df, cor=TRUE, scores=TRUE, formula=~0 + cond + batch))
    plot(another_pca)
    pca_hist = recordPlot()
    biplot(another_pca)
    pca_biplot = recordPlot()
    highest = NULL
    lowest = NULL
    for (pc in 1:length(colnames(another_pca$scores))) {
        tmphigh = another_pca$scores
        tmplow = another_pca$scores
        tmphigh = tmphigh[order(tmphigh[,pc], decreasing=TRUE),]
        tmphigh = head(tmphigh, n=20)
        tmplow = tmplow[order(tmplow[,pc], decreasing=FALSE),]
        tmplow = head(tmplow, n=20)
        high_column = paste0(signif(tmphigh[,pc], 4), ":", rownames(tmphigh))
        low_column = paste0(signif(tmplow[,pc], 4), ":", rownames(tmplow))
        highest = cbind(highest, high_column)
        lowest = cbind(lowest, low_column)
    }
    colnames(highest) = colnames(another_pca$scores)
    colnames(lowest) = colnames(another_pca$scores)
    ret_list = list(
        pca_hist=pca_hist,
        pca_biplot=pca_biplot,
        highest=highest,
        lowest=lowest)
    return(ret_list)
}

## EOF
