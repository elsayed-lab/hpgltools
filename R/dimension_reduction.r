#' Collect the r^2 values from a linear model fitting between a singular
#' value decomposition and factor.
#'
#' @param svd_v V' V = I portion of a fast.svd call.
#' @param fact Experimental factor from the original data.
#' @param type Make this categorical or continuous with factor/continuous.
#' @return The r^2 values of the linear model as a percentage.
#' @seealso \pkg{corpcor}
#'  \code{\link[corpcor]{fast.svd}}
#' @export
factor_rsquared <- function(svd_v, fact, type="factor") {
  if (type == "factor") {
    fact <- as.factor(fact)
  } else if (type == "numeric") {
    fact <- as.numeric(fact)
  } else {
    fact <- as.factor(as.numeric(fact))
  }
  ## FIXME! This is not the correct way to handle this
  if (length(levels(fact)) == length(fact)) {
    fact <- as.numeric(fact)
  }
  svd_lm <- try(stats::lm(svd_v ~ fact))
  if (class(svd_lm) == "try-error") {
    result <- 0
  } else {
    lm_summary <- stats::summary.lm(svd_lm)
    r_squared <- lm_summary[["r.squared"]]
    result <- round(r_squared * 100, 3)
  }
  return(result)
}

#' Gather information about principle components.
#'
#' Calculate some information useful for generating PCA plots.
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
#' @param expt_data  the data to analyze (usually exprs(somedataset)).
#' @param expt_design   a dataframe describing the experimental design, containing columns with
#'  useful information like the conditions, batches, number of cells, whatever...
#' @param expt_factors   a character list of experimental conditions to query
#'  for R^2 against the fast.svd of the data.
#' @param num_components   a number of principle components to compare the design factors against.
#'  If left null, it will query the same number of components as factors asked for.
#' @param plot_pcas   plot the set of PCA plots for every pair of PCs queried.
#' @param ...  Extra arguments for the pca plotter
#' @return a list of fun pca information:
#'  svd_u/d/v: The u/d/v parameters from fast.svd
#'  rsquared_table: A table of the rsquared values between each factor and principle component
#'  pca_variance: A table of the pca variances
#'  pca_data: Coordinates for a pca plot
#'  pca_cor: A table of the correlations between the factors and principle components
#'  anova_fstats: the sum of the residuals with the factor vs without (manually calculated)
#'  anova_f: The result from performing anova(withfactor, withoutfactor), the F slot
#'  anova_p: The p-value calculated from the anova() call
#'  anova_sums: The RSS value from the above anova() call
#'  cor_heatmap: A heatmap from recordPlot() describing pca_cor.
#' @seealso \pkg{corpcor} \pkg{stats}
#'  \code{\link[corpcor]{fast.svd}}, \code{\link[stats]{lm}}
#' @examples
#' \dontrun{
#'  pca_info = pca_information(exprs(some_expt$expressionset), some_design, "all")
#'  pca_info
#' }
#' @export
pca_information <- function(expt_data, expt_design=NULL, expt_factors=c("condition", "batch"),
                            num_components=NULL, plot_pcas=FALSE, ...) {
  ## Start out with some sanity tests
  colors_chosen <- NULL
  exprs_data <- NULL
  data_class <- class(expt_data)[1]
  if (data_class == "expt") {
    expt_design <- expt_data[["design"]]
    colors_chosen <- expt_data[["colors"]]
    exprs_data <- exprs(expt_data[["expressionset"]])
  } else if (data_class == "ExpressionSet") {
    exprs_data <- exprs(expt_data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    exprs_data <- as.matrix(expt_data)
  } else {
    stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Make sure colors get chosen.
  if (is.null(colors_chosen)) {
    colors_chosen <- as.numeric(as.factor(expt_design[["condition"]]))
    num_chosen <- max(3, length(levels(as.factor(colors_chosen))))
    colors_chosen <- RColorBrewer::brewer.pal(num_chosen, "Dark2")[colors_chosen]
    names(colors_chosen) <- rownames(expt_design)
  }

  ## Extract the various pieces of data we will use later.
  expt_means <- rowMeans(exprs_data)
  decomposed <- corpcor::fast.svd(exprs_data - expt_means)
  positives <- decomposed[["d"]]
  u <- decomposed[["u"]]
  v <- decomposed[["v"]]
  ## A neat idea from Kwame, rank order plot the U's in the svd version of:
  ## [Covariates] = [U][diagonal][V] for a given PC (usually/always PC1)
  ## The idea being: the resulting decreasing line should be either a slow even
  ## decrease if many genes are contributing to the given component
  ## Conversely, that line should drop suddenly if dominated by one/few genes.
  rownames(u) <- rownames(exprs_data)
  rownames(v) <- colnames(exprs_data)
  u_plot <- u_plot(u)

  ## Extract the variance for each of the included PCs
  ## Also extract the r^2 values by component.
  component_variance <- round((positives ^ 2) / sum(positives ^ 2) * 100, 3)
  cumulative_pc_variance <- cumsum(component_variance)
  if (is.null(expt_factors)) {
    expt_factors <- colnames(expt_design)
    expt_factors <- expt_factors[expt_factors != "sampleid"]
  } else if (expt_factors[1] == "all") {
    expt_factors <- colnames(expt_design)
    expt_factors <- expt_factors[expt_factors != "sampleid"]
  }
  component_rsquared_table <- data.frame(
    "prop_var" = component_variance,
    "cumulative_prop_var" = cumulative_pc_variance)
  ## Fill in the R squared table.
  for (component in expt_factors) {
    comp <- expt_design[[component]]
    if (is.null(comp)) {
      message("The given component is not in the design: ", comp)
      next
    }
    column <- apply(v, 2, factor_rsquared, fact=comp)
    component_rsquared_table[[component]] <- column
  }

  pca_variance <- round((positives ^ 2) / sum(positives ^ 2) * 100, 2)
  xl <- sprintf("PC1: %.2f%% variance", pca_variance[1])
  yl <- sprintf("PC2: %.2f%% variance", pca_variance[2])

  if (is.null(expt_design[["batch"]])) {
    expt_design[["batch"]] <- "undefined"
  }

  ## Create an initial dataframe to be used by ggplot for plotting the various PCs
  pca_data <- data.frame(
    "sampleid" = rownames(expt_design),
    "labels" = rownames(expt_design),
    "condition" = as.character(expt_design[["condition"]]),
    "colors" = colors_chosen,
    "batch" = as.character(expt_design[["batch"]]),
    "batch_int" = as.integer(as.factor(expt_design[["batch"]])),
    stringsAsFactors=FALSE)

  pc_df <- data.frame(
    "sampleid" = rownames(expt_design))
  rownames(pc_df) <- rownames(expt_design)

  if (is.null(num_components)) {
    num_components <- length(expt_factors)
  }
  max_components <- ncol(v)
  if (max_components < num_components) {
    message("The u and v components of SVD have only ", max_components,
            " columns, but the list of factors is ", num_components, " long.")
    message("Therefore, only searching for ", max_components, " PCs.")
    num_components <- max_components
  }

  ## Now fill in the pca_df with the data from the various PCs
  for (pc in 1:num_components) {
    name <- paste("PC", pc, sep="")
    ## v is a matrix, don't forget that.
    pca_data[[name]] <- v[, pc] ## note you _must_ not shortcut this with [[pc]]
    pc_df[[name]] <- v[, pc]
  }
  pc_df <- pc_df[, -1, drop=FALSE]
  ## Now that we have filled in a pca data frame, we may plot PCx vs PCy for all x,y.
  pca_plots <- list()
  if (isTRUE(plot_pcas)) {
    for (pc in 1:num_components) {
      next_pc <- pc + 1
      name <- paste("PC", pc, sep="")
      for (second_pc in next_pc:num_components) {
        if (pc < second_pc & second_pc <= num_components) {
          second_name <- paste("PC", second_pc, sep="")
          list_name <- paste(name, "_", second_name, sep="")
          ## Sometimes these plots fail because too many grid operations are happening.
          tmp_plot <- try(plot_pcs(pca_data,
                                   design=expt_design,
                                   variances=pca_variance,
                                   first=name,
                                   second=second_name))
          pca_plots[[list_name]] <- tmp_plot
        }
      }
    }
  }

  ## Now start filling in data which may be used for correlations/fstats/etc.
  factor_df <- data.frame(
    "sampleid" = rownames(expt_design))
  rownames(factor_df) <- rownames(expt_design)
  for (fact in expt_factors) {
    if (!is.null(expt_design[[fact]])) {
      factor_df[[fact]] <- as.numeric(as.factor(as.character(expt_design[, fact])))
    } else {
      message("The column ", fact, " seems to be missing from the design.")
      message("The available columns are: ", toString(colnames(expt_design)), ".")
    }
  }
  factor_df <- factor_df[, -1, drop=FALSE]

  ## Perform the correlations/fstats/anova here
  cor_df <- data.frame()
  anova_rss <- data.frame()
  anova_sums <- data.frame()
  anova_f <- data.frame()
  anova_p <- data.frame()
  anova_rss <- data.frame()
  anova_fstats <- data.frame()
  for (fact in expt_factors) {
    for (pc in 1:num_components) {
      factor_name <- names(factor_df[fact])
      pc_name <- names(pc_df[pc])
      tmp_df <- merge(factor_df, pc_df, by="row.names")
      rownames(tmp_df) <- tmp_df[, 1]
      tmp_df <- tmp_df[, -1, drop=FALSE]
      lmwithfactor_test <- try(stats::lm(formula=get(pc_name) ~ 1 + get(factor_name), data=tmp_df))
      lmwithoutfactor_test <- try(stats::lm(formula=get(pc_name) ~ 1, data=tmp_df))
      ## This fstat provides a metric of how much variance is removed by including this specific factor
      ## in the model vs not.  Therefore higher numbers tell us that adding that factor
      ## removed more variance and are more important.
      fstat <- sum(residuals(lmwithfactor_test) ^ 2) / sum(residuals(lmwithoutfactor_test) ^ 2)
      ##1.  Perform lm(pc ~ 1 + factor) which is fit1
      ##2.  Perform lm(pc ~ 1) which is fit2
      ##3.  The Fstat is then defined as (sum(residuals(fit1)^2) / sum(residuals(fit2)^2))
      ##4.  The resulting p-value is 1 - pf(Fstat, (n-(#levels in the factor)), (n-1))
      ##    n is the number of samples in the fit
      ##5.  Look at anova.test() to see if this provides similar/identical information
      another_fstat <- try(stats::anova(lmwithfactor_test, lmwithoutfactor_test), silent=TRUE)
      if (class(another_fstat)[1] == "try-error") {
        anova_sums[fact, pc] <- 0
        anova_f[fact, pc] <- 0
        anova_p[fact, pc] <- 0
        anova_rss[fact, pc] <- 0
      } else {
        anova_sums[fact, pc] <- another_fstat$S[2]
        anova_f[fact, pc] <- another_fstat$F[2]
        anova_p[fact, pc] <- another_fstat$P[2]
        anova_rss[fact, pc] <- another_fstat$RSS[1]
      }
      anova_fstats[fact, pc] <- fstat
      cor_test <- NULL
      tryCatch({
        cor_test <- cor.test(tmp_df[, factor_name], tmp_df[, pc_name], na.rm=TRUE)
      },
      error=function(cond) {
        message(paste("The correlation failed for ", factor_name, " and ", pc_name, ".", sep=""))
        cor_test <- 0
      },
      warning=function(cond) {
        message(paste("The standard deviation was 0 for ", factor_name, " and ", pc_name, ".", sep=""))
      },
      finally={
      }
      ) ## End of the tryCatch
      if (class(cor_test) == "try-error" | is.null(cor_test)) {
        cor_df[fact, pc] <- 0
      } else {
        cor_df[fact, pc] <- cor_test$estimate
      }
    }
  }
  ## Sanitize the resulting matrices and get them ready for plotting.
  rownames(cor_df) <- colnames(factor_df)
  colnames(cor_df) <- colnames(pc_df)
  colnames(anova_sums) <- colnames(pc_df)
  colnames(anova_f) <- colnames(pc_df)
  colnames(anova_p) <- colnames(pc_df)
  colnames(anova_rss) <- colnames(pc_df)
  colnames(anova_fstats) <- colnames(pc_df)
  ## Finally, plot them.
  silly_colors <- grDevices::colorRampPalette(c("purple", "black", "yellow"))(100)
  cor_df <- cor_df[complete.cases(cor_df), ]
  pc_factor_corheat <- heatmap.3(as.matrix(cor_df), scale="none", trace="none",
                                 linewidth=0.5, keysize=2, margins=c(8, 8),
                                 col=silly_colors, dendrogram="none", Rowv=FALSE,
                                 Colv=FALSE, main="cor(factor, PC)")
  pc_factor_corheat <- grDevices::recordPlot()
  anova_f_colors <- grDevices::colorRampPalette(c("blue", "black", "red"))(100)
  anova_f_heat <- heatmap.3(as.matrix(anova_f), scale="none", trace="none",
                            linewidth=0.5, keysize=2, margins=c(8, 8),
                            col=anova_f_colors, dendrogram = "none", Rowv=FALSE,
                            Colv=FALSE, main="anova fstats for (factor, PC)")
  anova_f_heat <- grDevices::recordPlot()
  anova_fstat_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  anova_fstat_heat <- heatmap.3(as.matrix(anova_fstats), scale="none", trace="none",
                                linewidth=0.5, keysize=2, margins=c(8, 8),
                                col=anova_fstat_colors, dendrogram="none", Rowv=FALSE,
                                Colv=FALSE, main="anova fstats for (factor, PC)")
  anova_fstat_heat <- grDevices::recordPlot()
  neglog_p <- -1 * log(as.matrix(anova_p) + 1)
  anova_neglogp_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  anova_neglogp_heat <- heatmap.3(as.matrix(neglog_p), scale="none", trace="none",
                                  linewidth=0.5, keysize=2, margins=c(8, 8),
                                  col=anova_f_colors, dendrogram="none", Rowv=FALSE,
                                  Colv=FALSE, main="-log(anova_p values)")
  anova_neglogp_heat <- grDevices::recordPlot()
  ## Another option: -log10 p-value of the ftest for this heatmap.
  ## covariate vs PC score
  ## Analagously: boxplot(PCn ~ batch)

  ## Finally, return the set of materials collected.
  pca_list <- list(
    "pc1_trend" = u_plot,
    "svd_d" = positives,
    "svd_u" = u,
    "svd_v" = v,
    "rsquared_table" = component_rsquared_table,
    "pca_variance" = pca_variance,
    "pca_data" = pca_data,
    "anova_fstats" = anova_fstats,
    "anova_sums" = anova_sums,
    "anova_f" = anova_f,
    "anova_p" = anova_p,
    "pca_cor" = cor_df,
    "cor_heatmap" = pc_factor_corheat,
    "anova_f_heatmap" = anova_f_heat,
    "anova_fstat_heatmap" = anova_fstat_heat,
    "anova_neglogp_heatmaph" = anova_neglogp_heat,
    "pca_plots" = pca_plots)
  return(pca_list)
}

#' Get the highest/lowest scoring genes for every principle component.
#'
#' This function uses princomp to acquire a principle component biplot
#' for some data and extracts a dataframe of the top n genes for each
#' component by score.
#'
#' @param expt Experiment to poke.
#' @param n Number of genes to extract.
#' @param cor Perform correlations?
#' @param vs Do a mean or median when getting ready to perform the pca?
#' @param logged  Check for the log state of the data and adjust as deemed necessary?
#' @return a list including the princomp biplot, histogram, and tables
#'  of top/bottom n scored genes with their scores by component.
#' @seealso \pkg{stats}
#'  \code{\link[stats]{princomp}}
#' @examples
#' \dontrun{
#'  information <- pca_highscores(df=df, conditions=cond, batches=bat)
#'  information$pca_bitplot  ## oo pretty
#' }
#' @export
pca_highscores <- function(expt, n=20, cor=TRUE, vs="means", logged=TRUE) {
  if (isTRUE(logged)) {
    if (expt[["state"]][["transform"]] == "raw") {
      expt <- sm(normalize_expt(expt, transform="log2"))
    }
  }

  data <- as.data.frame(exprs(expt))
  if (!is.null(vs)) {
    if (vs == "means") {
      data <- as.matrix(data) - rowMeans(as.matrix(data))
    } else if (vs == "medians") {
      data <- as.matrix(data) - rowMedians(as.matrix(data))
    }
  }
  another_pca <- try(princomp(x=data, cor=cor))
  plot(another_pca)
  pca_hist <- grDevices::recordPlot()
  biplot(another_pca)
  pca_biplot <- grDevices::recordPlot()
  highest <- NULL
  lowest <- NULL
  for (pc in 1:length(colnames(another_pca[["scores"]]))) {
    tmphigh <- another_pca[["scores"]]
    tmplow <- another_pca[["scores"]]
    tmphigh <- tmphigh[order(tmphigh[, pc], decreasing=TRUE), ]
    tmphigh <- head(tmphigh, n=n)
    tmplow <- tmplow[order(tmplow[, pc], decreasing=FALSE), ]
    tmplow <- head(tmplow, n=n)
    high_column <- paste0(signif(tmphigh[, pc], 4), ":", rownames(tmphigh))
    low_column <- paste0(signif(tmplow[, pc], 4), ":", rownames(tmplow))
    highest <- cbind(highest, high_column)
    lowest <- cbind(lowest, low_column)
  }
  colnames(highest) <- colnames(another_pca[["scores"]])
  colnames(lowest) <- colnames(another_pca[["scores"]])
  retlist <- list(
    "scores" = another_pca[["scores"]],
    "pca_hist" = pca_hist,
    "pca_biplot" = pca_biplot,
    "highest" = highest,
    "lowest" = lowest,
    "result" = another_pca)
  retlist[["score_heat"]] <- plot_disheat(expt_data=another_pca[["scores"]])
  return(retlist)
}

#' Compute variance of each principal component and how they correlate with batch and cond
#'
#' This was copy/pasted from cbcbSEQ
#' https://github.com/kokrah/cbcbSEQ/blob/master/R/explore.R
#'
#' @param v from makeSVD
#' @param d from makeSVD
#' @param condition factor describing experiment
#' @param batch factor describing batch
#' @return A dataframe containig variance, cum. variance, cond.R-sqrd, batch.R-sqrd
#' @seealso \code{\link{plot_pca}}
#' @export
pcRes <- function(v, d, condition=NULL, batch=NULL){
  pcVar <- round((d ^ 2) / sum(d ^ 2) * 100, 2)
  cumPcVar <- cumsum(pcVar)
  calculate_rsquared_condition <- function(data) {
    lm_result <- lm(data ~ condition)
  }
  if (!is.null(condition)) {
    cond.R2 <- function(y) {
      round(summary(lm(y ~ condition))$r.squared * 100, 2)
    }
    cond.R2 <- apply(v, 2, cond.R2)
  }
  if (!is.null(batch)) {
    batch.R2 <- function(y) {
      round(summary(lm(y ~ batch))$r.squared * 100, 2)
    }
    batch.R2 <- apply(v, 2, batch.R2)
  }
  if (is.null(condition) & is.null(batch)) {
    res <- data.frame("propVar" = pcVar,
                      "cumPropVar" = cumPcVar)
  }
  if (!is.null(batch) & is.null(condition)) {
    res <- data.frame("propVar" = pcVar,
                      "cumPropVar" = cumPcVar,
                      "batch.R2" = batch.R2)
  }
  if (!is.null(condition) & is.null(batch)) {
    res <- data.frame("propVar" = pcVar,
                      "cumPropVar" = cumPcVar,
                      "cond.R2" = cond.R2)
  }
  if (!is.null(condition) & !is.null(batch)) {
    res <- data.frame("propVar" = pcVar,
                      "cumPropVar" = cumPcVar,
                      "cond.R2" = cond.R2,
                      "batch.R2" = batch.R2)
  }
  return(res)
}

#' A quick and dirty PCA plotter of arbitrary components against one another.
#'
#' @param pca_data  Dataframe of principle components PC1 .. PCN with any other arbitrary information.
#' @param first   Principle component PCx to put on the x axis.
#' @param second   Principle component PCy to put on the y axis.
#' @param variances   List of the percent variance explained by each component.
#' @param design   Experimental design with condition batch factors.
#' @param plot_title   Title for the plot.
#' @param plot_labels   Parameter for the labels on the plot.
#' @param plot_size  Size of the dots on the plot
#' @param size_column  Experimental factor to use for sizing the glyphs
#' @param rug  Include the rugs on the sides of the plot?
#' @param cis  What (if any) confidence intervals to include.
#' @param ...  Extra arguments dropped into arglist
#' @return  gplot2 PCA plot
#' @seealso \pkg{ggplot2}
#'  \code{\link[directlabels]{geom_dl}}
#' @examples
#' \dontrun{
#'  pca_plot = plot_pcs(pca_data, first="PC2", second="PC4", design=expt$design)
#' }
#' @export
plot_pcs <- function(pca_data, first="PC1", second="PC2", variances=NULL,
                     design=NULL, plot_title=TRUE, plot_labels=NULL,
                     x_label=NULL, y_label=NULL,
                     plot_size=5, plot_alpha=NULL, size_column=NULL, rug=TRUE,
                     cis=c(0.95, 0.9), ...) {
  arglist <- list(...)
  batches <- pca_data[["batch"]]
  label_column <- "condition"
  if (!is.null(arglist[["label_column"]])) {
    label_column <- arglist[["label_column"]]
  }
  point_labels <- factor(pca_data[[label_column]])
  if (is.null(plot_title)) {
    plot_title <- paste(first, " vs. ", second, sep="")
  }
  num_batches <- length(unique(batches))
  if (!is.null(arglist[["base_size"]])) {
    base_size <<- arglist[["base_size"]]
  }
  label_size <- 4
  if (!is.null(arglist[["label_size"]])) {
    label_size <<- arglist[["label_size"]]
  }
  ci_group <- "condition"
  if (!is.null(arglist[["ci_group"]])) {
    ci_group <- arglist[["ci_group"]]
  }
  ci_fill <- "condition"
  if (!is.null(arglist[["ci_fill"]])) {
    ci_fill <- arglist[["ci_fill"]]
  }

  pca_plot <- NULL

  color_listing <- pca_data[, c("condition", "colors")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["colors"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  ## Ok, so this is shockingly difficult.  For <5 batch data I want properly colored points with black outlines
  ## The legend colors need to match, in addition, the legend needs to have the shapes noted.
  ## In order to do this, one must do, _in_order_:
  ## 1.  Set up the normal ggplot object
  ## 2.  Set up a geom_point with color _and_ fill as the proper color.
  ##     The color but _NOT_ fill is used to color the legend's copy of the glyph.
  ## 3.  Then set up a new geom_point with color=black _and_ show_guide=FALSE
  ## 4.  Then set scale_color_manual to the proper color_list
  ## 5.  Then set scale_fill_manual to the proper color_list
  ## 6.  Finally, set the shape manual with a guide_legend override

  ## Step 1

  if (is.null(plot_alpha)) {
    plot_alpha <- 1
  }

  pca_plot <- ggplot(data=as.data.frame(pca_data),
                     aes_string(
                       x="get(first)",
                       y="get(second)",
                       text="sampleid"))
  minimum_size <- 2
  maximum_size <- 2
  if (!is.null(size_column)) {
    maximum_size <- max(levels(pca_data[["size"]]))
  }

  if (is.null(size_column) & num_batches <= 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(size=plot_size, alpha=plot_alpha,
                          aes_string(shape="as.factor(batches)",
                                     colour="as.factor(condition)",
                                     fill="as.factor(condition)")) +
      ggplot2::geom_point(size=plot_size, alpha=plot_alpha, colour="black", show.legend=FALSE,
                          aes_string(shape="as.factor(batches)",
                                     fill="as.factor(condition)")) +
      ggplot2::scale_color_manual(name="Condition",
                                  guide="legend",
                                  values=color_list) +
      ggplot2::scale_fill_manual(name="Condition",
                                 guide="legend",
                                 values=color_list) +
      ggplot2::scale_shape_manual(
                 name="Batch",
                 labels=levels(as.factor(pca_data[["batch"]])),
                 guide=ggplot2::guide_legend(override.aes=list(size=plot_size, fill="grey")),
                 values=21:25)
  } else if (is.null(size_column) & num_batches > 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(size=plot_size, alpha=plot_alpha,
                          aes_string(shape="as.factor(batches)",
                                     fill="as.factor(condition)",
                                     colour="as.factor(condition)")) +
      ggplot2::scale_color_manual(name="Condition",
                                  guide="legend",
                                  values=color_list) +
      ggplot2::scale_fill_manual(name="Condition",
                                 guide=ggplot2::guide_legend(override.aes=list(size=plot_size)),
                                 values=color_list) +
      ggplot2::scale_shape_manual(name="Batch",
                                  labels=levels(as.factor(pca_data[["batch"]])),
                                  guide=ggplot2::guide_legend(overwrite.aes=list(size=plot_size)),
                                  values=1:num_batches)
  } else if (!is.null(size_column) & num_batches <= 5) {
    ## This will require the 6 steps above and one more
    pca_plot <- pca_plot +
      ggplot2::geom_point(alpha=plot_alpha,
                 aes_string(shape="as.factor(batches)",
                            size="size",
                            colour="as.factor(condition)",
                            fill="as.factor(condition)")) +
      ggplot2::geom_point(colour="black", alpha=plot_alpha, show.legend=FALSE,
                          aes_string(shape="as.factor(batches)",
                                     size="size",
                                     fill="as.factor(condition)")) +
      ggplot2::scale_color_manual(name="Condition",
                                  guide="legend",
                                  values=color_list) +
      ggplot2::scale_fill_manual(name="Condition",
                                 guide=ggplot2::guide_legend(override.aes=list(size=plot_size)),
                                 values=color_list) +
      ggplot2::scale_shape_manual(
                 name="Batch",
                 labels=levels(as.factor(pca_data[["batch"]])),
                 guide=ggplot2::guide_legend(override.aes=list(size=plot_size, fill="grey")),
                 values=21:25) +
      ggplot2::scale_size_manual(name=size_column,
                                 labels=levels(pca_data[[size_column]]),
                                 values=as.numeric(levels(pca_data[["size"]])))
  } else if (!is.null(size_column) & num_batches > 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(alpha=plot_alpha,
                          aes_string(shape="as.factor(batches)",
                                     colour="pca_data[['condition']]",
                                     size="size")) +
      ggplot2::scale_shape_manual(name="Batch",
                                  labels=levels(as.factor(pca_data[["batch"]])),
                                  guide=ggplot2::guide_legend(overwrite.aes=list(size=plot_size)),
                                  values=1:num_batches) +
      ggplot2::scale_color_identity(name="Condition",
                                    guide="legend",
                                    values=color_list) +
      ggplot2::scale_size_manual(name=size_column,
                                 labels=levels(pca_data[[size_column]]),
                                 values=as.numeric(levels(pca_data[["size"]])))
  } else {
    stop("This should be an impossible state.")
  }

  if (!is.null(x_label)) {
    pca_plot <- pca_plot +
      ggplot2::xlab(x_label)
    if (!is.null(y_label)) {
      pca_plot <- pca_plot +
        ggplot2::ylab(y_label)
    }
  } else if (!is.null(variances)) {
    x_var_num <- as.numeric(gsub("PC", "", first))
    y_var_num <- as.numeric(gsub("PC", "", second))
    x_label <- paste0("PC1", first, ": ", variances[[x_var_num]], "%  variance")
    y_label <- paste0("PC2", second, ": ", variances[[y_var_num]], "%  variance")
    pca_plot <- pca_plot +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label)
  }

  if (isTRUE(rug)) {
    pca_plot <- pca_plot + ggplot2::geom_rug(colour="gray50", alpha=0.7)
  }

  if (is.null(plot_labels)) {
    plot_labels <- "repel"
  }
  if (plot_labels == FALSE) {
    message("Not putting labels on the plot.")
  } else if (plot_labels == "normal") {
    pca_plot <- pca_plot +
      ggplot2::geom_text(aes_string(x="PC1", y="PC2", label="labels",
                                    angle=45, size=label_size, vjust=2))
  } else if (plot_labels == "repel") {
    pca_plot <- pca_plot +
      ggrepel::geom_text_repel(aes_string(label="labels"),
                               size=label_size, box.padding=ggplot2::unit(0.5, "lines"),
                               point.padding=ggplot2::unit(1.6, "lines"),
                               arrow=ggplot2::arrow(length=ggplot2::unit(0.01, "npc")))
  } else if (plot_labels == "dlsmart") {
    pca_plot <- pca_plot +
      directlabels::geom_dl(aes_string(label="labels"), method="smart.grid")
  } else {
    pca_plot <- pca_plot +
      directlabels::geom_dl(aes_string(label="labels"), method="first.qp")
  }

  ## Add a little check to only deal with the confidence-interval-able data.
  ci_keepers <- pca_data %>%
    dplyr::group_by(!!rlang::sym(ci_group)) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count > 3)
  if (nrow(ci_keepers) < 1) {
    cis <- NULL
  }
  if (!is.null(cis)) {
    ci_idx <- pca_data[[ci_group]] %in% ci_keepers[[ci_group]]
    ci_data <- pca_data[ci_idx, ]
    alpha <- 0
    for (ci in cis) {
      alpha <- alpha + 0.1
      pca_plot <- pca_plot +
        ggplot2::stat_ellipse(data=ci_data,
                              mapping=aes_string(group=ci_group, fill=ci_fill),
                              geom="polygon", type="t", level=ci, alpha=alpha)
    }
  }

  ## Set default font sizes and colors
  pca_plot <- pca_plot +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   legend.key.size=grid::unit(0.5, "cm"))

  return(pca_plot)
}

#' Make a ggplot PCA plot describing the samples' clustering.
#'
#' @param data  an expt set of samples.
#' @param design   a design matrix and.
#' @param plot_colors   a color scheme.
#' @param plot_title   a title for the plot.
#' @param plot_size   size for the glyphs on the plot.
#' @param plot_labels   add labels?  Also, what type?  FALSE, "default", or "fancy".
#' @param size_column use an experimental factor to size the glyphs of the plot
#' @param pc_method  how to extract the components? (svd
#' @param x_pc  Component to put on the x axis.
#' @param y_pc  Component to put on the y axis.
#' @param transpose  Perform pca on samples or genes?
#' @param num_pc  How many components to calculate, default to the number of rows in the metadata.
#' @param ...  Arguments passed through to the pca implementations and plotter.
#' @return a list containing the following (this is currently wrong)
#' \enumerate{
#'  \item  pca = the result of fast.svd()
#'  \item  plot = ggplot2 pca_plot describing the principle component analysis of the samples.
#'  \item  table = a table of the PCA plot data
#'  \item  res = a table of the PCA res data
#'  \item  variance = a table of the PCA plot variance
#' }
#' @seealso \pkg{directlabels}
#'  \code{\link[directlabels]{geom_dl}} \code{\link{plot_pcs}}
#' @examples
#' \dontrun{
#'  pca_plot <- plot_pca(expt=expt)
#'  pca_plot
#' }
#' @export
plot_pca <- function(data, design=NULL, plot_colors=NULL, plot_title=NULL,
                     plot_size=5, plot_alpha=NULL, plot_labels=NULL, size_column=NULL,
                     pc_method="fast_svd", x_pc=1, y_pc=2, pc_type="sample",
                     num_pc=NULL,
                     ...) {
  arglist <- list(...)
  plot_names <- arglist[["plot_names"]]
  ## Set default columns in the experimental design for condition and batch
  ## changing these may be used to query other experimental factors with pca.
  cond_column <- "condition"
  if (!is.null(arglist[["cond_column"]])) {
    cond_column <- arglist[["cond_column"]]
    message("Using ", cond_column, " as the condition column in the experimental design.")
  }
  batch_column <- "batch"
  if (!is.null(arglist[["batch_column"]])) {
    batch_column <- arglist[["batch_column"]]
    message("Using ", batch_column, " as the batch column in the experimental design.")
  }
  if (!is.null(arglist[["base_size"]])) {
    base_size <<- arglist[["base_size"]]
  }

  ## The following if() series is used to check the type of data provided and extract the available
  ## metadata from it.  Since I commonly use my ExpressionSet wrapper (expt), most of the material is
  ## specific to that.  However, the functions in this package should be smart enough to deal when
  ## that is not true.
  ## The primary things this particular function is seeking to acquire are: design, colors, counts.
  ## The only thing it absolutely requires to function is counts, it will make up the rest if it cannot
  ## find them.
  data_class <- class(data)[1]
  names <- NULL
  mtrx <- NULL
  if (data_class == "expt") {
    design <- pData(data)
    if (cond_column == "condition") {
      plot_colors <- data[["colors"]]
    } else {
      plot_colors <- NULL
    }
    plot_names <- data[["samplenames"]]
    mtrx <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    mtrx <- exprs(data)
    design <- pData(data)
  } else if (data_class == "list") {
    mtrx <- data[["count_table"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table element.")
    }
  } else if (data_class == "matrix" | data_class == "data.frame") {
    mtrx <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Small modification for reusing some of my very oldest experimental designs.
  if (is.null(design[["sampleid"]])) {
    design[["sampleid"]] <- rownames(design)
  }

  ## Check that the given design works with the data
  ## Prune the design if necessary
  ## Also take into account the fact that sometimes I change the case of hpgl<->HPGL
  given_samples <- tolower(colnames(mtrx))
  colnames(mtrx) <- given_samples
  ## I hate uppercase characters, I ADMIT IT.
  avail_samples <- tolower(rownames(design))
  rownames(design) <- avail_samples
  if (sum(given_samples %in% avail_samples) == length(given_samples)) {
    design <- design[given_samples, ]
  }

  ## If nothing has given this some colors for the plot, make them up now.
  if (is.null(plot_colors)) {
    plot_colors <- as.numeric(as.factor(design[[cond_column]]))
    plot_colors <- RColorBrewer::brewer.pal(12, "Dark2")[plot_colors]
  }

  ## Similarly, if there is no information which may be used as a design yet, make one up.
  if (is.null(design)) {
    message("No design was provided.  Making one with x conditions, 1 batch.")
    design <- cbind(plot_names, 1)
    design <- as.data.frame(design)
    design[["condition"]] <- as.numeric(design[["plot_labels"]])
    colnames(design) <- c("name", "batch", "condition")
    design <- design[, c("name", "condition", "batch")]
    plot_names <- design[["name"]]
  }

  ## Different folks like different labels.  I prefer hpglxxxx, but others have asked for
  ## condition_batch; this handles that as eloquently as I am able.
  label_list <- NULL
  if (is.null(arglist[["label_list"]]) & is.null(plot_names)) {
    label_list <- design[["sampleid"]]
  } else if (is.null(arglist[["label_list"]])) {
    label_list <- plot_names
  } else if (arglist[["label_list"]] == "concat") {
    label_list <- paste(design[[cond_column]], design[[batch_column]], sep="_")
  } else {
    label_list <- paste0(design[["sampleid"]], "_", design[[cond_column]])
  }

  ## I want to be able to seamlessly switch between PCs on samples vs. rows.
  if (pc_type != "sample") {
    mtrx <- t(mtrx)
  }
  mtrx <- as.matrix(mtrx)

  ## How many components should be calculated when that is possible to define?
  if (is.null(num_pc)) {
    num_pc <- nrow(design)
  }

  ## Pull out the batches and conditions used in this plot.
  ## Probably could have just used xxx[stuff, drop=TRUE]
  included_batches <- as.factor(as.character(design[[batch_column]]))
  included_conditions <- as.factor(as.character(design[[cond_column]]))

  pca_res <- NULL
  pca_result <- NULL
  switchret <- switch(
    pc_method,
    "fast_svd" = {
      svd_result <- corpcor::fast.svd(mtrx - rowMeans(mtrx))
      if (pc_type == "sample") {
        rownames(svd_result[["v"]]) <- rownames(design)
      }
      colnames(svd_result[["v"]]) <- paste0("PC", 1:ncol(svd_result[["v"]]))
      pc_table <- svd_result[["v"]]

      x_name <- paste0("PC", x_pc)
      y_name <- paste0("PC", y_pc)
      ## By a similar token, get the percentage of variance accounted for in each PC
      pca_variance <- round((svd_result[["d"]] ^ 2) / sum(svd_result[["d"]] ^ 2) * 100, 2)
      ## These will provide metrics on the x/y axes showing the amount of variance on those
      ## components of our plot.
      x_label <- sprintf("%s: %.2f%% variance", x_name, pca_variance[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, pca_variance[y_pc])
      ## Create a data frame with all the material of interest in the actual PCA plot

      ## Depending on how much batch/condition information is available, invoke pcRes() to get some idea of how
      ## much variance in a batch model is accounted for with each PC.
      if (pc_type == "sample") {
        if (length(levels(included_conditions)) == 1 & length(levels(included_batches)) == 1) {
          warning("There is only one condition and one batch, it is impossible to get meaningful pcRes information.")
        } else if (length(levels(included_conditions)) == 1) {
          warning("There is only one condition, but more than one batch.
Going to run pcRes with the batch information.")
          pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]], batch=design[, batch_column])
        } else if (length(levels(included_batches)) == 1) {
          message("There is just one batch in this data.")
          pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]],
                           condition=design[, cond_column])
        } else {
          pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]],
                           condition=design[, cond_column], batch=design[, batch_column])
        }
      } else {
        pca_res <- NULL
      }
    },
    "tsne" = {
      perplexity <- floor(ncol(mtrx) / 5)
      if (!is.null(arglist[["perplexity"]])) {
        perplexity <- arglist[["perplexity"]]
      }

      plotting_indexes <- 1:nrow(data)
      if (is.null(arglist[["chosen_features"]])) {
        variances <- matrixStats::rowVars(as.matrix(data))
        if (!is.null(arglist[["number_features"]])) {
          number_features <- min(number_features, nrow(data))
        } else {
          number_features <- nrow(data)
        }
        plotting_indexes <- order(variances, decreasing=TRUE)[1:number_features]
      }

      plotting_data <- data[plotting_indexes, ]
      ## This I do understand and think is cool
      ## Drop features with low variance
      min_variance <- 0.001
      if (!is.null(arglist[["min_variance"]])) {
        min_variance <- arglist[["min_variance"]]
      }
      keepers <- (matrixStats::rowVars(as.matrix(plotting_data)) >= min_variance)
      keepers[is.na(keepers)] <- FALSE ## Another nice idea
      plotting_data <- plotting_data[keepers, ]

      ## There is an interesting standardization idea in scater
      ## But I think I would prefer to have flexibility here
      ## exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))

      ## Spend a moment filling in the default values as specified in Rtsne.
      if (!is.null(arglist[["seed"]])) {
        set.seed(arglist[["seed"]])
      }
      theta <- 0.3
      if (!is.null(arglist[["theta"]])) {
        theta <- arglist[["theta"]]
      }
      iterations <- 1000
      if (!is.null(arglist[["iterations"]])) {
        iterations <- arglist[["iterations"]]
      }
      pca <- TRUE
      if (!is.null(arglist[["pca"]])) {
        pca <- arglist[["pca"]]
      }
      components <- 2
      if (!is.null(arglist[["components"]])) {
        components <- arglist[["components"]]
      }
      tsne_result <- Rtsne::Rtsne(plotting_data, check_duplicates=FALSE, dims=components,
                                  max_iter=iterations, pca=pca, theta=theta,
                                  perplexity=perplexity)
      pc_table <- as.data.frame(tsne_result[["Y"]])
      ## Changing these assignments because of my new attempts to use GSVA
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- paste0("PC", 1:ncol(pc_table))
      ##pc_table <- pc_table[, 1:components]
      pos_sing <- tsne_result[["costs"]]
      x_name <- paste0("Factor", x_pc)
      y_name <- paste0("Factor", y_pc)
    },
    "fast_ica" = {
      ## Fill in the defaults from the ica package.
      alg_type <- "parallel"  ## alg.typ is how they name this argument, which is just too
      ## weird looking for me to deal with.
      if (!is.null(arglist[["alg_type"]])) {
        alg_type <- arglist[["alg_type"]]
      }
      fun <- "logcosh"
      if (!is.null(arglist[["fun"]])) {
        fun <- arglist[["fun"]]
      }
      alpha <- 1.0
      if (!is.null(arglist[["alpha"]])) {
        alpha <- arglist[["alpha"]]
      }
      ica_method <- "C"
      if (!is.null(arglist[["ica_method"]])) {
        ica_method <- arglist[["ica_method"]]
      }
      row.norm <- FALSE
      if (!is.null(arglist[["row.norm"]])) {
        row.norm <- arglist[["row.norm"]]
      }
      maxit <- 200
      if (!is.null(arglist[["maxit"]])) {
        maxit <- arglist[["maxit"]]
      }
      tol <- 0.0001
      if (!is.null(arglist[["tol"]])) {
        tol <- arglist[["tol"]]
      }
      verbose <- FALSE
      if (!is.null(arglist[["verbose"]])) {
        verbose <- arglist[["verbose"]]
      }
      w.init <- NULL
      if (!is.null(arglist[["w.init"]])) {
        w.init <- arglist[["w.init"]]
      }
      ica_result <- fastICA::fastICA(mtrx, n.comp=num_pcs, alg.typ=alg_type, fun=fun,
                               alpha=alpha, method=ica_method, row.norm=row.norm,
                               maxit=maxit, tol=tol, verbose=verbose, w.init=w.init)
      pc_table <- ica_result[["S"]]
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- paste0("IC", 1:ncol(pc_table))
      x_label <- paste0("IC", x_pc)
      y_label <- paste0("IC", y_pc)
    },
    {
      if (pc_method == "bpca") {
        message("Just so you know, bpca is annoyingly slow.")
      }
      ## If none of the above were provided, then pick it up via the pcaMethods package.
      ## It handles the following: svd, ppca, bpca, svdi, nipals, nlpca.
      ## The following if statements pick up the default values.
      scale <- "none"
      if (!is.null(arglist[["scale"]])) {
        scale <- arglist[["scale"]]
      }
      center <- TRUE
      if (!is.null(arglist[["center"]])) {
        center <- arglist[["center"]]
      }
      completeObs <- TRUE
      if (!is.null(arglist[["completeObs"]])) {
        completeObs <- arglist[["completeObs"]]
      }
      subset <- NULL
      if (!is.null(arglist[["subset"]])) {
        subset <- arglist[["subset"]]
      }
      cv <- "none"
      if (!is.null(arglist[["cv"]])) {
        cv <- arglist[["cv"]]
      }
      eps <- 1e-12
      if (!is.null(arglist[["eps"]])) {
        eps <- arglist[["eps"]]
      }
      simple <- TRUE
      if (!is.null(arglist[["simple"]])) {
        simple <- arglist[["simple"]]
      }
      reverse <- FALSE
      if (!is.null(arglist[["reverse"]])) {
        reverse <- arglist[["reverse"]]
      }
      ready <- pcaMethods::prep(t(mtrx), scale=scale, center=center, eps=eps,
                                simple=simple, reverse=reverse)
      pcam_result <- pcaMethods::pca(ready, method=pc_method,
                                     nPcs=num_pc,
                                     scale=scale,
                                     center=center,
                                     completeObs=completeObs,
                                     subset=subset,
                                     cv=cv)

      pc_table <- as.data.frame(pcam_result@scores)
      pc_variance <- pcam_result@R2 * 100
      x_name <- paste0("PC", x_pc)
      y_name <- paste0("PC", y_pc)
      x_label <- sprintf("%s: %.2f%% variance", x_name, pc_variance[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, pc_variance[y_pc])
    }
)  ## End of the switch()
  if (pc_type == "sample") {
    comp_data <- data.frame(
      "sampleid" = as.character(design[["sampleid"]]),
      "condition" = design[[cond_column]],
      "batch" = design[[batch_column]],
      "batch_int" = as.integer(as.factor(design[[batch_column]])),
      "colors" = as.character(plot_colors),
      "labels" = label_list,
      stringsAsFactors=FALSE)
  } else {
    comp_data <- data.frame(
      "sampleid" = rownames(pc_table),
      "condition" = "a",
      "batch" = "a",
      "batch_int" = 1,
      "colors" = "black",
      "labels" = FALSE)
    plot_size <- 0.5
    plot_alpha <- 0.3
    plot_labels <- FALSE
  }
  comp_data[[x_name]] <- pc_table[, x_pc]
  comp_data[[y_name]] <- pc_table[, y_pc]

  if (!is.null(size_column)) {
    ## Adding a column with the same name as the size column from the experimental design
    ## and making sure it is a factor.
    comp_data[[size_column]] <- as.factor(design[, size_column])
    ## Just forcing the size to be numeric non-zero.
    comp_data[["size"]] <- as.factor(as.integer(design[[size_column]]) + 1)
  }

  if (isTRUE(plot_title)) {
    plot_title <- what_happened(expt=expt)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt=expt)
    plot_title <- paste0(plot_title, "; ", data_title)
  } else {
    ## Leave the title blank.
  }

  ## The plot_pcs() function gives a decent starting plot
  comp_plot <- plot_pcs(
    comp_data,
    first=x_name,
    second=y_name,
    design=design,
    plot_labels=plot_labels,
    x_label=x_label,
    y_label=y_label,
    plot_title=plot_title,
    plot_size=plot_size,
    size_column=size_column,
    plot_alpha=plot_alpha,
    ...)

  ## If plot_title is NULL, print nothing, if it is TRUE
  ## Then give some information about what happened to the data to make the plot.
  ## I tried foolishly to put this in plot_pcs(), but there is no way that receives
  ## my expt containing the normalization state of the data.
  if (isTRUE(plot_title)) {
    data_title <- what_happened(expt=expt)
    comp_plot <- comp_plot + ggplot2::ggtitle(data_title)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt=expt)
    plot_title <- paste0(plot_title, "; ", data_title)
    comp_plot <- comp_plot + ggplot2::ggtitle(plot_title)
  } else {
    ## Leave the title blank.
  }

  ## Finally, return a list of the interesting bits of what happened.
  pca_return <- list(
    "pca" = svd_result,
    "plot" = comp_plot,
    "table" = comp_data,
    "res" = pca_res,
    "variance" = pca_variance)
  return(pca_return)
}

## An alternate to plotting rank order of svd$u
## The plotted_u1s and such below
## y-axis is z(i), x-axis is i
## z(i) = cumulative sum of $u squared
## z = cumsum((svd$u ^ 2))

#' Plot the rank order svd$u elements to get a view of how much
#' the first genes contribute to the total variance by PC.
#'
#' @param plotted_us  a list of svd$u elements
#' @return a recordPlot() plot showing the first 3 PCs by rank-order svd$u.
#' @export
u_plot <- function(plotted_us) {
  plotted_us <- abs(plotted_us[, c(1, 2, 3)])
  plotted_u1s <- plotted_us[order(plotted_us[, 1], decreasing=TRUE), ]
  plotted_u2s <- plotted_us[order(plotted_us[, 2], decreasing=TRUE), ]
  plotted_u3s <- plotted_us[order(plotted_us[, 3], decreasing=TRUE), ]
  ## allS <- BiocGenerics::rank(allS, ties.method = "random")
  ## plotted_us$rank = rank(plotted_us[,1], ties.method="random")
  plotted_u1s <- cbind(plotted_u1s, rev(rank(plotted_u1s[, 1], ties.method="random")))
  plotted_u1s <- plotted_u1s[, c(1, 4)]
  colnames(plotted_u1s) <- c("PC1", "rank")
  plotted_u1s <- data.frame(plotted_u1s)
  plotted_u1s[["ID"]] <- as.character(rownames(plotted_u1s))
  plotted_u2s <- cbind(plotted_u2s, rev(rank(plotted_u2s[, 2], ties.method="random")))
  plotted_u2s <- plotted_u2s[, c(2, 4)]
  colnames(plotted_u2s) <- c("PC2", "rank")
  plotted_u2s <- data.frame(plotted_u2s)
  plotted_u2s[["ID"]] <- as.character(rownames(plotted_u2s))
  plotted_u3s <- cbind(plotted_u3s, rev(rank(plotted_u3s[, 3], ties.method="random")))
  plotted_u3s <- plotted_u3s[, c(3, 4)]
  colnames(plotted_u3s) <- c("PC3", "rank")
  plotted_u3s <- data.frame(plotted_u3s)
  plotted_u3s[["ID"]] <- as.character(rownames(plotted_u3s))
  plotted_us <- merge(plotted_u1s, plotted_u2s, by.x="rank", by.y="rank")
  plotted_us <- merge(plotted_us, plotted_u3s, by.x="rank", by.y="rank")
  colnames(plotted_us) <- c("rank", "PC1", "ID1", "PC2", "ID2", "PC3", "ID3")
  rm(plotted_u1s)
  rm(plotted_u2s)
  rm(plotted_u3s)
  ## top_threePC = head(plotted_us, n=20)
  plotted_us <- plotted_us[, c("PC1", "PC2", "PC3")]
  plotted_us[, "ID"] <- rownames(plotted_us)
  message("More shallow curves in these plots suggest more genes in this principle component.")
  plot(plotted_us)
  u_plot <- grDevices::recordPlot()
  return(u_plot)
}

## EOF
