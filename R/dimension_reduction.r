## dimension_reduction.r: Variants of PCA plotting.  The functions in this file
## seek to simplify performing the various dimension reduction methods and plot
## the results.

#' Collect the r^2 values from a linear model fitting between a singular
#' value decomposition and factor.
#'
#' @param datum Result from corpcor::fast.svd.
#' @param fact Experimental factor from the original data.
#' @param type Make this categorical or continuous with factor/continuous.
#' @return The r^2 values of the linear model as a percentage.
#' @seealso \pkg{corpcor}
#'  \code{\link[corpcor]{fast.svd}}
#' @export
factor_rsquared <- function(datum, fact, type = "factor") {
  if (type == "factor") {
    fact <- as.factor(fact)
  } else if (type == "numeric") {
    fact <- as.numeric(fact)
  } else {
    fact <- as.factor(as.numeric(fact))
  }
  ## FIXME! This is not the correct way to handle this
  if (length(levels(fact)) < 2) {
    svd_lm <- NULL
  } else {
    svd_lm <- stats::lm(datum ~ fact)
  }

  if (class(svd_lm)[1] == "try-error" || is.null(svd_lm)) {
    rsq_result <- NULL
  } else {
    lm_summaries <- summary(svd_lm)
    rsq_result <- c()
    for (i in 1:length(lm_summaries)) {
      rsq <- lm_summaries[[i]][["r.squared"]]
      rsq_result[i] <- rsq
    }
  }
  return(rsq_result)
}

#' Attempt to get residuals from tsne data
#'
#' I strongly suspect that this is not correct, but it is a start.
#'
#' @param svd_result The set of results from one of the many potential svd-ish methods.
#' @param design Experimental design from which to get experimental factors.
#' @param factors Set of experimental factors for which to calculate rsquared values.
#' @param res_slot Where is the res data in the svd result?
#' @param var_slot Where is the var data in the svd result?
get_res <- function(svd_result, design, factors = c("condition", "batch"),
                    res_slot = "v", var_slot = "d") {
  retlst <- list()
  vars <- svd_result[[var_slot]]
  rsquared_column <- round((vars ^ 2) / sum(vars ^ 2) * 100, 2)
  cumulative_sum_column <- cumsum(rsquared_column)
  datum <- svd_result[[res_slot]]

  res_df <- data.frame(
    "prop_var" = rsquared_column,
    "cum_prop_var" = cumulative_sum_column)

  for (j in 1:length(factors)) {
    factor <- factors[j]
    if (!is.null(design[[factor]])) {
      fact <- design[[factor]]
      column <- factor_rsquared(datum, fact)
      column_name <- glue::glue("{factor}_rsquared")
      res_df[[column_name]] <- column
    }
  }
  return(res_df)
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
#' @param expt Data to analyze (usually exprs(somedataset)).
#' @param expt_design Dataframe describing the experimental design, containing
#'   columns with useful information like the conditions, batches, number of
#'   cells, whatever...
#' @param expt_factors Character list of experimental conditions to query for
#'   R^2 against the fast.svd of the data.
#' @param num_components Number of principle components to compare the design
#'   factors against. If left null, it will query the same number of components
#'   as factors asked for.
#' @param plot_pcas Plot the set of PCA plots for every pair of PCs queried.
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
pca_information <- function(expt, expt_design = NULL, expt_factors = c("condition", "batch"),
                            num_components = NULL, plot_pcas = FALSE, ...) {
  ## Start out with some sanity tests
  colors_chosen <- NULL
  exprs_data <- NULL
  data_class <- class(expt)[1]
  if (data_class == "expt") {
    expt_design <- expt[["design"]]
    colors_chosen <- expt[["colors"]]
    exprs_data <- exprs(expt)
  } else if (data_class == "ExpressionSet") {
    exprs_data <- exprs(expt)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    exprs_data <- as.matrix(expt)
  } else {
    stop("This understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Make sure colors get chosen.
  if (is.null(colors_chosen)) {
    colors_chosen <- as.numeric(as.factor(expt_design[["condition"]]))
    num_chosen <- max(3, length(levels(as.factor(colors_chosen))))
    colors_chosen <- RColorBrewer::brewer.pal(num_chosen, "Dark2")[colors_chosen]
    names(colors_chosen) <- rownames(expt_design)
  }

  initial_pca <- plot_pca(expt)
  v <- initial_pca[["result"]][["v"]]
  u <- initial_pca[["result"]][["u"]]
  d <- initial_pca[["result"]][["d"]]
  u_plot <- u_plot(u)
  component_rsquared_table <- initial_pca[["residual_df"]]
  variance <- initial_pca[["prop_var"]]
  pca_data <- initial_pca[["table"]]

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
    name <- glue::glue("PC{pc}")
    ## v is a matrix, don't forget that.
    pca_data[[name]] <- v[, pc] ## note you _must_ not shortcut this with [[pc]]
  }
  ## Now that we have filled in a pca data frame, we may plot PCx vs PCy for all
  ## x,y.

  pca_plots <- list()
  if (isTRUE(plot_pcas)) {
    nminus_one <- num_components - 1
    for (pc in 1:nminus_one) {
      next_pc <- pc + 1
      name <- glue::glue("PC{pc}")
      for (second_pc in next_pc:num_components) {
        if (pc < second_pc & second_pc <= num_components) {
          second_name <- glue::glue("PC{second_pc}")
          list_name <- glue::glue("{name}_{second_name}")
          ## Sometimes these plots fail because too many grid operations are happening.
          tmp_plot <- try(plot_pcs(pca_data, design = expt_design, variances = variance,
                                   first = name, second = second_name))
          pca_plots[[list_name]] <- tmp_plot
        }
      }
    }
  }

  ## Now start filling in data which may be used for correlations/fstats/etc.
  factor_df <- data.frame(
    "sampleid" = tolower(rownames(expt_design)))
  rownames(factor_df) <- tolower(rownames(expt_design))
  for (fact in expt_factors) {
    if (!is.null(expt_design[[fact]])) {
      factor_df[[fact]] <- as.numeric(as.factor(as.character(expt_design[, fact])))
    } else {
      message("The column ", fact, " seems to be missing from the design.")
      message("The available columns are: ", toString(colnames(expt_design)), ".")
    }
  }
  factor_df <- factor_df[, -1, drop = FALSE]

  ## Perform the correlations/fstats/anova here
  cor_df <- data.frame()
  anova_rss <- data.frame()
  anova_sums <- data.frame()
  anova_f <- data.frame()
  anova_p <- data.frame()
  anova_rss <- data.frame()
  anova_fstats <- data.frame()
  for (f in 1:length(expt_factors)) {
    fact <- expt_factors[f]
    for (pc in 1:num_components) {
      pc_name <- glue::glue("pc_{pc}")
      tmp_df <- merge(factor_df, pca_data, by = "row.names")
      test_column <- glue::glue("{fact}.x")
      if (!is.null(tmp_df[[test_column]])) {
        col_idx <- which(colnames(tmp_df) == test_column)
        colnames(tmp_df)[col_idx] <- fact
      }
      rownames(tmp_df) <- tmp_df[["Row.names"]]
      tmp_df <- tmp_df[, -1, drop = FALSE]
      form <- as.formula(glue::glue("{pc_name} ~ 1 + {fact}"))
      lmwithfactor_test <- try(stats::lm(formula = form,
                                         data = tmp_df))
      form <- as.formula(glue::glue("{pc_name} ~ 1"))
      lmwithoutfactor_test <- try(stats::lm(formula = form, data = tmp_df))
      ## This fstat provides a metric of how much variance is removed by
      ## including this specific factor in the model vs not.  Therefore higher
      ## numbers tell us that adding that factor removed more variance and are more important.
      fstat <- sum(residuals(lmwithfactor_test) ^ 2) / sum(residuals(lmwithoutfactor_test) ^ 2)
      ## 1.  Perform lm(pc ~ 1 + factor) which is fit1
      ## 2.  Perform lm(pc ~ 1) which is fit2
      ## 3.  The Fstat is then defined as (sum(residuals(fit1)^2) / sum(residuals(fit2)^2))
      ## 4.  The resulting p-value is 1 - pf(Fstat, (n-(#levels in the factor)), (n-1))
      ##     n is the number of samples in the fit
      ## 5.  Look at anova.test() to see if this provides similar/identical information
      another_fstat <- try(stats::anova(lmwithfactor_test, lmwithoutfactor_test), silent = TRUE)
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
        cor_test <- cor.test(tmp_df[[fact]], tmp_df[[pc_name]], na.rm = TRUE)
      },
      error = function(cond) {
        message("The correlation failed for ", fact, " and ", pc_name, ".")
        cor_test <- 0
      },
      warning = function(cond) {
        message("The standard deviation was 0 for ", fact, " and ", pc_name, ".")
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
  colnames(cor_df) <- glue::glue("PC{1:ncol(cor_df)}")
  colnames(anova_sums) <- glue::glue("PC{1:ncol(anova_sums)}")
  colnames(anova_f) <- glue::glue("PC{1:ncol(anova_f)}")
  colnames(anova_p) <- glue::glue("PC{1:ncol(anova_p)}")
  colnames(anova_rss) <- glue::glue("PC{1:ncol(anova_rss)}")
  colnames(anova_fstats) <- glue::glue("PC{1:ncol(anova_fstats)}")
  ## Finally, plot them.
  silly_colors <- grDevices::colorRampPalette(c("purple", "black", "yellow"))(100)
  cor_df <- cor_df[complete.cases(cor_df), ]
  pc_factor_corheat <- heatmap.3(as.matrix(cor_df), scale = "none", trace = "none",
                                 linewidth = 0.5, keysize = 2, margins = c(8, 8),
                                 col = silly_colors, dendrogram = "none", Rowv = FALSE,
                                 Colv = FALSE, main = "cor(factor, PC)")
  pc_factor_corheat <- grDevices::recordPlot()
  anova_f_colors <- grDevices::colorRampPalette(c("blue", "black", "red"))(100)
  anova_f_heat <- heatmap.3(as.matrix(anova_f), scale = "none", trace = "none",
                            linewidth = 0.5, keysize = 2, margins = c(8, 8),
                            col = anova_f_colors, dendrogram = "none", Rowv = FALSE,
                            Colv = FALSE, main = "anova fstats for (factor, PC)")
  anova_f_heat <- grDevices::recordPlot()
  anova_fstat_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  anova_fstat_heat <- heatmap.3(as.matrix(anova_fstats), scale = "none", trace = "none",
                                linewidth = 0.5, keysize = 2, margins = c(8, 8),
                                col = anova_fstat_colors, dendrogram = "none", Rowv = FALSE,
                                Colv = FALSE, main = "anova fstats for (factor, PC)")
  anova_fstat_heat <- grDevices::recordPlot()
  ## I had this as log(anova_p + 1) !! I am a doofus; too many times I have been log2-ing counts.
  ## The messed up part is that I did not notice this for multiple years.
  neglog_p <- -1 * log(as.matrix(anova_p) + 0.00001)
  anova_neglogp_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  anova_neglogp_heat <- heatmap.3(as.matrix(neglog_p), scale = "none", trace = "none",
                                  linewidth = 0.5, keysize = 2, margins = c(8, 8),
                                  col = anova_f_colors, dendrogram = "none", Rowv = FALSE,
                                  Colv = FALSE, main = "-log(anova_p values)")
  anova_neglogp_heat <- grDevices::recordPlot()
  ## Another option: -log10 p-value of the ftest for this heatmap.
  ## covariate vs PC score
  ## Analagously: boxplot(PCn ~ batch)

  ## Finally, return the set of materials collected.
  pca_list <- list(
    "pc1_trend" = u_plot,
    "svd_d" = d,
    "svd_u" = u,
    "svd_v" = v,
    "rsquared_table" = component_rsquared_table,
    "variance" = variance,
    "pca_data" = pca_data,
    "anova_fstats" = anova_fstats,
    "anova_sums" = anova_sums,
    "anova_f" = anova_f,
    "anova_p" = anova_p,
    "pca_cor" = cor_df,
    "cor_heatmap" = pc_factor_corheat,
    "anova_f_heatmap" = anova_f_heat,
    "anova_fstat_heatmap" = anova_fstat_heat,
    "anova_neglogp_heatmap" = anova_neglogp_heat,
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
#'  information <- pca_highscores(df = df, conditions = cond, batches = bat)
#'  information$pca_bitplot  ## oo pretty
#' }
#' @export
pca_highscores <- function(expt, n = 20, cor = TRUE, vs = "means", logged = TRUE) {
  if (isTRUE(logged)) {
    if (expt[["state"]][["transform"]] == "raw") {
      expt <- sm(normalize_expt(expt, transform = "log2", filter = TRUE))
    }
  }

  data <- as.data.frame(exprs(expt))
  na_idx <- is.na(data)
  data[na_idx] <- 0
  if (!is.null(vs)) {
    if (vs == "means") {
      data <- as.matrix(data) - rowMeans(as.matrix(data))
    } else if (vs == "medians") {
      data <- as.matrix(data) - rowMedians(as.matrix(data))
    }
  }
  another_pca <- try(princomp(x = data, cor = cor))
  plot(another_pca)
  pca_hist <- grDevices::recordPlot()
  biplot(another_pca)
  pca_biplot <- grDevices::recordPlot()
  highest <- NULL
  lowest <- NULL
  for (pc in 1:length(colnames(another_pca[["scores"]]))) {
    tmphigh <- another_pca[["scores"]]
    tmplow <- another_pca[["scores"]]
    tmphigh <- tmphigh[order(tmphigh[, pc], decreasing = TRUE), ]
    tmphigh <- head(tmphigh, n = n)
    tmplow <- tmplow[order(tmplow[, pc], decreasing = FALSE), ]
    tmplow <- head(tmplow, n = n)
    high_column <- glue::glue("{signif(tmphigh[, pc], 4)}:{rownames(tmphigh)}")
    low_column <- glue::glue("{signif(tmplow[, pc], 4)}:{rownames(tmplow)}")
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
  retlist[["score_heat"]] <- plot_disheat(expt_data = another_pca[["scores"]])
  return(retlist)
}

#' Something silly for Najib.
#'
#' This will make him very happy, but I remain skeptical.
#'
#' @param pc_result  The result from plot_pca()
#' @param components  List of three axes by component.
#' @param file  File to write the created plotly object.
#' @export
plot_3d_pca <- function(pc_result, components = c(1,2,3), file = "3dpca.html") {
  image_dir <- dirname(as.character(file))
  if (!file.exists(image_dir)) {
    dir.create(image_dir, recursive = TRUE)
  }

  x_axis <- glue::glue("pc_{components[1]}")
  y_axis <- glue::glue("pc_{components[2]}")
  z_axis <- glue::glue("pc_{components[3]}")
  table <- pc_result[["table"]]
  color_levels <- levels(as.factor(table[["colors"]]))
  silly_plot <- plotly::plot_ly(table, x = as.formula(glue::glue("~{x_axis}")),
                                y = as.formula(glue::glue("~{y_axis}")), z = as.formula(glue::glue("~{z_axis}")),
                                color = as.formula("~condition"), colors = color_levels,
                                text=~paste0('condition: ', condition, ' batch: ', batch)) %>%
    plotly::add_markers()
  widget <- htmlwidgets::saveWidget(
                           plotly::as_widget(silly_plot), file = file, selfcontained = TRUE)
  retlist <- list(
    "plot" = silly_plot,
    "file" = file)
  return(retlist)
}

#' Make a PCA plot describing the samples' clustering.
#'
#' @param data an expt set of samples.
#' @param design a design matrix and.
#' @param plot_colors a color scheme.
#' @param plot_title a title for the plot.
#' @param plot_size size for the glyphs on the plot.
#' @param plot_alpha Add an alpha channel to the dots?
#' @param plot_labels add labels?  Also, what type?  FALSE, "default", or "fancy".
#' @param size_column use an experimental factor to size the glyphs of the plot
#' @param pc_method how to extract the components? (svd
#' @param x_pc Component to put on the x axis.
#' @param y_pc Component to put on the y axis.
#' @param num_pc How many components to calculate, default to the number of
#'   rows in the metadata.
#' @param expt_names Column or character list of preferred sample names.
#' @param label_chars Maximum number of characters before abbreviating sample names.
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
#'  pca_plot <- plot_pca(expt = expt)
#'  pca_plot
#' }
#' @export
plot_pca <- function(data, design = NULL, plot_colors = NULL, plot_title = NULL,
                     plot_size = 5, plot_alpha = NULL, plot_labels = NULL, size_column = NULL,
                     pc_method = "fast_svd", x_pc = 1, y_pc = 2,
                     num_pc = NULL, expt_names = NULL, label_chars = 10,
                     ...) {
  arglist <- list(...)
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

  if (!is.null(arglist[["transform"]]) || !is.null(arglist[["convert"]]) ||
      !is.null(arglist[["filter"]]) || !is.null(arglist[["norm"]]) ||
      !is.null(arglist[["batch"]])) {
    data <- normalize_expt(data, transform = arglist[["transform"]],
                           convert = arglist[["convert"]],
                           filter = arglist[["filter"]],
                           batch = arglist[["batch"]],
                           norm = arglist[["norm"]])
  }

  ## The following if() series is used to check the type of data provided and
  ## extract the available metadata from it.  Since I commonly use my
  ## ExpressionSet wrapper (expt), most of the material is specific to that.
  ## However, the functions in this package should be smart enough to deal when
  ## that is not true. The primary things this particular function is seeking to
  ## acquire are: design, colors, counts. The only thing it absolutely requires
  ## to function is counts, it will make up the rest if it cannot find them.
  data_class <- class(data)[1]
  mtrx <- NULL
  if (data_class == "expt") {
    design <- pData(data)
    if (cond_column == "condition") {
      plot_colors <- data[["colors"]]
    } else {
      plot_colors <- NULL
    }
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
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    mtrx <- as.data.frame(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Get rid of NAs, this is relevant because of recent changes in how I handle
  ## proteomic data with missing values.
  na_idx <- is.na(mtrx)
  mtrx[na_idx] <- 0

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
    design <- cbind(given_samples, 1)
    design <- as.data.frame(design)
    design[["condition"]] <- as.numeric(design[["plot_labels"]])
    colnames(design) <- c("name", "batch", "condition")
    design <- design[, c("name", "condition", "batch")]
    plot_names <- design[["name"]]
  }

  ## Different folks like different labels.  I prefer hpglxxxx, but others have asked for
  ## condition_batch; this handles that as eloquently as I am able.
  label_list <- NULL
  if (is.null(arglist[["label_list"]]) & is.null(expt_names)) {
    label_list <- design[["sampleid"]]
  } else if (class(expt_names) == "character" & length(expt_names) == 1) {
    label_list <- design[[expt_names]]
  } else if (is.null(arglist[["label_list"]])) {
    label_list <- given_samples
  } else if (arglist[["label_list"]] == "concat") {
    label_list <- paste(design[[cond_column]], design[[batch_column]], sep = "_")
  } else {
    label_list <- glue::glue("{design[['sampleid']]}_{design[[cond_column]]}")
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    label_list <- abbreviate(label_list, minlength = label_chars)
  }

  ## This line should be redundant
  mtrx <- as.matrix(mtrx)
  ## How many components should be calculated when that is possible to define?
  if (is.null(num_pc)) {
    num_pc <- nrow(design) - 1
  }

  ## Pull out the batches and conditions used in this plot.
  ## Probably could have just used xxx[stuff, drop = TRUE]
  included_batches <- as.factor(as.character(design[[batch_column]]))
  included_conditions <- as.factor(as.character(design[[cond_column]]))

  ## Expected things to retrieve from any dimension reduction method.
  x_label <- NULL
  y_label <- NULL
  residual_df <- NULL
  prop_lst <- NULL
  svd_result <- NULL
  pc_summary <- NULL

  switchret <- switch(
    pc_method,
    "fast_svd" = {
      svd_result <- corpcor::fast.svd(mtrx - rowMeans(mtrx))
      rownames(svd_result[["v"]]) <- rownames(design)
      colnames(svd_result[["v"]]) <- glue::glue("PC{1:ncol(svd_result[['v']])}")
      pc_table <- svd_result[["v"]]
      x_name <- glue::glue("PC{x_pc}")
      y_name <- glue::glue("PC{y_pc}")
      ## Depending on how much batch/condition information is available, invoke
      ## pcRes() to get some idea of how much variance in a batch model is
      ## accounted for with each PC.
      residual_df <- get_res(svd_result, design)
      prop_lst <- residual_df[["prop_var"]]
      ## get the percentage of variance accounted for in each PC
      x_label <- sprintf("%s: %.2f%% variance", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, prop_lst[y_pc])
    },
    "tsne" = {
      plotting_indexes <- 1:nrow(mtrx)
      if (is.null(arglist[["chosen_features"]])) {
        variances <- matrixStats::rowVars(as.matrix(mtrx))
        if (!is.null(arglist[["number_features"]])) {
          number_features <- min(number_features, nrow(mtrx))
        } else {
          number_features <- nrow(mtrx)
        }
        plotting_indexes <- order(variances, decreasing = TRUE)[1:number_features]
      }
      plotting_data <- mtrx[plotting_indexes, ]

      ## This I do understand and think is cool
      ## Drop features with low variance
      min_variance <- 0.001
      if (!is.null(arglist[["min_variance"]])) {
        min_variance <- arglist[["min_variance"]]
      }
      keepers <- (matrixStats::rowVars(as.matrix(plotting_data)) >= min_variance)
      keepers[is.na(keepers)] <- FALSE ## Another nice idea
      plotting_data <- plotting_data[keepers, ]

      perplexity <- 1
      plotting_data <- t(plotting_data)
      perplexity <- floor(nrow(plotting_data) / 5)
      if (!is.null(arglist[["perplexity"]])) {
        perplexity <- arglist[["perplexity"]]
      }
      if (perplexity <= 0) {
        warning("TSNE: Attempting to auto-detect perplexity failed, setting it to 1.")
        perplexity <- 1
      }

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
      svd_result <- Rtsne::Rtsne(plotting_data, check_duplicates = FALSE, dims = components,
                                 max_iter = iterations, pca = pca, theta = theta,
                                 perplexity = perplexity)
      pc_table <- as.data.frame(svd_result[["Y"]])
      ## Changing these assignments because of my new attempts to use GSVA
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- glue::glue("PC{1:ncol(pc_table)}")
      ##pc_table <- pc_table[, 1:components]
      pos_sing <- svd_result[["costs"]]
      x_name <- glue::glue("Factor{x_pc}")
      y_name <- glue::glue("Factor{y_pc}")

      ## Pull out the batches and conditions used in this plot.
      ## Probably could have just used xxx[stuff, drop = TRUE]
      included_batches <- as.factor(as.character(design[[batch_column]]))
      included_conditions <- as.factor(as.character(design[[cond_column]]))

      residual_df <- get_res(svd_result, design, res_slot = "Y", var_slot = "itercosts")
      prop_lst <- residual_df[["prop_var"]]
      if (x_pc > components | y_pc > components) {
        stop("The components plotted must be smaller than the number of components calculated.")
      }

      x_label <- sprintf("%s: %.2f tsne 'variance'", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f tsne 'variance'", y_name, prop_lst[y_pc])
    },
    "umap" = {
      message("Using uwot's implementation of umap.")
      plotting_data <- t(mtrx)
      pc_table <- as.data.frame(uwot::umap(X = plotting_data, n_components = num_pc, ...))
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- glue::glue("PC{1:ncol(pc_table)}")
      x_name <- glue::glue("Factor{x_pc}")
      y_name <- glue::glue("Factor{y_pc}")
      x_label <- x_name
      y_label <- y_name
      included_batch <- as.factor(as.character(design[[batch_column]]))
      included_conditions <- as.factor(as.character(design[[cond_column]]))
    },
    "uwot" = {
      plotting_data <- t(mtrx)
      pc_table <- as.data.frame(uwot::umap(X = plotting_data, n_components = num_pc, ...))
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- glue::glue("PC{1:ncol(pc_table)}")
      x_name <- glue::glue("Factor{x_pc}")
      y_name <- glue::glue("Factor{y_pc}")
      x_label <- x_name
      y_label <- y_name
      included_batch <- as.factor(as.character(design[[batch_column]]))
      included_conditions <- as.factor(as.character(design[[cond_column]]))
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
      svd_result <- fastICA::fastICA(t(mtrx), n.comp = num_pc, alg.typ = alg_type, fun = fun,
                                     alpha = alpha, method = ica_method, row.norm = row.norm,
                                     maxit = maxit, tol = tol, verbose = verbose, w.init = w.init)
      pc_table <- svd_result[["S"]]
      residual_df <- get_res(svd_result, design, res_slot = "S", var_slot = "W")
      prop_lst <- residual_df[[1]]
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- glue::glue("IC{1:ncol(pc_table)}")
      x_label <- glue::glue("IC{x_pc}")
      y_label <- glue::glue("IC{y_pc}")
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
      ready <- pcaMethods::prep(t(mtrx), scale = scale, center = center, eps = eps,
                                simple = simple, reverse = reverse)
      pca_result <- try(pcaMethods::pca(ready, method = pc_method, nPcs = num_pc, scale = scale,
                                        center = center, completeObs = completeObs,
                                        subset = subset, cv = cv))
      if (class(pca_result)[1] == "try-error") {
        message("pca invocation failed.  A likely reason if that too many PCs were requested.")
        message("Trying again with 1/2 the PCs.")
        num_pc <- floor(num_pc / 2)
        pca_result <- pcaMethods::pca(ready, method = pc_method, nPcs = num_pc, scale = scale,
                                      center = center, completeObs = completeObs,
                                      subset = subset, cv = cv)
      }
      ## This is mostly guessing on my part.
      svd_result <- list(
        "v" = pca_result@scores,
        "d" = pca_result@sDev)
      residual_df <- get_res(svd_result, design)
      pc_table <- as.data.frame(pca_result@scores)
      prop_lst <- pca_result@R2 * 100
      x_name <- glue::glue("PC{x_pc}")
      y_name <- glue::glue("PC{y_pc}")
      x_label <- sprintf("%s: %.2f%% variance", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, prop_lst[y_pc])
    })  ## End of the switch()

  comp_data <- data.frame(
    "sampleid" = as.character(design[["sampleid"]]),
    "condition" = design[[cond_column]],
    "batch" = design[[batch_column]],
    "batch_int" = as.integer(as.factor(design[[batch_column]])),
    "colors" = as.character(plot_colors),
    "labels" = label_list,
    stringsAsFactors = FALSE)

  rownames(comp_data) <- rownames(pc_table)
  comp_data[[x_name]] <- pc_table[, x_pc]
  comp_data[[y_name]] <- pc_table[, y_pc]
  tmp <- as.data.frame(pc_table)
  for (pc in 1:num_pc) {
    oldname <- glue::glue("PC{pc}")
    pc_name <- glue::glue("pc_{pc}")
    comp_data[[pc_name]] <- tmp[[oldname]]
  }

  if (!is.null(size_column)) {
    ## Adding a column with the same name as the size column from the experimental design
    ## and making sure it is a factor.
    if (is.null(arglist[["size_order"]])) {
      comp_data[[size_column]] <- factor(design[[size_column]])
    } else {
      comp_data[[size_column]] <- factor(design[[size_column]], levels = arglist[["size_order"]])
    }
    ## Just forcing the size to be numeric non-zero.
    comp_data[["size"]] <- as.factor(as.integer(comp_data[[size_column]]) + 1)
  }

  if (isTRUE(plot_title)) {
    plot_title <- what_happened(expt = data)
  }

  ## Perform a check of the PC table.
  if (sum(is.na(comp_data)) > 0) {
    message("Potentially check over the experimental design, there appear to be missing values.")
    warning("There are NA values in the component data.  This can lead to weird plotting errors.")
  }

  if (nrow(comp_data) > 100 & is.null(plot_labels)) {
    message("plot labels was not set and there are more than 100 samples, disabling it.")
    plot_labels <- FALSE
  }

  ## The plot_pcs() function gives a decent starting plot
  comp_plot <- plot_pcs(
    comp_data, first = x_name, second = y_name, design = design,
    plot_labels = plot_labels, x_label = x_label, y_label = y_label,
    plot_title = plot_title, plot_size = plot_size, size_column = size_column,
    plot_alpha = plot_alpha,
    ...)

  ## If plot_title is NULL, print nothing, if it is TRUE
  ## Then give some information about what happened to the data to make the plot.
  ## I tried foolishly to put this in plot_pcs(), but there is no way that receives
  ## my expt containing the normalization state of the data.
  if (isTRUE(plot_title)) {
    data_title <- what_happened(expt = data)
    comp_plot <- comp_plot + ggplot2::ggtitle(data_title)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt = data)
    plot_title <- glue::glue("{plot_title}; {data_title}")
    comp_plot <- comp_plot + ggplot2::ggtitle(plot_title)
  } else {
    ## Leave the title blank.
  }

  ## Finally, return a list of the interesting bits of what happened.
  pca_return <- list(
    "residual_df" = residual_df,
    "prop_var" = prop_lst,
    "plot" = comp_plot,
    "table" = comp_data,
    "result" = svd_result)
  return(pca_return)
}

#' Make a PC plot describing the gene' clustering.
#'
#' @param data an expt set of samples.
#' @param design a design matrix and.
#' @param plot_colors a color scheme.
#' @param plot_title a title for the plot.
#' @param plot_size size for the glyphs on the plot.
#' @param plot_alpha Add an alpha channel to the dots?
#' @param plot_labels add labels?  Also, what type?  FALSE, "default", or "fancy".
#' @param size_column use an experimental factor to size the glyphs of the plot
#' @param pc_method how to extract the components? (svd
#' @param x_pc Component to put on the x axis.
#' @param y_pc Component to put on the y axis.
#' @param label_column Which metadata column to use for labels.
#' @param num_pc How many components to calculate, default to the number of
#'   rows in the metadata.
#' @param expt_names Column or character list of preferred sample names.
#' @param label_chars Maximum number of characters before abbreviating sample names.
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
#'  pca_plot <- plot_pca(expt = expt)
#'  pca_plot
#' }
#' @export
plot_pca_genes <- function(data, design = NULL, plot_colors = NULL, plot_title = NULL,
                           plot_size = 2, plot_alpha = 0.4, plot_labels = FALSE, size_column = NULL,
                           pc_method = "fast_svd", x_pc = 1, y_pc = 2, label_column = "description",
                           num_pc = 2, expt_names = NULL, label_chars = 10,
                           ...) {
  arglist <- list(...)
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

  if (!is.null(arglist[["transform"]]) || !is.null(arglist[["convert"]]) ||
      !is.null(arglist[["filter"]]) || !is.null(arglist[["norm"]]) ||
      !is.null(arglist[["batch"]])) {
    data <- normalize_expt(data, transform = arglist[["transform"]],
                           convert = arglist[["convert"]],
                           filter = arglist[["filter"]],
                           batch = arglist[["batch"]],
                           norm = arglist[["norm"]])
  }

  ## The following if() series is used to check the type of data provided and
  ## extract the available metadata from it.  Since I commonly use my
  ## ExpressionSet wrapper (expt), most of the material is specific to that.
  ## However, the functions in this package should be smart enough to deal when
  ## that is not true. The primary things this particular function is seeking to
  ## acquire are: design, colors, counts. The only thing it absolutely requires
  ## to function is counts, it will make up the rest if it cannot find them.
  data_class <- class(data)[1]
  mtrx <- NULL
  if (data_class == "expt") {
    design <- pData(data)
    if (cond_column == "condition") {
      plot_colors <- data[["colors"]]
    } else {
      plot_colors <- NULL
    }
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
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    mtrx <- as.data.frame(data)
  } else {
    stop("This understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
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
    design <- cbind(given_samples, 1)
    design <- as.data.frame(design)
    design[["condition"]] <- as.numeric(design[["plot_labels"]])
    colnames(design) <- c("name", "batch", "condition")
    design <- design[, c("name", "condition", "batch")]
    plot_names <- design[["name"]]
  }

  ## Different folks like different labels.  I prefer hpglxxxx, but others have asked for
  ## condition_batch; this handles that as eloquently as I am able.
  label_list <- NULL
  if (is.null(arglist[["label_list"]]) & is.null(expt_names)) {
    label_list <- design[["sampleid"]]
  } else if (class(expt_names) == "character" & length(expt_names) == 1) {
    label_list <- design[[expt_names]]
  } else if (is.null(arglist[["label_list"]])) {
    label_list <- given_samples
  } else if (arglist[["label_list"]] == "concat") {
    label_list <- paste(design[[cond_column]], design[[batch_column]], sep = "_")
  } else {
    label_list <- glue::glue("{design[['sampleid']]}_{design[[cond_column]]}")
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    label_list <- abbreviate(label_list, minlength = label_chars)
  }

  mtrx <- as.matrix(t(mtrx))

  ## How many components should be calculated when that is possible to define?
  if (is.null(num_pc)) {
    num_pc <- nrow(design) - 1
  }

  ## Pull out the batches and conditions used in this plot.
  ## Probably could have just used xxx[stuff, drop = TRUE]
  included_batches <- as.factor(as.character(design[[batch_column]]))
  included_conditions <- as.factor(as.character(design[[cond_column]]))

  ## Expected things to retrieve from any dimension reduction method.
  x_label <- NULL
  y_label <- NULL
  residual_df <- NULL
  prop_lst <- NULL
  svd_result <- NULL
  pc_summary <- NULL

  switchret <- switch(
    pc_method,
    "fast_svd" = {
      svd_result <- corpcor::fast.svd(mtrx - rowMeans(mtrx))
      fdata_data <- fData(data)
      fdata_rows <- rownames(fdata_data)
      pc_rows <- rownames(pc_table)
      kept_rows <- fdata_rows %in% pc_rows
      kept_fdata <- fdata_data[kept_rows, ]
      rownames(svd_result[["v"]]) <- rownames(kept_fdata)
      colnames(svd_result[["v"]]) <- glue::glue("PC{1:ncol(svd_result[['v']])}")
      pc_table <- svd_result[["v"]]
      x_name <- glue::glue("PC{x_pc}")
      y_name <- glue::glue("PC{y_pc}")
      ## Depending on how much batch/condition information is available, invoke
      ## pcRes() to get some idea of how much variance in a batch model is
      ## accounted for with each PC.
      residual_df <- get_res(svd_result, kept_fdata)
      prop_lst <- residual_df[["prop_var"]]
      ## get the percentage of variance accounted for in each PC
      x_label <- sprintf("%s: %.2f%% variance", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, prop_lst[y_pc])
    },
    "tsne" = {
      plotting_indexes <- 1:nrow(mtrx)
      if (is.null(arglist[["chosen_features"]])) {
        variances <- matrixStats::rowVars(as.matrix(mtrx))
        if (!is.null(arglist[["number_features"]])) {
          number_features <- min(number_features, nrow(mtrx))
        } else {
          number_features <- nrow(mtrx)
        }
        plotting_indexes <- order(variances, decreasing = TRUE)[1:number_features]
      }
      plotting_data <- mtrx[plotting_indexes, ]

      ## This I do understand and think is cool
      ## Drop features with low variance
      min_variance <- 0.001
      if (!is.null(arglist[["min_variance"]])) {
        min_variance <- arglist[["min_variance"]]
      }
      keepers <- (matrixStats::rowVars(as.matrix(plotting_data)) >= min_variance)
      keepers[is.na(keepers)] <- FALSE ## Another nice idea
      plotting_data <- plotting_data[keepers, ]

      perplexity <- 1
      plotting_data <- t(plotting_data)
      perplexity <- floor(nrow(plotting_data) / 5)
      if (!is.null(arglist[["perplexity"]])) {
        perplexity <- arglist[["perplexity"]]
      }
      if (perplexity <= 0) {
        warning("TSNE: Attempting to auto-detect perplexity failed, setting it to 1.")
        perplexity <- 1
      }

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
      svd_result <- Rtsne::Rtsne(plotting_data, check_duplicates = FALSE, dims = components,
                                 max_iter = iterations, pca = pca, theta = theta,
                                 perplexity = perplexity)
      pc_table <- as.data.frame(svd_result[["Y"]])
      ## Changing these assignments because of my new attempts to use GSVA
      rownames(pc_table) <- rownames(plotting_data)
      colnames(pc_table) <- glue::glue("PC{1:ncol(pc_table)}")
      ##pc_table <- pc_table[, 1:components]
      pos_sing <- svd_result[["costs"]]
      x_name <- glue::glue("Factor{x_pc}")
      y_name <- glue::glue("Factor{y_pc}")

      ## Pull out the batches and conditions used in this plot.
      ## Probably could have just used xxx[stuff, drop = TRUE]
      included_batches <- as.factor(as.character(design[[batch_column]]))
      included_conditions <- as.factor(as.character(design[[cond_column]]))

      residual_df <- get_res(svd_result, fData(data), res_slot = "Y", var_slot = "itercosts")
      prop_lst <- residual_df[["prop_var"]]
      if (x_pc > components | y_pc > components) {
        stop("The components plotted must be smaller than the number of components calculated.")
      }

      x_label <- sprintf("%s: %.2f tsne 'variance'", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f tsne 'variance'", y_name, prop_lst[y_pc])
    },
    "umap" = {
      message("Not yet implemented.")
    },
    "uwot" = {
      plotting_data <- t(mtrx)
      pc_table <- as.data.frame(uwot::umap(X = plotting_data, n_components = num_pc, ...))
      rownames(pc_table) <- rownames(fData(data))
      colnames(pc_table) <- glue::glue("PC{1:ncol(pc_table)}")
      x_name <- glue::glue("Factor{x_pc}")
      y_name <- glue::glue("Factor{y_pc}")
      x_label <- x_name
      y_label <- y_name
      included_batch <- as.factor(as.character(design[[batch_column]]))
      included_conditions <- as.factor(as.character(design[[cond_column]]))
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
      svd_result <- fastICA::fastICA(t(mtrx), n.comp = num_pc, alg.typ = alg_type, fun = fun,
                                     alpha = alpha, method = ica_method, row.norm = row.norm,
                                     maxit = maxit, tol = tol, verbose = verbose, w.init = w.init)
      pc_table <- svd_result[["S"]]
      residual_df <- get_res(svd_result, design, res_slot = "S", var_slot = "W")
      prop_lst <- residual_df[[1]]
      rownames(pc_table) <- rownames(design)
      colnames(pc_table) <- glue::glue("IC{1:ncol(pc_table)}")
      x_label <- glue::glue("IC{x_pc}")
      y_label <- glue::glue("IC{y_pc}")
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
      ready <- pcaMethods::prep(t(mtrx), scale = scale, center = center, eps = eps,
                                simple = simple, reverse = reverse)
      pca_result <- try(pcaMethods::pca(ready, method = pc_method, nPcs = num_pc, scale = scale,
                                        center = center, completeObs = completeObs,
                                        subset = subset, cv = cv))
      if (class(pca_result)[1] == "try-error") {
        message("pca invocation failed.  A likely reason if that too many PCs were requested.")
        message("Trying again with 1/2 the PCs.")
        num_pc <- floor(num_pc / 2)
        pca_result <- pcaMethods::pca(ready, method = pc_method, nPcs = num_pc, scale = scale,
                                      center = center, completeObs = completeObs,
                                      subset = subset, cv = cv)
      }
      ## This is mostly guessing on my part.
      svd_result <- list(
        "v" = pca_result@scores,
        "d" = pca_result@sDev)
      residual_df <- get_res(svd_result, design)
      pc_table <- as.data.frame(pca_result@scores)
      prop_lst <- pca_result@R2 * 100
      x_name <- glue::glue("PC{x_pc}")
      y_name <- glue::glue("PC{y_pc}")
      x_label <- sprintf("%s: %.2f%% variance", x_name, prop_lst[x_pc])
      y_label <- sprintf("%s: %.2f%% variance", y_name, prop_lst[y_pc])
    })  ## End of the switch()

  ## An important caveat, some dimension reduction methods remove rows from the data.
  fdata_data <- fData(data)
  fdata_rows <- rownames(fdata_data)
  pc_rows <- rownames(pc_table)
  kept_rows <- fdata_rows %in% pc_rows
  kept_fdata <- fdata_data[kept_rows, ]

  comp_data <- data.frame(
    "sampleid" = rownames(pc_table),
    "condition" = "a",
    "batch" = "a",
    "batch_int" = 1,
    "colors" = "black",
    "text" = kept_fdata[[label_column]],
    "labels" = FALSE)

  comp_data[[x_name]] <- pc_table[, x_pc]
  comp_data[[y_name]] <- pc_table[, y_pc]
  tmp <- as.data.frame(pc_table)
  for (pc in 1:num_pc) {
    oldname <- glue::glue("PC{pc}")
    pc_name <- glue::glue("pc_{pc}")
    comp_data[[pc_name]] <- tmp[[oldname]]
  }
  if (!is.null(size_column)) {
    ## Adding a column with the same name as the size column from the experimental design
    ## and making sure it is a factor.
    comp_data[[size_column]] <- as.factor(design[, size_column])
    ## Just forcing the size to be numeric non-zero.
    comp_data[["size"]] <- as.factor(as.integer(design[[size_column]]) + 1)
  }

  if (isTRUE(plot_title)) {
    plot_title <- what_happened(expt = data)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt = data)
    plot_title <- glue::glue("{plot_title}; {data_title}")
  } else {
    ## Leave the title blank.
  }

  ## The plot_pcs() function gives a decent starting plot
  ## plot_labels = comp_data[["text"]] I think is incorrect, as the plot_labels
  ## parameters is intended to only be a single word defining how to place the
  ## labels.  I think this should instead be point_labels -- but if we set that, given
  ## the large number of dots, we will need to stop any attempted 'smart'
  ## placement of the labels, otherwise it will make the computer very sad.
  ##comp_plot <- plot_pcs(
  ##  comp_data, first = x_name, second = y_name, design = design,
  ##  plot_labels = comp_data[["text"]], x_label = x_label, y_label = y_label,
  ##  plot_title = plot_title, plot_size = plot_size, size_column = size_column,
  ##  plot_alpha = plot_alpha,
  ##  ...)
  comp_plot <- plot_pcs(
    comp_data, first = x_name, second = y_name, design = design,
    x_label = x_label, y_label = y_label,
    plot_title = plot_title, plot_size = plot_size, size_column = size_column,
    plot_alpha = plot_alpha,
    ...)

  ## If plot_title is NULL, print nothing, if it is TRUE
  ## Then give some information about what happened to the data to make the plot.
  ## I tried foolishly to put this in plot_pcs(), but there is no way that receives
  ## my expt containing the normalization state of the data.
  if (isTRUE(plot_title)) {
    data_title <- what_happened(expt = data)
    comp_plot <- comp_plot + ggplot2::ggtitle(data_title)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt = data)
    plot_title <- glue::glue("{plot_title}; {data_title}")
    comp_plot <- comp_plot + ggplot2::ggtitle(plot_title)
  } else {
    ## Leave the title blank.
  }

  ## Finally, return a list of the interesting bits of what happened.
  pca_return <- list(
    "residual_df" = residual_df,
    "prop_var" = prop_lst,
    "plot" = comp_plot,
    "table" = comp_data,
    "result" = svd_result)
  return(pca_return)
}

#' Print a plot of the top-n most PC loaded genes.
#'
#' Sometimes it is nice to know what is happening with the genes which have the
#' greatest effect on a given principal component.  This function provides that.
#'
#' @param expt Input expressionset.
#' @param genes How many genes to observe?
#' @param desired_pc Which component to examine?
#' @param which_scores Perhaps one wishes to see the least-important genes, if
#'   so set this to low.
#' @param ... Extra arguments passed, currently to nothing.
#' @return List containing an expressionset of the subset and a plot of their
#'   expression.
#' @export
plot_pcload <- function(expt, genes = 40, desired_pc = 1, which_scores = "high",
                        ...) {
  arglist <- list(...)
  scores <- pca_highscores(expt, n = genes)

  desired <- data.frame()
  if (which_scores == "high") {
    desired <- scores[["highest"]]
  } else if (which_scores == "low") {
    desired <- scores[["lowest"]]
  } else {
    stop("This only accepts high or low to extract PC scored genes.")
  }

  comp_genes <- desired[, desired_pc]
  comp_genes <- gsub(pattern = "^\\d+\\.\\d+:", replacement = "", x = comp_genes)
  comp_genes_subset <- sm(exclude_genes_expt(expt, ids = comp_genes, method = "keep"))
  samples <- plot_sample_heatmap(comp_genes_subset, row_label = NULL)
  sample_plot <- grDevices::recordPlot()

  retlist <- list(
    "comp_genes_expt" = comp_genes_subset,
    "plot" = sample_plot)
  return(retlist)
}

#' Plot principle components and make them pretty.
#'
#' All the various dimension reduction methods share some of their end-results
#' in common. Most notably a table of putative components which may be plotted
#' against one another so that one may stare at the screen and look for
#' clustering among the samples/genes/whatever.  This function attempts to make
#' that process as simple and pretty as possible.
#'
#' @param pca_data Dataframe of principle components PC1 .. PCN with any other
#'   arbitrary information.
#' @param first Principle component PCx to put on the x axis.
#' @param second Principle component PCy to put on the y axis.
#' @param variances List of the percent variance explained by each component.
#' @param design Experimental design with condition batch factors.
#' @param plot_title Title for the plot.
#' @param plot_labels Parameter for the labels on the plot.
#' @param x_label Label for the x-axis.
#' @param y_label Label for the y-axis.
#' @param plot_size Size of the dots on the plot
#' @param outlines Add a black outline to the plotted shapes?
#' @param plot_alpha Add an alpha channel to the dots?
#' @param size_column Experimental factor to use for sizing the glyphs
#' @param rug Include the rugs on the sides of the plot?
#' @param cis What (if any) confidence intervals to include.
#' @param ... Extra arguments dropped into arglist
#' @return  gplot2 PCA plot
#' @seealso \pkg{ggplot2}
#'  \code{\link[directlabels]{geom_dl}}
#' @examples
#' \dontrun{
#'  pca_plot = plot_pcs(pca_data, first = "PC2", second = "PC4", design = expt$design)
#' }
#' @export
plot_pcs <- function(pca_data, first = "PC1", second = "PC2", variances = NULL,
                     design = NULL, plot_title = TRUE, plot_labels = NULL,
                     x_label = NULL, y_label = NULL, plot_size = 5, outlines = TRUE,
                     plot_alpha = NULL, size_column = NULL, rug = TRUE,
                     cis = c(0.95, 0.9), ...) {
  arglist <- list(...)
  batches <- as.factor(pca_data[["batch"]])
  label_column <- "condition"
  if (!is.null(arglist[["label_column"]])) {
    label_column <- arglist[["label_column"]]
  }
  point_labels <- factor(pca_data[[label_column]])
  if (!is.null(arglist[["point_labels"]])) {
    point_labels <- arglist[["point_labels"]]
  }
  if (is.null(plot_title)) {
    plot_title <- paste(first, " vs. ", second, sep = "")
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
  plot_legend <- TRUE
  if (!is.null(arglist[["plot_legend"]])) {
    plot_legend <- arglist[["plot_legend"]]
  }

  pca_plot <- NULL
  color_listing <- pca_data[, c("condition", "colors")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["colors"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  ## Ok, so this is shockingly difficult.  For <5 batch data I want properly
  ## colored points with black outlines The legend colors need to match, in
  ## addition, the legend needs to have the shapes noted.
  ## In order to do this, one must do, _in_order_:
  ## 1.  Set up the normal ggplot object
  ## 2.  Set up a geom_point with color _and_ fill as the proper color.
  ##     The color but _NOT_ fill is used to color the legend's copy of the glyph.
  ## 3.  Then set up a new geom_point with color = black _and_ show_guide = FALSE
  ## 4.  Then set scale_color_manual to the proper color_list
  ## 5.  Then set scale_fill_manual to the proper color_list
  ## 6.  Finally, set the shape manual with a guide_legend override

  ## Step 1
  if (is.null(plot_alpha)) {
    plot_alpha <- 1
  }

  pca_data <- as.data.frame(pca_data)
  pca_data[["condition"]] <- as.factor(pca_data[["condition"]])
  pca_data[["batch"]] <- as.factor(pca_data[["batch"]])
  pca_plot <- ggplot(data = pca_data,
                     aes_string(x = "get(first)", y = "get(second)", text = "sampleid"))

  ## Add a little check to only deal with the confidence-interval-able data.
  count <- NULL
  ci_keepers <- pca_data %>%
    group_by(!!sym(ci_group)) %>%
    summarise(count = n()) %>%
    filter(count > 3)
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
        ggplot2::stat_ellipse(data = ci_data,
                              mapping = aes_string(group = ci_group, fill = ci_fill),
                              geom = "polygon", type = "t", level = ci, alpha = alpha)
    }
  }

  minimum_size <- 2
  maximum_size <- 2
  if (!is.null(size_column)) {
    maximum_size <- max(levels(pca_data[["size"]]))
  }

  if (is.null(size_column) & num_batches <= 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(size = plot_size, alpha = plot_alpha,
                          aes_string(shape = "batches",
                                     colour = "condition",
                                     fill = "condition"))
    if (isTRUE(outlines)) {
      pca_plot <- pca_plot +
        ggplot2::geom_point(size = plot_size, alpha = plot_alpha,
                            colour = "black",
                            ## show.legend = FALSE,
                            aes_string(shape = "batches",
                                       fill = "condition"))
    }
    pca_plot <- pca_plot +
      ggplot2::scale_color_manual(name = "Condition",
                                  ## guide = "legend",
                                  guide = "legend",
                                  values = color_list) +
      ggplot2::scale_fill_manual(name = "Condition",
                                 ## guide = "legend",
                                 guide = "none",
                                 values = color_list) +
      ggplot2::scale_shape_manual(
                 name = "Batch",
                 labels = levels(as.factor(pca_data[["batch"]])),
                 guide = ggplot2::guide_legend(override.aes = list(size = plot_size, fill = "grey")),
                 values = 21:25)
  } else if (is.null(size_column) & num_batches > 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(size = plot_size, alpha = plot_alpha,
                          aes_string(shape = "batches",
                                     fill = "condition",
                                     colour = "condition")) +
      ggplot2::scale_color_manual(name = "Condition",
                                  ##guide = ggplot2::guide_legend(override.aes = list(size = plot_size)),
                                  guide = "legend",
                                  values = color_list) +
      ggplot2::scale_fill_manual(name = "Condition",
                                 ##guide = ggplot2::guide_legend(override.aes = list(size = plot_size)),
                                 guide = "none",
                                 values = color_list) +
      ggplot2::scale_shape_manual(name = "Batch",
                                  labels = levels(as.factor(pca_data[["batch"]])),
                                  guide = ggplot2::guide_legend(overwrite.aes = list(size = plot_size)),
                                  values = 1:num_batches)
  } else if (!is.null(size_column) & num_batches <= 5) {
    ## This will require the 6 steps above and one more
    pca_plot <- ggplot(data = as.data.frame(pca_data),
                       aes_string(x = "get(first)", y = "get(second)", text = "sampleid",
                                  shape = "batches")) +
      ggplot2::geom_point(alpha = plot_alpha,
                          aes_string(shape = "batches",
                                     size = "size",
                                     colour = "condition",
                                     fill = "condition"))
    if (isTRUE(outlines)) {
      pca_plot <- pca_plot +
        ggplot2::geom_point(alpha = plot_alpha, colour = "black", show.legend = FALSE,
                            aes_string(size = "size", shape = "batches", fill = "condition"))

      ##size = plot_size, alpha = plot_alpha, colour = "black", show.legend = FALSE,
      ##aes_string(shape = "batches",
      ##fill = "condition"))
    }
    pca_plot <- pca_plot +
      ggplot2::geom_point(colour = "black", alpha = plot_alpha, show.legend = FALSE,
                          aes_string(shape = "batches",
                                     size = "size",
                                     fill = "condition")) +
      ggplot2::scale_color_manual(name = "Condition",
                                  guide = "legend",
                                  values = color_list) +
      ggplot2::scale_fill_manual(name = "Condition",
                                 guide = ggplot2::guide_legend(override.aes = list(size = plot_size)),
                                 values = color_list) +
      ggplot2::scale_shape_manual(
                 name = "Batch",
                 labels = levels(as.factor(pca_data[["batch"]])),
                 guide = ggplot2::guide_legend(override.aes = list(size = plot_size, fill = "grey")),
                 values = 21:25) +
      ggplot2::scale_size_manual(name = size_column,
                                 labels = levels(pca_data[[size_column]]),
                                 values = as.numeric(levels(pca_data[["size"]])))
  } else if (!is.null(size_column) & num_batches > 5) {
    pca_plot <- pca_plot +
      ggplot2::geom_point(alpha = plot_alpha,
                          aes_string(shape = "batches",
                                     colour = "pca_data[['condition']]",
                                     size = "size")) +
      ggplot2::scale_shape_manual(name = "Batch",
                                  labels = levels(as.factor(pca_data[["batch"]])),
                                  guide = ggplot2::guide_legend(overwrite.aes = list(size = plot_size)),
                                  values = 1:num_batches) +
      ggplot2::scale_color_manual(name = "Condition",
                                  guide = "legend",
                                  values = color_list) +
      ##ggplot2::scale_color_identity(name = "Condition",
      ##                              guide = "legend",
      ##                              values = color_list) +
      ggplot2::scale_size_manual(name = size_column,
                                 labels = levels(pca_data[[size_column]]),
                                 values = as.numeric(levels(pca_data[["size"]])))
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
    x_var_num <- as.numeric(gsub(pattern = "PC", replacement = "", x = first))
    y_var_num <- as.numeric(gsub(pattern = "PC", replacement = "", x = second))
    x_label <- glue::glue("PC1{first}: {variances[[x_var_num]]}% variance")
    y_label <- glue::glue("PC2{second}: {variances[[y_var_num]]}% variance")
    pca_plot <- pca_plot +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label)
  }

  if (isTRUE(rug)) {
    pca_plot <- pca_plot + ggplot2::geom_rug(colour = "gray50", alpha = 0.7)
  }

  if (is.null(plot_labels)) {
    plot_labels <- "repel"
  }
  if (isFALSE(plot_labels)) {
    message("Not putting labels on the PC plot.")
  } else if (plot_labels == "normal") {
    pca_plot <- pca_plot +
      ggplot2::geom_text(aes_string(x = "PC1", y = "PC2", label = "labels",
                                    angle = 45, size = label_size, vjust = 2))
  } else if (plot_labels == "repel") {
    pca_plot <- pca_plot +
      ggrepel::geom_text_repel(aes_string(label = "labels"),
                               size = label_size, box.padding = ggplot2::unit(0.5, "lines"),
                               point.padding = ggplot2::unit(1.6, "lines"),
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  } else if (plot_labels == "dlsmart") {
    pca_plot <- pca_plot +
      directlabels::geom_dl(aes_string(label = "labels"), method = "smart.grid")
  } else {
    pca_plot <- pca_plot +
      directlabels::geom_dl(aes_string(label = "labels"), method = "first.qp")
  }

  legend_position <- "right"
  if (isFALSE(plot_legend)) {
    legend_position <- "none"
  }

  ## Set default font sizes and colors
  pca_plot <- pca_plot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   legend.position = legend_position,
                   legend.key.size = grid::unit(0.5, "cm"))

  return(pca_plot)
}

#' Plot a PC plot with options suitable for ggplotly.
#'
#' @param data an expt set of samples.
#' @param design a design matrix and.
#' @param plot_colors a color scheme.
#' @param plot_title a title for the plot.
#' @param plot_size size for the glyphs on the plot.
#' @param plot_alpha Add an alpha channel to the dots?
#' @param plot_labels add labels?  Also, what type?  FALSE, "default", or "fancy".
#' @param size_column use an experimental factor to size the glyphs of the plot
#' @param pc_method how to extract the components? (svd
#' @param x_pc Component to put on the x axis.
#' @param y_pc Component to put on the y axis.
#' @param outlines Include black outlines around glyphs?
#' @param num_pc How many components to calculate, default to the number of
#'   rows in the metadata.
#' @param expt_names Column or character list of preferred sample names.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param tooltip Which columns to include in the tooltip.
#' @param ...  Arguments passed through to the pca implementations and plotter.
#' @return This passes directly to plot_pca(), so its returns should be
#'   applicable along with the result from ggplotly.
#' @export
plotly_pca <-  function(data, design = NULL, plot_colors = NULL, plot_title = NULL,
                        plot_size = 5, plot_alpha = NULL, plot_labels = NULL, size_column = NULL,
                        pc_method = "fast_svd", x_pc = 1, y_pc = 2, outlines = FALSE,
                        num_pc = NULL, expt_names = NULL, label_chars = 10,
                        tooltip = c("shape", "fill", "sampleid"),
                        ...) {
  pca_result <- plot_pca(data, design = design, plot_colors = plot_colors,
                         plot_title = plot_title, plot_size = plot_size,
                         plot_alpha = plot_alpha, plot_labels = plot_labels,
                         size_column = size_column, pc_method = pc_method,
                         x_pc = x_pc, y_pc = y_pc, outlines = outlines, num_pc = num_pc,
                         expt_names = expt_names, label_chars = label_chars, ...)
  plotly_result <- plotly::ggplotly(pca_result[["plot"]], tooltip = tooltip)
  retlist <- pca_result
  retlist[["plotly"]] <- plotly_result
  return(retlist)
}

#' Shortcut to plot_pca(pc_method = "tsne")
#'
#' @param ...  Arguments for plot_pca()
#' @export
plot_tsne <- function(...) {
  plot_pca(..., pc_method = "tsne")
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
  plotted_u1s <- plotted_us[order(plotted_us[, 1], decreasing = TRUE), ]
  plotted_u2s <- plotted_us[order(plotted_us[, 2], decreasing = TRUE), ]
  plotted_u3s <- plotted_us[order(plotted_us[, 3], decreasing = TRUE), ]
  ## allS <- BiocGenerics::rank(allS, ties.method = "random")
  ## plotted_us$rank = rank(plotted_us[,1], ties.method = "random")
  plotted_u1s <- cbind(plotted_u1s, rev(rank(plotted_u1s[, 1], ties.method = "random")))
  plotted_u1s <- plotted_u1s[, c(1, 4)]
  colnames(plotted_u1s) <- c("PC1", "rank")
  plotted_u1s <- data.frame(plotted_u1s)
  plotted_u1s[["ID"]] <- as.character(rownames(plotted_u1s))
  plotted_u2s <- cbind(plotted_u2s, rev(rank(plotted_u2s[, 2], ties.method = "random")))
  plotted_u2s <- plotted_u2s[, c(2, 4)]
  colnames(plotted_u2s) <- c("PC2", "rank")
  plotted_u2s <- data.frame(plotted_u2s)
  plotted_u2s[["ID"]] <- as.character(rownames(plotted_u2s))
  plotted_u3s <- cbind(plotted_u3s, rev(rank(plotted_u3s[, 3], ties.method = "random")))
  plotted_u3s <- plotted_u3s[, c(3, 4)]
  colnames(plotted_u3s) <- c("PC3", "rank")
  plotted_u3s <- data.frame(plotted_u3s)
  plotted_u3s[["ID"]] <- as.character(rownames(plotted_u3s))
  plotted_us <- merge(plotted_u1s, plotted_u2s, by.x = "rank", by.y = "rank")
  plotted_us <- merge(plotted_us, plotted_u3s, by.x = "rank", by.y = "rank")
  colnames(plotted_us) <- c("rank", "PC1", "ID1", "PC2", "ID2", "PC3", "ID3")
  rm(plotted_u1s)
  rm(plotted_u2s)
  rm(plotted_u3s)
  ## top_threePC = head(plotted_us, n = 20)
  plotted_us <- plotted_us[, c("PC1", "PC2", "PC3")]
  plotted_us[, "ID"] <- rownames(plotted_us)
  message("More shallow curves in these plots suggest more genes in this principle component.")
  plot(plotted_us)
  u_plot <- grDevices::recordPlot()
  return(u_plot)
}


## EOF
