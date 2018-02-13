#' Plot tnse data for the genes in a data set.
#'
#' @param data  Some data!
#' @param design a design!
#' @param plot_colors  Some colors!
#' @param seed for tsne
#' @param chosen_features  Use these genes for labeling
#' @param number_features  Somethingsomething
#' @param perplexity  for tsne
#' @param min_variance  Only include genes with more variance than this.
#' @param plot_title  A title!
#' @param components  How many components to plot.
#' @param iterations  Use x tsne iterations.
#' @param theta  Yay greek!
#' @param pca  Seed this with an initial pca plot?
#' @param component_x  Put which component on the x-axis?
#' @param component_y  And which component on the y-axis?
#' @param ...  Arglist arguments.
#' @export
plot_tsne_genes <- function(data, design=NULL, plot_colors=NULL, seed=1,
                            chosen_features=NULL, number_features=NULL,
                            perplexity=NULL, min_variance=0.01, plot_title=NULL,
                            components=2, iterations=1000, theta=0.3, pca=TRUE,
                            component_x=1, component_y=2,  ...) {
  ## I have been using hpgl_env for keeping aes() from getting contaminated.
  ## I think that this is no longer needed because I have been smater(sic) about how
  ## I invoke aes_string() and ggplot2()
  hpgl_env <- environment()
  arglist <- list(...)
  plot_names <- arglist[["plot_names"]]
  ## Set default columns in the experimental design for condition and batch
  ## changing these may be used to query other experimental factors with pca.
  cond_column <- "condition"
  if (!is.null(arglist[["cond_column"]])) {
    cond_column <- arglist[["cond_column"]]
    message(paste0("Using ", cond_column, " as the condition column in the experimental design."))
  }
  batch_column <- "batch"
  if (!is.null(arglist[["batch_column"]])) {
    batch_column <- arglist[["batch_column"]]
    message(paste0("Using ", batch_column, " as the batch column in the experimental design."))
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
  expt <- NULL
  if (data_class == "expt") {
    expt <- data
    design <- data[["design"]]
    if (cond_column == "condition") {
      plot_colors <- data[["colors"]]
    } else {
      plot_colors <- NULL
    }
    plot_names <- data[["samplenames"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "list") {
    data <- data[["count_table"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table element.")
    }
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Check that the given design works with the data
  ## Prune the design if necessary
  ## Also take into account the fact that sometimes I change the case of hpgl<->HPGL
  given_samples <- tolower(colnames(data))
  colnames(data) <- given_samples
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
  ## All of the above is logic stolen from plot_pca()
  ## I increasingly think I should just fold this into it.

  ## A bunch of the logic in this section is taken from scater.
  if (is.null(perplexity)) {
    perplexity <- floor(ncol(data) / 5)
  }

  ## I am curious to know why the order of the genes by variance is significant.
  plotting_indexes <- 1:nrow(data)
  if (is.null(chosen_features)) {
    variances <- matrixStats::rowVars(as.matrix(data))
    if (!is.null(number_features)) {
      number_features <- min(number_features, nrow(data))
    } else {
      number_features <- nrow(data)
    }
    plotting_indexes <- order(variances, decreasing=TRUE)[1:number_features]
  }

  plotting_data <- data[plotting_indexes, ]
  ## This I do understand and think is cool
  ## Drop features with low variance
  keepers <- (matrixStats::rowVars(as.matrix(plotting_data)) >= min_variance)
  keepers[is.na(keepers)] <- FALSE ## Another nice idea
  plotting_data <- plotting_data[keepers, ]

  ## There is an interesting standardization idea in scater
  ## But I think I would prefer to have flexibility here
  ## exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))
  if (!is.null(seed)) {
    set.seed(seed)
  }

  sne <- Rtsne::Rtsne(plotting_data,
                      check_duplicates=FALSE,
                      dims=components,
                      max_iter=iterations,
                      pca=pca,
                      theta=theta,
                      perplexity=perplexity)
  sne_df <- as.data.frame(sne[["Y"]])

  rownames(sne_df) <- rownames(plotting_data)
  sne_df <- sne_df[, 1:components]

  ## Pull out the batches and conditions used in this plot.
  ## Probably could have just used xxx[stuff, drop=TRUE]
  included_batches <- as.factor(as.character(design[[batch_column]]))
  included_conditions <- as.factor(as.character(design[[cond_column]]))


  tsne_data <- data.frame(
    "sampleid" = as.character(design[["sampleid"]]),
    "condition" = as.character(design[[cond_column]]),
    "batch" = as.character(design[[batch_column]]),
    "batch_int" = as.integer(as.factor(design[[batch_column]])),
    "colors" = as.character(plot_colors),
    "labels" = label_list)
  ##tsne_data[[compname_x]] <- sne_df[[paste0("V", component_x)]]
  ##tsne_data[[compname_y]] <- sne_df[[paste0("V", component_y)]]

  a_plot <- plot_scatter(sne_df)
  return(a_plot)
}

#' Make a ggplot TSNE plot describing the samples' clustering.
#'
#' @param data  an expt set of samples.
#' @param design   a design matrix and.
#' @param plot_colors   a color scheme.
#' @param seed  A seed for Rtsne
#' @param chosen_features  Use these features?
#' @param number_features  And this number.
#' @param perplexity I am perplexed.
#' @param min_variance  Only include genes with more than this variance.
#' @param plot_title   a title for the plot.
#' @param plot_size   size for the glyphs on the plot.
#' @param plot_labels   add labels?  Also, what type?  FALSE, "default", or "fancy".
#' @param size_column use an experimental factor to size the glyphs of the plot
#' @param components  Look for n components.
#' @param iterations  Perform n iterations of tsne.
#' @param theta  Yay greek!
#' @param pca  Seed with an initial pca plot?
#' @param component_x  Which component goes on the x-axis?
#' @param component_y  And which goes on the y-axis?
#' @param ...  arglist from elipsis!
#' @return a list containing the following:
#' \enumerate{
#'  \item  plot = a plot
#' }
#' @seealso \pkg{directlabels}
#'  \code{\link[directlabels]{geom_dl}} \code{\link{plot_pcs}}
#' @examples
#' \dontrun{
#'  tsne_plot <- plot_tsne(expt=expt)
#'  tsne_plot
#' }
#' @export
plot_tsne <- function(data, design=NULL, plot_colors=NULL, seed=1,
                      chosen_features=NULL, number_features=NULL,
                      plot_labels=NULL, perplexity=NULL,
                      min_variance=0.001, plot_title=NULL, plot_size=5,
                      size_column=NULL, components=2, iterations=1000,
                      theta=0.3, pca=TRUE, component_x=1, component_y=2,  ...) {
  ## I have been using hpgl_env for keeping aes() from getting contaminated.
  ## I think that this is no longer needed because I have been smater(sic) about how
  ## I invoke aes_string() and ggplot2()
  hpgl_env <- environment()
  arglist <- list(...)
  plot_names <- arglist[["plot_names"]]
  ## Set default columns in the experimental design for condition and batch
  ## changing these may be used to query other experimental factors with pca.
  cond_column <- "condition"
  if (!is.null(arglist[["cond_column"]])) {
    cond_column <- arglist[["cond_column"]]
    message(paste0("Using ", cond_column, " as the condition column in the experimental design."))
  }
  batch_column <- "batch"
  if (!is.null(arglist[["batch_column"]])) {
    batch_column <- arglist[["batch_column"]]
    message(paste0("Using ", batch_column, " as the batch column in the experimental design."))
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
  expt <- NULL
  if (data_class == "expt") {
    expt <- data
    design <- data[["design"]]
    if (cond_column == "condition") {
      plot_colors <- data[["colors"]]
    } else {
      plot_colors <- NULL
    }
    plot_names <- data[["samplenames"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "list") {
    data <- data[["count_table"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table element.")
    }
  } else if (data_class == "matrix" | data_class == "data.frame") {
    data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
  }

  ## Check that the given design works with the data
  ## Prune the design if necessary
  ## Also take into account the fact that sometimes I change the case of hpgl<->HPGL
  given_samples <- tolower(colnames(data))
  colnames(data) <- given_samples
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
    design <- cbind(plot_labels, 1)
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
  ## All of the above is logic stolen from plot_pca()
  ## I increasingly think I should just fold this into it.

  ## A bunch of the logic in this section is taken from scater.
  if (is.null(perplexity)) {
    perplexity <- floor(ncol(data) / 5)
  }

  ## I am curious to know why the order of the genes by variance is significant.
  plotting_indexes <- 1:nrow(data)
  if (is.null(chosen_features)) {
    variances <- matrixStats::rowVars(as.matrix(data))
    if (!is.null(number_features)) {
      number_features <- min(number_features, nrow(data))
    } else {
      number_features <- nrow(data)
    }
    plotting_indexes <- order(variances, decreasing=TRUE)[1:number_features]
  }

  plotting_data <- data[plotting_indexes, ]
  ## This I do understand and think is cool
  ## Drop features with low variance
  keepers <- (matrixStats::rowVars(as.matrix(plotting_data)) >= min_variance)
  keepers[is.na(keepers)] <- FALSE ## Another nice idea
  plotting_data <- plotting_data[keepers, ]

  ## There is an interesting standardization idea in scater
  ## But I think I would prefer to have flexibility here
  ## exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))
  if (!is.null(seed)) {
    set.seed(seed)
  }

  sne <- Rtsne::Rtsne(t(plotting_data),
                      check_duplicates=FALSE,
                      dims=components,
                      max_iter=iterations,
                      pca=pca,
                      theta=theta,
                      perplexity=perplexity)
  sne_df <- as.data.frame(sne[["Y"]])
  rownames(sne_df) <- rownames(design)
  sne_df <- sne_df[, 1:components]

  ## Pull out the batches and conditions used in this plot.
  ## Probably could have just used xxx[stuff, drop=TRUE]
  included_batches <- as.factor(as.character(design[[batch_column]]))
  included_conditions <- as.factor(as.character(design[[cond_column]]))
  if (length(levels(included_conditions)) == 1 & length(levels(included_batches)) == 1) {
    warning("There is only one condition and one batch, it is impossible to get meaningful pcRes information.")
  } else if (length(levels(included_conditions)) == 1) {
    warning("There is only one condition, but more than one batch.
Going to run pcRes with the batch information.")
    tsne_residuals <- tsne_res(Y=sne[["Y"]],
                               costs=sne[["costs"]],
                               batch=design[, batch_column])
  } else if (length(levels(included_batches)) == 1) {
    message("There is just one batch in this data.")
    tsne_residuals <- tsne_res(Y=sne[["Y"]],
                               costs=sne[["costs"]],
                               condition=design[, cond_column])
  } else {
    tsne_residuals <- tsne_res(Y=sne[["Y"]],
                               costs=sne[["costs"]],
                               condition=design[, cond_column],
                               batch=design[, batch_column])
  }

  if (component_x > components | component_y > components) {
    stop("The components plotted must be smaller than the number of components calculated.")
  }
  compname_x <- paste0("Comp", component_x)
  compname_y <- paste0("Comp", component_y)
  tsne_data <- data.frame(
    "sampleid" = as.character(design[["sampleid"]]),
    "condition" = as.character(design[[cond_column]]),
    "batch" = as.character(design[[batch_column]]),
    "batch_int" = as.integer(as.factor(design[[batch_column]])),
    "colors" = as.character(plot_colors),
    "labels" = label_list)
  tsne_data[[compname_x]] <- sne_df[[paste0("V", component_x)]]
  tsne_data[[compname_y]] <- sne_df[[paste0("V", component_y)]]

  a_plot <- plot_pcs(
    tsne_data,
    first=compname_x,
    second=compname_y,
    design=design,
    plot_labels=plot_labels,
    plot_size=plot_size,
    size_column=size_column, ...)

  ## components of our plot.
  tsne_variance <- round((sne[["costs"]] ^ 2) / sum(sne[["costs"]] ^ 2) * 100, 2)
  xl <- sprintf("Comp%s: %.2f%% rsquared", component_x, tsne_residuals[["cond.R2"]][[1]])
  yl <- sprintf("Comp%s: %.2f%% rsquared", component_y, tsne_residuals[["batch.R2"]][[1]])
  ## The following are some pretty-ifiers for the plot, they should be moved into plot_pcs
  a_plot <- a_plot +
    ggplot2::xlab(xl) +
    ggplot2::ylab(yl) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.key.size=grid::unit(0.5, "cm"))

  ## If plot_title is NULL, print nothing, if it is TRUE
  ## Then give some information about what happened to the data to make the plot.
  if (isTRUE(plot_title)) {
    data_title <- what_happened(expt=expt)
    a_plot <- a_plot + ggplot2::ggtitle(data_title)
  } else if (!is.null(plot_title)) {
    data_title <- what_happened(expt=expt)
    plot_title <- paste0(plot_title, "; ", data_title)
    a_plot <- a_plot + ggplot2::ggtitle(plot_title)
  } else {
    ## Leave the title blank.
  }

  retlist <- list(
    "tsne_data" = sne,
    "plot" = a_plot,
    "table" = tsne_data,
    "res" = tsne_residuals,
    "variance" = tsne_variance)
  return(retlist)
}

tsne_res <- function(Y, costs, condition=NULL, batch=NULL) {
  tsne_var <- round((costs ^ 2) / sum(costs ^ 2) * 100, 2)
  tsne_cum <- cumsum(tsne_var)
  conditional_rsquared <- NULL
  batch_rsquared <- NULL
  if(!is.null(condition)) {
    cond.R2 <- function(y) {
      round(summary(lm(y ~ condition))[["r.squared"]] * 100, 2)
    }
    conditional_rsquared <- apply(Y, 2, cond.R2)
  }
  if (!is.null(batch)) {
    batch.R2 <- function(y) {
      round(summary(lm(y ~ batch))[["r.squared"]] * 100, 2)
    }
    batch_rsquared <- apply(Y, 2, batch.R2)
  }

  result <- list(
    "propVar" = tsne_var,
    "cumPropVar" = tsne_cum,
    "cond.R2" = conditional_rsquared,
    "batch.R2" = batch_rsquared)
  return(result)
}

## EOF
