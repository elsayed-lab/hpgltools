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
#' @export
pcRes <- function(v, d, condition=NULL, batch=NULL){
  pcVar <- round((d ^ 2) / sum(d ^ 2) * 100, 2)
  cumPcVar <- cumsum(pcVar)
  calculate_rsquared_condition <- function(data) {
      lm_result <- lm(data ~ condition)
  }
  if(!is.null(condition)) {
      cond.R2 <- function(y) {
          round(summary(lm(y ~ condition))$r.squared * 100, 2)
      }
      cond.R2 <- apply(v, 2, cond.R2)
  }
  if(!is.null(batch)) {
      batch.R2 <- function(y) {
          round(summary(lm(y ~ batch))$r.squared * 100, 2)
      }
      batch.R2 <- apply(v, 2, batch.R2)
  }
  if(is.null(condition) & is.null(batch)){
      res <- data.frame("propVar"=pcVar,
                        "cumPropVar"=cumPcVar)
  }
  if(!is.null(batch) & is.null(condition)){
    res <- data.frame("propVar"=pcVar,
                      "cumPropVar"=cumPcVar,
                      "batch.R2"=batch.R2)
  }
  if(!is.null(condition) & is.null(batch)){
    res <- data.frame("propVar"=pcVar,
                      "cumPropVar"=cumPcVar,
                      "cond.R2"=cond.R2)
  }
  if(!is.null(condition) & !is.null(batch)){
    res <- data.frame("propVar"=pcVar,
                      "cumPropVar"=cumPcVar,
                      "cond.R2"=cond.R2,
                      "batch.R2"=batch.R2)
  }
  return(res)
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
#' @param ...  arglist from elipsis!
#' @return a list containing the following:
#'  \enumerate{
#'   \item  pca = the result of fast.svd()
#'   \item  plot = ggplot2 pca_plot describing the principle component analysis of the samples.
#'   \item  table = a table of the PCA plot data
#'   \item  res = a table of the PCA res data
#'   \item  variance = a table of the PCA plot variance
#'  }
#' @seealso \code{\link{makeSVD}},
#' \code{\link[directlabels]{geom_dl}} \code{\link{plot_pcs}}
#' @examples
#' \dontrun{
#'  pca_plot <- plot_pca(expt=expt)
#'  pca_plot
#' }
#' @export
plot_pca <- function(data, design=NULL, plot_colors=NULL, plot_labels=NULL,
<<<<<<< HEAD
                     plot_title=NULL, plot_size=5, size_column=NULL, ...) {
=======
                     plot_title=TRUE, plot_size=5, size_column=NULL, ...) {
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    ## I have been using hpgl_env for keeping aes() from getting contaminated.
    ## I think that this is no longer needed because I have been smater(sic) about how
    ## I invoke aes_string() and ggplot2()
    hpgl_env <- environment()
    arglist <- list(...)
    plot_names <- arglist[["plot_names"]]
<<<<<<< HEAD
=======

>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
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
        data <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
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
        colnames(design) <- c("name","batch","condition")
        design <- design[, c("name","condition","batch")]
        plot_names <- design[["name"]]
    }

<<<<<<< HEAD
    ## Different folks like different labels.  I prefer hpglxxxx, but others have asked for
    ## condition_batch; this handles that as eloquently as I am able.
=======
    ## Different folks like different labels.  I prefer hpglxxxx, but others have asked for condition_batch
    ## This handles that as eloquently as I am able.
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
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

<<<<<<< HEAD
    svd_result <- corpcor::fast.svd(as.matrix(data) - rowMeans(as.matrix(data)))
    v_vector <- svd_result[["v"]]
    rownames(v_vector) <- colnames(data)
    svd_result[["v"]] <- v_vector
    ##pca <- makeSVD(data)
    ## Pull out the batches and conditions used in this plot.
    ## Probably could have just used xxx[stuff, drop=TRUE]
=======
    ## I think the makeSVD function is kind of dumb.  It calls fast.svd() and gives me back the pieces.
    ## But I keep this function call as a reminder that Kwame (who wrote our initial pca plotter) was a nice guy.
    ## Also, keep in mind that this is a fast.svd of the (matrix - rowmeans(matrix))
    pca <- makeSVD(data)
    ## Pull out the batches and conditions used in this plot.  Probably could have just used xxx[stuff, drop=TRUE]
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    included_batches <- as.factor(as.character(design[[batch_column]]))
    included_conditions <- as.factor(as.character(design[[cond_column]]))

    ## Depending on how much batch/condition information is available, invoke pcRes() to get some idea of how
    ## much variance in a batch model is accounted for with each PC.
    if (length(levels(included_conditions)) == 1 & length(levels(included_batches)) == 1) {
        warning("There is only one condition and one batch, it is impossible to get meaningful pcRes information.")
    } else if (length(levels(included_conditions)) == 1) {
        warning("There is only one condition, but more than one batch.   Going to run pcRes with the batch information.")
        pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]], batch=design[, batch_column])
    } else if (length(levels(included_batches)) == 1) {
        message("There is just one batch in this data.")
<<<<<<< HEAD
        pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]],
                         condition=design[, cond_column])
=======
        pca_res <- pcRes(v=pca[["v"]], d=pca[["d"]], condition=design[, cond_column])
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    } else {
        pca_res <- pcRes(v=svd_result[["v"]], d=svd_result[["d"]],
                         condition=design[, cond_column], batch=design[, batch_column])
    }

    ## By a similar token, get the percentage of variance accounted for in each PC
<<<<<<< HEAD
    pca_variance <- round((svd_result[["d"]] ^ 2) / sum(svd_result[["d"]] ^ 2) * 100, 2)
    ## These will provide metrics on the x/y axes showing the amount of variance on those
    ## components of our plot.
=======
    pca_variance <- round((pca[["d"]] ^ 2) / sum(pca[["d"]] ^ 2) * 100, 2)
    ## These will provide metrics on the x/y axes showing the amount of variance on those components of our plot.
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    xl <- sprintf("PC1: %.2f%% variance", pca_variance[1])
    yl <- sprintf("PC2: %.2f%% variance", pca_variance[2])

    ## Create a data frame with all the material of interest in the actual PCA plot
    pca_data <- data.frame(
        "sampleid" = as.character(design[["sampleid"]]),
        "condition" = as.character(design[[cond_column]]),
        "batch" = as.character(design[[batch_column]]),
        "batch_int" = as.integer(as.factor(design[[batch_column]])),
        "PC1" = svd_result[["v"]][, 1],
        "PC2" = svd_result[["v"]][, 2],
        "colors" = as.character(plot_colors),
        "labels" = label_list)

    ## Add an optional column which may be used to change the glyph sizes in the plot
    if (!is.null(size_column)) {
        pca_data[[size_column]] <- as.integer(as.factor(design[, size_column]))
        pca_data[[size_column]] <- pca_data[[size_column]] + 1
    }

    pca_plot <- NULL
    ## The plot_pcs() function gives a decent starting plot
    pca_plot <- plot_pcs(pca_data, first="PC1", second="PC2",
                         design=design, plot_labels=plot_labels,
                         plot_size=plot_size, size_column=size_column, ...)

    ## The following are some pretty-ifiers for the plot, they should be moved into plot_pcs
    pca_plot <- pca_plot +
        ggplot2::xlab(xl) +
        ggplot2::ylab(yl) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key.size=grid::unit(0.5, "cm"))

    ## If plot_title is NULL, print nothing, if it is TRUE
    ## Then give some information about what happened to the data to make the plot.
    if (isTRUE(plot_title)) {
        data_title <- what_happened(expt=expt)
        pca_plot <- pca_plot + ggplot2::ggtitle(data_title)
    } else if (!is.null(plot_title)) {
        data_title <- what_happened(expt=expt)
        plot_title <- paste0(plot_title, "; ", data_title)
<<<<<<< HEAD
    } else {
        ## Leave the title blank.
=======
        pca_plot <- pca_plot + ggplot2::ggtitle(plot_title)
    } else {
        ## Leave the title blank
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    }

    ## Finally, return a list of the interesting bits of what happened.
    pca_return <- list(
        "pca" = svd_result,
        "plot" = pca_plot,
        "table" = pca_data,
        "res" = pca_res,
        "variance" = pca_variance)
    return(pca_return)
}

#' Collect the r^2 values from a linear model fitting between a singular
#' value decomposition and factor.
#'
#' @param svd_v V' V = I portion of a fast.svd call.
#' @param fact Experimental factor from the original data.
#' @param type Make this categorical or continuous with factor/continuous.
#' @return The r^2 values of the linear model as a percentage.
#' @seealso \code{\link[corpcor]{fast.svd}}
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

#' A quick and dirty PCA plotter of arbitrary components against one another.
#'
#' @param pca_data  a dataframe of principle components PC1 .. PCN with any other arbitrary information.
#' @param first   principle component PCx to put on the x axis.
#' @param second   principle component PCy to put on the y axis.
#' @param variances   a list of the percent variance explained by each component.
#' @param design   the experimental design with condition batch factors.
#' @param plot_title   a title for the plot.
#' @param plot_labels   a parameter for the labels on the plot.
#' @param plot_size  The size of the dots on the plot
#' @param size_column an experimental factor to use for sizing the glyphs
#' @param ... extra arguments dropped into arglist
#' @return a ggplot2 PCA plot
#' @seealso \pkg{ggplot2} \code{\link[directlabels]{geom_dl}}
#' @examples
#' \dontrun{
#'  pca_plot = plot_pcs(pca_data, first="PC2", second="PC4", design=expt$design)
#' }
#' @export
plot_pcs <- function(pca_data, first="PC1", second="PC2", variances=NULL,
<<<<<<< HEAD
                     design=NULL, plot_title=TRUE, plot_labels=NULL,
                     plot_size=5, size_column=NULL, ...) {
=======
                     design=NULL, plot_title=TRUE, plot_labels=NULL, plot_size=5, size_column=NULL, ...) {
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
    arglist <- list(...)
    hpgl_env <- environment()
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
    pca_plot <- NULL

    color_listing <- pca_data[, c("condition","colors")]
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
    pca_plot <- ggplot(data=as.data.frame(pca_data), aes_string(x="get(first)", y="get(second)"), environment=hpgl_env)

    if (is.null(size_column) & num_batches <= 5) {
        pca_plot <- pca_plot +
            ggplot2::geom_point(size=plot_size,
                                aes_string(shape="as.factor(batches)",
                                           colour="as.factor(condition)",
                                           fill="as.factor(condition)")) +
            ggplot2::geom_point(size=plot_size, colour="black", show.legend=FALSE,
                                aes_string(shape="as.factor(batches)",
                                           fill="as.factor(condition)")) +
            ggplot2::scale_color_manual(name="Condition",
                                        guide="legend",
                                        values=color_list) +
            ggplot2::scale_fill_manual(name="Condition",
                                       guide="legend",
                                       values=color_list) +
            ggplot2::scale_shape_manual(name="Batch",
                                        labels=levels(as.factor(pca_data[["batch"]])),
                                        guide=ggplot2::guide_legend(override.aes=list(size=plot_size, fill="grey")),
                                        values=21:25)
    } else if (is.null(size_column) & num_batches > 5) {
        pca_plot <- pca_plot +
            ggplot2::geom_point(size=plot_size,
                                aes_string(shape="as.factor(batches)",
                                           colour="as.factor(condition)")) +
            ggplot2::scale_color_manual(name="Condition",
                                        guide="legend",
                                        values=color_list) +
            ggplot2::scale_shape_manual(name="Batch",
                                        labels=levels(as.factor(pca_data[["batch"]])),
                                        guide=ggplot2::guide_legend(overwrite.aes=list(size=plot_size)),
                                        values=1:num_batches)
    } else if (!is.null(size_column) & num_batches <= 5) {
        ## This will require the 6 steps above and one more
        maxsize <- max(pca_data[[size_column]])
        pca_plot <- pca_plot +
            ggplot2::geom_point(ggplot2::aes_string(shape="as.factor(batches)",
                                                    size=size_column,
                                                    colour="as.factor(condition)",
                                                    fill="as.factor(condition)")) +
            ggplot2::geom_point(colour="black", show.legend=FALSE,
                                aes_string(shape="as.factor(batches)",
                                           size=size_column,
                                           fill="as.factor(condition)")) +
            ggplot2::scale_color_manual(name="Condition",
                                        guide="legend",
                                        values=color_list) +
            ggplot2::scale_fill_manual(name="Condition",
                                       guide=ggplot2::guide_legend(override.aes=list(size=plot_size)),
                                       values=color_list) +
            ggplot2::scale_shape_manual(name="Batch",
                                        labels=levels(as.factor(pca_data[["batch"]])),
                                        guide=ggplot2::guide_legend(override.aes=list(size=plot_size, fill="grey")),
                                        values=21:25) +
            ggplot2::scale_size(range=c(2,7))
    } else if (!is.null(size_column) & num_batches > 5) {
        maxsize <- max(pca_data[[size_column]])
        pca_plot <- pca_plot +
            ggplot2::geom_point(ggplot2::aes_string(shape="as.factor(batches)",
                                                    colour="pca_data[['condition']]",
                                                    size=size_column)) +
            ggplot2::scale_shape_manual(name="Batch",
                                        labels=levels(as.factor(pca_data[["batch"]])),
                                        guide=ggplot2::guide_legend(overwrite.aes=list(size=plot_size)),
                                        values=1:num_batches) +
            ggplot2::scale_color_identity(name="Condition",
                                          guide="legend",
                                          values=color_list) +
            ggplot2::scale_size(range=c(2,7))
    } else {
        stop("This should be an impossible state.")
    }

    if (!is.null(variances)) {
        x_var_num <- as.numeric(gsub("PC", "", first))
        y_var_num <- as.numeric(gsub("PC", "", second))
        x_label <- paste("PC", x_var_num, ": ", variances[[x_var_num]], "%  variance", sep="")
        y_label <- paste("PC", y_var_num, ": ", variances[[y_var_num]], "%  variance", sep="")
        pca_plot <- pca_plot + ggplot2::xlab(x_label) + ggplot2::ylab(y_label)
    }

    if (is.null(plot_labels)) {
        plot_labels <- "repel"
    }
    if (plot_labels == FALSE) {
        message("Not putting labels on the plot.")
    } else if (plot_labels == "normal") {
        pca_plot <- pca_plot +
            ggplot2::geom_text(ggplot2::aes_string(x="PC1", y="PC2", label="labels",
                                                   angle=45, size=4, vjust=2))
    } else if (plot_labels == "repel") {
        pca_plot <- pca_plot +
            ggrepel::geom_text_repel(ggplot2::aes_string(label="labels"),
                                     size=5, box.padding=ggplot2::unit(0.5, 'lines'),
                                     point.padding=ggplot2::unit(1.6, 'lines'),
                                     arrow=ggplot2::arrow(length=ggplot2::unit(0.01, 'npc')))
    } else if (plot_labels == "dlsmart") {
        pca_plot <- pca_plot +
            directlabels::geom_dl(ggplot2::aes_string(label="labels"), method="smart.grid")
    } else {
        pca_plot <- pca_plot +
            directlabels::geom_dl(ggplot2::aes_string(label="labels"), method="first.qp")
    }

    return(pca_plot)
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
    plotted_us <- abs(plotted_us[,c(1,2,3)])
    plotted_u1s <- plotted_us[order(plotted_us[,1], decreasing=TRUE),]
    plotted_u2s <- plotted_us[order(plotted_us[,2], decreasing=TRUE),]
    plotted_u3s <- plotted_us[order(plotted_us[,3], decreasing=TRUE),]
    ## allS <- BiocGenerics::rank(allS, ties.method = "random")
    ## plotted_us$rank = rank(plotted_us[,1], ties.method="random")
    plotted_u1s <- cbind(plotted_u1s, rev(rank(plotted_u1s[,1], ties.method="random")))
    plotted_u1s <- plotted_u1s[,c(1,4)]
    colnames(plotted_u1s) <- c("PC1","rank")
    plotted_u1s <- data.frame(plotted_u1s)
    plotted_u1s$ID <- as.character(rownames(plotted_u1s))
    plotted_u2s <- cbind(plotted_u2s, rev(rank(plotted_u2s[, 2], ties.method="random")))
    plotted_u2s <- plotted_u2s[, c(2,4)]
    colnames(plotted_u2s) <- c("PC2","rank")
    plotted_u2s <- data.frame(plotted_u2s)
    plotted_u2s$ID <- as.character(rownames(plotted_u2s))
    plotted_u3s <- cbind(plotted_u3s, rev(rank(plotted_u3s[, 3], ties.method="random")))
    plotted_u3s <- plotted_u3s[, c(3,4)]
    colnames(plotted_u3s) <- c("PC3","rank")
    plotted_u3s <- data.frame(plotted_u3s)
    plotted_u3s$ID <- as.character(rownames(plotted_u3s))
    plotted_us <- merge(plotted_u1s, plotted_u2s, by.x="rank", by.y="rank")
    plotted_us <- merge(plotted_us, plotted_u3s, by.x="rank", by.y="rank")
    colnames(plotted_us) <- c("rank","PC1","ID1","PC2","ID2","PC3","ID3")
    rm(plotted_u1s)
    rm(plotted_u2s)
    rm(plotted_u3s)
    ## top_threePC = head(plotted_us, n=20)
    plotted_us <- plotted_us[, c("PC1","PC2","PC3")]
    plotted_us[, "ID"] <- rownames(plotted_us)
    message("More shallow curves in these plots suggest more genes in this principle component.")
    plot(plotted_us)
    u_plot <- grDevices::recordPlot()
    return(u_plot)
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
#'   useful information like the conditions, batches, number of cells, whatever...
#' @param expt_factors   a character list of experimental conditions to query
#'   for R^2 against the fast.svd of the data.
#' @param num_components   a number of principle components to compare the design factors against.
#'   If left null, it will query the same number of components as factors asked for.
#' @param plot_pcas   plot the set of PCA plots for every pair of PCs queried.
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
#' @seealso \code{\link[corpcor]{fast.svd}}, \code{\link[stats]{lm}}
#' @examples
#' \dontrun{
#'  pca_info = pca_information(exprs(some_expt$expressionset), some_design, "all")
#'  pca_info
#' }
#' @export
pca_information <- function(expt_data, expt_design=NULL, expt_factors=c("condition","batch"),
                            num_components=NULL, plot_pcas=FALSE) {
    colors_chosen <- expt_data[["colors"]]
    data_class <- class(expt_data)[1]
    if (data_class == "expt") {
        expt_design <- expt_data[["design"]]
        expt_data <- Biobase::exprs(expt_data[["expressionset"]])
    } else if (data_class == "ExpressionSet") {
        expt_data <- Biobase::exprs(expt_data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        expt_data <- as.matrix(expt_data)
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    expt_data <- as.matrix(expt_data)
    expt_means <- rowMeans(expt_data)
    decomposed <- corpcor::fast.svd(expt_data - expt_means)
    positives <- decomposed[["d"]]
    u <- decomposed[["u"]]
    v <- decomposed[["v"]]
    ## A neat idea from Kwame, rank order plot the U's in the svd version of:
    ## [Covariates] = [U][diagonal][V] for a given PC (usually/always PC1)
    ## The idea being: the resulting decreasing line should be either a slow even
    ## decrease if many genes are contributing to the given component
    ## Conversely, that line should drop suddenly if dominated by one/few genes.
    rownames(u) <- rownames(expt_data)
    rownames(v) <- colnames(expt_data)
    u_plot <- u_plot(u)
    component_variance <- round((positives^2) / sum(positives^2) * 100, 3)
    cumulative_pc_variance <- cumsum(component_variance)
    ## Include in this table the fstatistic and pvalue described in rnaseq_bma.rmd
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
    for (component in expt_factors) {
        message(component)
        ##comp <- factor(as.character(expt_design[, component]), exclude=FALSE)
        comp <- expt_design[[component]]
        if (is.null(comp)) {
            message(paste0("The given component is not in the design: ", comp))
            next
        }
        column <- apply(v, 2, factor_rsquared, fact=comp)
        component_rsquared_table[[component]] <- column
    }

    pca_variance <- round((positives ^ 2) / sum(positives ^2) * 100, 2)
    xl <- sprintf("PC1: %.2f%% variance", pca_variance[1])
    ##print(xl)
    yl <- sprintf("PC2: %.2f%% variance", pca_variance[2])
    ##print(yl)

    pca_data <- data.frame(
        "sampleid" = rownames(expt_design),
        "labels" = rownames(expt_design),
        "condition" = as.character(expt_design[["condition"]]),
        "batch" = as.character(expt_design[["batch"]]),
        "batch_int" = as.integer(as.factor(expt_design[["batch"]])),
        "colors" = colors_chosen)
    pc_df <- data.frame(
        "sampleid" = rownames(expt_design))
    rownames(pc_df) <- rownames(expt_design)

    if (is.null(num_components)) {
        num_components <- length(expt_factors)
    }
    max_components <- ncol(v)
    if (max_components < num_components) {
        message(paste0("The u and v components of SVD have only ", max_components, " columns, but the list of factors is ", num_components, " long."))
        message(paste0("Therefore, only searching for ", max_components, " PCs."))
        num_components <- max_components
    }
    for (pc in 1:num_components) {
        name <- paste("PC", pc, sep="")
        pca_data[[name]] <- v[, pc] ## note you _must_ not shortcut this with [[pc]]
        pc_df[[name]] <- v[, pc]
    }
    pc_df <- pc_df[, -1, drop=FALSE]
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
                    tmp_plot <- try(print(plot_pcs(pca_data,
                                                   design=expt_design,
                                                   variances=pca_variance,
                                                   first=name,
                                                   second=second_name)))
                    pca_plots[[list_name]] <- tmp_plot
                }
            }
        }
    }
    factor_df <- data.frame(
        "sampleid" = rownames(expt_design))
    rownames(factor_df) <- rownames(expt_design)
    for (fact in expt_factors) {
        if (!is.null(expt_design[[fact]])) {
            factor_df[[fact]] <- as.numeric(as.factor(as.character(expt_design[, fact])))
        } else {
            message(paste0("The column ", fact, " seems to be missing from the design."))
            message(paste0("The available columns are: ", toString(colnames(expt_design)), "."))
        }
    }
    factor_df <- factor_df[, -1, drop=FALSE]
    ## fit_one = data.frame()
    ## fit_two = data.frame()
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
            ##4.  The resulting p-value is 1 - pf(Fstat, (n-(#levels in the factor)), (n-1))  ## n is the number of samples in the fit
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
            tryCatch(
                {
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
    rownames(cor_df) <- colnames(factor_df)
    colnames(cor_df) <- colnames(pc_df)
    colnames(anova_sums) <- colnames(pc_df)
    colnames(anova_f) <- colnames(pc_df)
    colnames(anova_p) <- colnames(pc_df)
    colnames(anova_rss) <- colnames(pc_df)
    colnames(anova_fstats) <- colnames(pc_df)
    cor_df <- as.matrix(cor_df)
    silly_colors <- grDevices::colorRampPalette(c("purple","black","yellow"))(100)
    cor_df <- cor_df[complete.cases(cor_df), ]
    pc_factor_corheat <- heatmap.3(cor_df, scale="none", trace="none", linewidth=0.5,
                                   keysize=2, margins=c(8,8), col=silly_colors,
                                   dendrogram="none", Rowv=FALSE, Colv=FALSE,
                                   main="cor(factor, PC)")
    pc_factor_corheat <- grDevices::recordPlot()
    anova_f_colors <- grDevices::colorRampPalette(c("blue","black","red"))(100)
    anova_f_heat <- heatmap.3(as.matrix(anova_f), scale="none", trace="none",
                              linewidth=0.5, keysize=2, margins=c(8,8), col=anova_f_colors,
                              dendrogram = "none", Rowv=FALSE, Colv=FALSE,
                              main="anova fstats for (factor, PC)")
    anova_f_heat <- grDevices::recordPlot()
    anova_fstat_colors <- grDevices::colorRampPalette(c("blue","white","red"))(100)
    anova_fstat_heat <- heatmap.3(as.matrix(anova_fstats), scale="none", trace="none", linewidth=0.5,
                                  keysize=2, margins=c(8,8), col=anova_fstat_colors, dendrogram="none",
                                  Rowv=FALSE, Colv=FALSE, main="anova fstats for (factor, PC)")
    anova_fstat_heat <- grDevices::recordPlot()
    neglog_p <- -1 * log(as.matrix(anova_p) + 1)
    anova_neglogp_colors <- grDevices::colorRampPalette(c("blue","white","red"))(100)
    anova_neglogp_heat <- heatmap.3(as.matrix(neglog_p), scale="none", trace="none", linewidth=0.5,
                                    keysize=2, margins=c(8,8), col=anova_f_colors, dendrogram="none",
                                    Rowv=FALSE, Colv=FALSE, main="-log(anova_p values)")
    anova_neglogp_heat <- grDevices::recordPlot()
    ## Another option: -log10 p-value of the ftest for this heatmap.
    ## covariate vs PC score
    ## Analagously: boxplot(PCn ~ batch)
    pca_list <- list(
        pc1_trend=u_plot,svd_d=positives, svd_u=u, svd_v=v, rsquared_table=component_rsquared_table,
        pca_variance=pca_variance, pca_data=pca_data, anova_fstats=anova_fstats, anova_sums=anova_sums,
        anova_f=anova_f, anova_p=anova_p, pca_cor=cor_df, cor_heatmap=pc_factor_corheat,
        anova_f_heatmap=anova_f_heat, anova_fstat_heatmap=anova_fstat_heat,
        anova_neglogp_heatmaph=anova_neglogp_heat, pca_plots=pca_plots
    )
    return(pca_list)
}

#' Get the highest/lowest scoring genes for every principle component.
#'
#' This function uses princomp to acquire a principle component biplot
#' for some data and extracts a dataframe of the top n genes for each
#' component by score.
#'
#' @param df   a dataframe of (pseudo)counts
#' @param conditions   a factor or character of conditions in the experiment.
#' @param batches   a factor or character of batches in the experiment.
#' @param n   the number of genes to extract.
#' @return a list including the princomp biplot, histogram, and tables
#' of top/bottom n scored genes with their scores by component.
#' @seealso \code{\link[stats]{princomp}}
#' @examples
#' \dontrun{
#'  information = pca_highscores(df=df, conditions=cond, batches=bat)
#'  information$pca_bitplot  ## oo pretty
#' }
#' @export
pca_highscores <- function(df=NULL, conditions=NULL, batches=NULL, n=20) {
    ## Another method of using PCA
    ## cond = as.factor(as.numeric(conditions))
    ## batch = as.factor(as.numeric(batches))
    another_pca <- try(stats::princomp(x=df, cor=TRUE, scores=TRUE, formula=~0 + cond + batch))
    plot(another_pca)
    pca_hist <- grDevices::recordPlot()
    biplot(another_pca)
    pca_biplot <- grDevices::recordPlot()
    highest <- NULL
    lowest <- NULL
    for (pc in 1:length(colnames(another_pca$scores))) {
        tmphigh <- another_pca$scores
        tmplow <- another_pca$scores
        tmphigh <- tmphigh[order(tmphigh[, pc], decreasing=TRUE),]
        tmphigh <- head(tmphigh, n=20)
        tmplow <- tmplow[order(tmplow[, pc], decreasing=FALSE),]
        tmplow <- head(tmplow, n=20)
        high_column <- paste0(signif(tmphigh[, pc], 4), ":", rownames(tmphigh))
        low_column <- paste0(signif(tmplow[, pc], 4), ":", rownames(tmplow))
        highest <- cbind(highest, high_column)
        lowest <- cbind(lowest, low_column)
    }
    colnames(highest) <- colnames(another_pca$scores)
    colnames(lowest) <- colnames(another_pca$scores)
    ret_list <- list(
        "pca_hist" = pca_hist,
        "pca_biplot" = pca_biplot,
        "highest" = highest,
        "lowest" = lowest)
    return(ret_list)
}

## EOF
