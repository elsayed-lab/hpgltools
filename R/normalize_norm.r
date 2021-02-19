## normalize_norm.r: Invoke the various normalization methods on expressionset
## data.  For my purposes I keep this separate from the various log
## transformations and scaling/conversion methods (cpm/rpkm) because ... well
## because they just feel different to me.

#' Perform a simple normalization of a count table.
#'
#' This provides shortcut interfaces for normalization functions from
#' deseq2/edger and friends.
#'
#' @param data Matrix of count data.
#' @param design Dataframe describing the experimental
#'  design. (conditions/batches/etc)
#' @param norm Normalization to perform:
#'  'sf|quant|qsmooth|tmm|upperquartile|tmm|rle' I keep wishy-washing on
#'  whether design is a required argument.
#' @param ... More arguments might be necessary.
#' @return Dataframe of normalized(counts)
#' @seealso \pkg{edgeR} \pkg{limma} \pkg{DESeq2}
#' @examples
#' \dontrun{
#'  norm_table = normalize_counts(count_table, design = design, norm='qsmooth')
#' }
#' @export
normalize_counts <- function(data, design = NULL, method = "raw", ...) {
  arglist <- list(...)
  if (!is.null(arglist[["norm"]])) {
    method <- arglist[["norm"]]
  }
  ## Note that checkUsage flagged my 'libsize = ' calls
  ## I set norm_libsize at the bottom of the function
  ## but perhaps instead I should be using these libsizes?
  data_class <- class(data)[1]
  norm_performed <- "raw"
  if (data_class == "expt") {
    design <- data[["design"]]
    count_table <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    count_table <- exprs(data)
  } else if (data_class == "list") {
    count_table <- data[["count_table"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table.")
    }
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit
    count_table <- as.data.frame(data)
  } else {
    stop(glue("You provided a class of type: {data_class}.
This works with: expt, ExpressionSet, data.frame, and matrices.
"))
  }

  switchret <- switch(
    method,
    "qshrink" = {
      count_table <- hpgl_qshrink(exprs = count_table, groups = design[["condition"]],
                                  plot = TRUE)
      norm_performed <- "qshrink"
    },
    "qshrink_median" = {
      count_table <- hpgl_qshrink(exprs = count_table, groups = design[["condition"]],
                                  plot = TRUE, refType = "median",
                                  groupLoc = "median", window = 50)
      norm_performed <- "qshrink_median"
    },
    "quant" = {
      ## Quantile normalization (Bolstad et al., 2003)
      count_rownames <- rownames(count_table)
      count_colnames <- colnames(count_table)
      count_table <- preprocessCore::normalize.quantiles(
                                       as.matrix(count_table))
      rownames(count_table) <- count_rownames
      colnames(count_table) <- count_colnames
      norm_performed <- "quant"
    },
    "quant_robust" = {
      count_rownames <- rownames(count_table)
      count_colnames <- colnames(count_table)
      ## 20181210 -- this gives a pthread_create() error 22.
      ##count_table <- preprocessCore::normalize.quantiles(
      ##                                 as.matrix(count_table), copy = TRUE)
      count_table <- preprocessCore::normalize.quantiles.robust(
                                       as.matrix(count_table))
      rownames(count_table) <- count_rownames
      colnames(count_table) <- count_colnames
      norm_performed <- "quant_robust"
    },
    "quantile" = {
      ## Quantile normalization (Bolstad et al., 2003)
      count_rownames <- rownames(count_table)
      count_colnames <- colnames(count_table)
      count_table <- preprocessCore::normalize.quantiles(
                                       as.matrix(count_table), copy = TRUE)
      rownames(count_table) <- count_rownames
      colnames(count_table) <- count_colnames
      norm_performed <- "quant"
    },
    "rle" = {
      ## Get the tmm normalization factors
      count_table <- edgeR::DGEList(counts = count_table)
      norms <- edgeR::calcNormFactors(count_table, method = "RLE")
      ## libsizes = count_table$samples$lib.size
      factors <- norms[["samples"]][["norm.factors"]]
      counts <- norms[["counts"]]
      tmm_counts <- counts / factors
      count_table <- as.matrix(tmm_counts)
      norm_performed <- "rle"
    },
    "sf" = {
      original_cols <- colnames(count_table)
      conds <- design[["conditions"]]
      if (is.null(conds)) {
        conds <- original_cols
      }
      ## cds <- DESeq::newCountDataSet(count_table, conditions = conds)
      factors <- BiocGenerics::estimateSizeFactors(count_table)
      count_table <- BiocGenerics::counts(factors, normalized = TRUE)
      norm_performed <- "sf"
    },
    "sf2" = {
      ## Size-factored normalization is a part of DESeq
      factors <- DESeq2::estimateSizeFactorsForMatrix(count_table)
      num_rows <- dim(count_table)[1]
      sf_counts <- count_table / do.call(rbind, rep(list(factors), num_rows))
      ##sf_counts = counts / (libsizes * factors)
      count_table <- as.matrix(sf_counts)
      norm_performed <- "sf2"
    },
    "tmm" = {
      ## TMM normalization is documented in edgeR
      ## Set up the edgeR data structure
      count_table <- edgeR::DGEList(counts = count_table)
      norms <- edgeR::calcNormFactors(count_table, method = "TMM")
      ## libsizes = count_table$samples$lib.size
      factors <- norms[["samples"]][["norm.factors"]]
      counts <- norms[["counts"]]
      tmm_counts <- counts / factors
      count_table <- as.matrix(tmm_counts)
      norm_performed <- "tmm"
    },
    "upperquartile" = {
      ## Get the tmm normalization factors
      count_table <- edgeR::DGEList(counts = count_table)
      norms <- edgeR::calcNormFactors(count_table, method = "upperquartile")
      ## libsizes = count_table$samples$lib.size
      factors <- norms[["samples"]][["norm.factors"]]
      counts <- norms[["counts"]]
      tmm_counts <- counts / factors
      count_table <- as.matrix(tmm_counts)
      norm_performed <- "upperquartile"
    },
    "vsd" = {
      original_cols <- colnames(count_table)
      conds <- design[["conditions"]]
      if (is.null(conds)) {
        conds <- design[["condition"]]
        if (is.null(conds)) {
          conds <- original_cols
        }
      }
      fit_type <- "parametric"
      if (!is.null(arglist[["fit_type"]])) {
        fit_type <- arglist[["fit_type"]]
      }
      tt <- sm(requireNamespace("locfit"))
      tt <- sm(requireNamespace("DESeq2"))
      cds <- DESeq2::DESeqDataSetFromMatrix(
                       countData = count_table, colData = design, design=~condition)
      cds <- DESeq2::estimateSizeFactors(cds)
      cds <- DESeq2::estimateDispersions(cds, fitType = fit_type)
      count_table <- DESeq2::getVarianceStabilizedData(cds)
      norm_performed <- "vsd"
    },
    {
      message("Did not recognize the normalization, leaving the table alone.
  Recognized normalizations include: 'qsmooth', 'sf', 'sf2', 'vsd', 'quant',
  'tmm', 'qsmooth_median', 'upperquartile', and 'rle.'")
      count_table <- as.matrix(count_table)
    }
  ) ## End of the switch statement.
  norm_libsize <- colSums(count_table, na.rm = TRUE)
  norm_counts <- list(count_table = count_table, libsize = norm_libsize,
                      norm_performed = norm_performed)
  return(norm_counts)
}

#' A hacked copy of Kwame's qsmooth/qstats code.
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data Count table to modify
#' @param groups Factor of the experimental conditions
#' @param refType Method for grouping conditions
#' @param groupLoc Method for grouping groups
#' @param window Window, for looking!
#' @param groupCol Column to define conditions
#' @param plot Plot the quantiles?
#' @param ... More options
#' @return New data frame of normalized counts
#' @seealso \pkg{qsmooth}
#' @examples
#' \dontrun{
#'  df <- hpgl_qshrink(data)
#' }
#' @export
hpgl_qshrink <- function(data = NULL, groups = NULL, refType = "mean",
                         groupLoc = "mean", window = 99,
                         groupCol = NULL, plot = TRUE, ...) {
  data <- as.matrix(data)
  if (is.null(groups)) {
    message("Groups were not provided.  Performing a simple quantile")
    message("normalization. This is probably not what you actually want!")
    count_rownames <- rownames(data)
    count_colnames <- colnames(data)
    normExprs <- preprocessCore::normalize.quantiles(as.matrix(data), copy = TRUE)
    rownames(normExprs) <- count_rownames
    colnames(normExprs) <- count_colnames
    return(normExprs)
  }
  res <- hpgl_qstats(exprs(data), groups, refType = refType,
                     groupLoc = groupLoc, window = window)
  QBETAS <- res[["QBETAS"]]
  Qref <- res[["Qref"]]
  X <- res[["model"]]
  w <- res[["smoothWeights"]]
  wQBETAS <- QBETAS * (1 - w)
  wQBETAS <- X %*% t(wQBETAS)
  wQref <- Qref * w
  wQref <- matrix(rep(1, nrow(X)), ncol = 1) %*% t(wQref)
  normExprs <- t(wQBETAS + wQref)
  RANKS <- t(matrixStats::colRanks(data, ties.method = "average"))
  for (k in 1:ncol(normExprs)) {
    x <- normExprs[, k]
    normExprs[, k] <- x[RANKS[, k]]
  }

  aveTies <- function (ranks, y) {
    tab <- table(ranks)
    sel <- tab > 1
    if (sum(sel) != 0) {
      ties <- as.numeric(names(tab[sel]))
      for (k in ties) {
        sel <- ranks==k
        y[sel] <- mean(y[sel])
      }
    }
    y
  }
  normExprs <- aveTies(RANKS, normExprs)

  rownames(normExprs) <- rownames(data)
  colnames(normExprs) <- colnames(data)
  if (plot) {
    oldpar <- par(mar = c(4, 4, 1.5, 0.5))
    lq <- length(Qref)
    u <- (1:lq - 0.5)/lq
    if (length(u) > 10000) {
      sel <- sample(1:lq, 10000)
      plot(u[sel], w[sel], pch = ".", main = "Quantile reference weights",
           xlab = "u (norm. gene ranks)", ylab = "Weight", ylim = c(0, 1), ...)
    } else {
      plot(u, w, pch = ".", main = "Quantile reference weights",
           xlab = "u (norm. gene ranks)", ylab = "Weight", ylim = c(0, 1), ...)
    }
    abline(h = 0.5, v = 0.5, col = "red", lty = 2)
    newpar <- par(oldpar)
  }
  return(normExprs)
}

#' A hacked copy of Kwame's qsmooth/qstats code.
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data Initial count data
#' @param groups Experimental conditions as a factor.
#' @param refType Method to separate groups, mean or median.
#' @param groupLoc I don't remember what this is for.
#' @param window Window for basking!
#' @return Some new data.
#' @seealso \pkg{matrixStats}
#' @examples
#' \dontrun{
#'  qstatted <- hpgl_qstats(data, conditions)
#' }
#' @export
hpgl_qstats <- function(data, groups, refType = "mean",
                        groupLoc = "mean", window = 99) {
  Q <- apply(data, 2, sort)
  if (refType == "median") {
    Qref <- matrixStats::rowMedians(Q)
  }
  if (refType == "mean") {
    Qref <- rowMeans(Q)
  }

  QBETAS <- c()
  SIGMA <- c()
  uGroups <- unique(groups)
  for (g in uGroups) {
    index <- (g == groups)
    if (sum(index) == 1) {
      message("There was only replicate of type: ", g)
      message("This will likely do terrible things to qsmooth.")
      QBETAS <- cbind(QBETAS, Q[, index])
      SIGMA <- cbind(SIGMA, 0)
    } else if (sum(index) > 1) {
      if (groupLoc == "mean") {
        QBETAS <- cbind(QBETAS, rowMeans(Q[, index]))
        SIGMA <- cbind(SIGMA, matrixStats::rowVars(Q[, g == groups]))
      } else if (groupLoc == "median") {
        QBETAS <- cbind(QBETAS, matrixStats::rowMedians(Q[, index]))
        SIGMA <- cbind(SIGMA, (matrixStats::rowMads(Q[, g == groups]))^2)
      }
    } else {
      warning(glue("There were 0 of type: {g}"))
    }
  }
  colnames(QBETAS) <- uGroups
  colnames(SIGMA) <- uGroups
  if (groupLoc == "mean") {
    TAU <- matrixStats::rowVars(QBETAS)
    SIGMA <- rowMeans(SIGMA)
  } else {
    ## median
    TAU <- matrixStats::rowMads(QBETAS)^2
    SIGMA <- matrixStats::rowMedians(SIGMA)
  }
  roughWeights <- SIGMA / (SIGMA + TAU)
  roughWeights[is.nan(roughWeights)] <- 0 ## is this backward?
  roughWeights[SIGMA < 10 ^ -6 & TAU < 10 ^ -6] <- 1
  smoothWeights <- stats::runmed(roughWeights, k = window, endrule = "constant")
  qstats_model <- model.matrix(~0 + factor(groups, levels = uGroups))
  qstats_result <- list(Q = Q, Qref = Qref, QBETAS = QBETAS, TAU = TAU,
                        SIGMA = SIGMA, roughWeights = roughWeights,
                        smoothWeights = smoothWeights, model = qstats_model)
  return(qstats_result)
}

## EOF
