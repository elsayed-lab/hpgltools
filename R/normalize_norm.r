## Time-stamp: <Tue May 10 14:07:48 2016 Ashton Trey Belew (abelew@gmail.com)>

#'   Perform a simple normalization of a count table
#'
#' @param data A matrix of count data
#' @param design  A dataframe describing the experimental design
#' (conditions/batches/etc)
#' @param norm  A normalization to perform:
#' 'sf|quant|qsmooth|tmm|upperquartile|tmm|rle'
#' I keep wishy-washing on whether design is a required argument.
#' @return dataframe of normalized(counts)
#' @seealso \pkg{edgeR} \pkg{limma} \pkg{DESeq2}
#' @examples
#' \dontrun{
#' norm_table = normalize_counts(count_table, design=design, norm='qsmooth')
#' }
#' @export
normalize_counts <- function(data, design=NULL, norm="raw", ...) {
    arglist <- list(...)
    ## Note that checkUsage flagged my 'libsize = ' calls
    ## I set norm_libsize at the bottom of the function
    ## but perhaps instead I should be using these libsizes?
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data[["design"]]
        count_table <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == "ExpressionSet") {
        count_table <- Biobase::exprs(data)
    } else if (data_class == "list") {
        count_table <- data[["count_table"]]
        if (is.null(data)) {
            stop("The list provided contains no count_table.")
        }
    } else if (data_class == "matrix" | data_class == "data.frame") {
        ## some functions prefer matrix, so I am keeping this explicit
        count_table <- as.data.frame(data)
    } else {
        stop(paste0("You provided a class of type: ", data_class, ".
This works with: expt, ExpressionSet, data.frame, and matrices.
"))
    }
    if (norm == "sf") {
        ## Size-factored normalization is a part of DESeq
        factors <- DESeq2::estimateSizeFactorsForMatrix(count_table)
        num_rows <- dim(count_table)[1]
        sf_counts <- count_table / do.call(rbind, rep(list(factors), num_rows))
        ##sf_counts = counts / (libsizes * factors)
        count_table <- as.matrix(sf_counts)
        norm_performed <- "sf"
    } else if (norm == "sf2") {
        original_cols <- colnames(count_table)
        conds <- design[["conditions"]]
        if (is.null(conds)) {
            conds <- original_cols
        }
        cds <- DESeq::newCountDataSet(count_table, conditions=conds)
        factors <- BiocGenerics::estimateSizeFactors(cds)
        count_table <- BiocGenerics::counts(factors, normalized=TRUE)
        norm_performed <- "sf2"
    } else if (norm == "vsd") {
        original_cols <- colnames(count_table)
        conds <- design[["conditions"]]
        if (is.null(conds)) {
            conds <- original_cols
        }
        cds <- DESeq::newCountDataSet(count_table, conditions=conds)
        factors <- BiocGenerics::estimateSizeFactors(cds)
        dispersions <- BiocGenerics::estimateDispersions(factors, method='blind')
        count_table <- DESeq::getVarianceStabilizedData(dispersions)
        norm_performed <- "vsd"
    } else if (norm == "quant") {
        # Quantile normalization (Bolstad et al., 2003)
        count_rownames <- rownames(count_table)
        count_colnames <- colnames(count_table)
        count_table <- preprocessCore::normalize.quantiles(as.matrix(count_table), copy=TRUE)
        rownames(count_table) <- count_rownames
        colnames(count_table) <- count_colnames
        norm_performed <- "quant"
    } else if (norm == "qsmooth") {
        count_table <- qsmooth::qsmooth(count_table, groups=design$condition, plot=TRUE)
        norm_performed <- "qsmooth"
    } else if (norm == "qshrink") {
        count_table <- hpgl_qshrink(exprs=count_table, groups=design$condition,
                                    plot=TRUE)
        norm_performed <- "qshrink"
    } else if (norm == "qshrink_median") {
        count_table <- hpgl_qshrink(exprs=count_table, groups=design$condition,
                                    plot=TRUE, refType="median",
                                    groupLoc="median", window=50)
        norm_performed <- "qshrink_median"
    } else if (norm == "tmm") {
        ## TMM normalization is documented in edgeR
        ## Set up the edgeR data structure
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="TMM")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "tmm"
    } else if (norm == "upperquartile") {
        ## Get the tmm normalization factors
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="upperquartile")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "upperquartile"
    } else if (norm == "rle") {
        ## Get the tmm normalization factors
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="RLE")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "rle"
    } else {
        message("Did not recognize the normalization, leaving the table alone.
  Recognized normalizations include: 'qsmooth', 'sf', 'sf2', 'vsd', 'quant',
  'tmm', 'qsmooth_median', 'upperquartile', and 'rle.'
")
        count_table <- as.matrix(count_table)
    }
    norm_libsize <- colSums(count_table)
    norm_counts <- list(count_table=count_table, libsize=norm_libsize, norm_performed=norm_performed)
    return(norm_counts)
}

#' A hacked copy of Kwame's qsmooth/qstats code
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data count table to modify
#' @param groups factor of the experimental conditions
#' @param refType method for grouping conditions
#' @param groupLoc method for grouping groups
#' @param window a window, for looking!
#' @param groupCol column to define conditions
#' @param plot plot the quantiles?
#' @param ... more options
#' @return data a new data frame of normalized counts
#' @seealso \pkg{qsmooth}
#' @examples
#' \dontrun{
#' df <- hpgl_qshrink(data)
#' }
#' @export
hpgl_qshrink <- function(data=NULL, groups=NULL, refType="mean",
                        groupLoc="mean", window=99,
                        groupCol=NULL, plot=TRUE, ...) {
    data <- as.matrix(data)
    if (is.null(groups)) {
        message("Groups were not provided.  Performing a simple quantile")
        message("normalization. This is probably not what you actually want!")
        count_rownames <- rownames(data)
        count_colnames <- colnames(data)
        normExprs <- preprocessCore::normalize.quantiles(as.matrix(data), copy=TRUE)
        rownames(normExprs) <- count_rownames
        colnames(normExprs) <- count_colnames
        return(normExprs)
    }
    res <- hpgl_qstats(Biobase::exprs, groups, refType=refType,
                       groupLoc=groupLoc, window=window)
    QBETAS <- res$QBETAS
    Qref <- res$Qref
    X <- res$model
    w <- res$smoothWeights
    wQBETAS <- QBETAS * (1 - w)
    wQBETAS <- X %*% t(wQBETAS)
    wQref <- Qref * w
    wQref <- matrix(rep(1, nrow(X)), ncol=1) %*% t(wQref)
    normExprs <- t(wQBETAS + wQref)
    RANKS <- t(matrixStats::colRanks(data, ties.method="average"))
    for (k in 1:ncol(normExprs)) {
        x <- normExprs[, k]
        normExprs[, k] <- x[RANKS[, k]]
    }

    aveTies = function (ranks, y) {
        tab = table(ranks)
        sel = tab > 1
        if (sum(sel) != 0) {
            ties = as.numeric(names(tab[sel]))
            for (k in ties) {
                sel = ranks==k
                y[sel] = mean(y[sel])
            }
        }
        y
    }
    normExprs <- aveTies(RANKS, normExprs)

    rownames(normExprs) <- rownames(data)
    colnames(normExprs) <- colnames(data)
    if (plot) {
        oldpar <- par(mar=c(4, 4, 1.5, 0.5))
        lq <- length(Qref)
        u <- (1:lq - 0.5)/lq
        if (length(u) > 10000) {
            sel <- sample(1:lq, 10000)
            plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
                 xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
            ## plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
            ##      xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1))
        } else {
            plot(u, w, pch=".", main="Quantile reference weights",
                 xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
            ## plot(u, w, pch=".", main="Quantile reference weights",
            ##      xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1))
        }
        abline(h=0.5, v=0.5, col="red", lty=2)
        newpar <- par(oldpar)
    }
    return(normExprs)
}

#' A hacked copy of Kwame's qsmooth/qstats code
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data the initial count data
#' @param groups the experimental conditions as a factor
#' @param refType  (or median) the method to separate groups
#' @param groupLoc   I don't remember
#' @param window window for basking
#' @return new data
#' @examples
#' \dontrun{
#' qstatted <- hpgl_qstats(data, conditions)
#' }
#' @export
hpgl_qstats <- function (data, groups, refType="mean",
                         groupLoc="mean", window=99) {
    ## require.auto("matrixStats")
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
            message(paste0("There was only replicate of type: ", g))
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
            warning(paste0("There were 0 of type: ", g))
        }
    }
    colnames(QBETAS) <- uGroups
    colnames(SIGMA) <- uGroups
    if (groupLoc == "mean") {
        TAU <- matrixStats::rowVars(QBETAS)
        SIGMA <- rowMeans(SIGMA)
    } else { ## median
        TAU <- matrixStats::rowMads(QBETAS)^2
        SIGMA <- matrixStats::rowMedians(SIGMA)
    }
    roughWeights <- SIGMA/(SIGMA + TAU)
    roughWeights[is.nan(roughWeights)] = 0 ## is this backward?
    roughWeights[SIGMA < 10^(-6) & TAU < 10^(-6)] = 1
    smoothWeights <- stats::runmed(roughWeights, k=window, endrule="constant")
    qstats_model <- model.matrix(~0 + factor(groups, levels=uGroups))
    qstats_result <- list(Q=Q, Qref=Qref, QBETAS=QBETAS, TAU=TAU,
                          SIGMA=SIGMA, roughWeights=roughWeights,
                          smoothWeights=smoothWeights, model=qstats_model)
    return(qstats_result)
}

