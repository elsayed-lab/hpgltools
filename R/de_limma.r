## Time-stamp: <Thu Apr 14 10:27:00 2016 Ashton Trey Belew (abelew@gmail.com)>

#'   Plot out 2 coefficients with respect to one another from limma
#'
#' It can be nice to see a plot of two coefficients from a limma comparison with respect to one another
#' This hopefully makes that easy.
#'
#' @param output the set of pairwise comparisons provided by limma_pairwise()
#' @param toptable  use this to get up/downs and color them on the scatter plot
#' @param x   the name or number of the first coefficient column to extract
#' @param y   the name or number of the second coefficient column to extract
#' @param gvis_filename   A filename for plotting gvis interactive graphs of the data.
#' @param gvis_trendline   add a trendline to the gvis plot?
#' @param z   how far from the median to color the plot red and green
#' @param tooltip_data   a dataframe of gene annotations to be used in the gvis plot
#' @param base_url   a basename for gvis plots
#' @param up_color   color for the ups
#' @param down_color   color for the downs
#' @param ...  more parameters to make you happy
#' @return a ggplot2 plot showing the relationship between the two coefficients
#' @seealso \link{hpgl_linear_scatter} \link{limma_pairwise}
#' @examples
#' \dontrun{
#'  pretty = coefficient_scatter(limma_data, x="wt", y="mut")
#' }
#' @export
limma_coefficient_scatter <- function(output, toptable=NULL, x=1, y=2, ##gvis_filename="limma_scatter.html",
                                      gvis_filename=NULL, gvis_trendline=TRUE, z=1.5,
                                      tooltip_data=NULL, base_url=NULL,
                                      up_color="#7B9F35", down_color="#DD0000", ...) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    arglist <- list(...)
    qlimit <- 0.1
    if (!is.null(arglist$qlimit)) {
        qlimit <- arglist$qlimit
    }
    message("This can do comparisons among the following columns in the limma result:")
    thenames <- colnames(output$pairwise_comparisons$coefficients)
    message(thenames)
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
    coefficients <- output$pairwise_comparisons$coefficients
    coefficients <- coefficients[, c(x,y)]
    maxvalue <- max(coefficients) + 1
    plot <- hpgl_linear_scatter(df=coefficients, loess=TRUE, gvis_filename=gvis_filename,
                                gvis_trendline=gvis_trendline, first=xname, second=yname,
                                tooltip_data=tooltip_data, base_url=base_url, pretty_colors=FALSE)
    plot$scatter <- plot$scatter +
        ggplot2::scale_x_continuous(limits=c(0, maxvalue)) +
        ggplot2::scale_y_continuous(limits=c(0, maxvalue))
    if (!is.null(toptable)) {
        theplot <- plot$scatter + ggplot2::theme_bw()
        sig <- limma_subset(toptable, z=z)
        sigup <- sig$up
        sigdown <- sig$down
        ## sigup <- subset(sigup, qvalue < 0.1)
        sigup <- sigup[sigup$qvalue <= qlimit, ]
        ## sigdown <- subset(sigdown, qvalue < 0.1)
        sigdown <- sigdown[sigdown$qvalue <= qlimit, ]
        up_index <- rownames(coefficients) %in% rownames(sigup)
        down_index <- rownames(coefficients) %in% rownames(sigdown)
        up_df <- as.data.frame(coefficients[up_index, ])
        down_df <- as.data.frame(coefficients[down_index, ])
        colnames(up_df) <- c("first","second")
        colnames(down_df) <- c("first","second")
        theplot <- theplot +
            ggplot2::geom_point(data=up_df, colour=up_color) +
            ggplot2::geom_point(data=down_df, colour=down_color)
        plot$scatter <- theplot
    }
    plot$df <- coefficients
    return(plot)
}

#'   A slight modification of limma's voom() function.
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
#' @param model   an experimental model defining batches/conditions/etc
#' @param libsize   the size of the libraries (usually provided by
#' edgeR).
#' @param stupid   whether or not to cheat when the resulting matrix is not solvable.
#' @param logged   whether the input data is known to be logged.
#' @param converted   whether the input data is known to be cpm converted.
#' @return an EList containing the following information:
#'   E = The normalized data
#'   weights = The weights of said data
#'   design = The resulting design
#'   lib.size = The size in pseudocounts of the library
#'   plot = A ggplot of the mean/variance trend with a blue loess fit and red trend fit
#' @seealso \link[limma]{voom} \link[cbcbSEQ]{voomMod} \link[limma]{lmFit}
#' @examples
#' \dontrun{
#'  funkytown = hpgl_voom(samples, model)
#' }
#' @export
hpgl_voom <- function(dataframe, model=NULL, libsize=NULL, stupid=FALSE, logged=FALSE, converted=FALSE) {
    out <- list()
    if (is.null(libsize)) {
        libsize <- colSums(dataframe, na.rm=TRUE)
    }
    if (converted == 'cpm') {
        converted <- TRUE
    }
    if (!isTRUE(converted)) {
        message("The voom input was not cpm, converting now.")
        posed <- t(dataframe + 0.5)
        dataframe <- t(posed / (libsize + 1) * 1e+06)
        ##y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1000000)) ## from voom()
    }
    if (logged == 'log2') {
        logged <- TRUE
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
        dataframe <- log2(dataframe)
    }
    dataframe <- as.matrix(dataframe)

    if (is.null(model)) {
        model <- matrix(1, ncol(dataframe), 1)
        rownames(model) <- colnames(dataframe)
        colnames(model) <- "GrandMean"
    }
    linear_fit <- limma::lmFit(dataframe, model, method="ls")
    if (is.null(linear_fit$Amean)) {
        linear_fit$Amean <- rowMeans(dataframe, na.rm=TRUE)
    }
    sx <- linear_fit$Amean + mean(log2(libsize + 1)) - log2(1e+06)
    sy <- sqrt(linear_fit$sigma)
    if (is.na(sum(sy))) { ## 1 replicate
        return(NULL)
    }
    allzero <- rowSums(dataframe) == 0
    stupid_NAs <- is.na(sx)
    sx <- sx[!stupid_NAs]
    stupid_NAs <- is.na(sy)
    sy <- sy[!stupid_NAs]
    if (any(allzero == TRUE, na.rm=TRUE)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    fitted <- gplots::lowess(sx, sy, f=0.5)
    f <- stats::approxfun(fitted, rule=2)
    mean_var_df <- data.frame(mean=sx, var=sy)
    mean_var_plot <- ggplot2::ggplot(mean_var_df, ggplot2::aes_string(x="mean", y="var")) +
        ggplot2::geom_point() +
        ggplot2::xlab("Log2(count size + 0.5)") +
        ggplot2::ylab("Square root of the standard deviation.") +
        ## stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE, show_guide=FALSE) +
        ggplot2::stat_density2d(geom="tile", ggplot2::aes_string(fill="..density..^0.25"),
                                contour=FALSE, show.legend=FALSE) +
        ggplot2::scale_fill_gradientn(colours=grDevices::colorRampPalette(c("white","black"))(256)) +
        ggplot2::geom_smooth(method="loess") +
        ggplot2::stat_function(fun=f, colour="red") +
        ggplot2::theme(legend.position="none")
    if (is.null(linear_fit$rank)) {
        message("Some samples cannot be balanced across the experimental design.")
        if (isTRUE(stupid)) {
            ## I think this is telling me I have confounded data, and so
            ## for those replicates I will have no usable coefficients, so
            ## I say set them to 1 and leave them alone.
            linear_fit$coefficients[is.na(linear_fit$coef)] <- 1
            fitted.values <- linear_fit$coef %*% t(linear_fit$design)
        }
    } else if (linear_fit$rank < ncol(linear_fit$design)) {
        j <- linear_fit$pivot[1:linear_fit$rank]
        fitted.values <- linear_fit$coef[, j, drop=FALSE] %*% t(linear_fit$design[, j, drop=FALSE])
    } else {
        fitted.values <- linear_fit$coef %*% t(linear_fit$design)
    }
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (libsize + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1 / f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    rownames(w) <- rownames(dataframe)
    colnames(w) <- colnames(dataframe)
    out$E <- dataframe
    out$weights <- w
    out$design <- model
    out$lib.size <- libsize
    out$plot <- mean_var_plot
    new("EList", out)
}

#' Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using voom/limma.
#'
#' @param input  a dataframe/vector or expt class containing count tables, normalization state, etc.
#' @param conditions   a factor of conditions in the experiment
#' @param batches   a factor of batches in the experiment
#' @param model_cond   include condition in the model?
#' @param model_batch   include batch in the model? This is hopefully TRUE.
#' @param model_intercept   perform a cell-means or intercept model?
#' A little more difficult for me to understand.  I have tested and get the same answer either way.
#' @param extra_contrasts   some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param alt_model   a separate model matrix instead of the normal condition/batch.
#' @param libsize   I've recently figured out that libsize is far more important than I previously realized.  Play with it here.
#' @param annot_df   data frame for annotations
#' @param ... The elipsis parameter is fed to write_limma() at the end.
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
#' @seealso \link{write_limma}
#' @examples
#' \dontrun{
#' pretend = balanced_pairwise(data, conditions, batches)
#' }
#' @export
limma_pairwise <- function(input, conditions=NULL, batches=NULL, model_cond=TRUE,
                           model_batch=TRUE, model_intercept=FALSE, extra_contrasts=NULL,
                           alt_model=NULL, libsize=NULL, annot_df=NULL, ...) {
    arglist <- list(...)
    message("Starting limma pairwise comparison.")
    input_class <- class(input)[1]
    if (input_class == 'expt') {
        conditions <- input[["conditions"]]
        batches <- input[["batches"]]
        data <- Biobase::exprs(input[["expressionset"]])
        if (is.null(libsize)) {
            message("libsize was not specified, this parameter has profound effects on limma's result.")
            if (!is.null(input[["best_libsize"]])) {
                message("Using the libsize from expt$best_libsize.")
                ## libsize = expt$norm_libsize
                libsize <- input[["best_libsize"]]
            } else if (!is.null(input[["normalized"]][["intermediate_counts"]][["normalization"]][["libsize"]])) {
                libsize <- colSums(data)
            } else {
                message("Using the libsize from expt$normalized$intermediate_counts$normalization$libsize")
                libsize <- input[["normalized"]][["intermediate_counts"]][["normalization"]][["libsize"]]
            }
        } else {
            message("libsize was specified.  This parameter has profound effects on limma's result.")
        }
    } else {  ## Not an expt class, data frame or matrix
        data <- as.data.frame(input)
    }
    if (is.null(libsize)) {
        libsize <- colSums(data)
    }
    condition_table <- table(conditions)
    batch_table <- table(batches)
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    ## Make a model matrix which will have one entry for each of these condition/batches
    cond_model <- stats::model.matrix(~ 0 + conditions)  ## I am not putting a try() on this, because if it fails, then we are effed.
    batch_model <- try(stats::model.matrix(~ 0 + batches), silent=TRUE)
    condbatch_model <- try(stats::model.matrix(~ 0 + conditions + batches), silent=TRUE)
    batch_int_model <- try(stats::model.matrix(~ batches), silent=TRUE)
    cond_int_model <- try(stats::model.matrix(~ conditions), silent=TRUE)
    condbatch_int_model <- try(stats::model.matrix(~ conditions + batches), silent=TRUE)
    fun_model <- NULL
    fun_int_model <- NULL
    if (isTRUE(model_cond) & isTRUE(model_batch)) {
        if (class(condbatch_model) == 'try-error') {
            message("The condition+batch model failed.  Does your experimental design support both condition and batch?")
            message("Using only a conditional model.")
            fun_model <- cond_model
            fun_int_model <- cond_int_model
        } else {
            fun_model <- condbatch_model
            fun_int_model <- condbatch_int_model
        }
    } else if (class(model_batch) == 'matrix' | class(model_batch) == 'numeric') {
        message("Limma: Including multiple sv batch estimates from sva/ruv/pca in the limma model.")
        fun_model <- stats::model.matrix(~ 0 + conditions + model_batch)
        fun_int_model <- stats::model.matrix(~ conditions + model_batch)
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
    tmpnames <- gsub("-", "", tmpnames)
    tmpnames <- gsub("+", "", tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    colnames(fun_model) <- tmpnames
    fun_voom <- NULL
    message("Limma step 1/6: choosing model.")
    ## voom() it, taking into account whether the data has been log2 transformed.
    loggedp <- input[["state"]][["transform"]]
    if (is.null(loggedp)) {
        message("I don't know if this data is logged, testing if it is integer.")
        if (is.integer(data)) {
            loggedp <- FALSE
        } else {
            loggedp <- TRUE
        }
    } else {
        if (loggedp == "raw") {
            loggedp <- FALSE
        } else {
            loggedp <- TRUE
        }
    }
    convertedp = input[["state"]][["conversion"]]
    if (is.null(convertedp)) {
        message("I cannot determine if this data has been converted, assuming no.")
        convertedp <- FALSE
    } else {
        if (convertedp == "raw") {
            convertedp <- FALSE
        } else {
            convertedp <- TRUE
        }
    }
    ##fun_voom = voom(data, fun_model)
    ##fun_voom = hpgl_voom(data, fun_model, libsize=libsize)
    ##fun_voom = voomMod(data, fun_model, lib.size=libsize)
    message("Limma step 2/6: running voom")
    fun_voom <- hpgl_voom(data, fun_model, libsize=libsize, logged=loggedp, converted=convertedp)
    one_replicate <- FALSE
    if (is.null(fun_voom)) {
        message("voom returned null, I am not sure what will happen.")
        one_replicate <- TRUE
        fun_voom <- data
        fun_design <- NULL
    } else {
        fun_design <- fun_voom[["design"]]
    }

    ## Extract the design created by voom()
    ## This is interesting because each column of the design will have a prefix string 'macb' before the
    ## condition/batch string, so for the case of clbr_tryp_batch_C it will look like: macbclbr_tryp_batch_C
    ## This will be important in 17 lines from now.
    ## Do the lmFit() using this model
    message("limma step 3/6: running lmFit")
    fun_fit <- limma::lmFit(fun_voom, fun_model)
    ##fun_fit = lmFit(fun_voom)
    ## The following three tables are used to quantify the relative contribution of each batch to the sample condition.
    message("Limma step 4/6: making and fitting contrasts.")
    if (isTRUE(model_intercept)) {
        contrasts <- "intercept"
        identities <- NULL
        contrast_string <- NULL
        all_pairwise <- NULL
        all_pairwise_fits <- fun_fit
    } else {
        contrasts <- make_pairwise_contrasts(fun_model, conditions, extra_contrasts=extra_contrasts)
        all_pairwise_contrasts <- contrasts[["all_pairwise_contrasts"]]
        identities <- contrasts[["identities"]]
        contrast_string <- contrasts[["contrast_string"]]
        all_pairwise <- contrasts[["all_pairwise"]]
        ## Once all that is done, perform the fit
        ## This will first provide the relative abundances of each condition
        ## followed by the set of all pairwise comparisons.
        all_pairwise_fits <- limma::contrasts.fit(fun_fit, all_pairwise_contrasts)
    }
    all_tables <- NULL
    message("Limma step 5/6: Running eBayes and topTable.")
    if (isTRUE(one_replicate)) {
        all_pairwise_comparisons <- all_pairwise_fits[["coefficients"]]
    } else {
        all_pairwise_comparisons <- limma::eBayes(all_pairwise_fits)
        all_tables <- try(limma::topTable(all_pairwise_comparisons, number=nrow(all_pairwise_comparisons)))
    }
    message("Limma step 6/6: Writing limma outputs.")
    if (isTRUE(model_intercept)) {
        limma_result <- all_tables
    } else {
        limma_result <- try(write_limma(all_pairwise_comparisons, excel=FALSE))
    }
    result <- list(
        "all_pairwise" = all_pairwise,
        "all_tables" = limma_result,
        "batches" = batches,
        "batches_table" = batch_table,
        "conditions" = conditions,
        "conditions_table" = condition_table,
        "contrast_string" = contrast_string,
        "fit" = fun_fit,
        "identities" = identities,
        "input_data" = data,
        "model" = fun_model,
        "pairwise_fits" = all_pairwise_fits,
        "pairwise_comparisons" = all_pairwise_comparisons,
        "single_table" = all_tables,
        "voom_design" = fun_design,
        "voom_result" = fun_voom)
    return(result)
}

#' Plot arbitrary data from limma
#'
#' @param all_pairwise_result  the result from calling balanced_pairwise()
#' @param first_table   the first table from all_pairwise_result$limma_result to look at (may be a name or number)
#' @param first_column   the name of the column to plot from the first table
#' @param second_table   the second table inside all_pairwise_result$limma_result (name or number)
#' @param second_column  a column to compare against
#' @param type A type of scatter plot (linear model, distance, vanilla)
#' @param ... so that you may feed it the gvis/tooltip information to make clicky graphs if so desired.
#' @return a hpgl_linear_scatter() set of plots comparing the chosen columns
#' If you forget to specify tables to compare, it will try the first vs the second.
#' @seealso \link{hpgl_linear_scatter} \link[limma]{topTable}
#' @examples
#' \dontrun{
#' compare_logFC = limma_scatter(all_pairwise, first_table="wild_type", second_column="mutant",
#'                               first_table="AveExpr", second_column="AveExpr")
#' compare_B = limma_scatter(all_pairwise, first_column="B", second_column="B")
#' }
#' @export
limma_scatter <- function(all_pairwise_result, first_table=1, first_column="logFC",
                         second_table=2, second_column="logFC", type="linear_scatter", ...) {
    tables <- all_pairwise_result$all_tables
    if (is.numeric(first_table)) {
        x_name <- paste(names(tables)[first_table], first_column, sep=":")
    }
    if (is.numeric(second_table)) {
        y_name <- paste(names(tables)[second_table], second_column, sep=":")
    }

    ## This section is a little bit paranoid
    ## I want to make absolutely certain that I am adding only the
    ## two columns I care about and that nothing gets reordered
    ## As a result I am explicitly pulling a single column, setting
    ## the names, then pulling the second column, then cbind()ing them.
    x_name <- paste(first_table, first_column, sep=":")
    y_name <- paste(second_table, second_column, sep=":")
    df <- data.frame(x=tables[[first_table]][[first_column]])
    rownames(df) <- rownames(tables[[first_table]])
    second_column_list <- tables[[second_table]][[second_column]]
    names(second_column_list) <- rownames(tables[[second_table]])
    df <- cbind(df, second_column_list)
    colnames(df) <- c(x_name, y_name)
    plots <- NULL
    if (type == "linear_scatter") {
        plots <- hpgl_linear_scatter(df, loess=TRUE, ...)
    } else if (type == "dist_scatter") {
        plots <- hpgl_dist_scatter(df, ...)
    } else {
        plots <- hpgl_scatter(df, ...)
    }
    plots[['dataframe']] <- df
    return(plots)
}

#'  A quick and dirty way to pull the top/bottom genes from toptable()
#'
#' If neither n nor z is provided, it assumes you want 1.5 z-scores from the median.
#'
#' @param table  the original data from limma
#' @param n   a number of genes to keep
#' @param z   a number of z-scores from the mean
#' @return a dataframe subset from toptable
#' @seealso \pkg{limma}
#' @examples
#' \dontrun{
#'  subset = limma_subset(df, n=400)
#'  subset = limma_subset(df, z=1.5)
#' }
#' @export
limma_subset <- function(table, n=NULL, z=NULL) {
    if (is.null(n) & is.null(z)) {
        z <- 1.5
    }
    if (is.null(n)) {
        out_summary <- summary(table$logFC)
        out_mad <- stats::mad(table$logFC, na.rm=TRUE)
        up_median_dist <- out_summary["Median"] + (out_mad * z)
        down_median_dist <- out_summary["Median"] - (out_mad * z)

        up_genes <- table[ which(table$logFC >= up_median_dist), ]
        ## up_genes = subset(table, logFC >= up_median_dist)
        down_genes <- table[ which(table$logFC <= down_median_dist), ]
        ## down_genes = subset(table, logFC <= down_median_dist)
    } else if (is.null(z)) {
        upranked <- table[ order(table$logFC, decreasing=TRUE),]
        up_genes <- head(upranked, n=n)
        down_genes <- tail(upranked, n=n)
    }
    ret_list <- list(up=up_genes, down=down_genes)
    return(ret_list)
}

#'  Perform a simple experimental/control comparison
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
#'  Currently this assumes that a variant of toptable was used which
#'  gives adjusted p-values.  This is not always the case and I should
#'  check for that, but I have not yet.
#'
#' @param subset  an experimental subset with two conditions to compare.
#' @param workbook   an excel workbook to which to write.
#' @param sheet   an excel worksheet to which to write.
#' @param basename   a url to which to send click evens in clicky volcano/ma plots.
#' @param batch   whether or not to include batch in limma's model.
#' @param combat   whether or not to use combatMod().
#' @param combat_noscale   whether or not to include combat_noscale (makes combat a little less heavy-handed).
#' @param pvalue_cutoff   p-value definition of 'significant.'
#' @param logfc_cutoff   fold-change cutoff of significance. 0.6 on the low end and therefore 1.6 on the high.
#' @param tooltip_data   text descriptions of genes if one wants google graphs.
#' @param verbose   be verbose?
#' @param ... more parameters!
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
#' @seealso \link{hpgl_gvis_ma_plot} \link[limma]{toptable}
#' \link[limma]{voom} \link[cbcbSEQ]{voomMod} \link{hpgl_voom}
#' \link[limma]{lmFit} \link[limma]{makeContrasts} \link[limma]{contrasts.fit}
#' @examples
#' \dontrun{
#'  model = model.matrix(~ 0 + subset$conditions)
#'  simple_comparison(subset, model)
#' }
#' @export
simple_comparison <- function(subset, workbook="simple_comparison.xls", sheet="simple_comparison",
                              basename=NA, batch=TRUE, combat=FALSE, combat_noscale=TRUE,
                              pvalue_cutoff=0.05, logfc_cutoff=0.6, tooltip_data=NULL,
                              verbose=FALSE, ...) {
    condition_model <- stats::model.matrix(~ 0 + subset$condition)
    if (length(levels(subset$batch)) == 1) {
        message("There is only one batch! I can only include condition in the model.")
        condbatch_model <- stats::model.matrix(~ 0 + subset$condition)
    } else {
        condbatch_model <- stats::model.matrix(~ 0 + subset$condition + subset$batch)
    }
    if (isTRUE(batch)) {
        model <- condbatch_model
    } else {
        model <- condition_model
    }
    expt_data <- as.data.frame(Biobase::exprs(subset$expressionset))
    if (combat) {
#        expt_data = ComBat(expt_data, subset$batches, condition_model)
        expt_data <- cbcbSEQ::combatMod(expt_data, subset$batches, subset$conditions)
    }
    expt_voom <- hpgltools::hpgl_voom(expt_data, model, libsize=subset$original_libsize,
                                      logged=subset$transform, converted=subset$convert)
    lf <- limma::lmFit(expt_voom)
    colnames(lf$coefficients)
    coefficient_scatter <- hpgltools::hpgl_linear_scatter(lf$coefficients)
    colnames(lf$design)[1] <- "changed"
    colnames(lf$coefficients)[1] <- "changed"
    colnames(lf$design)[2] <- "control"
    colnames(lf$coefficients)[2] <- "control"
    ## Now make sure there are no weird characters in the column names...
    if (length(colnames(lf$design)) >= 3) {
        for (counter in 3:length(colnames(lf$design))) {
            oldname <- colnames(lf$design)[counter]
            newname <- gsub("\\$","_", oldname, perl=TRUE)
            colnames(lf$design)[counter] <- newname
            colnames(lf$coefficients)[counter] <- newname
        }
    }
    contrast_matrix <- limma::makeContrasts(changed_v_control="changed-control", levels=lf$design)
    ## contrast_matrix = limma::makeContrasts(changed_v_control=changed-control, levels=lf$design)
    cond_contrasts <- limma::contrasts.fit(lf, contrast_matrix)
    hist_df <- data.frame(values=cond_contrasts$coefficients)
    contrast_histogram <- hpgl_histogram(hist_df)
    hist_df <- data.frame(values=cond_contrasts$Amean)
    amean_histogram <- hpgltools::hpgl_histogram(hist_df, fillcolor="pink", color="red")
    coef_amean_cor <- stats::cor.test(cond_contrasts$coefficients, cond_contrasts$Amean, exact=FALSE)
    cond_comparison <- limma::eBayes(cond_contrasts)
    hist_df <- data.frame(values=cond_comparison$p.value)
    pvalue_histogram <- hpgl_histogram(hist_df, fillcolor="lightblue", color="blue")
    cond_table <- limma::topTable(cond_comparison, number=nrow(expt_voom$E),
                                  coef="changed_v_control", sort.by="logFC")
    if (!is.na(basename)) {
        vol_gvis_filename <- paste(basename, "volplot.html", sep="_")
        a_volcano_plot <- hpgl_volcano_plot(cond_table, gvis_filename=vol_gvis_filename,
                                            tooltip_data=tooltip_data)
    } else {
        a_volcano_plot <- hpgl_volcano_plot(cond_table)
    }
    if (!is.na(basename)) {
        ma_gvis_filename <- paste(basename, "maplot.html", sep="_")
        an_ma_plot <- hpgl_ma_plot(expt_voom$E, cond_table, gvis_filename=ma_gvis_filename,
                                   tooltip_data=tooltip_data)
    } else {
        an_ma_plot <- hpgltools::hpgl_ma_plot(expt_voom$E, cond_table)
    }
    write_xls(cond_table, sheet, file=workbook, rowname="row.names")
    ## upsignificant_table = subset(cond_table, logFC >=  logfc_cutoff)
    upsignificant_table <- cond_table[ which(cond_table$logFC >= logfc_cutoff), ]
    ## downsignificant_table = subset(cond_table, logFC <= (-1 * logfc_cutoff))
    downsignificant_table <- cond_table[ which(cond_table$logFC <= (-1 * logfc_cutoff)), ]
    ## psignificant_table = subset(cond_table, adj.P.Val <= pvalue_cutoff)
    ## psignificant_table = subset(cond_table, P.Value <= pvalue_cutoff)
    psignificant_table <- cond_table[ which(cond_table$P.Value <= pvalue_cutoff), ]

    if (isTRUE(verbose)) {
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
    return_info <- list(
        amean_histogram=amean_histogram, coef_amean_cor=coef_amean_cor,
        coefficient_scatter=coefficient_scatter$scatter,
        coefficient_x=coefficient_scatter$x_histogram,
        coefficient_y=coefficient_scatter$y_histogram,
        coefficient_both=coefficient_scatter$both_histogram,
        coefficient_lm=coefficient_scatter$lm_model,
        coefficient_lmsummary=coefficient_scatter$lm_summary,
        coefficient_weights=coefficient_scatter$lm_weights,
        comparisons=cond_comparison, contrasts=cond_contrasts,
        contrast_histogram=contrast_histogram,
        downsignificant=downsignificant_table,
        fit=lf, ma_plot=an_ma_plot, psignificant=psignificant_table,
        pvalue_histogram=pvalue_histogram, table=cond_table,
        upsignificant=upsignificant_table,
        volcano_plot=a_volcano_plot, voom_data=expt_voom,
        voom_plot=expt_voom$plot)
    return(return_info)
}

#' Writes out the results of a limma search using toptable()
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data  the output from eBayes()
#' @param adjust   the pvalue adjustment chosen.
#' @param n   the number of entries to report, 0 says do them all.
#' @param coef   which coefficients/contrasts to report, NULL says do them all.
#' @param workbook   an excel filename into which to write the data
#' @param excel   write an excel workbook?
#' @param csv   write out csv files of the tables?
#' @param annot_df   an optional data frame including annotation information to include with the tables.
#' @return a list of data frames comprising the toptable output for each coefficient,
#'    I also added a qvalue entry to these toptable() outputs.
#' @seealso \link[limma]{toptable} \link{write_xls}
#' @examples
#' \dontrun{
#'  finished_comparison = eBayes(limma_output)
#'  data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
#' }
#' @export
write_limma <- function(data, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls",
                       excel=FALSE, csv=FALSE, annot_df=NULL) {
    testdir <- dirname(workbook)
    if (n == 0) {
        n <- dim(data$coefficients)[1]
    }
    if (is.null(coef)) {
        coef <- colnames(data$contrasts)
    } else {
        coef <- as.character(coef)
    }
    return_data <- list()
    end <- length(coef)
    for (c in 1:end) {
        comparison <- coef[c]
        message(paste0("limma step 6/6: ", c, "/", end, ": Printing table: ", comparison, "."))
        data_table <- limma::topTable(data, adjust=adjust, n=n, coef=comparison)
        ## Reformat the numbers so they are not so obnoxious
##        data_table$logFC <- refnum(data_table$logFC, sci=FALSE)
##        data_table$AveExpr <- refnum(data_table$AveExpr, sci=FALSE)
##        data_table$t <- refnum(data_table$t, sci=FALSE)
##        data_table$P.Value <- refnum(data_table$P.Value)
##        data_table$adj.P.Val <- refnum(data_table$adj.P.Val)
##        data_table$B <- refnum(data_table$B, sci=FALSE)
        data_table$logFC <- signif(x=as.numeric(data_table$logFC), digits=4)
        data_table$AveExpr <- signif(x=as.numeric(data_table$AveExpr), digits=4)
        data_table$t <- signif(x=as.numeric(data_table$t), digits=4)
        data_table$P.Value <- signif(x=as.numeric(data_table$P.Value), digits=4)
        data_table$adj.P.Val <- signif(x=as.numeric(data_table$adj.P.Val), digits=4)
        data_table$B <- signif(x=as.numeric(data_table$B), digits=4)
        data_table$qvalue <- tryCatch(
        {
            ##as.numeric(format(signif(
            ##    suppressWarnings(qvalue::qvalue(
            ##        as.numeric(data_table$P.Value), robust=TRUE))$qvalues, 4),
            ##    scientific=TRUE))
            ttmp <- as.numeric(data_table$P.Value)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
##            ttmp <- signif(ttmp, 4)
##            ttmp <- format(ttmp, scientific=TRUE)
##            ttmp
        },
        error=function(cond) {
            message(paste("The qvalue estimation failed for ", comparison, ".", sep=""))
            return(1)
        },
        ##warning=function(cond) {
        ##    message("There was a warning?")
        ##    message(cond)
        ##    return(1)
        ##},
        finally={
        })
        if (!is.null(annot_df)) {
            data_table <- merge(data_table, annot_df, by.x="row.names", by.y="row.names")
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
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwritefile=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename <- gsub(".xls$", "", workbook)
            csv_filename <- paste0(csv_filename, "_", comparison, ".csv")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] <- data_table
    }
    return(return_data)
}

## EOF
