#' all_pairwise(): Wrap up limma/DESeq2/EdgeR pairwise analyses in one call.
#'
#' @param expt an expt class containing count tables, normalization state, etc.
#' @param conditions a factor of conditions in the experiment
#' @param batches a factor of batches in the experiment
#' @param extra_contrasts some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_minus_b_ctrla = (C-B)-A, e_minus_d_ctrla = (E-D)-A,
#'  de_minus_cb = (E-D)-(C-B),"
#' @param model_cond Include condition in the model?  This should pretty much always be true.
#' @param model_batch Include batch in the model? FALSE by default, but hopefully true often.
#' @param model_intercept Perform a cell-means or intercept model?  FALSE by default because I understand subtraction math better.
#'   But I have tested and get the same answer either way.
#' @param libsize I've recently figured out that libsize is far more important than I previously realized.  Play with it here.
#' @param ... The elipsis parameter is fed to write_limma() at the end.
#'
#' @return A list of limma, deseq, edger results.
#' @export
#' @examples
#' ## finished_comparison = eBayes(limma_output)
#' ## data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
all_pairwise = function(expt=NULL, data=NULL, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, model_intercept=FALSE, extra_contrasts=NULL, alt_model=NULL, libsize=NULL) {
    
    limma_result = limma_pairwise(expt=expt, data=data, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize)
    deseq_result = deseq2_pairwise(expt=expt, data=data, conditions=conditions, batches=batches) ## The rest of the arguments should be added back sooner than later.
    edger_result = edger_pairwise(expt=expt, data=data, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize)

    result_comparison = compare_tables(limma=limma_result, deseq=deseq_result, edger=edger_result)
    ret = list(limma=limma_result, deseq=deseq_result, edger=edger_result, comparison=result_comparison)
    return(ret)
}

#' write_limma(): Writes out the results of a limma search using toptable()
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data The output from eBayes()
#' @param adjust The pvalue adjustment chosen (fdr by default)
#' @param n The number of entries to report, defaults to 0, which says do them all
#' @param coef which coefficients/contrasts to report, NULL says do them all
#' @param workbook an excel filename into which to write the data, used for csv files too.
#' @param excel T/F whether or not to write an excel workbook (useful if they are too big)
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
write_limma = function(data=NULL, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls", excel=FALSE, csv=TRUE, annotation=NULL) {
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
    for (c in 1:length(coef)) {
        comparison = coef[c]
        message(paste(c, ": Printing table: ", comparison, sep=""))        
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
                message(paste("Creating directory: ", testdir, " for writing excel/csv data.", sep=""))
            }
        }
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwrite=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename = gsub(".xls$", "", workbook)            
            csv_filename = paste(csv_filename, "_", comparison, ".csv", sep="")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] = data_table
    }
    return(return_data)
}

#' write_edger(): Writes out the results of an EdgeR search using topTags()
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data The output from eBayes()
#' @param adjust The pvalue adjustment chosen (fdr by default)
#' @param n The number of entries to report, defaults to 0, which says do them all
#' @param coef which coefficients/contrasts to report, NULL says do them all
#' @param workbook an excel filename into which to write the data, used for csv files too.
#' @param excel T/F whether or not to write an excel workbook (useful if they are too big)
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
write_edger = function(data=NULL, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls", excel=FALSE, csv=TRUE, annotation=NULL) {
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
    for (c in 1:length(coef)) {
        comparison = coef[c]
        message(paste("Printing table: ", comparison, sep=""))        
        data_table = topTable(data, adjust=adjust, n=n, coef=comparison)

        data_table$qvalue = tryCatch(
            {
                format(signif(qvalue(data_table$P.Value, gui=FALSE)$qvalues, 4), scientific=TRUE)
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
        if (!is.null(annotation)) {
            data_table = merge(data_table, annotation, by.x="row.names", by.y="row.names")
            ###data_table = data_table[-1]
        }
        ## This write_xls runs out of memory annoyingly often
        if (isTRUE(excel) | isTRUE(csv)) {
            if (!file.exists(testdir)) {
                dir.create(testdir)
                message(paste("Creating directory: ", testdir, " for writing excel/csv data.", sep=""))
            }
        }
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwrite=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename = gsub(".xls$", "", workbook)            
            csv_filename = paste(csv_filename, "_", comparison, ".csv", sep="")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] = data_table
    }
    return(return_data)
}

#' write_deseq(): Writes out the results of a DESeq2 search using results()
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data The output from eBayes()
#' @param adjust The pvalue adjustment chosen (fdr by default)
#' @param n The number of entries to report, defaults to 0, which says do them all
#' @param coef which coefficients/contrasts to report, NULL says do them all
#' @param workbook an excel filename into which to write the data, used for csv files too.
#' @param excel T/F whether or not to write an excel workbook (useful if they are too big)
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
write_deseq2 = function(data=NULL, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls", excel=FALSE, csv=TRUE, annotation=NULL) {
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
    for (c in 1:length(coef)) {
        comparison = coef[c]
        message(paste("Printing table: ", comparison, sep=""))        
        data_table = topTable(data, adjust=adjust, n=n, coef=comparison)

        data_table$qvalue = tryCatch(
            {
                format(signif(qvalue(data_table$P.Value, gui=FALSE)$qvalues, 4), scientific=TRUE)
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
        if (!is.null(annotation)) {
            data_table = merge(data_table, annotation, by.x="row.names", by.y="row.names")
            ###data_table = data_table[-1]
        }
        ## This write_xls runs out of memory annoyingly often
        if (isTRUE(excel) | isTRUE(csv)) {
            if (!file.exists(testdir)) {
                dir.create(testdir)
                message(paste("Creating directory: ", testdir, " for writing excel/csv data.", sep=""))
            }
        }
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwrite=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename = gsub(".xls$", "", workbook)            
            csv_filename = paste(csv_filename, "_", comparison, ".csv", sep="")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] = data_table
    }
    return(return_data)
}

#' make_SVD() is a function scabbed from Hector and Kwame's cbcbSEQ
#' It just does fast.svd of a matrix against its rowMeans().
#'
#' @param data A data frame to decompose
#'
#' @return a list containing the s,v,u from fast.svd
#' @seealso \code{\link{fast.svd}}
#'
#' @export
#' @examples
#' ## svd = makeSVD(data)
makeSVD = function (x) {
    x = as.matrix(x)
    s = fast.svd(x - rowMeans(x))
    v = s$v
    rownames(v) = colnames(x)
    s = list(v=v, u=s$u, d=s$d)
    return(s)
}

## Some code to more strongly remove batch effects
remove_batch_effect = function(normalized_counts, model) {
    ## model = model.matrix(~ condition + batch)
    voomed = hpgl_voom(normalized_counts, model)
    voomed_fit = lmFit(voomed)
    modified_model = model
    modified_model = modified_model[,grep("batch", colnames(modified_model))] = 0 ## Drop batch from the model
    new_data = tcrossprod(voomed_fit$coefficient, modified_model) + residuals(voomed_fit, normalized_counts)
    return(new_data)
}

#' hpgl_voom():  A slight modification of limma's voom() function.
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
#' @param model an experimental model defining batches/conditions/etc
#' @param libsize the size of the libraries (usually provided by
#' edgeR).  NULL by default.
#' @param stupid An TRUE/FALSE of whether or not to cheat when the
#' resulting matrix is not solvable. FALSE by default.
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
hpgl_voom = function(dataframe, model, libsize=NULL, stupid=FALSE, logged=FALSE, converted=FALSE) {
    out = list()
    if (is.null(libsize)) {
        libsize = colSums(dataframe, na.rm=TRUE)
    }
    if (!isTRUE(converted)) {
        message("The voom input was not cpm, converting now.")
        posed = t(dataframe + 0.5)
        dataframe = t(posed/(libsize + 1) * 1e+06)
        ##y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1000000)) ## from voom()
    }
    if (!isTRUE(logged)) {
        message("The voom input was not log2, transforming now.")
        dataframe = log2(dataframe)
    }
    dataframe = as.matrix(dataframe)
    linear_fit = limma::lmFit(dataframe, model, method="ls")
    if (is.null(linear_fit$Amean)) {
        linear_fit$Amean = rowMeans(dataframe, na.rm=TRUE)
    }
    sx = linear_fit$Amean + mean(log2(libsize + 1)) - log2(1e+06)
    sy = sqrt(linear_fit$sigma)
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


#' limma_subset():  A quick and dirty way to pull the top/bottom genes from toptable()
#'
#' @param df The original data from limma
#' @param n A number of genes to keep
#' @param z A number of z-scores from the mean
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



#' limma_pairwise():  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using voom/limma.
#'
#' @param expt an expt class containing count tables, normalization state, etc.
#' @param conditions a factor of conditions in the experiment
#' @param batches a factor of batches in the experiment
#' @param extra_contrasts some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_minus_b_ctrla = (C-B)-A, e_minus_d_ctrla = (E-D)-A,
#'  de_minus_cb = (E-D)-(C-B),"
#' @param model_cond Include condition in the model?  This should pretty much always be true.
#' @param model_batch Include batch in the model? FALSE by default, but hopefully true often.
#' @param model_intercept Perform a cell-means or intercept model?  FALSE by default because I understand subtraction math better.
#'   But I have tested and get the same answer either way.
#' @param libsize I've recently figured out that libsize is far more important than I previously realized.  Play with it here.
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
limma_pairwise = function(expt=NULL, data=NULL, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, model_intercept=FALSE, extra_contrasts=NULL, alt_model=NULL, libsize=NULL) {
    if (is.null(expt) & is.null(data)) {
        stop("This requires either an expt or data/condition/batches")
    } else if (!is.null(expt)) {
        conditions = expt$conditions
        batches = expt$batches
        data = exprs(expt$expressionset)
        if (is.null(libsize)) {
            message("libsize was not specified, this parameter has profound effects on limma's result.")
            if (!is.null(expt$best_libsize)) {
                message("Using the libsize from expt$best_libsize.")
                ## libsize = expt$norm_libsize
                libsize = expt$best_libsize
            } else {
                message("Using the libsize from expt$normalized$normalized_counts.")                
                libsize = expt$normalized$normalized_counts$libsize
            }
        } else {
            message("libsize was specified.  This parameter has profound effects on limma's result.")
        }
    } ## End if expt is defined.
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
    logged = expt$transform
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
    converted = expt$convert
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
    
    ## Extract the design created by voom()
    ## This is interesting because each column of the design will have a prefix string 'macb' before the
    ## condition/batch string, so for the case of clbr_tryp_batch_C it will look like: macbclbr_tryp_batch_C
    ## This will be important in 17 lines from now.
    fun_design = fun_voom$design    
    ## Do the lmFit() using this model    
    ##fun_fit = lmFit(fun_voom, fun_model)
    fun_fit = lmFit(fun_voom)
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
    all_pairwise_comparisons = eBayes(all_pairwise_fits)
    all_tables = try(topTable(all_pairwise_comparisons, number=nrow(all_pairwise_comparisons)))
    if (isTRUE(model_intercept)) {
        limma_result = all_tables
    } else {
        limma_result = write_limma(all_pairwise_comparisons, excel=FALSE)
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

coefficient_scatter = function(limma_output, x=NULL, y=NULL, gvis_filename="limma_scatter.html", gvis_trendline=TRUE, tooltip_data=NULL) {
    ##  If taking a limma_pairwise output, then this lives in
    ##  output$pairwise_comparisons$coefficients
    print("This can do comparisons among the following columns in the limma result:")
    print(colnames(limma_output$pairwise_comparisons$coefficients))
    coefficients = limma_output$pairwise_comparisons$coefficients
    coefficients = coefficients[,c(x,y)]
    plot = hpgl_linear_scatter(df=coefficients, loess=TRUE, gvis_filename=gvis_filename, gvis_trendline=gvis_trendline)
    return(plot)
}

#' edger_pairwise():  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using EdgeR.
#'
#' @param expt a expt class containing data, normalization state, etc.
#' @param conditions a factor of conditions in the experiment
#' @param batches a factor of batches in the experiment
#' @param model_cond Include condition in the experimental model?  This is pretty much always true.
#' @param model_batch Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept Use cell means or intercept? (I default to the former, but they work out the same)
#' @param extra_contrasts some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_minus_b_ctrla = (C-B)-A, e_minus_d_ctrla = (E-D)-A,
#'  de_minus_cb = (E-D)-(C-B),"
#' @param ... The elipsis parameter is fed to write_limma() at the end.
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
edger_pairwise = function(expt=NULL, data=NULL, conditions=NULL, batches=NULL, model_cond=TRUE, model_batch=FALSE, model_intercept=FALSE, alt_model=NULL, extra_contrasts=NULL, ...) {
    if (is.null(expt) & is.null(data)) {
        stop("This requires either an expt or data+conditions+batches")
    } else if (!is.null(expt)) {
        conditions = expt$conditions
        batches = expt$batches
        data = as.data.frame(exprs(expt$expressionset))
        ## As I understand it, edgeR fits a binomial distribution
        ## and expects data as floating point counts,
        ##not a log2 transformation.
        if (!is.null(expt$transform)) {
            if (expt$transform == "log2") {
                ##data = (2^data) - 1
                data = expt$normalized$normalized_counts$count_table
            }
        }        
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
    message("Estimating the common dispersion.")
    disp_norm = estimateCommonDisp(norm)
    message("Estimating dispersion across genes.")
    tagdisp_norm = estimateTagwiseDisp(disp_norm)
    message("Estimating GLM Common dispersion.")
    glm_norm = estimateGLMCommonDisp(tagdisp_norm, fun_model)
    message("Estimating GLM Trended dispersion.")
    glm_trended = estimateGLMTrendedDisp(glm_norm, fun_model)
    message("Estimating GLM Tagged dispersion.")
    glm_tagged = estimateGLMTagwiseDisp(glm_trended, fun_model)
    cond_fit = edgeR::glmFit(glm_tagged, design=fun_model)

    apc = make_pairwise_contrasts(fun_model, conditions, do_identities=FALSE)
    ## This is pretty weird because glmLRT only seems to take up to 7 contrasts at a time...
    contrast_list = list()
    result_list = list()
    lrt_list = list()
    single_contrasts = list()
    sc = vector("list", length(apc$names))
    for (con in 1:length(apc$names)) {
        name = apc$names[[con]]
        message(paste0(con, ": Performing ", name, " contrast.")) ## correct
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
                format(signif(qvalue(res$PValue, robust=TRUE)$qvalues, 4), scientific=TRUE)
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

    final = list(
        contrasts=apc,
        lrt=lrt_list,
        contrast_list=contrast_list,
        all_tables=result_list)

    return(final)
}

#' deseq2_pairwise():  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions using DESeq2.
#'
#' @param expt a expt class containing data, normalization state, etc.
#' @param conditions a factor of conditions in the experiment
#' @param batches a factor of batches in the experiment
#' @param model_cond Include condition in the experimental model?  This is pretty much always true.
#' @param model_batch Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept Use cell means or intercept? (I default to the former, but they work out the same)
#' @param ... The elipsis parameter is fed to write_limma() at the end.
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
#' ## pretend = edger_pairwise(data, conditions, batches)
deseq2_pairwise = function(expt=NULL, data=NULL, conditions=NULL, batches=NULL) {
    if (is.null(expt) & is.null(data)) {
        stop("This requires either an expt or data+conditions+batches")
    } else if (!is.null(expt)) {
        conditions = expt$conditions
        batches = expt$batches
        data = as.data.frame(exprs(expt$expressionset))
        ## As I understand it, DESeq2 (and edgeR) fits a binomial distribution
        ## and expects data as floating point counts,
        ## not a log2 transformation.
        if (!is.null(expt$norm)) {
            if (expt$norm != "raw") {
                data = exprs(expt$original_expressionset)
            } else if (!is.null(expt$transform)) {
                if (expt$transform == "log2") {
                    ##data = (2^data) - 1
                    data = expt$normalized$normalized_counts$count_table
                }
            }
        }
    } ## expt is not null.
    condition_table = table(conditions)
    batch_table = table(batches)
    conditions = as.factor(conditions)
    batches = as.factor(batches)
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    cond_model = model.matrix(~ 0 + conditions)
    tmpnames = colnames(cond_model)
    tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    tmpnames = gsub("conditions", "", tmpnames)
    colnames(cond_model) = tmpnames

    ## An interesting note about the use of formulae in DESeq:
    ## "you should put the variable of interest at the end of the formula and make sure the control level is the first level."
    ## Thus, all these formulae should have condition(s) at the end.
    summarized = DESeqDataSetFromMatrix(countData=data, colData=pData(expt$expressionset), design=~0+condition)
    ## If making a model ~0 + condition -- then must set betaPrior=FALSE
    dataset = DESeqDataSet(se=summarized, design=~ 0 + condition)
    deseq_run = DESeq(dataset, betaPrior=FALSE)
    ## Set contrast= for each pairwise comparison here!

    denominators = list()
    numerators = list()
    result_list = list()
    result_mle_list = list()    
    condition_list = resultsNames(deseq_run)
    for (c in 1:(length(condition_list) - 1)) {
        denominator = names(condition_table[c])
        nextc = c + 1
        for (d in nextc:length(condition_list)) {
            numerator = names(condition_table[d])            
            result = as.data.frame(results(deseq_run, contrast=c("condition", numerator, denominator), format="DataFrame"))
            result = result[order(result$log2FoldChange),]
            colnames(result) = c("baseMean","logFC", "lfcSE","stat","P.Value","adj.P.Val")
            result[is.na(result$P.Value), "P.Value"] = 1 ## Some p-values come out as NA
            result[is.na(result$adj.P.Val), "adj.P.Val"] = 1 ## Some p-values come out as NA            
            result$qvalue = tryCatch(
                {
                    format(signif(qvalue(result$P.Value, robust=TRUE)$qvalues, 4), scientific=TRUE)
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
            result$P.Value = format(signif(result$P.Value, 4), scientific=TRUE)
            result$adj.P.Val = format(signif(result$adj.P.Val, 4), scientific=TRUE)
            result_name = paste0(numerator, "_minus_", denominator)
            denominators[[result_name]] = denominator
            numerators[[result_name]] = numerator
            result_list[[result_name]] = result
        }
    }
    
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
    
    ret_list = list(
        run=deseq_run,
        denominators=denominators,
        numerators=numerators,
        conditions=condition_list,
        all_tables=result_list        
    )
    return(ret_list)
}

compare_tables = function(limma=NULL, deseq=NULL, edger=NULL) {
    ## Fill each column/row of these with the correlation between tools for one contrast performed
    if (class(limma) == "list") { ## Then this was fed the raw output from limma_pairwise, lets assume the same is true for deseq/edger too and pull out the result tables.
        limma = limma$all_tables
        deseq = deseq$all_tables
        edger = edger$all_tables
    }
        
    len = length(names(deseq))
    limma_vs_edger = list()
    limma_vs_deseq = list()
    edger_vs_deseq = list()
    
    cc = 0
    last = length(names(deseq))
    for (comp in names(deseq)) {  ## assume all three have the same names() -- note that limma has more than the other two though
        cc = cc + 1
        message(paste0(cc, "/", len, ": Comparing analyses for: ", comp))
        l = data.frame(limma[[comp]])
        e = data.frame(edger[[comp]])
        d = data.frame(deseq[[comp]])
        le = merge(l, e, by.x="row.names", by.y="row.names")
        le = le[,c("logFC.x","logFC.y")]
        lec = cor.test(le[,1], le[,2])$estimate
        ld = merge(l, d, by.x="row.names", by.y="row.names")
        ld = ld[,c("logFC.x","logFC.y")]            
        ldc = cor.test(ld[,1], ld[,2])$estimate
        ed = merge(e, d, by.x="row.names", by.y="row.names")
        ed = ed[,c("logFC.x","logFC.y")]
        edc = cor.test(ed[,1], ed[,2])$estimate
        limma_vs_edger[[comp]] = lec
        limma_vs_deseq[[comp]] = ldc
        edger_vs_deseq[[comp]] = edc
    } ## End loop
    names(limma_vs_edger) = names(deseq)
    names(limma_vs_deseq) = names(deseq)
    names(edger_vs_deseq) = names(deseq)

    tt = rbind(as.numeric(limma_vs_edger), as.numeric(limma_vs_deseq))
    tt = rbind(tt, as.numeric(edger_vs_deseq))
    tt = as.matrix(tt)
    rownames(tt) = c("le", "ld", "ed")
    colnames(tt) = names(deseq)
    heat_colors = colorRampPalette(c("white","black"))
    sillytime = heatmap.3(tt, scale="none", trace="none", linewidth=0.5, keysize=2, margins=c(8,8), col=heat_colors, dendrogram="none", Rowv=FALSE, Colv=FALSE, main="Compare DE tools")
    heat = recordPlot()
    ret = list(limma_vs_edger=limma_vs_edger, limma_vs_deseq=limma_vs_deseq, edger_vs_deseq, heat=heat)
    return(ret)
}

#' make_pairwise_contrasts(): Run makeContrasts() with all pairwise comparisons.
#' 
#' @param model A model describing the conditions/batches/etc in the experiment
#' @param conditions A factor of conditions in the experiment
#' @param do_identities Whether or not to include all the identity strings.
#' Limma can handle this, edgeR cannot.  True by default.
#' @param do_pairwise Whether or not to include all the pairwise strings.
#' This shouldn't need to be set to FALSE, but just in case.
#' @param extra_contrasts An optional string of extra contrasts to include.
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
    ## makeContrasts(bob=bob, alice=alice, bob_minus_alice=(bob)-(alice), levels=design)
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
            minus_string = paste(d_name, "_minus_", c_name, sep="")
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

#' Plot arbitrary data from limma
#'
#' @param all_pairwise_result The result from calling balanced_pairwise()
#' @param first_name A table inside all_pairwise_result$limma_result
#' @param first_type A column within the chosen table
#' @param second_name Another table inside all_pairwise_result$limma_result
#' @param second_type A column to compare against
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
    tables = all_pairwise_result$limma_result
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

#' simple_comparison():  Perform a simple experimental/control comparison
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
#' @param subset an experimental subset with two conditions to compare.
#' @param model a model describing the experiment.
#' @param workbook an excel workbook to which to write. simple_comparison.xls by default.
#' @param worksheet an excel worksheet to which to write. simple_comparison by default.
#' @param basename a url to which to send click evens in clicky volcano/ma plots.  NA by default.
#' @param batch whether or not to include batch in limma's model.  TRUE by default.
#' @param combat whether or not to use combatMod().  FALSE by default.
#' @param combat_noscale whether or not to include combat_noscale (makes combat a little less heavy-handed).  TRUE by default.
#' @param pvalue_cutoff p-value definition of 'significant.'  0.05 by default.
#' @param logfc_cutoff fold-change cutoff of significance. 0.6 (and therefor 1.6) by default.
#' @param tooltip_data Text descriptions of genes if one wants google graphs.
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
        print("Hey! There is only one batch! I am only including condition in the model!")
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
    colnames(expt_data)
    if (combat) {
#        expt_data = ComBat(expt_data, subset$batches, condition_model)
        expt_data = cbcbSEQ::combatMod(expt_data, subset$batches, subset$conditions)
    }
    expt_voom = hpgltools::hpgl_voom(expt_data, model)
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

#' make_exampledata():  A small hack of limma's exampleData()
#' function to allow for arbitrary data set sizes.
#'
#' @param ngenes how many genes in the fictional data set.  1000 by default.
#' @param columns how many samples in this data set.  5 by default.
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
