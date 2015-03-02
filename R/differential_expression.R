## Time-stamp: <Sun Mar  1 13:15:26 2015 Ashton Trey Belew (abelew@gmail.com)>
## differential_expression.R contains functions useful for differential expression tasks.

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
        message(paste("Printing table: ", comparison, sep=""))        
        data_table = topTable(data, adjust=adjust, n=n, coef=comparison)

        data_table$qvalue = tryCatch(
            {
                qvalue(data_table$P.Value, gui=FALSE)$qvalues
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
hpgl_voom = function(dataframe, model, libsize=NULL, stupid=FALSE) {
    out = list()
    if (is.null(libsize)) {
        libsize = colSums(dataframe, na.rm=TRUE)
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

#'
#' 
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



#' balanced_pairwise():  Set up a model matrix and set of contrasts to do
#' a pairwise comparison of all conditions/batches.  In this case, there
#' must be a balanced set of batches for each condition.
#'
#' @param conditions a factor of conditions in the experiment
#' @param batches a factor of batches in the experiment
#' @param extra_contrasts some extra contrasts to add to the list
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_minus_b_ctrla = (C-B)-A, e_minus_d_ctrla = (E-D)-A,
#'  de_minus_cb = (E-D)-(C-B),"
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
#' ## pretend = balanced_pairwise(data, conditions, batches)
balanced_pairwise = function(data, conditions, batches, extra_contrasts=NULL, ...) {
    condition_table = table(conditions)
    batch_table = table(batches)
    conditions = as.factor(conditions)
    batches = as.factor(batches)
    ## Make a model matrix which will have one entry for each of these condition/batches
    fun_model = model.matrix(~0 + conditions + batches)
    tmpnames = colnames(fun_model)
    tmpnames = gsub("data[[:punct:]]", "", tmpnames)
    tmpnames = gsub("conditions", "", tmpnames)
    colnames(fun_model) = tmpnames
    ## voom() it
    fun_voom = hpgl_voom(data, fun_model)
    ## Extract the design created by voom()
    ## This is interesting because each column of the design will have a prefix string 'macb' before the
    ## condition/batch string, so for the case of clbr_tryp_batch_C it will look like: macbclbr_tryp_batch_C
    ## This will be important in 17 lines from now.
    fun_design = fun_voom$design    
    ## Do the lmFit() using this model    
    fun_fit = lmFit(fun_voom, fun_model)
    ## The following three tables are used to quantify the relative contribution of each batch to the sample condition.
    identities = list()
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
            exprs_string = paste(minus_string, " = ", d_name, "-", c_name, ",", sep="")
            all_pairwise[minus_string] = exprs_string
        }
    }
    ## At this point, I have strings which represent the definition of every
    ## sample condition as well as strings which represent every possible
    ## B-A where B comes somewhere after A in the model matrix.
    ## The goal now is to create the variables in the R environment
    ## and add them to makeContrasts()
    eval_strings = append(identities, all_pairwise)
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
    contrast_string = paste(contrast_string, "levels=fun_design)")
    eval(parse(text=contrast_string))
    ## I like to change the column names of the contrasts because by default
    ## they are kind of obnoxious and too long to type
    
    if (!is.null(extra_contrasts)) {
        eval_names = append(eval_names, extra_eval_names)
    }
    colnames(all_pairwise_contrasts) = eval_names
    ## Once all that is done, perform the fit
    ## This will first provide the relative abundances of each condition
    ## followed by the set of all pairwise comparisons.
    all_pairwise_fits = contrasts.fit(fun_fit, all_pairwise_contrasts)
    all_pairwise_comparisons = eBayes(all_pairwise_fits)
    limma_result = write_limma(all_pairwise_comparisons, excel=FALSE)
    result = list(
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
        limma_result=limma_result)
    return(result)
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
#' @param logfc_cutoff fold-change cutoff of significance. 0.6 (and therefore 1.6) by default.
#' @param tooltip_data Text descriptions of genes if one wants google graphs.
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
