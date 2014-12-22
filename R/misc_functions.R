## Time-stamp: "Wed Dec 17 15:24:47 2014 Ashton Trey Belew (abelew@gmail.com)"
## my_functions.R contains a host of (hopefully) useful functions when
## dealing with genomic data in R.

### This file is intended to keep a set of functions which will prove
### elsewhere in the analysis of ribosome profiling data.  These are
### mostly stolen from other people's code, notably Ramzi Temanni,
### Keith Hughitt, Laura Dillon, Yuan Li, and Kwame.

#' Wrap cor() to include robust correlations
#'
#' @param df a data frame to test
#' @param method Correlation method to use.  Defaults to pearson.
#' Includes pearson, spearman, kendal, robust.
#' 
#' @return  correlation some fun correlation statistics
#' @seealso \code{\link{cor}}, \code{\link{cov}}, \code{\link{covRob}}
#' 
#' @export
#' @examples
#' ## my_cor(df=df)
#' ## my_cor(df=df, method="robust")
my_cor = function(df=NULL, method="pearson", ...) {
    if (method == "robust") {
        robust_cov = robust::covRob(df, corr=TRUE)
        correlation = robust_cov$cov
    } else {
        correlation = stats::cor(df, method=method, ...)
    }
    return(correlation)
}

#' Write limma
#'
#' @param data The output from eBayes()
#' @param adjust The pvalue adjustment chosen (fdr by default)
#' @param n The number of entries to report, defaults to 0, which says do them all
#' @param coef which coefficients/contrasts to report, NULL says do them all
#' @param workbook an excel filename into which to write the data, used for csv files too.
#' @param excel T/F whether or not to write an excel workbook (useful if they are too big)
#'
#' @return a list of data frames comprising the toptable output for each coefficient,
#'    I also added a qvalue entry.
#' @seealso \code{\link{toptable}}. \code{\link{write_xls}}
#'
#' @export
#' @examples
#' ## finished_comparison = eBayes(limma_output)
#' ## data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
write_limma = function(data=NULL, adjust="fdr", n=0, coef=NULL, workbook="excel/limma.xls", excel=TRUE) {
    csv_basename = gsub("xls$", "", workbook)
    testdir = dirname(workbook)
    if (!file.exists(testdir)) {
        dir.create(testdir)
        print(paste("Creating directory: ", testdir, " for writing excel/csv data.", sep=""))
    }
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
        print(paste("Printing table: ", comparison, sep=""))        
        data_table = topTable(data, adjust=adjust, n=n, coef=comparison)
        data_table$qvalue = qvalue(data_table$P.Value)$qvalues
        ## This write_xls runs out of memory annoyingly often
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwrite=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        csv_filename = paste(csv_basename, "_", comparison, ".csv", sep="")
        write.csv(data_table, file=csv_filename)
        return_data[[comparison]] = data_table
    }
    return(return_data)
}

    
#' Write a dataframe to an excel spreadsheet sheet.
#'
#' @param data a dataframe of information
#' @param sheet the name of an excel sheet in a workbook.
#' @param file an excel workbook to which to write.  Defaults to "excel/workbook.xls"
#' @param rowname include rownames?  Defalts to no.
#' 
#' @return NULL, on the say it creates a workbook if necessary,
#' creates a sheet, and writes the data to it.
#' 
#' @seealso \code{\link{loadWorkbook}}, \code{\link{createSheet}},
#' \code{\link{writeWorksheet}}, \code{\link{saveWorkbook}}
#' 
#' @export
#' @examples
#' ## write_xls(dataframe, "my_data")
#' ## Sometimes it is a good idea to go in and delete the workbook and
#' ## re-create it if this is used heavily, because it will get crufty.
write_xls = function(data=NULL, sheet="first", file="excel/workbook.xls", rowname="rownames", overwrite=FALSE) {
    ## Testing parameters
    ##data = epi_cl14_table
    ##sheet = "epi_cl14"
    ##file = "excel/abundances.xls"
    ##rowname = NA
    ## End testing parameters
    if (is.na(file.info(file)$size)) {
        xls = XLConnect::loadWorkbook(file, create=TRUE)
    } else if (file.info(file)$size == 0) {
        file.remove(file)
        xls = XLConnect::loadWorkbook(file, create=TRUE)        
    } else {
        xls = XLConnect::loadWorkbook(file)
    }

    write_sheet = 0
    if (isTRUE(existsSheet(xls, sheet))) {
        if (isTRUE(overwrite)) {
            write_sheet = 1
        } else {
            write_sheet = 0
        }
    } else {
        write_sheet = 1
    }

    if (write_sheet == 1) {
        XLConnect::createSheet(xls, name=sheet)
        if (is.na(rowname)) {
            XLConnect::writeWorksheet(xls, data, sheet=sheet)
        } else {
            XLConnect::writeWorksheet(xls, data, sheet=sheet, rowname=rowname)
        }
        XLConnect::saveWorkbook(xls)
    }
}

#' Compile a table including annotation information and differential expression
#'
#' @param datTable df of differential expression data (toptable)
#' @param anot.file annotation file (gff)
#' @param source.format the format of the annotation
#' information. Defaults to "gff".
#' @return resanot a df of the DE/annotation information.
#' 
#' @seealso \code{\link{readGff3}}
#' 
#' @export
annotate_data = function (datTable, anot.file, source.format="gff",
    SelFeature="gene", fields=c("ID","Name","description"), cleandescription=FALSE) {
    if (source.format=="gff") {
        genomeIntervals::gff = readGff3(anot.file,isRightOpen=FALSE)
        anotattribs = unique(as.character(annotation(gff)$type))
        if ((!SelFeature%in%anotattribs) & length(intersect(fields,anotattribs) > 0 )) {
            stop("SelFeature not found or none of the element of fields found in gff column 9 !")
        } else {
            tmp = gff[annotation(gff)$type==SelFeature]
            myanot = as.data.frame(cbind(as.character(annotation(tmp)$seq_name),tmp[,1:2],getGffAttribute(tmp,fields)))
            colnames(myanot)[1:3] = c("chr","start","end")
            if (cleandescription) {
                for (i in fields) {
                    myanot = ascii2char(myanot,i)
                }
            }
            ##apply( myanot$description, 1, function(x) gsub("%2C",",",x) )
            ##myanot=getGffAttribute(gff[annotation(gff)$type==SelFeature],fields)
            resanot = cbind(datTable,myanot[match(tolower(datTable$id),tolower(myanot$ID)),])
            return(resanot)
        }  
    }
    ## anotate from txt tab delimited txt file 
    ## key column name in annotation file must be "ID"
    if (source.format=="txt") {
        myanot = read.delim(anot.file,header=T)
        resanot=cbind(datTable,myanot[match(tolower(datTable$id),tolower(myanot$ID)),])
        return(resanot)
    }
}

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
#' edgeR).  By default this is NULL.
#' 
#' @return an EList containing the modified data, weights, design,
#' libsize, and mean-variance plot.
#' 
#' @seealso \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{lmFit}}
#' 
#' @export
#' @examples
#' ## funkytown = my_voom(samples, model)
my_voom = function(dataframe, model, libsize=NULL, stupid=FALSE) {
    ## Test options
###    dataframe = expt_data
###    model = condbatch_model
###    libsize = NULL
    ##
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
##        theme_bw()
    if (is.null(linear_fit$rank)) {
        print("Some samples cannot be balanced across the experimental design.")
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

#' A stupid distance function of a point against two axes
#'
#' @param firstterm the x-values of the points
#' @param secondterm the y-values of the points
#' @param firstaxis the x-value of the vertical axis
#' @param secondaxis the y-value of the second axis
#'
#' @return dataframe of the distances
#' @export
sillydist = function(firstterm, secondterm, firstaxis, secondaxis) {
    dataframe = data.frame(firstterm, secondterm)
    dataframe$x = (abs(dataframe[,1]) - abs(firstaxis)) / abs(firstaxis)
    dataframe$y = abs((dataframe[,2] - secondaxis) / secondaxis)
    dataframe$x = abs(dataframe[,1] / max(dataframe$x))
    dataframe$y = abs(dataframe[,2] / max(dataframe$y))
    dataframe$dist = abs(dataframe$x * dataframe$y)
    dataframe$dist = dataframe$dist / max(dataframe$dist)
    return(dataframe)
}

#' Perform a simple experimental/control comparison
#'
#' @param subset an experimental subset with two conditions to compare
#' @param model a model describing the experiment
#' @param workbook an excel workbook to write data to, defaults to simple_comparison.xls
#' @param worksheet an excel worksheet to write data to, also defaults to simple_comparison
#' @param basename a url to which to send click evens in clicky volcano/ma plots
#' @return lots of text, some plots, and a list which includes:
#' the contrast model, the outputs from voom/lmFit/eBayes/toptable
#' histograms of the pvalues and mean values, and the intermediate plots
#' 
#' @seealso \code{\link{my_gvis_ma_plot}}, \code{\link{toptable}},
#' \code{\link{voom}}, \code{\link{voomMod}}, \code{\link{my_voom}},
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
simple_comparison = function(subset, workbook="simple_comparison.xls", sheet="simple_comparison", basename=NA, batch=TRUE, combat=FALSE, combat_noscale=FALSE, pvalue_cutoff=0.05, logfc_cutoff=0.6, tooltip_data=NULL, verbose=FALSE, ...) {
    ### Testing parameters
    ##subset = all_norm_expt_5448
    ##basename="html/t0t1"
    ##tooltip_data=tooltip_data_5448
    ### End Testing parameters
    condition_model = stats::model.matrix(~ 0 + subset$condition)
    condbatch_model = stats::model.matrix(~ 0 + subset$condition + subset$batch)
    if (batch) {
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
    expt_voom = myr::my_voom(expt_data, model)
    lf = limma::lmFit(expt_voom)
    colnames(lf$coefficients)
    coefficient_scatter = myr::my_linear_scatter(lf$coefficients)
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
    contrast_histogram = myr::my_histogram(hist_df)
    hist_df = data.frame(values=cond_contrasts$Amean)
    amean_histogram = myr::my_histogram(hist_df, fillcolor="pink", color="red")
    coef_amean_cor = cor.test(cond_contrasts$coefficients, cond_contrasts$Amean, exact=FALSE)
    cond_comparison = limma::eBayes(cond_contrasts)
    hist_df = data.frame(values=cond_comparison$p.value)
    pvalue_histogram = myr::my_histogram(hist_df, fillcolor="lightblue", color="blue")
    cond_table = limma::topTable(cond_comparison, number=nrow(expt_voom$E), coef="changed_v_control", sort.by="logFC")
    if (!is.na(basename)) {
        vol_gvis_filename = paste(basename, "volplot.html", sep="_")
        a_volcano_plot = myr::my_volcano_plot(cond_table, gvis_filename=vol_gvis_filename, tooltip_data=tooltip_data)
    } else {
        a_volcano_plot = myr::my_volcano_plot(cond_table)
    }
    if (!is.na(basename)) {
        ma_gvis_filename=paste(basename, "maplot.html", sep="_")
        an_ma_plot = myr::my_ma_plot(expt_voom$E, cond_table, gvis_filename=ma_gvis_filename, tooltip_data=tooltip_data)
    } else {
        an_ma_plot = myr::my_ma_plot(expt_voom$E, cond_table)
    }
    myr::write_xls(cond_table, sheet, file=workbook, rowname="row.names")
    upsignificant_table = subset(cond_table, logFC >=  logfc_cutoff)
    downsignificant_table = subset(cond_table, logFC <= (-1 * logfc_cutoff))
#    psignificant_table = subset(cond_table, adj.P.Val <= pvalue_cutoff)
    psignificant_table = subset(cond_table, P.Value <= pvalue_cutoff)    

    if (verbose) {
        print("The model looks like:")
        print(model)
        print("The mean:variance trend follows")
        plot(expt_voom$plot)
        print("Drawing a scatterplot of the genes.")
        print("The following statistics describe the relationship between:")
        print(coefficient_scatter$scatter)
        print(paste("Setting the column:", colnames(lf$design)[2], "to control"))
        print(paste("Setting the column:", colnames(lf$design)[1], "to changed"))
        print("Performing contrasts of the experimental - control.")        
        print("Taking a histogram of the subtraction values.")
        print(contrast_histogram)
        print("Taking a histogram of the mean values across samples.")
        print("The subtraction values should not be related to the mean values.")
        print(coef_amean_cor)
        print("Making a table of the data including p-values and F-statistics.")
        print("Taking a histogram of the p-values.")
        print(pvalue_histogram)
        print("Printing a volcano plot of this data.")
        print("Printing an maplot of this data.")
        print(paste("Writing excel sheet:", sheet))               
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

Beta.NA = function(y,X) {
    des = X[!is.na(y),]
    y1 = y[!is.na(y)]
    B = solve(t(des)%*%des)%*%t(des)%*%y1
    B
}

#' limma's exampleData() function for arbitrary data set sizes
#'
#' @param ngenes how many genes in the fictional data set (1000 by default)
#' @param columns how many samples in this data set (5)
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
