## Time-stamp: <Mon Nov 23 10:07:00 2015 Ashton Trey Belew (abelew@gmail.com)>

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

#' Beta.NA: Perform a quick solve to gather residuals etc
#' This was provided by Kwame for something which I don't remember a loong time ago.
Beta.NA = function(y,X) {
    des = X[!is.na(y),]
    y1 = y[!is.na(y)]
    B = solve(t(des)%*%des)%*%t(des)%*%y1
    return(B)
}

#' get_genelengths()  Grab gene lengths from a gff file.
#'
#' @param gff  a gff file with (hopefully) IDs and widths
#' @param type default='gene'  the annotation type to use.
#' @param key default='ID'  the identifier in the 10th column of the gff file to use.
#'
#' This function attempts to be robust to the differences in output from importing gff2/gff3 files.  But it certainly isn't perfect.
#'
#' @return  a data frame of gene IDs and widths.
#' @export
#' @seealso \code{\link{import.gff3}}, \code{\link{import.gff}}, \code{\link{import.gff2}}
#'
#' @examples
#' ## tt = hpgltools:::get_genelengths('reference/fun.gff.gz')
#' ## head(tt)
#' ##          ID width
#' ##1   YAL069W   312
#' ##2   YAL069W   315
#' ##3   YAL069W     3
#' ##4 YAL068W-A   252
#' ##5 YAL068W-A   255
#' ##6 YAL068W-A     3
get_genelengths = function(gff, type="gene", key='ID') {
    ret = gff2df(gff)
    ret = ret[ret$type == type,]
    ret = ret[,c(key,"width")]
    colnames(ret) = c("ID","width")
    if (dim(genelengths)[1] == 0) {
        stop(paste0("No genelengths were found.  Perhaps you are using the wrong 'type' or 'key' arguments, type is: ", type, ", key is: ", key))
    }
    return(ret)
}


#' sum_exons()  Given a data frame of exon counts and annotation information, sum the exons.
#'
#' @param data  a count table by exon
#' @param gff default=NULL  a gff filename
#' @param annotdf default=NULL  a dataframe of annotations (probably from gff2df)
#' @param parent default='Parent'  a column from the annotations with the gene names
#' @param child default='row.names'  a column from the annotations with the exon names
#'
#' This function will merge a count table to an annotation table by the child column.
#' It will then sum all rows of exons by parent gene and sum the widths of the exons.
#' Finally it will return a list containing a df of gene lengths and summed counts.
#'
#' @return  a list of 2 data frames.
#' @export
sum_exons = function(data, gff=NULL, annotdf=NULL, parent='Parent', child='row.names') {
    if (is.null(annotdf) & is.null(gff)) {
        stop("I need either a df with parents, children, and widths; or a gff filename.")
    } else if (is.null(annotdf)) {
        annotdf = gff2df(gff)
    }

    tmp_data = merge(data, annotdf, by=child)
    rownames(tmp_data) = tmp_data$Row.names
    tmp_data = tmp_data[-1]
    ## Start out by summing the gene widths
    column = aggregate(tmp_data[,"width"], by=list(Parent=tmp_data[,parent]), FUN=sum)
    new_data = data.frame(column$x)
    rownames(new_data) = column$Parent
    colnames(new_data) = c("width")

    for (c in 1:length(colnames(data))) {
        column_name = colnames(data)[[c]]
        column = aggregate(tmp_data[,column_name], by=list(Parent=tmp_data[,parent]), FUN=sum)
        rownames(column) = column$Parent
        new_data = cbind(new_data, column$x)
    } ## End for loop
    width_df = data.frame(new_data$width)
    rownames(width_df) = rownames(new_data)
    colnames(width_df) = c("width")
    new_data = new_data[-1]
    colnames(new_data) = colnames(data)
    rownames(new_data) = rownames(column)
    ret = list(width=width_df, counts=new_data)
    return(ret)
}

#' make_report()  Make a knitr report with some defaults set
#'
#' @param type default='pdf'  html/pdf/fancy html reports?
#'
#' @return a dated report file
make_report = function(name="report", type='pdf') {
    opts_knit$set(progress=FALSE, verbose=FALSE, error=FALSE, fig.width=7, fig.height=7)
    theme_set(theme_bw(base_size=10))
    options(java.parameters="-Xmx8g")
    set.seed(1)
    output_date = format(Sys.time(), "%Y%m%d-%H%M")
    input_filename = name
    ## In case I add .rmd on the end.
    input_filename = gsub("\\.rmd", "", input_filename, perl=TRUE)
    input_filename = paste0(input_filename, ".rmd")
    if (type == 'html') {
        output_filename = paste0(name, "-", output_date, ".html")
        output_format = 'html_document'
        render(output_filename, output_format)
    } else if (type == 'pdf') {
        output_filename = paste0(name, "-", output_date, ".pdf")
        output_format = 'pdf_document'
    } else {
        output_filename = paste0(name, "-", output_date, ".html")
        output_format = 'knitrBootstrap::bootstrap_document'
    }
    message(paste0("About to run: render(input=", input_filename, ", output_file=", output_filename, " and output_format=", output_format))
    result = try(render(input=input_filename, output_file=output_filename, output_format=output_format), silent=TRUE)
    return(result)
}


#' gff2df()  Try to make import.gff a little more robust
#'
#' @param gff  a gff filename
#'
#' This function wraps import.gff/import.gff3/import.gff2 calls in try()
#' Because sometimes those functions fail in unpredictable ways.
#'
#' @export
#' @return  a df!
gff2df = function(gff) {
    ret = NULL
    annotations = try(import.gff3(gff), silent=TRUE)
    if (class(annotations) == 'try-error') {
        annotations = try(import.gff2(gff), silent=TRUE)
        if (class(annotations) == 'try-error') {
            stop("Could not extract the widths from the gff file.")
        } else {
            ret = annotations
        }
    } else {
        ret = annotations
    }
    ## The call to as.data.frame must be specified with the GenomicRanges namespace, otherwise one gets an error about
    ## no method to coerce an S4 class to a vector.
    ret = GenomicRanges::as.data.frame(ret)
    return(ret)
}

#' gff2irange()  Try to make import.gff a little more robust
#'
#' @param gff  a gff filename
#'
#' This function wraps import.gff/import.gff3/import.gff2 calls in try()
#' Because sometimes those functions fail in unpredictable ways.
#'
#' @export
#' @return  an iranges! (useful for getSeq())
gff2irange = function(gff) {
    ret = NULL
    annotations = try(import.gff3(gff), silent=TRUE)
    if (class(annotations) == 'try-error') {
        annotations = try(import.gff2(gff), silent=TRUE)
        if (class(annotations) == 'try-error') {
            stop("Could not extract the widths from the gff file.")
        } else {
            ret = annotations
        }
    } else {
        ret = annotations
    }
    ## The call to as.data.frame must be specified with the GenomicRanges namespace, otherwise one gets an error about
    ## no method to coerce an S4 class to a vector.
     return(ret)
}

#' hpgl_cor()  Wrap cor() to include robust correlations.
#'
#' @param df  a data frame to test.
#' @param method default='pearson'  correlation method to use. Includes pearson, spearman, kendal, robust.
#' @param ...  other options to pass to stats::cor()
#'
#' @return  correlation some fun correlation statistics
#' @seealso \code{\link{cor}}, \code{\link{cov}}, \code{\link{covRob}}
#'
#' @export
#' @examples
#' ## hpgl_cor(df=df)
#' ## hpgl_cor(df=df, method="robust")
hpgl_cor = function(df, method="pearson", ...) {
    if (method == "robust") {
        robust_cov = robust::covRob(df, corr=TRUE)
        correlation = robust_cov$cov
    } else {
        correlation = stats::cor(df, method=method, ...)
    }
    return(correlation)
}

#' make_tooltips()  Create a simple df from gff which contains tooltip usable information for gVis graphs.
#'
#' @param gff or annotations: Either a gff file or annotation data frame (which likely came from a gff file.)
#'
#' @return a df of tooltip information
make_tooltips = function(annotations=NULL, gff=NULL, desc_col='description') {
    if (is.null(annotations) & is.null(gff)) {
        stop("I need either a data frame or gff file.")
    } else {
        if (!is.null(annotations)) {
            tooltip_data = annotations[,c("ID", desc_col)]
        } else {
            ret = NULL
            annotations = try(import.gff3(gff), silent=TRUE)
            if (class(annotations) == 'try-error') {
                annotations = try(import.gff2(gff), silent=TRUE)
                if (class(annotations) == 'try-error') {
                    stop("Could not extract the widths from the gff file.")
                } else {
                    ret = annotations
                }
            } else {
                ret = annotations
            }
            ## The call to as.data.frame must be specified with the GenomicRanges namespace, otherwise one gets an error about
            ## no method to coerce an S4 class to a vector.
            tooltip_data = GenomicRanges::as.data.frame(ret)
        }
    }
    tooltip_data$tooltip = ""
    if (is.null(tooltip_data[[desc_col]])) {
        stop("I need a name!")
    } else {
        tooltip_data$tooltip = paste0(tooltip_data$ID, ': ', tooltip_data[[desc_col]])
    }
    tooltip_data$tooltip = gsub("\\+", " ", tooltip_data$tooltip)
    tooltip_data$tooltip = gsub(": $", "", tooltip_data$tooltip)
    tooltip_data$tooltip = gsub("^: ", "", tooltip_data$tooltip)
    rownames(tooltip_data) = make.names(tooltip_data$ID, unique=TRUE)
    colnames(tooltip_data) = c("1.tooltip")
    return(tooltip_data)
}

#' pattern_count_genome()  Find how many times a given pattern occurs in every gene of a genome.
#'
#' @param fasta  a fasta genome
#' @param gff default=NULL  an optional gff of annotations (if not provided it will just ask the whole genome.
#' @param pattern default='TA'  what pattern to search for?  This was used for tnseq and TA is the mariner insertion point.
#' @param key default='locus_tag'  what type of entry of the gff file to key from?
#'
#' @return num_pattern a data frame of names and numbers.
#' @export
#' @seealso \code{\link{PDict}} \code{\link{FaFile}}
#' @examples
#' ## num_pattern = pattern_count_genome('mgas_5005.fasta', 'mgas_5005.gff')
pattern_count_genome = function(fasta, gff=NULL, pattern='TA', type='gene', key='locus_tag') {
    rawseq = FaFile(fasta)
    if (is.null(gff)) {
        entry_sequences = rawseq
    } else {
        entries = import.gff3(gff, asRangedData=FALSE)
        type_entries = subset(entries, type==type)
        names(type_entries) = rownames(type_entries)
        entry_sequences = getSeq(rawseq, type_entries)
        names(entry_sequences) = entry_sequences[[,key]]
    }
    dict = PDict(pattern, max.mismatch=0)
    result = vcountPDict(dict, entry_sequences)
    num_pattern = data.frame(name=names(entry_sequences), num=as.data.frame(t(result)))
    return(num_pattern)
}

#' sillydist()  A stupid distance function of a point against two axes.
#'
#' @param firstterm  the x-values of the points.
#' @param secondterm  the y-values of the points.
#' @param firstaxis  the x-value of the vertical axis.
#' @param secondaxis  the y-value of the second axis.
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

#' write_xls()  Write a dataframe to an excel spreadsheet sheet.
#'
#' @param data  a dataframe of information.
#' @param sheet default='first'  the name of an excel sheet in a workbook.
#' @param file default='excel/workbook.xls'  an excel workbook to which to write.
#' @param rowname default='rownames'  what will the rownames be?
#' @param overwritefile default=FALSE  overwrite the xls file with this new data, or use the original?
#' @param overwritesheet default=TRUE  overwrite the xls sheet with this new data?  (if true it will make a backup sheet .bak).
#'
#' @return NULL, on the say it creates a workbook if necessary,
#' creates a sheet, and writes the data to it.
#'
#' @seealso \code{\link{loadWorkbook}}, \code{\link{createSheet}},
#' \code{\link{writeWorksheet}}, \code{\link{saveWorkbook}}
#'
#' @export
#' @examples
#' ## write_xls(dataframe, "hpgl_data")
#' ## Sometimes it is a good idea to go in and delete the workbook and
#' ## re-create it if this is used heavily, because it will get crufty.
write_xls = function(data, sheet="first", file="excel/workbook", rowname="rownames", overwritefile=FALSE, overwritesheet=TRUE, dated=TRUE, suffix=".xls") {
    excel_dir = dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    file = gsub(pattern="\\.xls.", replacement="", file, perl=TRUE)
    filename = NULL
    if (isTRUE(dated)) {
        timestamp = format(Sys.time(), "%Y%m%d%H")
        filename = paste0(file, "-", timestamp, suffix)
    } else {
        filename = paste0(file, suffix)
    }

    if (file.exists(filename)) {
        if (isTRUE(overwritefile)) {
            backup_file(filename)
        }
    }
    xls = loadWorkbook(filename, create=TRUE)

    if (isTRUE(overwritesheet)) {
        newname = paste0(sheet, '.bak')
        if (existsSheet(xls, newname)) {
            removeSheet(xls, sheet=newname)
        }
        if (existsSheet(xls, sheet)) {
            renameSheet(xls, sheet=sheet, newName=newname)
        }
    }

    createSheet(xls, name=sheet)
    if (is.na(rowname)) {
        writeWorksheet(xls, data, sheet=sheet)
    } else {
        writeWorksheet(xls, data, sheet=sheet, rowname=rowname)
    }
    saveWorkbook(xls)
}

#' backup_file()  Make a backup of an existing file with n revisions, like VMS!
#'
#' @param file  the file to backup.
#' @param backups default=10  how many revisions?
backup_file = function(backup_file, backups=10) {
    if (file.exists(backup_file)) {
        for (i in backups:01) {
            j = i + 1
            i = sprintf("%02d", i)
            j = sprintf("%02d", j)
            test = paste0(backup_file, ".", i)
            new = paste0(backup_file, ".", j)
            if (file.exists(test)) {
                file.rename(test, new)
            }
        }
        newfile = paste0(backup_file, ".", i)
        message(paste0("Renaming ", backup_file, " to ", newfile, "."))
        file.rename(backup_file, newfile)
    } else {
        message("The file does not yet exist.")
    }
}

loadme = function(dir="savefiles") {
    savefile = paste0(getwd(), "/", dir, "/RData.rda.xz")
    message(paste0("Loading the savefile: ", savefile))
    load_string <- paste0("load(", savefile, ", envir=globalenv())")
    eval(parse(text=load_string))
}

#' saveme()  Make a backup rdata file for future reference
#'
#' @param dir  the directory to save the Rdata file.
#' @param backups default=10  how many revisions?
saveme = function(directory="savefiles", backups=4) {
    environment()
    if (!file.exists(directory)) {
        dir.create(directory)
    }
    savefile = paste0(getwd(), "/", directory, "/RData.rda.xz")
    message(paste0("The savefile is: ", savefile))
    backup_file(savefile, backups=backups)
    ## The following save strings work:
    ## save_string <- paste0("save(list=ls(all.names=TRUE, envir=globalenv()), envir=globalenv(), file='", savefile, "')")
    ## save_string <- paste0("con <- base::pipe(paste0('pigz -p8 > ", savefile, "'), 'wb');\n save(list=ls(all.names=TRUE, envir=globalenv(), envir=globalenv(), file=con);\n close(con)")
    save_string <- paste0("con <- base::pipe(paste0('pxz -T4 > ", savefile, "'), 'wb');\n save(list=ls(all.names=TRUE, envir=globalenv()), envir=globalenv(), file=con, compress=FALSE);\n close(con)")
    message(paste0("The save string is: ", save_string))
    eval(parse(text=save_string))
}

## EOF
