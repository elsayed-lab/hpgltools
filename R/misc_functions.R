## Time-stamp: <Wed Sep 16 11:29:09 2015 Ashton Trey Belew (abelew@gmail.com)>

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
    ret = ret[ret$type == type,]
    ret = ret[,c(key,"width")]
    colnames(ret) = c("ID","width")
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

make_tooltips = function(annotations=NULL, gff=NULL) {
    if (is.null(annotations) & is.null(gff)) {
        stop("I need either a data frame or gff file.")
    } else {
        if (!is.null(annotations)) {
            tooltip_data = annotations[,c("ID","Name","locus_tag")]
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
    if (is.null(tooltip_data$locus_tag) & is.null(tooltip_data$Name) & is.null(tooltip_data$ID)) {
        stop("I need a name!")
    } else {
        if (is.null(tooltip_data$locus_tag)) {
            tooltip_data$locus_tag = ""
        } else {
            tooltip_data$tooltip = paste(tooltip_data$tooltip, tooltip_data$locus_tag, sep=": ")
        } 
        if (is.null(tooltip_data$Name)) {
            tooltip_data$Name = ""
        } else {
            tooltip_data$tooltip = paste(tooltip_data$tooltip, tooltip_data$Name, sep=": ")
        }
        if (is.null(tooltip_data$ID)) {
            tooltip_data$ID = ""
        } else {
            tooltip_data$tooltip = paste(tooltip_data$tooltip, tooltip_data$ID, sep=": ")
        }
    }
    tooltip_data$tooltip = gsub("\\+", " ", tooltip_data$tooltip)
    tooltip_data$tooltip = gsub(": $", "", tooltip_data$tooltip)
    tooltip_data$tooltip = gsub("^: ", "", tooltip_data$tooltip)
    rownames(tooltip_data) = make.names(tooltip_data$ID, unique=TRUE)
    tooltip_data = tooltip_data[,c("ID","Name","locus_tag", "1.tooltip")]
    tooltip_data = tooltip_data[-1]
    tooltip_data = tooltip_data[-1]
    tooltip_data = tooltip_data[-1]
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
write_xls = function(data, sheet="first", file="excel/workbook.xls", rowname="rownames", overwritefile=FALSE, overwritesheet=TRUE) {
    excel_dir = dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    if (file.exists(file)) {
        if (isTRUE(overwritefile)) {
            backup_file(file)
        }
    }
    xls = loadWorkbook(file, create=TRUE)

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
backup_file = function(file, backups=10) {
    for (i in backups:01) {
        j = i + 1
        i = sprintf("%02d", i)
        j = sprintf("%02d", j)
        test = paste0(file, ".", i)
        new = paste0(file, ".", j)
        if (file.exists(test)) {
            file.rename(test, new)
        }
    }
    file.rename(file, paste0(file, ".", i))
}

## EOF
