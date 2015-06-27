## Time-stamp: <Tue Jun 23 13:13:50 2015 Ashton Trey Belew (abelew@gmail.com)>

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

#' Beta.NA: Perform a quick solve to gather residuals etc
#' This was provided by Kwame for something which I don't remember a loong time ago.
Beta.NA = function(y,X) {
    des = X[!is.na(y),]
    y1 = y[!is.na(y)]
    B = solve(t(des)%*%des)%*%t(des)%*%y1
    return(B)
}

#' Grab gene lengths from a gff file
#'
#' @param gff a gff file with (hopefully) IDs and widths
#' @param type The annotation type to use (gene)
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
get_genelengths = function(gff, type="gene") {
    ret = NULL
    annotations = try(data.frame(import.gff3(gff)), silent=TRUE)
    if (class(annotations) == 'try-error') {
        annotations = try(BiocGenerics:::as.data.frame(import.gff2(gff)), silent=TRUE)
        if (class(annotations) == 'try-error') {
            stop("Could not extract the widths from the gff file.")
        } else {
            ret = annotations
        }
    } else {
        ret = annotations
    }
    ##ret = subset(ret, type==type)
    ret = ret[ret$type == type,]
    ret = ret[,c("ID","width")]
    return(ret)
}

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
#' ## hpgl_cor(df=df)
#' ## hpgl_cor(df=df, method="robust")
hpgl_cor = function(df=NULL, method="pearson", ...) {
    if (method == "robust") {
        robust_cov = robust::covRob(df, corr=TRUE)
        correlation = robust_cov$cov
    } else {
        correlation = stats::cor(df, method=method, ...)
    }
    return(correlation)
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
#' ## write_xls(dataframe, "hpgl_data")
#' ## Sometimes it is a good idea to go in and delete the workbook and
#' ## re-create it if this is used heavily, because it will get crufty.
write_xls = function(data=NULL, sheet="first", file="excel/workbook.xls", rowname="rownames", overwrite=FALSE) {
    excel_dir = dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }
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

## EOF

