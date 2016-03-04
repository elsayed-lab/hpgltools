## Time-stamp: <Thu Feb 25 14:31:11 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Beta.NA: Perform a quick solve to gather residuals etc
#' This was provided by Kwame for something which I don't remember a loong time ago.
#'
#' @param y  a y
#' @param X  a x
Beta.NA <- function(y,X) {
    des <- X[!is.na(y),]
    y1 <- y[!is.na(y)]
    B <- solve(t(des)%*%des)%*%t(des)%*%y1
    return(B)
}

#' Grab gene lengths from a gff file.
#'
#' This function attempts to be robust to the differences in output from importing gff2/gff3 files.  But it certainly isn't perfect.
#'
#' @param gff  a gff file with (hopefully) IDs and widths
#' @param type   the annotation type to use.
#' @param key   the identifier in the 10th column of the gff file to use.
#' @return  a data frame of gene IDs and widths.
#' @seealso \pkg{rtracklayer} \link[rtracklayer]{import.gff}
#' @examples
#' \dontrun{
#'  tt = get_genelengths('reference/fun.gff.gz')
#'  head(tt)
#' #          ID width
#' #1   YAL069W   312
#' #2   YAL069W   315
#' #3   YAL069W     3
#' #4 YAL068W-A   252
#' #5 YAL068W-A   255
#' #6 YAL068W-A     3
#' }
#' @export
get_genelengths <- function(gff, type="gene", key='ID') {
    ret <- gff2df(gff)
    ret <- ret[ret$type == type,]
    ret <- ret[,c(key,"width")]
    colnames(ret) <- c("ID","width")
    if (dim(ret)[1] == 0) {
        stop(paste0("No genelengths were found.  Perhaps you are using the wrong 'type' or 'key' arguments, type is: ", type, ", key is: ", key))
    }
    return(ret)
}

#' Given a data frame of exon counts and annotation information, sum the exons.
#'
#' This function will merge a count table to an annotation table by the child column.
#' It will then sum all rows of exons by parent gene and sum the widths of the exons.
#' Finally it will return a list containing a df of gene lengths and summed counts.
#'
#' @param data  a count table by exon
#' @param gff   a gff filename
#' @param annotdf   a dataframe of annotations (probably from gff2df)
#' @param parent   a column from the annotations with the gene names
#' @param child   a column from the annotations with the exon names
#'
#' @return  a list of 2 data frames, counts and lengths by summed exons
#' @seealso \pkg{rtracklayer}
#' @examples
#' \dontrun{
#' summed <- sum_exons(counts, gff='reference/xenopus_laevis.gff.xz')
#' }
#' @export
sum_exons <- function(data, gff=NULL, annotdf=NULL, parent='Parent', child='row.names') {
    if (is.null(annotdf) & is.null(gff)) {
        stop("I need either a df with parents, children, and widths; or a gff filename.")
    } else if (is.null(annotdf)) {
        annotdf <- gff2df(gff)
    }

    tmp_data <- merge(data, annotdf, by=child)
    rownames(tmp_data) <- tmp_data$Row.names
    tmp_data <- tmp_data[-1]
    ## Start out by summing the gene widths
    column <- aggregate(tmp_data[,"width"], by=list(Parent=tmp_data[,parent]), FUN=sum)
    new_data <- data.frame(column$x)
    rownames(new_data) <- column$Parent
    colnames(new_data) <- c("width")

    for (c in 1:length(colnames(data))) {
        column_name <- colnames(data)[[c]]
        column <- aggregate(tmp_data[,column_name], by=list(Parent=tmp_data[,parent]), FUN=sum)
        rownames(column) <- column$Parent
        new_data <- cbind(new_data, column$x)
    } ## End for loop
    width_df <- data.frame(new_data$width)
    rownames(width_df) <- rownames(new_data)
    colnames(width_df) <- c("width")
    new_data <- new_data[-1]
    colnames(new_data) <- colnames(data)
    rownames(new_data) <- rownames(column)
    ret <- list(width=width_df, counts=new_data)
    return(ret)
}

#' Make a knitr report with some defaults set
#'
#' @param name   Name the document!
#' @param type   html/pdf/fancy html reports?
#' @return a dated report file
#' @seealso \pkg{knitr} \pkg{rmarkdown} \pkg{knitrBootstrap}
#' @export
make_report <- function(name="report", type='pdf') {
    knitr::opts_knit$set(progress=FALSE, verbose=FALSE, error=FALSE, fig.width=7, fig.height=7)
    ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
    set.seed(1)
    output_date <- format(Sys.time(), "%Y%m%d-%H%M")
    input_filename <- name
    ## In case I add .rmd on the end.
    input_filename <- gsub("\\.rmd", "", input_filename, perl=TRUE)
    input_filename <- paste0(input_filename, ".rmd")
    if (type == 'html') {
        output_filename <- paste0(name, "-", output_date, ".html")
        output_format <- 'html_document'
        rmarkdown::render(output_filename, output_format)
    } else if (type == 'pdf') {
        output_filename <- paste0(name, "-", output_date, ".pdf")
        output_format <- 'pdf_document'
    } else {
        output_filename <- paste0(name, "-", output_date, ".html")
        output_format <- 'knitrBootstrap::bootstrap_document'
    }
    message(paste0("About to run: render(input=", input_filename, ", output_file=", output_filename, " and output_format=", output_format))
    result <- try(rmarkdown::render(
        input=input_filename,
        output_file=output_filename,
        output_format=output_format), silent=TRUE)
    return(result)
}

#'  Implement the arescan function in R
#'
#' This function was taken almost verbatim from AREScore() in SeqTools
#' Available at: https://github.com/lianos/seqtools.git
#' At least on my computer I could not make that implementation work
#' So I rewrapped its apply() calls and am now hoping to extend its logic
#' a little to make it more sensitive and get rid of some of the spurious
#' parameters or at least make them more transparent.
#'
#' Note that I did this two months ago and haven't touched it since...
#'
#' @param x  A DNA/RNA StringSet containing the UTR sequences of interest
#' @param basal   I dunno.
#' @param overlapping default=1.5
#' @param d1.3 default=0.75  These parameter names are so stupid, lets be realistic
#' @param d4.6 default=0.4
#' @param d7.9 default=0.2
#' @param within.AU default=0.3
#' @param aub.min.length default=10
#' @param aub.p.to.start default=0.8
#' @param aub.p.to.end default=0.55
#' @return a DataFrame of scores
#' @seealso \pkg{IRanges} \pkg{Biostrings}
#' @examples
#' \dontrun{
#' ## Extract all the genes from my genome, pull a static region 120nt following the stop
#' ## and test them for potential ARE sequences.
#' ## FIXME: There may be an error in this example, another version I have handles the +/- strand
#' ## genes separately, I need to return to this and check if it is providing the 5' UTR for 1/2
#' ## the genome, which would be unfortunate -- but the logic for testing remains the same.
#' are_candidates <- hpgl_arescore(genome)
#' utr_genes <- subset(lmajor_annotations, type == 'gene')
#' threep <- GenomicRanges::GRanges(seqnames=Rle(utr_genes[,1]),
#'                                ranges=IRanges(utr_genes[,3], end=(utr_genes[,3] + 120)), strand=Rle(utr_genes[,5]),
#'                                name=Rle(utr_genes[,10]))
#' threep_seqstrings <- Biostrings::getSeq(lm, threep)
#' are_test <- hpgltools:::hpgl_arescore(x=threep_seqstrings)
#' are_genes <- rownames(are_test[ which(are_test$score > 0), ])
#' }
#' @export
hpgl_arescore <- function (x, basal=1, overlapping=1.5, d1.3=0.75, d4.6=0.4,
                           d7.9=0.2, within.AU=0.3, aub.min.length=10, aub.p.to.start=0.8,
                           aub.p.to.end=0.55) {
    ## The seqtools package I am using is called in R 'SeqTools' (note the capital S T)
    ## However, the repository I want for it is 'seqtools'
    ## Ergo my stupid require.auto() will be confused by definition because it assumes equivalent names
    ##if (isTRUE('SeqTools' %in% .packages(all.available=TRUE))) {
    ##    library('SeqTools')
    ##} else {
    ##    require.auto("lianos/seqtools/R/pkg")
    ##    library('SeqTools')
    ##}
    xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
    if (xtype == "DNA") {
        pentamer <- "ATTTA"
        overmer <- "ATTTATTTA"
    } else {
        pentamer <- "AUUUA"
        overmer <- "AUUUAUUUA"
    }
    x <- as(x, "DNAStringSet")
    pmatches <- Biostrings::vmatchPattern(pentamer, x)
    omatches <- Biostrings::vmatchPattern(overmer, x)
    basal.score <- S4Vectors::elementLengths(pmatches) * basal
    over.score <- S4Vectors::elementLengths(omatches) * overlapping
    no.cluster <- data.frame(d1.3 = 0, d4.6 = 0, d7.9 = 0)
    clust <- lapply(pmatches, function(m) {
        if (length(m) < 2) {
            return(no.cluster)
        }
        wg <- BiocGenerics::width(IRanges::gaps(m))
        data.frame(d1.3=sum(wg <= 3), d4.6=sum(wg >= 4 & wg <= 6), d7.9=sum(wg >= 7 & wg <= 9))
    })
    clust <- do.call(rbind, clust)
    dscores <- clust$d1.3 * d1.3 + clust$d4.6 * d4.6 + clust$d7.9 *  d7.9
    ## require.auto("Biostrings")
    au.blocks <- my_identifyAUBlocks(x, aub.min.length, aub.p.to.start, aub.p.to.end)
    aub.score <- sum(IRanges::countOverlaps(pmatches, au.blocks) * within.AU)
    score <- basal.score + over.score + dscores + aub.score
    ans <- S4Vectors::DataFrame(score=score,
                                n.pentamer=S4Vectors::elementLengths(pmatches),
                                n.overmer=S4Vectors::elementLengths(omatches),
                                au.blocks=au.blocks,
                                n.au.blocks=S4Vectors::elementLengths(au.blocks))
    cbind(ans, S4Vectors::DataFrame(clust))
}

#' copy/paste the function from SeqTools
#' and find where it falls on its ass.
#'
#' Yeah, I do not remember what I changed in this function.
#'
#' @param x  A sequence object
#' @param min.length   I dunno.
#' @param p.to.start   the p to start of course
#' @param p.to.end   and the p to end
#'
#' @return a list of IRanges which contain a bunch of As and Us.
my_identifyAUBlocks <- function (x, min.length=20, p.to.start=0.8, p.to.end=0.55) {
    xtype = match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
    stopifnot(S4Vectors::isSingleNumber(min.length) && min.length >= 5 &&  min.length <= 50)
    stopifnot(S4Vectors::isSingleNumber(p.to.start) && p.to.start >= 0.5 && p.to.start <= 0.95)
    stopifnot(S4Vectors::isSingleNumber(p.to.end) && p.to.end >= 0.2 && p.to.end <= 0.7)
    stopifnot(p.to.start > p.to.end)
    if (xtype == "DNA") {
        AU <- "AT"
    } else {
        AU <- "AU"
    }
    y <- as(x, sprintf("%sStringSet", xtype))

    widths <- BiocGenerics::width(x)
    fun <- function(i) {
        one_seq <- x[[i]]
        au <- Biostrings::letterFrequencyInSlidingView(one_seq, min.length, AU, as.prob=TRUE)
        if (is.null(au) | nrow(au) == 0) {
                return(IRanges::IRanges())
            }
        au <- as.numeric(au)
        can.start <- au >= p.to.start
        can.end <- au <= p.to.end
        posts <- .Call("find_au_start_end", au, p.to.start, p.to.end, PACKAGE="SeqTools")
        blocks <- IRanges::IRanges(posts$start, posts$end + min.length -  1L)
        stats::end(blocks) <- ifelse(stats::end(blocks) > widths[i], widths[i], stats::end(blocks))
        IRanges::reduce(blocks)
    }
    au.blocks = lapply(1:length(x), fun)
    ret <- IRanges::IRangesList(au.blocks)
    return(ret)
}

#' Try to make import.gff a little more robust
#' I acquire (hopefully) valid gff3 files from various sources:
#' yeastgenome.org, microbesonline, tritrypdb, ucsc, ncbi.
#' To my eyes, they all look like reasonably good gff3 files, but
#' some of them must be loaded with import.gff2, import.gff3, etc.
#' That is super annoying.
#' Also, I pretty much always just do as.data.frame() when I get something
#' valid from rtracklayer, so this does that for me, I have another function
#' which returns the iranges etc.
#'
#' This function wraps import.gff/import.gff3/import.gff2 calls in try()
#' Because sometimes those functions fail in unpredictable ways.
#'
#' @param gff  a gff filename
#' @param type   subset the gff file for entries of a specific type
#' @return  a df!
#' @seealso \pkg{rtracklayer} \link[rtracklayer]{import.gff} \link[rtracklayer]{import.gff2} \link[rtracklayer]{import.gff3}
#' @examples
#' \dontrun{
#' funkytown <- gff2df('reference/gff/saccharomyces_cerevsiae.gff.xz')
#' }
#' @export
gff2df <- function(gff, type=NULL) {
    ret <- NULL
    gff_test <- grepl("\\.gff", gff)
    gtf_test <- grepl("\\.gtf", gff)
    annotations <- NULL
    if (isTRUE(gtf_test)) {  ## Start with an attempted import of gtf files.
        ret <- try(rtracklayer::import.gff(gff, format="gtf"), silent=TRUE)
    } else {
        annotations <- try(rtracklayer::import.gff3(gff), silent=TRUE)
        if (class(annotations) == 'try-error') {
            annotations <- try(rtracklayer::import.gff2(gff), silent=TRUE)
            if (class(annotations) == 'try-error') {
                stop("Could not extract the widths from the gff file.")
            } else {
                ret <- annotations
            }
        } else {
            ret <- annotations
        }
    } ## End else this is not a gtf file
    ## The call to as.data.frame must be specified with the GenomicRanges namespace, otherwise one gets an error about
    ## no method to coerce an S4 class to a vector.
    ret <- GenomicRanges::as.data.frame(ret)
    if (!is.null(type)) {
        index <- ret[, "type"] == type
        ret <- ret[index, ]
    }
    return(ret)
}

#' Try to make import.gff a little more robust
#'
#' @param gff  a gff filename
#' @param type   a subset to extract
#'
#' Essentially gff2df() above, but returns data suitable for getSet()
#'
#' @return  an iranges! (useful for getSeq())
#' @seealso \pkg{rtracklayer} \link{gff2df} \link[Biostrings]{getSeq}
#' @examples
#' \dontrun{
#' library(BSgenome.Tcruzi.clbrener.all)
#' tc_clb_all <- BSgenome.Tcruzi.clbrener.all
#' cds_ranges <- gff2irange('reference/gff/tcruzi_clbrener.gff.xz', type='CDS')
#' cds_sequences <- Biostrings::getSeq(tc_clb_all, cds_ranges)
#' }
#' @export
gff2irange <- function(gff, type=NULL) {
    ret <- NULL
    annotations <- try(rtracklayer::import.gff3(gff), silent=TRUE)
    if (class(annotations) == 'try-error') {
        annotations <- try(rtracklayer::import.gff2(gff), silent=TRUE)
        if (class(annotations) == 'try-error') {
            stop("Could not extract the widths from the gff file.")
        } else {
            ret <- annotations
        }
    } else {
        ret <- annotations
    }
    ## The call to as.data.frame must be specified with the GenomicRanges namespace, otherwise one gets an error about
    ## no method to coerce an S4 class to a vector.
    if (!is.null(type)) {
        index <- ret[, "type"] == type
        ret <- ret[index, ]
    }
    return(ret)
}

#' Wrap cor() to include robust correlations.
#'
#' @param df  a data frame to test.
#' @param method correlation method to use. Includes pearson, spearman, kendal, robust.
#' @param ...  other options to pass to stats::cor()
#' @return  correlation some fun correlation statistics
#' @seealso \pkg{robust} \link{cor} \link{cov} \link[robust]{covRob}
#' @examples
#' \dontrun{
#' hpgl_cor(df=df)
#' hpgl_cor(df=df, method="robust")
#' }
#' @export
hpgl_cor <- function(df, method="pearson", ...) {
    if (method == "robust") {
        robust_cov <- robust::covRob(df, corr=TRUE)
        correlation <- robust_cov$cov
    } else {
        correlation <- stats::cor(df, method=method, ...)
    }
    return(correlation)
}

#'   Create a simple df from gff which contains tooltip usable
#' information for gVis graphs. The tooltip column is also a handy proxy for
#' anontations information when it would otherwise be too troublesome.
#'
#' @param annotations  Either a gff file or annotation data frame (which likely came from a gff file.)
#' @param desc_col   a column from a gff file to grab the data from
#' @return a df of tooltip information or name of a gff file
#' @seealso \pkg{googleVis} \link{gff2df}
#' @examples
#' \dontrun{
#' tooltips <- make_tooltips('reference/gff/saccharomyces_cerevisiae.gff.gz')
#' }
#' @export
make_tooltips <- function(annotations, desc_col='description') {
    tooltip_data <- NULL
    if (class(annotations) == 'character') {
        tooltip_data <- gff2df(annotations)
    } else if (class(annotations) == 'data.frame') {
        tooltip_data <- annotations
    } else {
        stop("This requires either a filename or data frame.")
    }
    tooltip_data <- tooltip_data[,c("ID", desc_col)]
    tooltip_data$tooltip <- ""
    if (is.null(tooltip_data[[desc_col]])) {
        stop("I need a name!")
    } else {
        tooltip_data$tooltip <- paste0(tooltip_data$ID, ': ', tooltip_data[[desc_col]])
    }
    tooltip_data$tooltip <- gsub("\\+", " ", tooltip_data$tooltip)
    tooltip_data$tooltip <- gsub(": $", "", tooltip_data$tooltip)
    tooltip_data$tooltip <- gsub("^: ", "", tooltip_data$tooltip)
    rownames(tooltip_data) <- make.names(tooltip_data$ID, unique=TRUE)
    tooltip_data <- tooltip_data[-1]
    colnames(tooltip_data) <- c("short", "1.tooltip")
    tooltip_data <- tooltip_data[-1]
    return(tooltip_data)
}

#' Find how many times a given pattern occurs in every gene of a genome.
#'
#' @param fasta  a fasta genome
#' @param gff   an optional gff of annotations (if not provided it will just ask the whole genome.
#' @param pattern   what pattern to search for?  This was used for tnseq and TA is the mariner insertion point.
#' @param type  the column to get frmo the gff file
#' @param key   what type of entry of the gff file to key from?
#' @return num_pattern a data frame of names and numbers.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools} \link[Biostrings]{PDict} \link[Rsamtools]{FaFile}
#' @examples
#' \dontrun{
#' num_pattern = pattern_count_genome('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @export
pattern_count_genome <- function(fasta, gff=NULL, pattern='TA', type='gene', key='locus_tag') {
    rawseq <- Rsamtools::FaFile(fasta)
    if (is.null(gff)) {
        entry_sequences <- rawseq
    } else {
        entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
        ## type_entries <- subset(entries, type==type)
        type_entries <- entries[entries$type == type, ]
        names(type_entries) <- rownames(type_entries)
        entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
        names(entry_sequences) <- entry_sequences[[, key]]
    }
    dict <- Biostrings::PDict(pattern, max.mismatch=0)
    result <- Biostrings::vcountPDict(dict, entry_sequences)
    num_pattern <- data.frame("name"=names(entry_sequences),
                              "num"=as.data.frame(t(result)))
    return(num_pattern)
}

#' Gather some simple sequence attributes
#'
#' @param fasta  a fasta genome
#' @param gff   an optional gff of annotations (if not provided it will just ask the whole genome.
#' @param pattern   what pattern to search for?  This was used for tnseq and TA is the mariner insertion point.
#' @param type  the column to get frmo the gff file
#' @param key   what type of entry of the gff file to key from?
#' @return num_pattern a data frame of names and numbers.
#' @seealso \pkg{Biostrings} \pkg{Rsamtools} \link[Biostrings]{letterFrequency} \link[Rsamtools]{FaFile}
#' @examples
#' \dontrun{
#' num_pattern = sequence_attributes('mgas_5005.fasta', 'mgas_5005.gff')
#' }
#' @export
sequence_attributes <- function(fasta, gff=NULL, type='gene', key='locus_tag') {
    rawseq <- Rsamtools::FaFile(fasta)
    if (is.null(gff)) {
        entry_sequences <- rawseq
    } else {
        entries <- rtracklayer::import.gff3(gff, asRangedData=FALSE)
        ## type_entries <- subset(entries, type==type)
        type_entries <- entries[entries$type == type, ]
        ##names(type_entries) <- rownames(type_entries)
        entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
        names(entry_sequences) <- type_entries$Name
    }
    attribs <- data.frame(
        gc = Biostrings::letterFrequency(entry_sequences, "CG", as.prob=TRUE),
        at = Biostrings::letterFrequency(entry_sequences, "AT", as.prob=TRUE),
        gt = Biostrings::letterFrequency(entry_sequences, "GT", as.prob=TRUE),
        ac = Biostrings::letterFrequency(entry_sequences, "AC", as.prob=TRUE))
    rownames(attribs) <- type_entries$locus_tag
    colnames(attribs) <- c("gc","at","gt","ac")
    return(attribs)
}

#'   A stupid distance function of a point against two axes.
#'
#' @param firstterm  the x-values of the points.
#' @param secondterm  the y-values of the points.
#' @param firstaxis   the x-value of the vertical axis.
#' @param secondaxis   the y-value of the second axis.
#' @return dataframe of the distances
#' This just takes the abs(distances) of each point to the axes,
#' normalizes them against the largest point on the axes, multiplies
#' the result, and normalizes against the max of all points.
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#' mydist <- sillydist(df[,1], df[,2], first_median, second_median)
#' first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"), environment=hpgl_env) +
#'   ggplot2::xlab(paste("Expression of", df_x_axis)) +
#'   ggplot2::ylab(paste("Expression of", df_y_axis)) +
#'   ggplot2::geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
#'   ggplot2::geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
#'   ggplot2::geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
#'   ggplot2::geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
#'   ggplot2::geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
#'   ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
#'   ggplot2::geom_point(colour=grDevices::hsv(mydist$dist, 1, mydist$dist), alpha=0.6, size=size) +
#'   ggplot2::theme(legend.position="none")
#' first_vs_second  ## dots get colored according to how far they are from the medians
#' ## replace first_median, second_median with 0,0 for the axes
#' }
#' @export
sillydist <- function(firstterm, secondterm, firstaxis=0, secondaxis=0) {
    dataframe <- data.frame(firstterm, secondterm)
    dataframe$x <- (abs(dataframe[,1]) - abs(firstaxis)) / abs(firstaxis)
    dataframe$y <- abs((dataframe[,2] - secondaxis) / secondaxis)
    dataframe$x <- abs(dataframe[,1] / max(dataframe$x))
    dataframe$y <- abs(dataframe[,2] / max(dataframe$y))
    dataframe$dist <- abs(dataframe$x * dataframe$y)
    dataframe$dist <- dataframe$dist / max(dataframe$dist)
    return(dataframe)
}

#' Write a dataframe to an excel spreadsheet sheet.
#' I like to give folks data in any format they prefer, even though I sort
#' of hate excel.  Most people I work with use it, so therefore I do too.
#' This function has been through many iterations, first using XLConnect,
#' then xlsx, and now openxlsx.  Hopefully this will not change again.
#'
#' @param data  A data frame to print
#' @param sheet   Name of the sheet to write
#' @param file   The filename for the workbook.
#' @param overwrite_file   required for XLConnect, still used but perhaps not needed.
#' @param newsheet  same, but makes sure we don't overwrite an existing sheet
#' @param overwrite_sheet  yeah, I need to prune these options
#' @param dated   Append a date to the excel filename?
#' @param first_two_widths   I add long titles to the tops of the sheets
#'   setting this makes sure that those columns are not too wide
#' @param start_row   The first row of the sheet to write
#' @param start_col   The first column to write
#' @param ...  the set of arguments given to for openxlsx
#' @return a list containing the sheet and workbook written as well as the bottom-right coordinates of the last
#'   row/column written of the table.
#' @seealso \pkg{openxlsx} \link[openxlsx]{writeDataTable}
#' @examples
#' \dontrun{
#'  xls_coords <- write_xls(dataframe, sheet="hpgl_data")
#'  xls_coords <- write_xls(another_df, sheet="hpgl_data", start_row=xls_coords$end_col)
#' }
#' @export
write_xls <- function(data, sheet="first", file="excel/workbook.xlsx", overwrite_file=TRUE, newsheet=FALSE,
                      overwrite_sheet=TRUE, dated=TRUE, first_two_widths=c("30","60"),
                      start_row=1, start_col=1, ...) {
    arglist <- list(...)
    excel_dir <- dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    file <- gsub(pattern="\\.xlsx", replacement="", file, perl=TRUE)
    file <- gsub(pattern="\\.xls", replacement="", file, perl=TRUE)
    filename <- NULL
    if (isTRUE(dated)) {
        timestamp <- format(Sys.time(), "%Y%m%d%H")
        filename <- paste0(file, "-", timestamp, ".xlsx")
    } else {
        filename <- paste0(file, ".xlsx")
    }
    if (file.exists(filename)) {
        if (!isTRUE(newsheet) | !isTRUE(overwrite_file)) {
            backup_file(filename)
        }
    }

    if (class(data) == 'matrix') {
        data <- as.data.frame(data)
    }
    if (file.exists(filename)) {
        wb <- openxlsx::loadWorkbook(filename)
    } else {
        wb <- openxlsx::createWorkbook(creator="atb")
    }
    newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
    if (class(newsheet) == 'try-error') {
        ## assume for the moment that this is because it already exists.
        ## For the moment, leave it alone though and see what happens
        replace_sheet <- FALSE
        if (isTRUE(replace_sheet)) {
            openxlsx::removeWorksheet(wb=wb, sheet=sheet)
            try(openxlsx::addWorksheet(wb, sheetName=sheet))
        }
    }
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    new_row <- start_row
    new_col <- start_col
    ##print(paste0("GOT HERE openxlswrite, title? ", arglist$title))
    if (!is.null(arglist$title)) {
        openxlsx::writeData(wb, sheet, x=arglist$title, startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
    }

    ## I might have run into a bug in openxlsx, in WorkbookClass.R there is a call to is.nan() for a data.frame
    ## and it appears to me to be called oddly and causing problems
    ## I hacked the writeDataTable() function in openxlsx and sent a bug report.
    ## Another way to trip this up is for a column in the table to be of class 'list'
    for (col in colnames(data)) {
        ## data[[col]] <- as.character(data[[col]])
        ## print(paste0("TESTME: ", class(data[[col]])))
        if (class(data[[col]]) == 'list' | class(data[[col]]) == 'vector' | class(data[[col]]) == 'factor' | class(data[[col]]) == 'AsIs') {
            message(paste0("Converted ", col, " to characters."))
            data[[col]] <- as.character(data[[col]])
        }
    }
    openxlsx::writeDataTable(wb, sheet, x=data, tableStyle="TableStyleMedium9",
                             startRow=new_row, rowNames=TRUE, startCol=new_col)
    new_row <- new_row + nrow(data) + 2
    ## Going to make an assumption about columns 1,2
    ## Maybe make this a parameter? nah for now at least
    openxlsx::setColWidths(wb, sheet=sheet, widths=first_two_widths, cols=c(1,2))
    openxlsx::setColWidths(wb, sheet=sheet, widths="auto", cols=3:ncol(data))
    openxlsx::saveWorkbook(wb, filename, overwrite=overwrite_sheet)
    end_col <- ncol(data) + 1
    ret <- list(workbook=wb, file=filename, end_row=new_row, end_col=end_col)
    return(ret)
}

openxlsx_add_plot <- function(wb, plot) {
    ## Not implemented yet, my thought was to do a test for the worksheet
    ## if it is there, place the plot intelligently
    ## if not, create the worksheet and place the plot
    ## then make sure the workbook is saveable
}

#' Make a backup of an existing file with n revisions, like VMS!
#'
#' @param backup_file  the file to backup.
#' @param backups how many revisions?
backup_file <- function(backup_file, backups=4) {
    if (file.exists(backup_file)) {
        for (i in backups:01) {
            j <- i + 1
            i <- sprintf("%02d", i)
            j <- sprintf("%02d", j)
            test <- paste0(backup_file, ".", i)
            new <- paste0(backup_file, ".", j)
            if (file.exists(test)) {
                file.rename(test, new)
            }
        }
        newfile <- paste0(backup_file, ".", i)
        message(paste0("Renaming ", backup_file, " to ", newfile, "."))
        file.copy(backup_file, newfile)
    } else {
        message("The file does not yet exist.")
    }
}

#' Load a backup rdata file
#'
#' I often use R over a sshfs connection, sometimes with significant latency, and
#' I want to be able to save/load my R sessions relatively quickly.
#' Thus this function uses my backup directory to load its R environment.
#'
#' @param dir   the directory containing the RData.rda.xz file.
#' @return a bigger global environment
#' @seealso \link{load} \link{save}
#' @examples
#' \dontrun{
#' loadme()
#' saveme()
#' }
#' @export
loadme <- function(dir="savefiles") {
    savefile <- paste0(getwd(), "/", dir, "/RData.rda.xz")
    message(paste0("Loading the savefile: ", savefile))
    load_string <- paste0("load('", savefile, "', envir=globalenv())")
    message(paste0("Command run: ", load_string))
    eval(parse(text=load_string))
}

#' Make a backup rdata file for future reference
#'
#' I often use R over a sshfs connection, sometimes with significant latency, and
#' I want to be able to save/load my R sessions relatively quickly.
#' Thus this function uses pxz to compress the R session maximally and relatively fast.
#' This assumes you have pxz installed and >= 4 CPUs.
#'
#' @param directory  the directory to save the Rdata file.
#' @param backups   how many revisions?
#' @return the command used to save the global environment
#' @seealso \link{save} \link{pipe}
#' @examples
#' \dontrun{
#' saveme()
#' }
#' @export
saveme <- function(directory="savefiles", backups=4) {
    environment()
    if (!file.exists(directory)) {
        dir.create(directory)
    }
    savefile <- paste0(getwd(), "/", directory, "/RData.rda.xz")
    message(paste0("The savefile is: ", savefile))
    backup_file(savefile, backups=backups)
    ## The following save strings work:
    ## save_string <- paste0("save(list=ls(all.names=TRUE, envir=globalenv()), envir=globalenv(), file='", savefile, "')")
    ## save_string <- paste0("con <- base::pipe(paste0('pigz -p8 > ", savefile, "'), 'wb');\n save(list=ls(all.names=TRUE, envir=globalenv(), envir=globalenv(), file=con);\n close(con)")
    save_string <- paste0("con <- base::pipe(paste0('pxz -T4 > ", savefile, "'), 'wb');\n save(list=ls(all.names=TRUE, envir=globalenv()), envir=globalenv(), file=con, compress=FALSE);\n close(con)")
    message(paste0("The save string is: ", save_string))
    eval(parse(text=save_string))
}

#' Print a model as y = mx + b just like in grade school!
ymxb_print <- function(model) {
    intercept <- round(coefficients(model)[1], 2)
    x_name <- names(coefficients(model)[-1])
    slope <- round(coefficients(model)[-1], 2)
    ret <- paste0("y = ", slope, "*", x_name, " + ", intercept)
    message(ret)
    return(ret)
}


## EOF
