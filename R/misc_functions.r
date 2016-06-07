#' png() shortcut
#'
#' I hate remembering my options for png()
#'
#' @param file a filename to write
#' @return a png with height=width=9 inches and a high resolution
#' @export
pp <- function(file) {
    png(filename=file, width=9, height=9, units="in", res=180)
}

#' We want to catch *and* save both errors and warnings, and in the case of
#' a warning, also keep the computed result.
#'
#' This was taken from:http://r.789695.n4.nabble.com/How-to-catch-both-warnings-and-errors-td3073597.html
#' and http://tolstoy.newcastle.edu.au/R/help/04/06/0217.html
#'
#' @title tryCatch both warnings and errors
#' @param expr an expression to try
#' @return a list with 'value' and 'warning', where
#'  'value' may be an error caught.
#' @author Martin Maechler
tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
    }
    ret <- list(
        "value" = withCallingHandlers(tryCatch(expr,
                                               error = function(e) e),
                                      warning = w.handler),
        "warning" = W)
    return(ret)
}

#' Silence, peasant!
#'
#' Some libraries/functions just won't shut up.  Ergo, silence, peasant!
#' This function uses 2 invocations of capture.output and a try(silent=TRUE) to capture the strings
#' of the outputs from the given expression in 'output', and the messages in 'message'.  The result
#' of the expression goes into 'result.'  If there is an error in the expression, it is returned as
#' a try-error object which may therefore be inspected as needed.
#'
#' @param code Some code to shut up.
#' @return List of the output log, message log, and result of the expression.
#' @export
s_p <- function(code) {
    warnings <- NULL
    output_log <- NULL
    message_log <- NULL
    result <- NULL
    output_log <- capture.output(type="output", {
        message_log <- capture.output(type="message", {
            result <- try(code)
        })
    })
    retlist <- list(
        "output" = output_log,
        "message" = message_log,
        "warnings" = warnings,
        "result" = result)
    return(retlist)
}

sq <- function(code, ps=NULL) {
    warnings <- NULL
    output_log <- NULL
    message_log <- NULL
    tryCatch(result <- { output_log <<- capture.output(type="output", { code }) },
                       error = function(e) {
                           call <- conditionCall(e)
                           if (!is.null(call)) {
                               if (identical(call[[1L]], quote(doTryCatch))) {
                                   call <- sys.call(-4L)
                               }
                               dcall <- deparse(call)[1L]
                               prefix <- paste("Error in", dcall, ": ")
                               LONG <- 75L
                               msg <- conditionMessage(e)
                               sm <- strsplit(msg, "\n")[[1L]]
                               w <- 14L + nchar(dcall, type = "w") + nchar(sm[1L], type = "w")
                               if (is.na(w)) {
                                   w <- 14L + nchar(dcall, type = "b") + nchar(sm[1L], type = "b")
                               }
                               if (w > LONG) {
                                   prefix <- paste0(prefix, "\n  ")
                               }
                           } else {
                               prefix <- "Error : "
                               msg <- paste0(prefix, conditionMessage(e), "\n")
                               .Internal(seterrmessage(msg[1L]))
                               if (!silent && identical(getOption("show.error.messages"), TRUE)) {
                                   cat(msg, file = stderr())
                                   .Internal(printDeferredWarnings())
                               }
                               invisible(structure(msg, class = "try-error", condition = e))
                           }
                       },
                       warning = function(w) { message_log <- capture.output(type="message", { w } )},
                       finally = { ps })
    retlist <- list(
        "output" = output_log,
        "message" = message_log,
        "warnings" = warnings,
        "result" = result)
    return(retlist)
}

sr <- function(code) {
    warnings <- NULL
    output_log <- NULL
    message_log <- NULL
    result <- {
        output_log <- capture.output(type="output", {
            message_log <- capture.output(type="message", {
                code
            })
        })
    }
    retlist <- list(
        "output" = output_log,
        "message" = message_log,
        "result" = result)
}

#' Grab gene lengths from a gff file.
#'
#' This function attempts to be robust to the differences in output from importing gff2/gff3 files.
#' But it certainly isn't perfect.
#'
#' @param gff Gff file with (hopefully) IDs and widths.
#' @param type Annotation type to use (3rd column).
#' @param key Identifier in the 10th column of the gff file to use.
#' @return Data frame of gene IDs and widths.
#' @seealso \pkg{rtracklayer} \link[rtracklayer]{import.gff}
#' @examples
#' \dontrun{
#'  tt = get_genelengths('reference/fun.gff.gz')
#'  head(tt)
#' ##           ID width
#' ## 1   YAL069W   312
#' ## 2   YAL069W   315
#' ## 3   YAL069W     3
#' ## 4 YAL068W-A   252
#' ## 5 YAL068W-A   255
#' ## 6 YAL068W-A     3
#' }
#' @export
get_genelengths <- function(gff, type="gene", key="ID") {
    ret <- gff2df(gff)
    ret <- ret[ret[["type"]] == type, ]
    ret <- ret[, c(key, "width")]
    colnames(ret) <- c("ID", "width")
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
#' @param data Count tables of exons.
#' @param gff Gff filename.
#' @param annotdf Dataframe of annotations (probably from gff2df).
#' @param parent Column from the annotations with the gene names.
#' @param child Column from the annotations with the exon names.
#'
#' @return List of 2 data frames, counts and lengths by summed exons.
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
    rownames(tmp_data) <- tmp_data[["Row.names"]]
    tmp_data <- tmp_data[-1]
    ## Start out by summing the gene widths
    column <- aggregate(tmp_data[, "width"], by=list(Parent=tmp_data[, parent]), FUN=sum)
    new_data <- data.frame(column[["x"]])
    rownames(new_data) <- column[["Parent"]]
    colnames(new_data) <- c("width")

    for (c in 1:length(colnames(data))) {
        column_name <- colnames(data)[[c]]
        column <- aggregate(tmp_data[, column_name], by=list(Parent=tmp_data[, parent]), FUN=sum)
        rownames(column) <- column[["Parent"]]
        new_data <- cbind(new_data, column[["x"]])
    } ## End for loop

    width_df <- data.frame(new_data[["width"]])
    rownames(width_df) <- rownames(new_data)
    colnames(width_df) <- c("width")
    new_data <- new_data[-1]
    colnames(new_data) <- colnames(data)
    rownames(new_data) <- rownames(column)
    ret <- list(
        "width" = width_df,
        "counts" = new_data)
    return(ret)
}

#' Make a knitr report with some defaults set a priori.
#'
#' I keep forgetting to set appropriate options for knitr.  This tries to set them.
#'
#' @param name Name the document!
#' @param type Html or pdf reports?
#' @return Dated report file.
#' @seealso \pkg{knitr} \pkg{rmarkdown} \pkg{knitrBootstrap}
#' @export
make_report <- function(name="report", type='pdf') {
    knitr::opts_knit$set(
        progress = TRUE,
        verbose = TRUE,
        width = 90,
        echo = TRUE)
    knitr::opts_chunk$set(
        error = TRUE,
        fig.width = 8,
        fig.height = 8,
        dpi = 96)
    options(digits = 4,
            stringsAsFactors = FALSE,
            knitr.duplicate.label = "allow")
    ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
    set.seed(1)
    output_date <- format(Sys.time(), "%Y%m%d-%H%M")
    input_filename <- name
    ## In case I add .rmd on the end.
    input_filename <- gsub("\\.rmd", "", input_filename, perl=TRUE)
    input_filename <- gsub("\\.Rmd", "", input_filename, perl=TRUE)
    input_filename <- paste0(input_filename, ".Rmd")
    if (type == 'html') {
        output_filename <- paste0(name, "-", output_date, ".html")
        output_format <- 'html_document'
        rmarkdown::render(output_filename, output_format)
    } else {
        output_filename <- paste0(name, "-", output_date, ".pdf")
        output_format <- 'pdf_document'
    }
    message(paste0("About to run: render(input=", input_filename, ", output_file=",
                   output_filename, " and output_format=", output_format))
    result <- try(rmarkdown::render(
        "input" = input_filename,
        "output_file" = output_filename,
        "output_format" = output_format), silent=TRUE)

    return(result)
}

#' Implement the arescan function in R
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
#' @param x DNA/RNA StringSet containing the UTR sequences of interest
#' @param basal I dunno.
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
#'                                ranges=IRanges(utr_genes[,3], end=(utr_genes[,3] + 120)),
#'                                strand=Rle(utr_genes[,5]),
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

#' copy/paste the function from SeqTools and figure out where it falls on its ass.
#'
#' Yeah, I do not remember what I changed in this function.
#'
#' @param x Sequence object
#' @param min.length I dunno.
#' @param p.to.start P to start of course
#' @param p.to.end The p to end -- wtf who makes names like this?
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

#' Extract annotation information from a gff file into a df
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff files from various
#' sources: yeastgenome.org, microbesonline, tritrypdb, ucsc, ncbi. To my eyes, they all look like
#' reasonably good gff3 files, but some of them must be loaded with import.gff2, import.gff3, etc.
#' That is super annoying. Also, I pretty much always just do as.data.frame() when I get something
#' valid from rtracklayer, so this does that for me, I have another function which returns the
#' iranges etc.  This function wraps import.gff/import.gff3/import.gff2 calls in try() because
#' sometimes those functions fail in unpredictable ways.
#'
#' @param gff Gff filename.
#' @param type Subset the gff file for entries of a specific type.
#' @return Dataframe of the annotation information found in the gff file.
#' @seealso \pkg{rtracklayer} \link[rtracklayer]{import.gff} \link[rtracklayer]{import.gff2} \link[rtracklayer]{import.gff3}
#' @examples
#' \dontrun{
#' funkytown <- gff2df('reference/gff/saccharomyces_cerevsiae.gff.xz')
#' }
#' @export
gff2df <- function(gff, type=NULL) {
    ret <- NULL
    annotations <- NULL
    success <- FALSE
    ## First try gtf annotations

    attempts <- c("rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=TRUE)",
                  "rtracklayer::import.gff3(gff, sequenceRegionsAsSeqinfo=FALSE)",
                  "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=TRUE)",
                  "rtracklayer::import.gff2(gff, sequenceRegionsAsSeqinfo=TRUE)",
                  "rtracklayer::import.gff(gff, format='gtf')",
                  "rtracklayer::import.gff(gff)"
                  )
    for (att in 1:length(attempts)) {
        message(paste0("Trying attempt: ", attempts[[att]]))
        attempt <- attempts[[att]]
        eval_string <- paste0("annotations <- try(", attempt, ")")
        eval(parse(text=eval_string))
        if (class(annotations) == "try-error") {
            success <- FALSE
        } else {
            success <- TRUE
            message("Had a successful gff import with ", attempt)
            break
        }
    }
    ret <- NULL
    if (class(annotations)[[1]] == "GRanges") {
        ret <- GenomicRanges::as.data.frame(annotations)
    } else {
        stop("Unable to load gff file.")
    }
    if (!is.null(type)) {
        index <- ret[, "type"] == type
        ret <- ret[index, ]
    }
    message(paste0("Returning a df with ", ncol(ret), " columns and ", nrow(ret), " rows."))
    return(ret)
}

#' Extract annotation information from a gff file into an irange object.
#'
#' Try to make import.gff a little more robust; I acquire (hopefully) valid gff files from various
#' sources: yeastgenome.org, microbesonline, tritrypdb, ucsc, ncbi. To my eyes, they all look like
#' reasonably good gff3 files, but some of them must be loaded with import.gff2, import.gff3, etc.
#' That is super annoying. Also, I pretty much always just do as.data.frame() when I get something
#' valid from rtracklayer, so this does that for me, I have another function which returns the
#' iranges etc.  This function wraps import.gff/import.gff3/import.gff2 calls in try() because
#' sometimes those functions fail in unpredictable ways.
#'
#' This is essentially gff2df(), but returns data suitable for getSet()
#'
#' @param gff Gff filename.
#' @param type Subset to extract.
#'
#' @return Iranges! (useful for getSeq().)
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
#' Take covRob's robust correlation coefficient and add it to the set of correlations available when
#' one calls cor().
#'
#' @param df Data frame to test.
#' @param method Correlation method to use. Includes pearson, spearman, kendal, robust.
#' @param ... Other options to pass to stats::cor().
#' @return Some fun correlation statistics.
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

#' Create a simple df from a gff which contains tooltips usable in googleVis graphs.
#'
#' The tooltip column is a handy proxy for more thorough anontations information when it would
#' otherwise be too troublesome to acquire.
#'
#' @param annotations Either a gff file or annotation data frame (which likely came from a gff file.).
#' @param desc_col Gff column from which to gather data.
#' @param type Gff type to use as the master key.
#' @param id_col Which annotation column to cross reference against?
#' @param ... Extra arguments dropped into arglist.
#' @return Df of tooltip information or name of a gff file.
#' @seealso \pkg{googleVis} \link{gff2df}
#' @examples
#' \dontrun{
#' tooltips <- make_tooltips('reference/gff/saccharomyces_cerevisiae.gff.gz')
#' }
#' @export
make_tooltips <- function(annotations, desc_col='description', type="gene", id_col="ID", ...) {
    arglist <- list(...)
    tooltip_data <- NULL
    if (class(annotations) == 'character') {
        tooltip_data <- gff2df(gff=annotations, type=type)
    } else if (class(annotations) == 'data.frame') {
        tooltip_data <- annotations
    } else {
        stop("This requires either a filename or data frame.")
    }
    if (is.null(tooltip_data[[id_col]])) {
        tooltip_data[["ID"]] <- rownames(tooltip_data)
    }

    ## Attempt to use multiple columns if a c() was given
    tooltip_data[["tooltip"]] <- tooltip_data[[id_col]]
    for (col in desc_col) {
        if (is.null(tooltip_data[[col]])) {
            message(paste0("The column ", col, " is null, not using it."))
        } else {
            tooltip_data[["tooltip"]] <- paste0(tooltip_data[["tooltip"]], ": ", tooltip_data[[col]])
        }
    }
    tooltip_data[["tooltip"]] <- gsub("\\+", " ", tooltip_data[["tooltip"]])
    tooltip_data[["tooltip"]] <- gsub(": $", "", tooltip_data[["tooltip"]])
    tooltip_data[["tooltip"]] <- gsub("^: ", "", tooltip_data[["tooltip"]])
    rownames(tooltip_data) <- make.names(tooltip_data[[id_col]], unique=TRUE)
    ## Now remove extraneous columns
    tooltip_data <- tooltip_data[, c("tooltip"), drop=FALSE]
    colnames(tooltip_data) <- c("1.tooltip")
    return(tooltip_data)
}

#' Find how many times a given pattern occurs in every gene of a genome.
#'
#' There are times when knowing how many times a given string appears in a genome/CDS is helpful.
#' This function provides that information and is primarily used by cp_seq_m().
#'
#' @param fasta Genome sequence.
#' @param gff Gff of annotation information from which to acquire CDS (if not provided it will just
#'     query the entire genome).
#' @param pattern What to search for? This was used for tnseq and TA is the mariner insertion point.
#' @param type Column to use in the gff file.
#' @param key What type of entry of the gff file to key from?
#' @return Data frame of gene names and number of times the pattern appears/gene.
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
    num_pattern <- data.frame(
        "name" = names(entry_sequences),
        "num" = as.data.frame(t(result)))
    return(num_pattern)
}

#' Gather some simple sequence attributes.
#'
#' This extends the logic of the pattern searching in pattern_count_genome() to search on some other
#' attributes.
#'
#' @param fasta Genome encoded as a fasta file.
#' @param gff Optional gff of annotations (if not provided it will just ask the whole genome).
#' @param type Column of the gff file to use.
#' @param key What type of entry of the gff file to key from?
#' @return List of data frames containing gc/at/gt/ac contents.
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
        type_entries <- entries[entries[["type"]] == type, ]
        ##names(type_entries) <- rownames(type_entries)
        entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
        names(entry_sequences) <- type_entries[["Name"]]
    }
    attribs <- data.frame(
        "gc" = Biostrings::letterFrequency(entry_sequences, "CG", as.prob=TRUE),
        "at" = Biostrings::letterFrequency(entry_sequences, "AT", as.prob=TRUE),
        "gt" = Biostrings::letterFrequency(entry_sequences, "GT", as.prob=TRUE),
        "ac" = Biostrings::letterFrequency(entry_sequences, "AC", as.prob=TRUE))
    rownames(attribs) <- type_entries[["locus_tag"]]
    colnames(attribs) <- c("gc","at","gt","ac")
    return(attribs)
}

#' Calculate a simplistic distance function of a point against two axes.
#'
#' Sillydist provides a distance of any point vs. the axes of a plot.
#' This just takes the abs(distances) of each point to the axes,
#' normalizes them against the largest point on the axes, multiplies
#' the result, and normalizes against the max of all point.
#'
#' @param firstterm X-values of the points.
#' @param secondterm Y-values of the points.
#' @param firstaxis X-value of the vertical axis.
#' @param secondaxis Y-value of the second axis.
#' @return Dataframe of the distances.
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#' mydist <- sillydist(df[,1], df[,2], first_median, second_median)
#' first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"),
#'                                    environment=hpgl_env) +
#'   ggplot2::xlab(paste("Expression of", df_x_axis)) +
#'   ggplot2::ylab(paste("Expression of", df_y_axis)) +
#'   ggplot2::geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
#'   ggplot2::geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
#'   ggplot2::geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
#'   ggplot2::geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
#'   ggplot2::geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
#'   ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
#'   ggplot2::geom_point(colour=grDevices::hsv(mydist$dist, 1, mydist$dist),
#'                       alpha=0.6, size=size) +
#'   ggplot2::theme(legend.position="none")
#' first_vs_second  ## dots get colored according to how far they are from the medians
#' ## replace first_median, second_median with 0,0 for the axes
#' }
#' @export
sillydist <- function(firstterm, secondterm, firstaxis=0, secondaxis=0) {
    dataframe <- data.frame(firstterm, secondterm)
    dataframe[["x"]] <- (abs(dataframe[, 1]) - abs(firstaxis)) / abs(firstaxis)
    dataframe[["y"]] <- abs((dataframe[, 2] - secondaxis) / secondaxis)
    dataframe[["x"]] <- abs(dataframe[, 1] / max(dataframe[["x"]]))
    dataframe[["y"]] <- abs(dataframe[, 2] / max(dataframe[["y"]]))
    dataframe[["dist"]] <- abs(dataframe$x * dataframe[["y"]])
    dataframe[["dist"]] <- dataframe$dist / max(dataframe[["dist"]])
    return(dataframe)
}

#' Write a dataframe to an excel spreadsheet sheet.
#'
#' I like to give folks data in any format they prefer, even though I sort
#' of hate excel.  Most people I work with use it, so therefore I do too.
#' This function has been through many iterations, first using XLConnect,
#' then xlsx, and now openxlsx.  Hopefully this will not change again.
#'
#' @param data Data frame to print.
#' @param wb Workbook to which to write.
#' @param sheet Name of the sheet to write.
#' @param start_row First row of the sheet to write. Useful if writing multiple tables.
#' @param start_col First column to write.
#' @param ...  Set of extra arguments given to openxlsx.
#' @return List containing the sheet and workbook written as well as the bottom-right coordinates of
#'     the last row/column written to the worksheet.
#' @seealso \pkg{openxlsx} \link[openxlsx]{writeDataTable}
#' @examples
#' \dontrun{
#'  xls_coords <- write_xls(dataframe, sheet="hpgl_data")
#'  xls_coords <- write_xls(another_df, sheet="hpgl_data", start_row=xls_coords$end_col)
#' }
#' @export
write_xls <- function(data, wb=NULL, sheet="first",
                      start_row=1, start_col=1, ...) {
    arglist <- list(...)
    if (class(data) == 'matrix') {
        data <- as.data.frame(data)
    }
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    } else if (class(wb)[[1]] == "list") { ## In case the return from write_xls() was passed to write_xls()
        wb <- wb[["workbook"]]
    } else if (class(wb)[[1]] != "Workbook") {
        stop("A workbook was passed to this, but the format is not understood.")
    }

    newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
    if (class(newsheet) == 'try-error') {
        if (grepl(pattern="too long", x=newsheet[1])) {
            message("The sheet name was too long for Excel, truncating it by removing vowels.")
            sheet <- gsub(pattern="a|e|i|o|u", replacement="", x=sheet)
            newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet), silent=TRUE)
            if (class(newsheet) == 'try-error') {
                message("hmph, still too long, truncating to 30 characters.")
                sheet <- substr(sheet, start=1, stop=30)
                newsheet <- try(openxlsx::addWorksheet(wb, sheetName=sheet))
            }
        } else {
            message(paste0("The sheet ", sheet, " already exists, it will get overwritten"))
        }
    }
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    new_row <- start_row
    new_col <- start_col
    ##print(paste0("GOT HERE openxlswrite, title? ", arglist$title))
    if (!is.null(arglist$title)) {
        openxlsx::writeData(wb, sheet, x=arglist[["title"]], startRow=new_row)
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
        if (class(data[[col]]) == 'list' | class(data[[col]]) == 'vector' |
            class(data[[col]]) == 'factor' | class(data[[col]]) == 'AsIs') {
            message(paste0("Converted ", col, " to characters."))
            data[[col]] <- as.character(data[[col]])
        }
    }
    openxlsx::writeDataTable(wb, sheet, x=data, tableStyle="TableStyleMedium9",
                             startRow=new_row, rowNames=TRUE, startCol=new_col)
    new_row <- new_row + nrow(data) + 2
    openxlsx::setColWidths(wb, sheet=sheet, widths="auto", cols=1:ncol(data))
    end_col <- ncol(data) + 1
    ret <- list(
        "workbook" = wb,
        "end_row" = new_row,
        "end_col" = end_col)
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
#' Sometimes I just want to kick myself for overwriting important files and then I remember using
#' VMS and wish modern computers were a little more like it.
#'
#' @param backup_file Filename to backup.
#' @param backups How many revisions?
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
#' I often use R over a sshfs connection, sometimes with significant latency, and I want to be able
#' to save/load my R sessions relatively quickly. Thus this function uses my backup directory to
#' load its R environment.
#'
#' @param dir Directory containing the RData.rda.xz file.
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
#' @param directory Directory to save the Rdata file.
#' @param backups How many revisions?
#' @return Command string used to save the global environment.
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
#'
#' Because, why not!?
#'
#' @param model Model to print from glm/lm/robustbase.
#' @return a string representation of that model.
ymxb_print <- function(model) {
    intercept <- round(coefficients(model)[1], 2)
    x_name <- names(coefficients(model)[-1])
    slope <- round(coefficients(model)[-1], 2)
    ret <- paste0("y = ", slope, "*", x_name, " + ", intercept)
    message(ret)
    return(ret)
}

## EOF
