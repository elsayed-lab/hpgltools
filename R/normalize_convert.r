## Time-stamp: <Sat May 14 13:41:08 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Perform a cpm/rpkm/whatever transformation of a count table.
#'
#' I should probably tell it to also handle a simple df/vector/list of gene lengths, but I
#' haven't. cp_seq_m is a cpm conversion of the data followed by a rp-ish conversion which
#' normalizes by the number of the given oligo.  By default this oligo is 'TA' because it was used
#' for tnseq which should be normalized by the number of possible transposition sites by mariner.
#' It could, however, be used to normalize by the number of methionines, for example -- if one
#' wanted to do such a thing.
#'
#' @param data Matrix of count data.
#' @param convert Type of conversion to perform: edgecpm/cpm/rpkm/cp_seq_m.
#' @param annotations Set of gff annotations are needed if using rpkm so we can get gene lengths.
#' @param fasta Fasta for rpkmish normalization.
#' @param pattern Cp_seq_m counts require a pattern to search on.
#' @param entry_type Used when reading a gff to acquire gene lengths.
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return Dataframe of cpm/rpkm/whatever(counts)
#' @seealso \pkg{edgeR} \pkg{Biobase} \code{\link[edgeR]{cpm}}
#' @examples
#' \dontrun{
#'  converted_table = convert_counts(count_table, convert='edgecpm')
#' }
#' @export
convert_counts <- function(data, convert="raw", annotations=NULL, fasta=NULL, pattern='TA', entry_type='gene', ...) {
    arglist <- list(...)
    data_class <- class(data)[1]
    if (data_class == 'expt') {
        count_table <- Biobase::exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        count_table <- Biobase::exprs(data)
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        count_table <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    if (convert == "edgecpm") {
        count_table <- edgeR::cpm(count_table)
    } else if (convert == "cpm") {
        lib_size <- colSums(count_table)
        ## count_table = t(t((count_table$counts + 0.5) / (lib_size + 1)) * 1e+06)
        transposed <- t(count_table + 0.5)
        cp_counts <- transposed / (lib_size + 1)
        cpm_counts <- t(cp_counts * 1e+06)
        count_table <- cpm_counts
    } else if (convert == "rpkm") {
        if (is.null(annotations)) {
            stop("RPKM conversion requires gene lengths.")
        }
        count_table <- hpgl_rpkm(count_table, annotations=annotations)
    } else if (convert == "cp_seq_m") {
        counts <- edgeR::cpm(count_table)
        ## count_table = divide_seq(counts, ...)
        count_table <- divide_seq(counts, fasta=fasta, gff=annotations, pattern=pattern, entry_type=entry_type)
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Express a data frame of counts as reads per pattern per
#' million(library).
#'
#' This uses a sequence pattern rather than length to normalize sequence.  It is essentially rpkm
#' but fancy pants.
#'
#' @param counts Read count matrix.
#' @param pattern Pattern to search against.  Defaults to 'TA'.
#' @param fasta Fasta genome to search.
#' @param gff Gff set of annotations to define start/ends of genes.
#' @param entry_type Type of gff entry to search against.  Defaults to 'gene'.
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return The 'RPseqM' counts
#' @seealso \code{\link[Rsamtools]{FaFile}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#' cptam <- divide_seq(cont_table, fasta="mgas_5005.fasta.xz", gff="mgas_5005.gff.xz")
#' }
#' @export
divide_seq <- function(counts, pattern="TA", fasta="testme.fasta", gff="testme.gff",
                       entry_type="gene", ...) {
    arglist <- list(...)
    if (!file.exists(fasta)) {
        compressed_fasta <- paste0(fasta, '.xz')
        system(paste0("xz -d ", compressed_fasta))
    }
    raw_seq <- try(Rsamtools::FaFile(fasta))
    if (class(raw_seq)[1] == 'try-error') {
        stop(paste0("There was a problem reading: ", fasta))
    }
    gff_entries <- gff2irange(gff)
    ## print(head(gff_entries))
    ##    cds_entries = subset(gff_entries, type==entry_type)
    found_entries <- (gff_entries$type == entry_type)
    if (sum(found_entries) == 0) {
        message(paste0("There were no found entries of type: ", entry_type, "."))
        message("Going to try locus_tag, and failing that, mRNA.")
        locus_entries <- (gff_entries$type == 'locus_tag')
        mrna_entries <- (gff_entries$type == 'mRNA')
        if ((sum(locus_entries) > sum(mrna_entries)) & sum(locus_entries) > 100) {
            found_entries <- locus_entries
        } else if (sum(mrna_entries) > 100) {
            found_entries <- mrna_entries
        } else {
            stop("Unable to find any entries of type locus_tag nor mrna.")
        }
    }
    ##cds_entries = subset(gff_entries, type==entry_type)
    cds_entries <- gff_entries[found_entries, ]
    names(cds_entries) <- make.names(cds_entries$locus_tag, unique=TRUE)
    cds_seq <- Biostrings::getSeq(raw_seq, cds_entries)
    names(cds_seq) <- cds_entries$locus_tag
    dict <- Biostrings::PDict(pattern, max.mismatch=0)
    result <- Biostrings::vcountPDict(dict, cds_seq)
    num_tas <- data.frame(name=names(cds_seq), tas=as.data.frame(t(result)))
    rownames(num_tas) <- make.names(num_tas$name, unique=TRUE)
    colnames(num_tas) <- c("name","TAs")
    num_tas$TAs <- num_tas$TAs + 1
    factor <- median(num_tas$TAs)
    num_tas$TAs <- num_tas$TAs / factor
    merged_tas <- merge(counts, num_tas, by="row.names", all.x=TRUE)
    rownames(merged_tas) <- merged_tas$Row.names
    merged_tas <- merged_tas[-1]
    ##merged_tas <- subset(merged_tas, select=-c("name"))  ## Two different ways of removing columns...
    merged_tas <- merged_tas[, -which(colnames(merged_tas) %in% c("name"))]
    merged_tas <- merged_tas / merged_tas$TAs
    merged_tas <- merged_tas[, !(colnames(merged_tas) %in% c("TAs"))]  ## Here is another!
    return(merged_tas)
}

#' Converts count matrix to log2 counts-per-million reads.
#'
#' Based on the method used by limma as described in the Law et al. (2014) voom
#' paper.
#'
#' @param counts Read count matrix.
#' @param lib.size Library size.
#' @return log2-CPM read count matrix.
#' @seealso \pkg{edgeR}
#' @examples
#' \dontrun{
#' l2cpm <- hpgl_log2cpm(counts)
#' }
#' @export
hpgl_log2cpm <- function(counts, lib.size=NULL) {
    if (is.null(lib.size)) {
        lib.size <- colSums(counts)
    }
    transposed_adjust <- t(counts + 0.5)
    cpm <- (transposed_adjust / (lib.size + 1)) * 1e+06
    l2cpm <- t(log2(cpm))
    return(l2cpm)
}

#' Reads/(kilobase(gene) * million reads)
#'
#' Express a data frame of counts as reads per kilobase(gene) per million(library). This function
#' wraps EdgeR's rpkm in an attempt to make sure that the required gene lengths get sent along.
#'
#' @param df Data frame of counts, alternately an edgeR DGEList.
#' @param annotations Contains gene lengths, defaulting to 'gene_annotations'.
#' @return Data frame of counts expressed as rpkm.
#' @seealso \pkg{edgeR} and \code{\link[edgeR]{cpm}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#' rpkm_df = hpgl_rpkm(df, annotations=gene_annotations)
#' }
#' @export
hpgl_rpkm <- function(df, annotations=get0('gene_annotations')) {
    if (class(df) == "edgeR") {
        df <- df$counts
    }
    df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ])
    if (dim(df_in)[1] == 0) {
        message("When the annotations and df were checked against each other
  the result was null.  Perhaps your annotation or df's rownames are not set?
  Going to attempt to use the column 'ID'.
")
        rownames(annotations) = make.names(annotations$ID, unique=TRUE)
        df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ])
        if (dim(df_in)[1] == 0) {
            stop("The ID column failed too.")
        }
    }
    colnames(df_in) <- colnames(df)
    merged_annotations <- merge(df, annotations, by="row.names")
    rownames(merged_annotations) <- merged_annotations[,"Row.names"]
    ##rownames(df_in) = merged_annotations[,"Row.names"]
    ## Sometimes I am stupid and call it length...
    lenvec <- NULL
    if (is.null(merged_annotations$width)) {
        lenvec <- as.vector(merged_annotations$length)
    } else {
        lenvec <- as.vector(merged_annotations$width)
    }
    names(lenvec) <- rownames(merged_annotations)
    requireNamespace("edgeR")
    rpkm_df <- edgeR::rpkm(df_in, gene.length=lenvec)
    colnames(rpkm_df) <- colnames(df)
    return(rpkm_df)
}

## EOF
