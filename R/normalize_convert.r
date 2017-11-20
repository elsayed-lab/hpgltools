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
#' @param ... Options I might pass from other functions are dropped into arglist, used by rpkm (gene
#'  lengths) and divide_seq (genome, pattern to match, and annotation type).
#' @return Dataframe of cpm/rpkm/whatever(counts)
#' @seealso \pkg{edgeR} \pkg{Biobase}
#'  \code{\link[edgeR]{cpm}}
#' @examples
#' \dontrun{
#'  converted_table = convert_counts(count_table, convert='cbcbcpm')
#' }
#' @export
convert_counts <- function(data, convert="raw", ...) {
  arglist <- list(...)
  data_class <- class(data)[1]
  annotations <- arglist[["annotations"]]
  if (data_class == "expt") {
    annotations <- fData(data)
    count_table <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    annotations <- fData(data)
    count_table <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    count_table <- as.data.frame(data)
  } else {
    stop("This function currently only types: expt, ExpressionSet, data.frame, and matrix.")
  }

  switchret <- switch(
    convert,
    "cpm" = {
      count_table <- edgeR::cpm(count_table)
    },
    "cbcbcpm" = {
      lib_size <- colSums(count_table)
      ## count_table = t(t((count_table$counts + 0.5) / (lib_size + 1)) * 1e+06)
      transposed <- t(count_table + 0.5)
      cp_counts <- transposed / (lib_size + 1)
      cpm_counts <- t(cp_counts * 1e+06)
      count_table <- cpm_counts
    },
    "rpkm" = {
      count_table <- hpgl_rpkm(count_table, annotations=annotations, ...)
    },
    "cp_seq_m" = {
      counts <- edgeR::cpm(count_table)
      count_table <- divide_seq(counts, annotations=annotations, ...)
      ##count_table <- divide_seq(counts, annotations=annotations, genome=genome)
    },
    {
      message(paste0("Not sure what to do with the method: ", convert))
    }
  ) ## End of the switch

  libsize <- colSums(count_table)
  counts <- list(
    "count_table" = count_table,
    "libsize" = libsize)
  return(counts)
}

#' Express a data frame of counts as reads per pattern per million.
#'
#' This uses a sequence pattern rather than length to normalize sequence.
#' It is essentially fancy pants rpkm.
#'
#' @param counts Read count matrix.
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return The RPseqM counts
#' @seealso \pkg{edgeR} \pkg{Rsamtools}
#'  \code{\link[Rsamtools]{FaFile}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#'  cptam <- divide_seq(cont_table, fasta="mgas_5005.fasta.xz", gff="mgas_5005.gff.xz")
#' }
#' @export
divide_seq <- function(counts, ...) {
  arglist <- list(...)
  annotations <- arglist[["annotations"]]
  genome <- arglist[["genome"]]
  pattern <- arglist[["pattern"]]

  ## The annotation data needs to have columns 'start', 'end', and 'chromosome'
  ## If they are named something else, I need to rename those columns.
  start_column <- "start"
  if (!is.null(arglist[["start_column"]])) {
    start_column <- arglist[["start_column"]]
    if (start_column != "start") {
      hit_idx <- colnames(annotations) == start_column
      colnames(annotations)[hit_idx] <- "start"
    }
  }
  end_column <- "end"
  if (!is.null(arglist[["end_column"]])) {
    end_column <- arglist[["end_column"]]
    if (end_column != "end") {
      hit_idx <- colnames(annotations) == end_column
      colnames(annotations)[hit_idx] <- "end"
    }
  }
  chromosome_column <- "chromosome"
  if (!is.null(arglist[["chromosome_column"]])) {
    chromosome_column <- arglist[["chromosome_column"]]
    if (chromosome_column != "chromosome") {
      hit_idx <- colnames(annotations) == chromosome_column
      colnames(annotations)[hit_idx] <- "chromosome"
    }
  }

  if (is.null(pattern)) {
    pattern <- "TA"
  }
  message(paste0("Using pattern: ", pattern, " instead of length for an rpkm-ish normalization."))

  compression <- NULL
  genome_class <- class(genome)[1]
  if (genome_class == "character") {
    ## This is presumably a fasta file, then.
    ## Sadly as of the last time I checked, FaFile doesn't handle compressed fasta
    if (grepl(pattern="gz$", x=genome)) {
      compression <- "gzip"
      system(paste0("gunzip ", genome))
    } else if (grepl(pattern="xz$", x=genome)) {
      compression <- "xz"
      system(paste0("xz -d ", genome))
    }
    raw_seq <- try(Rsamtools::FaFile(genome))
    if (class(raw_seq)[1] == "try-error") {
      stop(paste0("There was a problem reading: ", genome))
    }
    system(paste0(compression, " ", sub("^([^.]*).*", "\\1", genome)))
  } else if (genome_class == "BSgenome") {
    raw_seq <- genome
  } else {
    stop("Need a genome to search.")
  }

  ## Test that the annotations and genome have the same seqnames
  genome_seqnames <- sort(levels(as.factor(GenomicRanges::seqnames(genome))))
  annotation_seqnames <- sort(levels(as.factor(annotations[["chromosome"]])))
  hits <- sum(annotation_seqnames %in% genome_seqnames)
  if (hits == 0) {
    ## These are mislabeled (it seems the most common error is a chromosome names 'chr4' vs. '4'
    annotations[["chromosome"]] <- paste0("chr", annotations[["chromosome"]])
  } else if (hits < length(annotation_seqnames)) {
    warning("Not all the annotation sequences were found, this will probably end badly.")
  }

  annotation_class <- class(annotations)[1]
  annotation_entries <- NULL
  if (annotation_class == "character") {
    ## This is presumably a gff file, then
    annotation_entries <- gff2irange(annotations, ...)
  } else if (annotation_class == "data.frame") {
    colnames(annotations) <- tolower(colnames(annotations))
    annotations <- annotations[complete.cases(annotations), ]
    numberp <- sum(grepl(pattern="1", x=annotations[["strand"]]))
    if (numberp > 0) {
      annotations[["strand"]] <- as.numeric(annotations[["strand"]])
      annotations[["strand"]] <- ifelse(annotations[["strand"]] > 0, "+", "-")
    }
    ## Remove entries in annotations with start==NA
    na_idx <- is.na(sm(as.numeric(annotations[["start"]])))
    annotations <- annotations[!na_idx, ]
    annotation_entries <- GenomicRanges::makeGRangesFromDataFrame(annotations)
  } else if (annotation_class == "Granges") {
    annotations <- as.data.frame(annotations, stringsAsFactors=FALSE)
    colnames(annotations) <- tolower(colnames(annotations))
    annotation_entries <- annotations
  } else if (annotation_class == "orgDb") {
    ## TODO: Extract the annotation data frame
  } else {
    stop("Need some annotation information.")
  }

  cds_seq <- Biostrings::getSeq(raw_seq, annotation_entries)
  ## names(cds_seq) <- annotation_entries[[entry_type]]
  dict <- Biostrings::PDict(pattern, max.mismatch=0)
  result <- Biostrings::vcountPDict(dict, cds_seq)
  num_tas <- data.frame(name=names(cds_seq), tas=as.data.frame(t(result)))
  rownames(num_tas) <- make.names(num_tas[["name"]], unique=TRUE)
  colnames(num_tas) <- c("name", "pattern")
  num_tas[["pattern"]] <- num_tas[["pattern"]] + 1  ## No division by 0
  factor <- median(num_tas[["pattern"]])
  num_tas[["pattern"]] <- num_tas[["pattern"]] / factor
  merged_tas <- merge(counts, num_tas, by="row.names", all.x=TRUE)
  rownames(merged_tas) <- merged_tas[["Row.names"]]
  merged_tas <- merged_tas[-1]
  merged_tas <- merged_tas[, -which(colnames(merged_tas) %in% c("name"))]
  merged_tas <- merged_tas / merged_tas[["pattern"]]
  ##merged_tas <- subset(merged_tas, select=-c("name"))  ## Two different ways of removing columns...
  merged_tas <- merged_tas[, !(colnames(merged_tas) %in% c("pattern"))]  ## Here is another!
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
#'  l2cpm <- hpgl_log2cpm(counts)
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
#' @param ... extra options including annotations for defining gene lengths.
#' @return Data frame of counts expressed as rpkm.
#' @seealso \pkg{edgeR}
#'  \code{\link[edgeR]{cpm}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#'  rpkm_df = hpgl_rpkm(df, annotations=gene_annotations)
#' }
#' @export
hpgl_rpkm <- function(df, ...) {
  arglist <- list(...)
  annotations <- arglist[["annotations"]]
  ## holy crapola I wrote this when I had no clue what I was doing.
  if (class(df) == "edgeR") {
    df <- df[["counts"]]
  }
  df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ],
                         stringsAsFactors=FALSE)
  if (dim(df_in)[1] == 0) {
    message("When the annotations and df were checked against each other
  the result was null.  Perhaps your annotation or df's rownames are not set?
  Going to attempt to use the column 'ID'.
")
    rownames(annotations) <- make.names(annotations[["ID"]], unique=TRUE)
    df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ],
                           stringsAsFactors=FALSE)
    if (dim(df_in)[1] == 0) {
      stop("The ID column failed too.")
    }
  }
  colnames(df_in) <- colnames(df)
  df_in[["temporary_id_number"]] <- 1:nrow(df_in)
  merged_annotations <- merge(df_in, annotations, by="row.names", all.x=TRUE)
  rownames(merged_annotations) <- merged_annotations[, "Row.names"]
  merged_annotations <- merged_annotations[-1]
  merged_annotations <- merged_annotations[order(merged_annotations[["temporary_id_number"]]), ]
  merged_counts <- merged_annotations[, colnames(merged_annotations) %in% colnames(df) ]
  merged_annot <- merged_annotations[, colnames(merged_annotations) %in% colnames(annotations) ]

  ##rownames(df_in) = merged_annotations[,"Row.names"]
  ## Sometimes I am stupid and call it length...
  lenvec <- NULL
  if (!is.null(arglist[["column"]])) {
    lenvec <- as.vector(as.numeric(merged_annot[[arglist[["column"]]]]))
  } else if (is.null(merged_annot[["width"]])) {
    lenvec <- as.vector(as.numeric(merged_annot[["length"]]))
  } else {
    lenvec <- as.vector(as.numeric(merged_annot[["width"]]))
  }
  names(lenvec) <- rownames(merged_annot)
  tt <- sm(requireNamespace("edgeR"))
  rpkm_df <- edgeR::rpkm(as.matrix(merged_counts), gene.length=lenvec)
  colnames(rpkm_df) <- colnames(df)
  return(rpkm_df)
}

## EOF
