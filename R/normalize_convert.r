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
    if (is.null(annotations)) {
      annotations <- fData(data)
    }
    count_table <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    if (is.null(annotations)) {
      annotations <- fData(data)
    }
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
      ## count_table <- divide_seq(counts, annotations=annotations, genome=genome)
    },
    {
      message("Not sure what to do with the method: ", convert)
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

  if (is.null(pattern)) {
    pattern <- "TA"
  }
  message("Using pattern: ", pattern, " instead of length for an rpkm-ish normalization.")

  compression <- NULL
  genome_class <- class(genome)[1]
  raw_seq <- NULL
  genome_granges <- NULL
  if (genome_class == "character") {
    ## This is presumably a fasta file, then.
    ## Sadly as of the last time I checked, FaFile doesn't handle compressed fasta
    compression <- NULL
    if (grepl(pattern="gz$", x=genome)) {
      compression <- "gzip"
      system(paste0("gunzip ", genome))
    } else if (grepl(pattern="xz$", x=genome)) {
      compression <- "xz"
      system(paste0("xz -d ", genome))
    }

    raw_seq <- try(Rsamtools::FaFile(genome))
    genome_granges <- Rsamtools::scanFaIndex(raw_seq)

    if (class(raw_seq)[1] == "try-error") {
      stop(paste0("There was a problem reading: ", genome))
    }

    ## Recompress the file if it was compressed before.
    if (!is.null(compression)) {
      system(paste0(compression, " ", sub("^([^.]*).*", "\\1", genome)))
    }

  } else if (genome_class == "BSgenome") {
    raw_seq <- genome
  } else if (genome_class == "GRanges") {
    raw_seq <- genome
    genome_granges <- genome
  } else {
    stop("Need a genome to search.")
  }

  ## First make sure we have a valid annotation_df
  ## Then, using it, make annotation_gr.
  annotation_class <- class(annotations)[1]
  annotation_df <- data.frame()
  annotation_gr <- NULL
  if (annotation_class == "character") {
    ## This is presumably a gff file, then
    annotation_df <- load_gff_annotations(annotations)
  } else if (annotation_class == "data.frame") {
    annotation_df <- annotations
  } else if (annotation_class == "Granges") {
    annotation_df <- as.data.frame(annotations, stringsAsFactors=FALSE)
    annotation_gr <- annotations
  } else if (annotation_class == "orgDb") {
    ## TODO: Extract the annotation data frame
  } else {
    stop("Need some annotation information.")
  }

  ## Now we should have annotation_df and maybe annotation_gr.
  ## Spend this time sanitizing annotation_df
  ## The annotation data needs to have columns 'start', 'end', and 'chromosome'
  ## If they are named something else, I need to rename those columns.
  colnames(annotation_df) <- tolower(colnames(annotation_df))
  start_column <- "start"
  if (!is.null(arglist[["start_column"]])) {
    start_column <- arglist[["start_column"]]
    if (start_column != "start") {
      hit_idx <- colnames(annotation_df) == start_column
      colnames(annotation_df)[hit_idx] <- "start"
    }
  }
  end_column <- "end"
  if (!is.null(arglist[["end_column"]])) {
    end_column <- arglist[["end_column"]]
    if (end_column != "end") {
      hit_idx <- colnames(annotation_df) == end_column
      colnames(annotation_df)[hit_idx] <- "end"
    }
  }
  chromosome_column <- "chromosome"
  if (!is.null(arglist[["chromosome_column"]])) {
    chromosome_column <- arglist[["chromosome_column"]]
    if (chromosome_column != "chromosome") {
      hit_idx <- colnames(annotation_df) == chromosome_column
      colnames(annotation_df)[hit_idx] <- "chromosome"
    }
  }
  if (is.null(annotation_df[["chromosome"]]) & !is.null(annotation_df[["seqnames"]])) {
    annotation_df[["chromosome"]] <- annotation_df[["seqnames"]]
  }

  numberp <- sum(grepl(pattern="1", x=annotation_df[["strand"]]))
  if (numberp > 0) {
    annotation_df[["strand"]] <- as.numeric(annotation_df[["strand"]])
    annotation_df[["strand"]] <- ifelse(annotation_df[["strand"]] > 0, "+", "-")
  }
  ## Remove entries in annotations with start==NA
  na_idx <- is.na(sm(as.numeric(annotation_df[["start"]])))
  annotation_df <- annotation_df[!na_idx, ]

  ## We should have a sanitized annotation_df now.
  if (is.null(annotation_gr)) {
    annot_df <- annotation_df
    annotation_gr <- GenomicRanges::makeGRangesFromDataFrame(
                                           annotation_df,
                                           seqnames.field="chromosome")
  }

  ## Test that the annotations and genome have the same seqnames
  genome_seqnames <- NULL
  if (is.null(genome_granges)) {
    ## This should be true if the provided genome is a BSGenome.
    genome_seqnames <- sort(levels(as.factor(GenomicRanges::seqnames(raw_seq))))
  } else {
    genome_seqnames <- sort(levels(as.factor(GenomicRanges::seqnames(genome_granges))))
  }
  annotation_seqnames <- sort(levels(as.factor(annotation_df[["chromosome"]])))
  hits <- sum(annotation_seqnames %in% genome_seqnames)
  if (hits == 0) {
    ## These are mislabeled (it seems the most common error is a chromosome names 'chr4' vs. '4'
    new_levels <- c(GenomeInfoDb::seqlevels(annotation_gr),
                    paste0("chr", unique(GenomicRanges::seqnames(annotation_gr))))
    GenomeInfoDb::seqlevels(annotation_gr) <- new_levels
    GenomicRanges::seqnames(annotation_gr) <- factor(paste0(
                     "chr", GenomicRanges::seqnames(annotation_gr)), levels=new_levels)
  } else if (hits < length(annotation_seqnames)) {
    warning("Not all the annotation sequences were found, this will probably end badly.")
  }

  cds_seq <- Biostrings::getSeq(raw_seq, annotation_gr)
  names(cds_seq) <- rownames(annotation_df)
  ##names(cds_seq) <- annotation_entries[[entry_type]]
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
  merged_tas <- merged_tas[, -1]
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
#' @param count_table Data frame of counts, alternately an edgeR DGEList.
#' @param ... extra options including annotations for defining gene lengths.
#' @return Data frame of counts expressed as rpkm.
#' @seealso \pkg{edgeR}
#'  \code{\link[edgeR]{cpm}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#'  rpkm_df = hpgl_rpkm(df, annotations=gene_annotations)
#' }
#' @export
hpgl_rpkm <- function(count_table, ...) {
  arglist <- list(...)
  annotations <- arglist[["annotations"]]
  ## holy crapola I wrote this when I had no clue what I was doing.
  if (class(count_table) == "edgeR") {
    count_table <- count_table[["counts"]]
  }
  count_table_in <- as.data.frame(count_table[rownames(count_table) %in% rownames(annotations), ],
                         stringsAsFactors=FALSE)
  if (dim(count_table_in)[1] == 0) {
    message("When the annotations and count_table were checked against each other
  the result was null.  Perhaps your annotation or count_table's rownames are not set?
  Going to attempt to use the column 'ID'.
")
    rownames(annotations) <- make.names(annotations[["ID"]], unique=TRUE)
    count_table_in <- as.data.frame(count_table[rownames(count_table) %in% rownames(annotations), ],
                           stringsAsFactors=FALSE)
    if (dim(count_table_in)[1] == 0) {
      stop("The ID column failed too.")
    }
  }
  colnames(count_table_in) <- colnames(count_table)
  count_table_in[["temporary_id_number"]] <- 1:nrow(count_table_in)
  merged_annotations <- merge(count_table_in, annotations, by="row.names", all.x=TRUE)
  rownames(merged_annotations) <- merged_annotations[, "Row.names"]
  merged_annotations <- merged_annotations[-1]
  merged_annotations <- merged_annotations[order(merged_annotations[["temporary_id_number"]]), ]
  merged_counts <- merged_annotations[, colnames(merged_annotations) %in% colnames(count_table)]
  merged_annot <- merged_annotations[, colnames(merged_annotations) %in% colnames(annotations)]

  ##rownames(count_table_in) = merged_annotations[,"Row.names"]
  ## Sometimes I am stupid and call it length...
  lenvec <- NULL
  if (is.null(arglist[["column"]]) &
      is.null(merged_annot[["length"]]) &
      is.null(merged_annot[["width"]])) {
    message("There appears to be no gene length annotation data, here are the possible columns:")
    message(toString(colnames(annotations)))
    message("If one is appropriate, redo the function call with: column='good column'")
    stop("There appears to be no annotation data providing gene length.")
  }

  chosen_column <- "width"
  if (!is.null(arglist[["column"]])) {
    chosen_column <- arglist[["column"]]
  } else if (is.null(merged_annot[["width"]])) {
    chosen_column <- "length"
  }
  ## Keep in mind that I set missing material to 'undefined'
  ## So lets set those to NA now.
  undef_idx <- merged_annot[[chosen_column]] == "undefined"
  merged_annot[undef_idx, chosen_column] <- NA
  lenvec <- as.vector(as.numeric(merged_annot[[chosen_column]]))
  names(lenvec) <- rownames(merged_annot)
  tt <- sm(requireNamespace("edgeR"))
  rpkm_count_table <- edgeR::rpkm(as.matrix(merged_counts), gene.length=lenvec)
  colnames(rpkm_count_table) <- colnames(count_table)
  return(rpkm_count_table)
}

## EOF
