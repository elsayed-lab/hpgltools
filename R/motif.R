## motif.r: Try running some motif finding/analysis tools.

#' run the rGADEM suite
#'
#' This should provide a set of rGADEM results given an input file of sequences
#' and a genome.
#'
#' @param inputfile Fasta or bed file containing sequences to search.
#' @param genome BSgenome to read.
#' @param p pvalue cutoff
#' @param e evalue cutoff
#' @param ... Parameters for plotting the gadem result.
#' @return A list containing slots for plots, the stdout output from gadem, the
#'  gadem result, set of occurences of motif, and the returned set of motifs.
#' @seealso [IRanges] [Biostrings] [rGADEM]
#' @export
simple_gadem <- function(inputfile, genome = "BSgenome.Hsapiens.UCSC.hs19",
                         p = 0.1, e = 0.0, ...) {
  arglist <- list(...)
  ext <- tools::file_ext(inputfile)
  sequences <- NULL
  if (ext == "bed") {
    bed <- read.table(inputfile, header = FALSE, sep = "\t")
    bed <- data.frame("chr" = as.factor(bed[, 1]),
                      "start" = as.numeric(bed[, 2]),
                      "end" = as.numeric(bed[, 3]))
    ## RangedData has been deprecated in favor of GRanges
    ## The following is unlikely to work, but since I never use this
    ## I am not certain how much I care at the moment.
    rg_bed <- GenomicRanges::GRanges(start = bed[, 2], end = bed[, 3])
    sequences <- GenomicRanges::GRangesList(rg_bed, space = bed[, 1])
    ## Here are the previous invocations:
    ## rg_bed <- IRanges::IRanges(start = bed[, 2], end = bed[, 3])
    ## sequences <- IRanges::RangedData(rg_bed, space = bed[, 1])
  } else if (ext == "fasta") {
    sequences <- Biostrings::readDNAStringSet(inputfile, "fasta")
  } else {
    stop("Unable to interpret files of type: ", ext)
  }

  gadem_out <- rGADEM::GADEM(sequences, verbose = 1,
                             pValue = p, eValue = e, genome = genome)
  gadem_occurences <- rGADEM::nOccurrences(gadem_out)
  gadem_consensus <- rGADEM::consensus(gadem_out)
  gadem_motifs <- rGADEM::getPWM(gadem_out)
  plots <- list()
  pvals <- list()
  motif_names <- names(gadem_motifs)
  for (m in seq_along(gadem_motifs)) {
    name <- motif_names[m]
    motif_matrix <- gadem_motifs[[name]]
    pwm <- seqLogo::makePWM(as.matrix(motif_matrix))
    plotted <- try(seqLogo::seqLogo(pwm))
    if (class(plotted)[1] != "try-error") {
      plot <- recordPlot()
    } else {
      plot <- NULL
    }
    plots[[name]] <- plot
    pvals[[name]] <- gadem_out[[m]]@alignList[[1]]@pval
  }
  retlist <- list(
      "plots" = plots,
      "result" = gadem_out,
      "occurences" = gadem_occurences,
      "consensus" = gadem_consensus,
      "motifs" = gadem_motifs,
      "pvals" = pvals)
  return(retlist)
}

#' Run motifRG on a fasta file.
#'
#' @param input_fasta Input file.
#' @param control_fasta control file.
#' @param maximum 3
#' @param title Output image title.
#' @param prefix Prefix for the output files.
#' @param genome Package containing the full genome.
#' @seealso [motifRG]
simple_motifRG <- function(input_fasta, control_fasta, maximum = 3,
                           title = "Motifs of XXX", prefix = "motif",
                           genome = "BSgenome.Hsapiens.UCSC.hg19") {
  motifs <- motifRG::findMotifFasta(input_fasta, control_fasta,
                                    both.strand = TRUE,
                                    enriched = TRUE, mask = FALSE, start.width = 4,
                                    min.cutoff = 10, min.frac = 0.001, max.motif = 10,
                                    max.width = 30, discretize = FALSE)
  table <- motifRG::motifLatexTable(main = title, motifs, prefix = prefix)
  return(table)
}

#' Extract sequence flanking a set of annotations (generally coding sequences)
#'
#' Given a set of annotations and genome, one might want to get the set of
#' adjacent sequences.
#'
#' @param bsgenome Genome sequence
#' @param annotation Set of annotations
#' @param distance How far from each annotation is desired?
#' @param type What type of annotation is desired?
#' @param prefix Provide a prefix to the names to distinguish them from the
#'  existing annotations.
#' @return List of sequences before and after each sequence.
#' @seealso [load_gff_annotations()] [GenomicRanges] [IRanges]
flanking_sequence <- function(bsgenome, annotation, distance = 200,
                              type = "gene", prefix = "") {
  if (class(annotation) == "character") {
    ## Assume it is a filename to a gff file
    annotations <- load_gff_annotations(annotation, type = type)
    name_key <- "gene_id"
  } else if (class(annotation) == "data.frame") {
    annotations <- annotation
    name_key <- "tx_name"
  } else {
    ## Then assume it is a GenomicRanges from a TxDb or somesuch
    annotations <- as.data.frame(annotation)
    name_key <- "tx_name"
  }
  seqinfo <- as.data.frame(bsgenome@seqinfo)
  annotations <- merge(annotations, seqinfo, by.x = "seqnames", by.y = "row.names", all.x = TRUE)

  before <- GenomicRanges::GRanges(
                               seqnames = S4Vectors::Rle(annotations[, "seqnames"]),
                               ranges = IRanges::IRanges(ifelse(annotations[, "start"] <= distance,
                                                                1,
                                                                annotations[, "start"] - distance),
                                                         end=(annotations[, "start"] + 2)),
                               strand = S4Vectors::Rle(annotations[, "strand"]),
                               name = S4Vectors::Rle(annotations[, name_key]))

  after <- GenomicRanges::GRanges(
                              seqnames = S4Vectors::Rle(annotations[, "seqnames"]),
                              ranges = IRanges::IRanges(
                                                    annotations[, "end"],
                                                    end = ifelse(annotations[, "seqlengths"] <=
                                                                 (annotations[, "end"] + distance),
                                                                 annotations[, "seqlengths"],
                                                                 annotations[, "end"] + distance)),
                              strand = S4Vectors::Rle(annotations[, "strand"]),
                              name = S4Vectors::Rle(annotations[, name_key]))
  before_seq <- Biostrings::getSeq(bsgenome, before)
  after_seq <- Biostrings::getSeq(bsgenome, after)
  retlist <- list("before" = before_seq,
                  "after" = after_seq)
  return(retlist)
}

## EOF
