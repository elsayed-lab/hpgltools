#' run the rGADEM suite
#'
#' This should provide a set of rGADEM results given an input file of sequences
#' and a genome.
#'
#' @param inputfile Fasta or bed file containing sequences to search.
#' @param genome BSgenome to read.
#' @param ... Parameters for plotting the gadem result.
#' @return A list containing slots for plots, the stdout output from gadem, the
#'   gadem result, set of occurences of motif, and the returned set of motifs.
#' @export
simple_gadem <- function(inputfile, genome="BSgenome.Hsapiens.UCSC.hs19", ...) {
  arglist <- list(...)
  ext <- tools::file_ext(inputfile)
  sequences <- NULL
  if (ext == "bed") {
    bed <- read.table(inputfile, header=FALSE, sep="\t")
    bed <- data.frame("chr" = as.factor(bed[, 1]),
                      "start" = as.numeric(bed[, 2]),
                      "end" = as.numeric(bed[, 3]))
    rg_bed <- IRanges::IRanges(start=bed[, 2], end=bed[, 3])
    sequences <- IRanges::RangedData(rg_bed, space=bed[, 1])
  } else if (ext == "fasta") {
    sequences <- Biostrings::readDNAStringSet(inputfile, "fasta")
  } else {
    stop("Unable to interpret files of type: ", ext)
  }

  gadem_and_prints <- sm(rGADEM::GADEM(sequences, verbose=1,
                                       pValue=0.05, eValue=10, genome=genome))
  gadem_result <- gadem_and_prints[["result"]]
  gadem_occurences <- rGADEM::nOccurrences(gadem_result)
  gadem_consensus <- rGADEM::consensus(gadem_occurences)
  gadem_motifs <- rGADEM::getPWM(gadem_consensus)
  plots <- list()
  count <- 1
  for (motif_name in names(gadem_motifs)) {
    motif_matrix <- gadem_motifs[[motif_name]]
    pwm <- seqLogo::makePWM(matrix)
    ## Stolen from seqLogo with a minor change for RNA
    seqLogo::seqLogo(pwm, ...)
    plot <- recordPlot()
    count <- count + 1
    plots[[count]] <- plot
  }
  retlist <- list(
    "plots" = plots,
    "log" = gadem_and_prints[["output"]],
    "result" = gadem_result,
    "occurences" = gadem_occurences,
    "consensus" = gadem_consensus,
    "motifs" = gadem_motifs)
  return(retlist)
}

simple_motifRG <- function(input_fasta, control_fasta, maximum=3,
                           title="Motifs of XXX", prefix="motif",
                           genome="BSgenome.Hsapiens.UCSC.hg19") {
  motifs <- motifRG::findMotifFasta(input_fasta, control_fasta,
                                    both.strand=TRUE,
                                    enriched=TRUE, mask=FALSE, start.width=4,
                                    min.cutoff=10, min.frac=0.001, max.motif=10,
                                    max.width=30, discretize=FALSE)
  table <- motifRG::motifLatexTable(main=title, motifs, prefix=prefix)
  return(table)
}

#' Extract sequence flanking a set of annotations (generally coding sequences)
#'
#' Given a set of annotations and genome, one might want to get the set of
#' adjacent sequences.
#'
#' @param bsgenome  Genome sequence
#' @param annotation  Set of annotations
#' @param distance  How far from each annotation is desired?
#' @param type  What type of annotation is desired?
#' @param prefix  Provide a prefix to the names to distinguish them from the
#'   existing annotations.
#' @return  A list of sequences before and after each sequence.
flanking_sequence <- function(bsgenome, annotation, distance=200,
                              type="gene", prefix="") {
  if (class(annotation) == "character") {
    ## Assume it is a filename to a gff file
    annotations <- load_gff_annotations(annotation, type=type)
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
  annotations <- merge(annotations, seqinfo, by.x="seqnames", by.y="row.names", all.x=TRUE)

  before <- GenomicRanges::GRanges(
                             seqnames=S4Vectors::Rle(annotations[, "seqnames"]),
                             ranges=IRanges::IRanges(ifelse(annotations[, "start"] <= distance,
                                                            1,
                                                            annotations[, "start"] - distance),
                                                     end=(annotations[, "start"] + 2)),
                             strand=S4Vectors::Rle(annotations[, "strand"]),
                             name=S4Vectors::Rle(annotations[, name_key]))

  after <- GenomicRanges::GRanges(
                            seqnames=S4Vectors::Rle(annotations[, "seqnames"]),
                            ranges=IRanges::IRanges(
                                              annotations[, "end"],
                                              end=ifelse(annotations[, "seqlengths"] <=
                                                         (annotations[, "end"] + distance),
                                                         annotations[, "seqlengths"],
                                                         annotations[, "end"] + distance)),
                            strand=S4Vectors::Rle(annotations[, "strand"]),
                            name=S4Vectors::Rle(annotations[, name_key]))
  before_seq <- Biostrings::getSeq(bsgenome, before)
  after_seq <- Biostrings::getSeq(bsgenome, after)
  retlist <- list("before" = before_seq,
                  "after" = after_seq)
  return(retlist)
}

## EOF
