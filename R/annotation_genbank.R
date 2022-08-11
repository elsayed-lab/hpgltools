## annotation_genbank.r: A group of functions to extract annotations from
## GenBank. Probably the second most likely annotation source is NCBI, so this
## file seeks to simplify extracting annotations from genbank flat files and/or
## the NCBI web interface.

#' Given a genbank accession, make a txDb object along with sequences, etc.
#'
#' Let us admit it, sometimes biomart is a pain.  It also does not have easily
#' accessible data for microbes.  Genbank does!
#'
#' Tested in test_40ann_biomartgenbank.R and test_70expt_spyogenes.R
#' This primarily sets some defaults for the genbankr service in order to
#' facilitate downloading genomes from genbank and dumping them into a local
#' txdb instance.
#'
#' @param accession Accession to download and import.
#' @param file Use a file instead of downloading the accession?
#' @param reread Re-read (download) the file from genbank.
#' @param savetxdb Attempt saving a txdb object?
#' @return List containing a txDb, sequences, and some other stuff which I
#'  haven't yet finalized.
#' @seealso [Biostrings] [GenomicFeatures] [genbankr::import()] [genbankr::readGenBank()]
#' @examples
#'  sagalacticae_genbank_annot <- load_genbank_annotations(accession = "AE009948")
#'  dim(as.data.frame(sagalacticae_genbank_annot$cds))
#' @export
load_genbank_annotations <- function(accession = "AE009949", file = NULL, sequence = TRUE,
                                     reread = TRUE, savetxdb = FALSE) {
  gbk <- NULL
  if (!is.null(file)) {
    gbk <- genbankr::import(file)
  } else {
    input_file <- glue("{accession}.gb")
    if (!isTRUE(reread) & file.exists(input_file)) {
      gbk <- genbankr::import(input_file)
      ## The file exists, read it
    } else {
      ## Note that different versions of genbankr require somewhat different
      ## methods of querying here.
      mesg("Downloading accession from genbank.")
      gba <- genbankr::GBAccession(accession)
      mesg("Reading genbank file.")
      gbk <- genbankr::readGenBank(gba, partial = TRUE, verbose = TRUE, ret.seq = sequence)
    }
  }
  gbr <- try(genbankr::makeTxDbFromGenBank(gbk))
  seq <- genbankr::getSeq(gbk)
  others <- genbankr::otherFeatures(gbk)
  genes <- GenomicFeatures::genes(gbk)
  exons <- GenomicFeatures::exons(gbk)
  cds <- GenomicFeatures::cds(gbk)
  ret <- list(
      "others" = others,
      "exons" = exons,
      "cds" = cds,
      "genes" = genes,
      "txdb" = gbr,
      "seq" = seq)
  return(ret)
}

#' A genbank accession downloader scurrilously stolen from ape.
#'
#' This takes and downloads genbank accessions.
#'
#' Tested in test_40ann_biomartgenbank.R
#' In this function I stole the same functionality from the ape package and set
#' a few defaults so that it hopefully fails less often.
#'
#' @param accessions An accession -- actually a set of them.
#' @param write Write the files?  Otherwise return a list of the strings
#' @return A list containing the number of files downloaded and the character
#'   strings acquired.
#' @seealso [ape]
#' @examples
#'  written <- download_gbk(accessions = "AE009949")
#'  written$written_file
#' @author The ape authors with some modifications by atb.
#' @export
download_gbk <- function(accessions = "AE009949", write = TRUE) {
  num_acc <- length(accessions)
  nrequest <- num_acc %/% 400 + as.logical(num_acc %% 400)
  downloaded <- character(0)
  num_downloaded <- 0
  strings <- list()
  written_file <- NULL
  for (i in seq_len(nrequest)) {
    a <- (i - 1) * 400 + 1
    b <- 400 * i
    if (i == nrequest) {
      b <- num_acc
    }
    accession <- accessions[i]

    url <- paste0(
        "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
        paste(accessions[a:b], collapse = ","), "&rettype=gb&retmode=text&report=gbwithparts")

    dl_file <- glue("{accession}.gb")
    data <- try(download.file(url = url, destfile = dl_file, method = "wget", quiet = TRUE))
    scanned <- NULL
    if (class(data) != "try-error") {
      scanned <- try(scan(file = dl_file, what = "", sep = "\n", quiet = TRUE))
    }
    if (class(scanned) != "try-error") {
      downloaded <- c(downloaded, scanned)
      num_downloaded <- num_downloaded + 1
    }
    strings[[accession]] <- downloaded
    written_file <- glue("{accessions[a]}.gb")
    if (isTRUE(write)) {
      file_connection <- file(written_file)
      writeLines(downloaded, file_connection)
      close(file_connection)
    }
  } ## End of for loop
  retlist <- list(
      "written_file" = written_file,
      "num_success" = num_downloaded,
      "strings" = strings)
  return(retlist)
}

## EOF
