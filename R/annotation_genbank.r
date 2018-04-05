#' Given a genbank accession, make a txDb object along with sequences, etc.
#'
#' Let us admit it, sometimes biomart is a pain.  It also does not have easily accessible data for
#' microbes.  Genbank does!
#'
#' Tested in test_40ann_biomartgenbank.R and test_70expt_spyogenes.R
#' This just sets some defaults for the genbankr service in order to facilitate downloading
#' genomes and such from genbank and dumping them into a local txdb instance.
#'
#' @param accession Accession to download and import
#' @param reread  Re-read (download) the file from genbank
#' @param savetxdb  Save a txdb package from this? FIXME THIS DOES NOT WORK.
#' @return List containing a txDb, sequences, and some other stuff which I haven't yet finalized.
#' @seealso \pkg{genbankr} \pkg{rentrez}
#'  \code{\link[genbankr]{import}}
#' @examples
#' \dontrun{
#'  txdb_result <- gbk2txdb(accession="AE009948", savetxdb=TRUE)
#' }
#' @export
load_genbank_annotations <- function(accession="AE009949", reread=TRUE, savetxdb=FALSE) {
  gbk <- NULL
  input_file <- paste0(accession, ".gb")
  if (!isTRUE(reread) & file.exists(input_file)) {
    gbk <- genbankr::import(input_file)
    ## The file exists, read it
  } else {
    gba <- genbankr::GBAccession(accession)
    gbk <- genbankr::readGenBank(gba, partial=TRUE, verbose=TRUE)
    ## Something here is fubar, it falls down with:
    ## Error in split.default(text, fldnames) :
    ## group length is 0 but data length > 0
    ## gbk <- genbankr::readGenBank(file=input_file, partial=TRUE, verbose=TRUE)
    ## gbk <- parseGenBank(input_file, partial=TRUE, verbose=TRUE, ret.anno=TRUE, ret.seq=TRUE)
  }
  gbr <- genbankr::makeTxDbFromGenBank(gbk)
  seq <- Biostrings::getSeq(gbk)
  others <- genbankr::otherFeatures(gbk)
  genes <- GenomicFeatures::genes(gbk)
  exons <- GenomicFeatures::exons(gbk)
  ## intergenic <- genbankr:::intergenic(gbk)
  cds <- GenomicFeatures::cds(gbk)
  if (isTRUE(savetxdb)) {
    message("The genbankr txdb objects are incomplete.  This does not work.")
  }
  ret <- list(
    "others" = others,
    "exons" = exons,
    "cds" = cds,
  ##  "intergenic" = intergenic,
    "genes" = genes,
    "txdb" = gbr,
    "seq" = seq)
  return(ret)
}

load_uniprot_annotations <- function(id=NULL, species="Mycobacterium tuberculosis",
                                     keytype="GI_NUMBER*", chosen_columns=NULL) {
  if (is.null(id)) {
    result_df <- UniProt.ws::availableUniprotSpecies(pattern=species)
    if (nrow(result_df) == 0) {
      message("Unable to find an id corresponding to the provided species.")
      return(data.frame())
    } else if (nrow(result_df) > 1) {
      message("More than 1 species was returned, they follow; the first was chosen arbitrarily.")
      print(result_df)
      id <- result_df[1, 1]
    } else {
      message(paste0("Found 1 species, using its ID: ", result_df[1, 1], "."))
      id <- result_df[1, 1]
    }
  }

  mtb_data <- UniProt.ws::UniProt.ws(as.numeric(id))
  ## keytypes which return something useful: EGGNOG, EMBL/GENBANK/DDBJ, ENSEMBL GENOMES,
  ## ENSEMBL_GENOMES PROTEIN, ENSEMBL_GENOMES TRANSCRIPT, GI_NUMBER*
  ## yeah there are more, but I am tired of waiting for this stupid thing.
  possible_keys <- AnnotationDbi::keys(x=mtb_data, keytype=keytype)
  possible_columns <- AnnotationDbi::columns(x=mtb_data)
  if (is.null(chosen_columns)) {
    chosen_columns <- c("ENTREZ_GENE", "GO", "INTERPRO", "PATHWAY", "LENGTH", "EGGNOG")
  }
  ret <- sm(select(x=mtb_data, keytype=keytype, columns=chosen_columns, keys=possible_keys))
  return(ret)
}

#' Extract some useful information from a gbk imported as a txDb.
#'
#' Maybe this should get pulled into the previous function?
#'
#' Tested in test_40ann_biomartgenbank.R
#' This function should provide a quick reminder of how to use the AnnotationDbi select function
#' if it does nothing else.  It also (hopefully helpfully) returns a granges object containing
#' the essential information one might want for printing out a gff or whatever.
#'
#' @param gbr TxDb object to poke at.
#' @return Granges data
#' @seealso \pkg{AnnotationDbi} \pkg{GenomeInfoDb} \pkg{GenomicFeatures}
#'  \code{\link[AnnotationDbi]{select}}
#' @examples
#' \dontrun{
#'  annotations <- gbk_annotations("saureus_txdb")
#' }
#' @export
gbk_annotations <- function(gbr) {
  ## chromosomes <- GenomeInfoDb::seqlevels(gbr)
  genes <- AnnotationDbi::keys(gbr)
  ## keytypes <- AnnotationDbi::keytypes(gbr)
  ## columns <- AnnotationDbi::columns(gbr)
  lengths <- sm(AnnotationDbi::select(gbr,
                                      ## columns=columns,
                                      columns=c("CDSNAME", "CDSCHROM", "CDSEND", "CDSSTART",
                                                "CDSSTRAND", "CDSID", "TXNAME"),
                                      keys=genes, keytype="GENEID"))
  ## keys=genes, keytype=keytypes))
  lengths[["length"]] <- abs(lengths[["CDSSTART"]] - lengths[["CDSEND"]])
  granges <- GenomicFeatures::transcripts(gbr)
  return(granges)
}

#' A genbank accession downloader scurrilously stolen from ape.
#'
#' This takes and downloads genbank accessions.
#'
#' Tested in test_40ann_biomartgenbank.R
#' In this function I stole the same functionality from the ape package and set a few defaults
#' so that it hopefully fails less often.
#'
#' @param accessions An accession -- actually a set of them.
#' @param write  Write the files?  Otherwise return a list of the strings
#' @return A list containing the number of files downloaded and the character strings acquired.
#' @seealso \pkg{ape}
#' @examples
#' \dontrun{
#'  gbk_file <- download_gbk(accessions=c("AE009949","AE009948"))
#' }
#' @export
download_gbk <- function(accessions="AE009949", write=TRUE) {
  num_acc <- length(accessions)
  nrequest <- num_acc %/% 400 + as.logical(num_acc %% 400)
  downloaded <- character(0)
  num_downloaded <- 0
  strings <- list()
  for (i in 1:nrequest) {
    a <- (i - 1) * 400 + 1
    b <- 400 * i
    if (i == nrequest) {
      b <- num_acc
    }

    url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                  paste(accessions[a:b], collapse = ","), "&rettype=gb&retmode=text&report=gbwithparts")

    dl_file <- paste0(accessions[1], ".gb")
    data <- try(download.file(url=url, destfile=dl_file, method="wget", quiet=TRUE))
    scanned <- NULL
    if (class(data) != "try-error") {
      scanned <- try(scan(file=dl_file, what="", sep="\n", quiet=TRUE))
    }
    if (class(scanned) != "try-error") {
      downloaded <- c(downloaded, scanned)
      num_downloaded <- num_downloaded + 1
    }
    strings[i] <- downloaded
    if (isTRUE(write)) {
      file_connection <- file(paste0(accessions[a], ".gb"))
      writeLines(downloaded, file_connection)
      close(file_connection)
    }
  } ## End of for loop
  retlist <- list(
    "num_success" = num_downloaded,
    "strings" = strings)
  return(retlist)
}

## EOF
