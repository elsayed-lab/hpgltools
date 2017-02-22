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
#' @param savetxdb  Save a txdb package from this? FIXME THIS DOES NOT WORK.
#' @return List containing a txDb, sequences, and some other stuff which I haven't yet finalized.
#' @seealso \pkg{genbankr} \pkg{rentrez}
#'  \code{\link[genbankr]{import}}
#' @examples
#' \dontrun{
#'  txdb_result <- gbk2txdb(accession="AE009948", savetxdb=TRUE)
#' }
#' @export
gbk2txdb <- function(accession="AE009949", savetxdb=FALSE) {
    gbk <- NULL
    input_file <- paste0(accession, ".gb")
    if (file.exists(input_file)) {
        gbk <- genbankr::import(input_file)
        ## The file exists, read it
    } else {
        tt <- sm(require.auto("rentrez"))
        gba <- genbankr::GBAccession(accession)
        gbk <- genbankr::readGenBank(gba, partial=FALSE)
    }
    gbr <- genbankr::makeTxDbFromGenBank(gbk)
    seq <- genbankr::getSeq(gbk)
    others <- genbankr::otherFeatures(gbk)
    genes <- genbankr::genes(gbk)
    exons <- genbankr::exons(gbk)
    intergenic <- genbankr::intergenic(gbk)
    cds <- genbankr::cds(gbk)
    if (isTRUE(savetxdb)) {
        message("The genbankr txdb objects are incomplete.  This does not work.")
    }
    ret <- list(
        "others" = others,
        "exons" = exons,
        "cds" = cds,
        "intergenic" = intergenic,
        "genes" = genes,
        "txdb" = gbr,
        "seq" = seq)
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
    chromosomes <- GenomeInfoDb::seqlevels(gbr)
    genes <- AnnotationDbi::keys(gbr)
    keytypes <- AnnotationDbi::keytypes(gbr)
    columns <- AnnotationDbi::columns(gbr)
    lengths <- sm(AnnotationDbi::select(gbr,
                                        columns=c("CDSNAME", "CDSCHROM", "CDSEND", "CDSSTART",
                                                  "CDSSTRAND", "CDSID", "TXNAME"),
                                        keys=genes, keytype="GENEID"))
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
    N <- length(accessions)
    nrequest <- N %/% 400 + as.logical(N%%400)
    downloaded <- character(0)
    num_downloaded <- 0
    strings <- list()
    for (i in 1:nrequest) {
        a <- (i - 1) * 400 + 1
        b <- 400 * i
        if (i == nrequest) {
            b <- N
        }

        URL <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                     paste(accessions[a:b], collapse = ","), "&rettype=gb&retmode=text&report=gbwithparts")

        dl_file <- paste0(accessions[1], ".gb")
        data <- try(download.file(url=URL, destfile=dl_file, method="wget", quiet=TRUE))
        scanned <- try(scan(file=dl_file, what="", sep="\n", quiet=TRUE))
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
