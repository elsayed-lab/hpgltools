#' Skip the db and download all the text annotations for a given species.
#'
#' The microbesonline publicly available mysqldb is rather more complex than I
#' prefer.  This skips that process and just grabs a tsv copy of everything and
#' loads it into a dataframe.  I have not yet figured out how to so-easily query
#' microbesonline for species IDs, thus one will have to manually query the
#' database to find species of interest.
#'
#' Tested in test_70expt_spyogenes.R
#' There is so much awesome information in microbesonline, but damn is it annoying to download.
#' This function makes that rather easier, or so I hope at least.
#'
#' @param id Microbesonline ID to query.
#' @return Dataframe containing the annotation information.
#' @seealso \pkg{RCurl}
#'  \code{\link[RCurl]{getURL}}
#' @examples
#' \dontrun{
#'  annotations <- get_microbesonline_annotation(ids=c("160490","160491"))
#' }
#' @author atb
#' @export
load_microbesonline_annotations <- function(id="160490") {
  prelude_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id)
  result <- xml2::read_html(prelude_url)
  titles <- rvest::html_nodes(result, "title")
  species <- (titles %>% rvest::html_text())[1]
  message("The species being downloaded is: ", species)
  url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
  ##string <- RCurl::getURL(url)
  ##con <- textConnection(string)
  ##  data <- readr::read_table(con, sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
  data <- readr::read_tsv(url)
  return(data)
}

#' Download the various file formats from microbesoline.
#'
#' Microbesonline provides an interesting set of file formats to download.  Each
#' format proves useful under one condition or another, ergo this defaults to
#' iterating through them all and getting every file.
#'
#' @param id  Species ID to query.
#' @param type  File type(s) to download, if left null it will grab the genbank,
#'   tab, protein fasta, transcript fasta, and genome.
#' @return  List describing the files downloaded and their locations.
#' @author atb
download_microbesonline_files <- function(id="160490", type=NULL) {
  retlist <- list()
  gbk <- FALSE
  if (type == "gbk") {
    gbk <- TRUE
  }
  tab <- FALSE
  if (type == "tab") {
    tab <- TRUE
  }
  prot <- FALSE
  if (type == "prot") {
    prot <- TRUE
  }
  tx <- FALSE
  if (type == "tx") {
    tx <- TRUE
  }
  genome <- FALSE
  if (type == "genome") {
    genome <- TRUE
  }
  if (is.null(type)) {
    gbk <- TRUE
    tab <- TRUE
    prot <- TRUE
    tx <- TRUE
    genome <- TRUE
  }
  if (length(type) > 1) {
    for (t in type) {
      download_microbesonline_files(id=id, type=t)
    }
  }

  prelude_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id)
  result <- xml2::read_html(prelude_url)
  titles <- rvest::html_nodes(result, "title")
  species <- (titles %>% rvest::html_text())[1]

  if (isTRUE(gbk)) {
    gbk_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=gbk")
    gbk_file <- paste0(id, ".gbk")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", gbk_file, ".")
    gbk_downloaded <- download.file(gbk_url, gbk_file, quiet=TRUE)
    retlist[["gbk"]] <- gbk_file
  }

  if (isTRUE(tab)) {
    tab_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
    tab_file <- paste0(id, ".tab")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", tab_file, ".")
    tab_downloaded <- download.file(tab_url, tab_file, quiet=TRUE)
    retlist[["tab"]] <- tab_file
  }

  if (isTRUE(prot)) {
    prot_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=proteomes")
    prot_file <- paste0(id, "_proteome.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", prot_file, ".")
    prot_downloaded <- download.file(prot_url, prot_file)
    retlist[["prot"]] <- prot_file
  }

  if (isTRUE(tx)) {
    tx_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=transcriptomes")
    tx_file <- paste0(id, "_tx.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", tx_file, ".")
    tx_downloaded <- download.file(tx_url, tx_file)
    retlist[["tx"]] <- tx_file
  }

  if (isTRUE(genome)) {
    genome_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=genomes")
    genome_file <- paste0(id, "_genome.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", genome_file, ".")
    genome_downloaded <- download.file(genome_url, genome_file)
    retlist[["genome"]] <- genome_file
  }
  return(retlist)
}

#' Extract the set of GO categories by microbesonline locus
#'
#' The microbesonline is such a fantastic resource, it is a bit of a shame that
#' it is such a pain to query.
#'
#' Tested in test_42ann_microbes.R
#' I am not 100% certain that this is giving me the full correct set of gene
#' ontology accessions. At the very least, it does return a large number of
#' them, which is a start.
#'
#' @param id Which species to query.
#' @param id_column  This no longer uses MySQL, so which column from the html
#'   table to pull?
#' @param data_column  Similar to above, there are lots of places from which one might
#'   extract the data.
#' @param name  Allowing for non-specific searches by species name.
#' @return data frame of GO terms from www.microbesonline.org
#' @examples
#' \dontrun{
#'  go_df <- get_loci_go(id="160490")
#' }
#' @author atb
#' @export
load_microbesonline_go <- function(id="160490", id_column="name", data_column="GO", name=NULL) {
  chosen <- id
  table <- download_microbesonline_files(id=id, type="tab")
  table_df <- readr::read_tsv(file=table[["tab"]])
  go_df <- table_df[, c(id_column, data_column)] %>%
    tidyr::separate_rows(data_column, sep=",")
  keep_idx <- go_df[[data_column]] != ""
  go_df <- go_df[keep_idx, ]
  return(go_df)
}

## EOF
