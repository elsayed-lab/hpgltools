#' Skip the db and download all the text annotations for a given species.
#'
#' Like I said, the microbesonline mysqldb is rather more complex than I prefer.  This shortcuts
#' that process and just grabs a tsv copy of everything and loads it into a dataframe.
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
#' @export
load_microbesonline_annotations <- function(id="160490") {
  prelude_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id)
  result <- xml2::read_html(prelude_url)
  titles <- rvest::html_nodes(result, "title")
  species <- (titles %>% rvest::html_text())[1]
  message("The species being downloaded is: ", species)
  url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
  string <- RCurl::getURL(url)
  con <- textConnection(string)
  data <- readr::read_csv(con, sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
  return(data)
}

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
    gbk_downloaded <- download.file(gbk_url, gbk_file)
    retlist[["gbk"]] <- gbk_file
  }

  if (isTRUE(tab)) {
    tab_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
    tab_file <- paste0(id, ".tab")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", tab_file, ".")
    tab_downloaded <- download.file(tab_url, tab_file)
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
#' The microbesonline is such a fantastic resource, it is a bit of a shame that it is such a pain
#' to query.
#'
#' Tested in test_42ann_microbes.R
#' I am not 100% certain that this is giving me the full correct set of gene ontology accessions.
#' At the very least, it does return a large number of them, which is a start.
#'
#' @param id Which species to query.
#' @param id_column  This no longer uses MySQL, so which column from the html table to pull?
#' @param data_column  Similar to above, there are lots of places from which one might
#'   extract the data.
#' @param name  Allowing for non-specific searches by species name.
#' @return data frame of GO terms from www.microbesonline.org
#' @examples
#' \dontrun{
#'  go_df <- get_loci_go(id="160490")
#' }
#' @export
load_microbesonline_go <- function(id="160490", id_column="name", data_column="GO", name=NULL) {
  chosen <- id
  table <- download_microbesonline_files(id=id, type="tab")
  table_df <- readr::read_csv(file=table[["tab"]], header=TRUE, sep="\t")
  go_df <- table_df[, c(id_column, data_column)] %>%
    tidyr::separate_rows(data_column, sep=",")
  keep_idx <- go_df[[data_column]] != ""
  go_df <- go_df[keep_idx, ]
  return(go_df)
}

#' Extract the set of KEGG categories by microbesonline locus
#'
#' The microbesonline is such a fantastic resource, it is a bit of a shame that it is such a pain
#' to query.
#'
#' Tested in test_42ann_microbes.R
#' I am not 100% certain that this is giving me the full correct set of gene ontology accessions.
#' At the very least, it does return a large number of them, which is a start.
#'
#' @param id Which species to query.
#' @return data frame of GO terms from pub.microbesonline.org
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'  go_df <- get_loci_go(id="160490")
#' }
#' @export
load_microbesonline_kegg <- function(id="160490") {
  chosen <- id

  tt <- sm(requireNamespace("RMySQL"))
  db_driver <- DBI::dbDriver("MySQL")
  connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                               host="pub.microbesonline.org", dbname="genomics")

  query <- paste0("SELECT * FROM KEGG2Taxonomy WHERE
     taxonomyId = '", id, "'")

  ## Adding suppressWarnings to stop stupidly unhelpful 'Unsigned INTEGER in col 1 imported as numeric'
  result <- suppressWarnings(DBI::dbSendQuery(connection, query))
  result_df <- DBI::fetch(result, n=-1)
  result_df <- unique(result_df)
  clear <- DBI::dbClearResult(result)
  disconnect <- DBI::dbDisconnect(connection)
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  org <- result_df[1, 1] ## Grab the identifier
  message("The abbreviation for ", id, " is ", org, ".")
  genepaths <- load_kegg_annotations(abbreviation=org)
  return(genepaths)
}

#' Get the description of a microbesonline genomics table
#'
#' This at least in theory is only used by get_microbesonline,  but if one needs a quick and dirty SQL query
#' it might prove useful.
#'
#' @param table  Choose a table to query.
#' @return Data frame describing the relevant table
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'  description <- mdesc_table(table="Locus2Go")
#' }
mdesc_table <- function(table="Locus2Go") {
  db_driver <- DBI::dbDriver("MySQL")
  connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                               host="pub.microbesonline.org", dbname="genomics")
  query <- paste0("DESCRIBE ", table)
  result <- DBI::dbSendQuery(connection, query)
  result_df <- DBI::fetch(result, n=-1)
  clear <- DBI::dbClearResult(result)
  disconnect <- DBI::dbDisconnect(connection)
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  return(result_df)
}

## EOF
