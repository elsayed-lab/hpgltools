#' Download the various file formats from microbesoline.
#'
#' Microbesonline provides an interesting set of file formats to download.  Each
#' format proves useful under one condition or another, ergo this defaults to
#' iterating through them all and getting every file.
#'
#' @param id Species ID to query.
#' @param type File type(s) to download, if left null it will grab the genbank,
#'  tab, protein fasta, transcript fasta, and genome.
#' @return List describing the files downloaded and their locations.
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

  prelude_url <- glue("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id}")
  result <- xml2::read_html(prelude_url)
  titles <- rvest::html_nodes(result, "title")
  species <- (titles %>% rvest::html_text())[1]

  if (isTRUE(gbk)) {
    gbk_url <- glue("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=gbk")
    gbk_file <- glue("{id}.gbk")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", gbk_file, ".")
    gbk_downloaded <- download.file(gbk_url, gbk_file, quiet=TRUE)
    retlist[["gbk"]] <- gbk_file
  }

  if (isTRUE(tab)) {
    tab_url <- glue("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=tab")
    tab_file <- glue("{id}.tab")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", tab_file, ".")
    tab_downloaded <- download.file(tab_url, tab_file, quiet=TRUE)
    retlist[["tab"]] <- tab_file
  }

  if (isTRUE(prot)) {
    prot_url <- glue(
      "http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=proteomes")
    prot_file <- glue("{id}_proteome.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", prot_file, ".")
    prot_downloaded <- download.file(prot_url, prot_file, quiet=TRUE)
    retlist[["prot"]] <- prot_file
  }

  if (isTRUE(tx)) {
    tx_url <- glue(
      "http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=transcriptomes")
    tx_file <- glue("{id}_tx.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", tx_file, ".")
    tx_downloaded <- download.file(tx_url, tx_file, quiet=TRUE)
    retlist[["tx"]] <- tx_file
  }

  if (isTRUE(genome)) {
    genome_url <- glue(
      "http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=genomes")
    genome_file <- glue("{id}_genome.fasta")
    message("The species being downloaded is: ", species,
            " and is being downloaded as ", genome_file, ".")
    genome_downloaded <- download.file(genome_url, genome_file, quiet=TRUE)
    retlist[["genome"]] <- genome_file
  }
  return(retlist)
}

#' Extract microbesonline taxon IDs without having to click on the weird boxes
#' at the top of the website.
#'
#' This should simplify getting material from microbesonline.
#'
#' @param species String to search the set of microbesonline taxa.
#' @return NULL or 1 or more taxon ids.
#' @examples
#'  coli_taxids <- get_microbesonline_taxid(species="coli S88")
#'  head(coli_taxids)
#' @export
get_microbesonline_taxid <- function(species="Acyrthosiphon pisum virus") {
  id_url <- "http://microbesonline.org/cgi-bin/fetchGenome2.cgi?taxId=g1&byFavorites=1"
  result <- xml2::read_html(id_url)
  id_nodes <- result %>%
    rvest::html_nodes("#GenomeList")
  . <- id_nodes  ## Shush, R CMD check
  id_links <- id_nodes %>%
    rvest::html_nodes("td:nth-child(1) a") %>%
    rvest::html_attr("href") %>%
    gsub(pattern="^.*?tId=([[:digit:]]+)$", replacement="\\1", x=.)
  id_table <- id_nodes %>%
    rvest::html_table(header=TRUE, fill=TRUE)
  id_df <- id_table[[1]]  ## Grab the first (only) element
  id_df[["tax_id"]] <- id_links
  foundp <- grep(pattern=species, x=id_df[["Genome"]])
  ret <- NULL
  if (length(foundp) == 0) {
    message("Did not find any entries with species: ", species)
  } else if (length(foundp) == 1) {
    result <- id_df[foundp, ]
    message("Found 1 entry.")
    message(result)
    ret <- result[["tax_id"]]
  } else {
    result <- id_df[foundp, ]
    message("Found multiple entries.")
    message(result)
    ret <- result[["tax_id"]]
    names(ret) <- result[["Genome"]]
  }
  return(ret)
}

#' Skip the db and download all the text annotations for a given species.
#'
#' The microbesonline publicly available mysqldb is rather more complex than I
#' prefer.  This skips that process and just grabs a tsv copy of everything and
#' loads it into a dataframe.  I have not yet figured out how to so-easily query
#' microbesonline for species IDs, thus one will have to manually query the
#' database to find species of interest.
#'
#' Tested in test_70expt_spyogenes.R
#' There is so much awesome information in microbesonline, but damn is it
#' annoying to download. This function makes that rather easier, or so I hope at
#' least.
#'
#' @param species Microbesonline species.
#' @param id Microbesonline ID to query.
#' @return Dataframe containing the annotation information.
#' @seealso \pkg{rvest}
#' @examples
#'  pa14_microbesonline_annot <- load_microbesonline_annotations(species="PA14")
#'  colnames(pa14_microbesonline_annot)
#' @export
load_microbesonline_annotations <- function(species=NULL, id=NULL) {
  if (is.null(id) & is.null(species)) {
    stop("This needs either a species or taxon id.")
  } else if (is.null(id)) {
    id <- get_microbesonline_taxid(species)
    if (is.null(id)) {
      stop("The species name did not find a taxon id.")
    } else if (length(id) > 1) {
      warning("There are more than 1 taxon id matches, arbitrarily choosing the first.")
      id <- id[1]
    }
  }
  prelude_url <- paste0("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id)
  result <- xml2::read_html(prelude_url)
  titles <- rvest::html_nodes(result, "title")
  species <- (titles %>% rvest::html_text())[1]
  message("The species being downloaded is: ", species)
  url <- glue::glue("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId={id};export=tab")
  message("Downloading: ", url)
  data <- sm(readr::read_tsv(url))
  return(data)
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
#' @param species Microbesonline species.
#' @param id Which species to query.
#' @param table_df Pre-existing data frame of annotations containing GO stuff.
#' @param id_column This no longer uses MySQL, so which column from the html
#'   table to pull?
#' @param data_column Similar to above, there are lots of places from which one might
#'   extract the data.
#' @param name Allowing for non-specific searches by species name.
#' @return data frame of GO terms from www.microbesonline.org
#' @examples
#'  pa14_microbesonline_go <- load_microbesonline_go(species="PA14")
#'  head(pa14_microbesonline_go)
#' @export
load_microbesonline_go <- function(id=NULL, species=NULL, table_df=NULL, id_column="name",
                                   data_column="GO", name=NULL) {
  if (is.null(id) & is.null(species)) {
    stop("This needs either a species or taxon id.")
  } else if (is.null(id)) {
    id <- get_microbesonline_taxid(species)
    if (is.null(id)) {
      stop("The species name did not find a taxon id.")
    } else if (length(id) > 1) {
      warning("There are more than 1 taxon id matches, arbitrarily choosing the first.")
      id <- id[1]
    }
  }
  chosen <- id
  if (is.null(table_df)) {
    table <- download_microbesonline_files(id=id, type="tab")
    table_df <- sm(readr::read_tsv(file=table[["tab"]]))
  }
  if (! id_column %in% colnames(table_df)) {
    message(id_column, " was not found in the table, here are the available columns: ",
            toString(colnames(table_df)))
    ## message(head(as.data.frame(table_df), n=2))
    stop()
  }
  go_df <- table_df[, c(id_column, data_column)] %>%
    tidyr::separate_rows(data_column, sep=",")
  keep_idx <- go_df[[data_column]] != ""
  go_df <- go_df[keep_idx, ]
  return(go_df)
}

## EOF
