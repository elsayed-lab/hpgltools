#' Use the publicly available microbesonline mysql instance to get species ids.
#'
#' The microbesonline mysql instance is more complex than I like.  Their id system is reminiscent of
#' KEGG's and similarly annoying.  Though I haven't figured out how the tables interact, a query to
#' get ids is simple enough.
#'
#' Tested in test_42ann_microbes.R
#' This function sets the defaults required for getting a quick and dirty connection to the public
#' microbesonline database and returning the ids associated with a given name.
#'
#' @param name Text string containing some part of the species name of interest.
#' @param exact Use an exact species name?
#' @return Dataframe of ids and names.
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'
#' }
#' @export
get_microbesonline_ids <- function(name="Escherichia", exact=FALSE) {
  tt <- sm(requireNamespace("RMySQL"))
  db_driver <- DBI::dbDriver("MySQL")
  connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                               host="pub.microbesonline.org", dbname="genomics")
  query <- "SELECT taxonomyId, shortName FROM Taxonomy WHERE shortName "
  if (isTRUE(exact)) {
    query <- paste0(query, "= '", name, "'")
  } else {
    query <- paste0(query, "like '%", name, "%'")
  }
  message(query)
  result <- DBI::dbSendQuery(connection, query)
  result_df <- DBI::fetch(result, n=-1)
  clear <- try(DBI::dbClearResult(result))
  disconnect <- try(DBI::dbDisconnect(connection))
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  return(result_df)
}

#' Use the publicly available microbesonline mysql instance to get species name(s).
#'
#' The microbesonline mysql instance is more complex than I like.  Their id system is reminiscent of
#' KEGG's and similarly annoying.  Though I haven't figured out how the tables interact, a query to
#' get ids is simple enough.
#'
#' Tested in test_42ann_microbesonline.R
#' This is essentially covered in get_micrboesonline_ids(), but this works too.
#'
#' @param id Text string containing some part of the species name of interest.
#' @return Dataframe of ids and names.
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'  names <- get_microbesonline_name(id=316385)
#' }
#' @export
get_microbesonline_name <- function(id=316385, name=NULL) {
  chosen <- id
  if (!is.null(name)) {
    chosen <- get_microbesonline_ids(name=name)
    message(paste0("Retrieved ", nrow(chosen), " IDs, choosing the first, which is for ", chosen[1, 2], "."))
    chosen <- chosen[1,1]
  }
  tt <- sm(requireNamespace("RMySQL"))
  db_driver <- DBI::dbDriver("MySQL")
  connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                               host="pub.microbesonline.org", dbname="genomics")
  query <- paste0("SELECT shortName FROM Taxonomy WHERE taxonomyId = '", chosen, "'")
  result <- DBI::dbSendQuery(connection, query)
  result_df <- DBI::fetch(result, n=-1)
  clear <- DBI::dbClearResult(result)
  disconnect <- DBI::dbDisconnect(connection)
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  return(result_df)
}

#' Skip the db and download all the text annotations for a given species.
#'
#' Like I said, the microbesonline mysqldb is rather more complex than I prefer.  This shortcuts
#' that process and just grabs a tsv copy of everything and loads it into a dataframe.
#'
#' Tested in test_70expt_spyogenes.R
#' There is so much awesome information in microbesonline, but damn is it annoying to download.
#' This function makes that rather easier, or so I hope at least.
#'
#' @param ids List of ids to query.
#' @param species Species name(s) to use instead.
#' @return List of dataframes with the annotation information.
#' @seealso \pkg{RCurl}
#'  \code{\link[RCurl]{getURL}}
#' @examples
#' \dontrun{
#'  annotations <- get_microbesonline_annotation(ids=c("160490","160491"))
#' }
#' @export
load_microbesonline_annotations <- function(ids="160490", name=NULL) {
  retlist <- list()
  id_list <- list()
  if (is.null(ids) & is.null(name)) {
    stop("Either ids or species must be defined.")
  } else if (is.null(ids)) {
    for (spec in name) {
      id_df <- get_microbesonline_ids(spec)
      id_shortlist <- as.list(id_df[[1]])
      names(id_shortlist) <- id_df[[2]]
      id_list <- append(x=id_list, values=id_shortlist)
    }
  } else {
    for (id in ids) {
      idl <- as.list(id)
      current_name <- get_microbesonline_name(id)
      names(idl) <- current_name[1, 1]
      id_list <- append(x=id_list, values=idl)
    }
  }

  retlist <- list()
  print(id_list)
  for (t in 1:length(id_list)) {
    name <- names(id_list)[[t]]
    message(paste0("Querying microbesonline for: ", name, "."))
    id <- id_list[[t]]
    url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
    string <- RCurl::getURL(url)
    con <- textConnection(string)
    data <- read.csv(con, sep="\t", header=TRUE, row.names=NULL)
    retlist[[name]] <- data
  }
  return(retlist)
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
  message(query)
  result <- DBI::dbSendQuery(connection, query)
  result_df <- DBI::fetch(result, n=-1)
  clear <- DBI::dbClearResult(result)
  disconnect <- DBI::dbDisconnect(connection)
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  return(result_df)
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
#' @return data frame of GO terms from pub.microbesonline.org
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'  go_df <- get_loci_go(id="160490")
#' }
#' @export
load_microbesonline_go <- function(id="160490", name_type="ncbi_tag", name=NULL) {
  chosen <- id
  if (!is.null(name)) {
    chosen <- get_microbesonline_ids(name=name)
    message(paste0("Retrieved ", nrow(chosen), " IDs, choosing the first, which is for ", chosen[1, 2], "."))
    chosen <- chosen[1, 1]
  }
  tt <- sm(requireNamespace("RMySQL"))
  db_driver <- DBI::dbDriver("MySQL")
  connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                               host="pub.microbesonline.org", dbname="genomics")

  ## Here is the result of asking 'select * FROM SynonymType':
  ##| type | description               |
  ##|    0 | gene name                 |
  ##|    1 | NCBI locus tag            |
  ##|    2 | GI                        |
  ##|    4 | alternative locus tag     |
  ##|    3 | NCBI accession number     |
  ##|    5 | NCBI GeneID               |
  ##|    6 | deprecated NCBI locus tag |
  ##|    7 | other synonyms            |
  ##|    8 | IMG gene_oid              |
  type_number <- 0
  if (name_type == "deprecated") {
    type_number <- 6
  } else if (name_type == "ncbi_tag") {
    type_number <- 1
  } else if (name_type == "gi") {
    type_number <- 2
  } else if (name_type == "alternative") {
    type_number <- 4
  } else if (name_type == "ncbi_accession") {
    type_number <- 3
  } else if (name_type == "ncbi_geneid") {
    type_number <- 5
  } else if (name_type == "other") {
    type_number <- 7
  } else if (name_type == "img") {
    type_number <- 8
  }

  query <- paste0("SELECT B.name, L.locusId, G.goID, T.acc_synonym, A.acc FROM
     genomics.term A, genomics.Scaffold S, genomics.term_synonym T,
       genomics.Locus L, genomics.Locus2Go G, genomics.Synonym B
     where S.TaxonomyId = '", chosen, "' and S.isGenomic = 1 and S.scaffoldId = L.scaffoldId
       and G.locusId = L.locusId and B.locusId = L.locusId and B.type = '",
     type_number, "' and T.term_id = G.goID and A.id = G.goID")

  ## Adding suppressWarnings to stop stupidly unhelpful 'Unsigned INTEGER in col 1 imported as numeric'
  result <- suppressWarnings(DBI::dbSendQuery(connection, query))  
  result_df <- DBI::fetch(result, n=-1)
  result_df <- unique(result_df)
  clear <- DBI::dbClearResult(result)
  disconnect <- DBI::dbDisconnect(connection)
  if (class(clear) == "try-error" | class(disconnect) == "try-error") {
    warning("Did not disconnect cleanly.")
  }
  return(result_df)
}

load_microbesonline_kegg <- function(id="160490", name=NULL) {
  chosen <- id
  if (!is.null(name)) {
    chosen <- get_microbesonline_ids(name=name)
    message(paste0("Retrieved ", nrow(chosen), " IDs, choosing the first, which is for ", chosen[1, 2], "."))
    chosen <- chosen[1,1]
  }

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
  message(paste0("The abbreviation for ", taxonid, " is ", org, "."))
  genepaths <- get_kegg_genepaths(abbreviation=org)
  return(genepaths)
}

## EOF
