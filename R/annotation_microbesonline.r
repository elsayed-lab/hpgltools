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
#'  microbes_ids <- get_microbesonline_ids(name="Streptococcus")
#' }
#' @export
get_microbesonline_ids <- function(name="Escherichia", exact=FALSE) {
    requireNamespace("RMySQL")
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
    clear <- DBI::dbClearResult(result)
    disconnect <- DBI::dbDisconnect(connection)
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
get_microbesonline_name <- function(id=316385) {
    requireNamespace("RMySQL")
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                                 host="pub.microbesonline.org", dbname="genomics")
    query <- paste0("SELECT shortName FROM Taxonomy WHERE taxonomyId = '", id, "'")
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    clear <- DBI::dbClearResult(result)
    disconnect <- DBI::dbDisconnect(connection)
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
get_microbesonline_annotation <- function(ids="160490", species=NULL) {
    retlist <- list()
    id_list <- list()
    if (is.null(ids) & is.null(species)) {
        stop("Either ids or species must be defined.")
    } else if (is.null(ids)) {
        for (spec in species) {
            id_df <- get_microbesonline_ids(spec)
            id_shortlist <- as.list(id_df[[1]])
            names(id_shortlist) <- id_df[[2]]
            id_list <- append(x=id_list, values=id_df[[1]])
        }
    } else {
        for (id in ids) {
            idl <- as.list(id)
            current_name <- get_microbesonline_name(id)
            names(idl) <- current_name[1,1]
            id_list <- append(x=id_list, values=idl)
        }
    }

    retlist <- list()
    print(id_list)
    for (t in 1:length(id_list)) {
        name <- names(id_list)[[t]]
        id <- id_list[[t]]
        url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
        string <- RCurl::getURL(url)
        con <- textConnection(string)
        data <- read.csv(con, sep="\t", header=TRUE)
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
#' @param taxonid Which species to query.
#' @return data frame of GO terms from pub.microbesonline.org
#' @seealso \pkg{DBI}
#'  \code{\link[DBI]{dbSendQuery}} \code{\link[DBI]{fetch}}
#' @examples
#' \dontrun{
#'  go_df <- get_loci_go(taxonid="160490")
#' }
#' @export
get_loci_go <- function(taxonid="160490") {
    requireNamespace("RMySQL")
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest",
                                 host="pub.microbesonline.org", dbname="genomics")
    ## cheese and crackers their database is entirely too complex and poorly documented.
    query <- paste0("SELECT L.locusId, G.goID, T.acc_synonym, A.acc FROM
     genomics.term A, genomics.Scaffold S, genomics.term_synonym T, genomics.Locus L, genomics.Locus2Go G
     where S.TaxonomyId = '", taxonid, "' and S.isGenomic = 1 and S.scaffoldId = L.scaffoldId
       and G.locusId = L.locusId and T.term_id = G.goID and A.id = G.goID")
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    clear <- DBI::dbClearResult(result)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

## EOF
