#' Use the publicly available microbesonline mysql instance to get species ids.
#'
#' The microbesonline mysql instance is more complex than I like.  Their id system is reminiscent of
#' KEGG's and similarly annoying.  Though I haven't figured out how the tables interact, a query to
#' get ids is simple enough.
#'
#' @param name Text string containing some part of the species name of interest.
#' @param exact Use an exact species name?
#' @return Dataframe of ids and names.
#' @export
get_microbesonline_ids <- function(name, exact=FALSE) {
    requireNamespace("RMySQL")
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
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
#' @param id Text string containing some part of the species name of interest.
#' @return Dataframe of ids and names.
#' @export
get_microbesonline_name <- function(id) {
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
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
#' @param ids List of ids to query.
#' @param species Species name(s) to use instead.
#' @return List of dataframes with the annotation information.
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
mdesc_table <- function(table="Locus2GO") {
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
    query <- paste0("DESCRIBE TABLE ", table)
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    clear <- DBI::dbClearResult(result)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

#' Extract the set of GO categories by microbesonline locus
#'
#' The microbesonline is such a fantastic resource, it is a bit of a shame that it is such a pain to query.
#'
#' @param taxonid Which species to query.
get_loci_go <- function(taxonid="160490") {
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
    ## cheese and crackers their database is entirely too complex and poorly documented.
    query <- paste0("SELECT L.locusId, G.goID, T.acc_synonym FROM genomics.Scaffold S, genomics.term_synonym T, genomics.Locus L, genomics.Locus2Go G where S.TaxonomyId = '",
                    taxonid, "' and S.isGenomic = 1 and S.scaffoldId = L.scaffoldId  and G.locusId = L.locusId and T.term_id = G.goID")
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    clear <- DBI::dbClearResult(result)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

## EOF
