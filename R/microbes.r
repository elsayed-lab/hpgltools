get_microbesonline_ids <- function(name, exact=FALSE) {
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
    return(result_df)
}

get_microbesonline_annotation <- function(id="160490") {
    id="160490"
    url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
    string <- RCurl::getURL(url)
    con <- textConnection(string)
    data <- read.csv(con, sep="\t", header=TRUE)
    return(data)
}
