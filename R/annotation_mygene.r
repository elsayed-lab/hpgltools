#' Use mygene's queryMany to translate gene ID types
#'
#' Juggling between entrez, ensembl, etc can be quite a hassel.  This hopes to make it easier.
#'
#' Tested in test_40ann_biomart.R
#' This function really just sets a couple of hopefully helpful defaults.  When I first attempted
#' to use queryMany, it seemed to need much more intervention than it does now.  But at the least
#' this function should provide a reminder of this relatively fast and useful ID translation service.
#'
#' @param queries Gene IDs to translate.
#' @param email  E-mail address to send with the querymany query.
#' @param fields Set of fields to request, pass null for all.
#' @param species Human readable species for translation (Eg. 'human' instead of 'hsapiens'.)
#' @return Df of translated IDs/accessions
#' @seealso \pkg{mygene}
#'  \code{\link[mygene]{queryMany}}
#' @examples
#' \dontrun{
#'  data <- translate_ids_querymany(genes)
#' }
#' @export
translate_ids_querymany <- function(queries,
                                    ##from="ensembl",
                                    email="abelew@gmail.com",
                                    ##fields=c("uniprot", "ensembl.gene", "entrezgene", "go"),
                                    fields="all",
                                    species="all") {
    ##from_field <- from
    ##if (!is.null(from)) {
    ##    if (from == "ensembl") {
    ##        from_field <- "ensembl.gene"
    ##    } else if (from == "entrez") {
    ##        from_field <- "entrezgene"
    ##    }
    ##}

    ##one_way <- sm(mygene::queryMany(queries, scopes=from_field,
    ##                                fields=fields, species=species, returnall=TRUE))
    responses <- list()
    for (id in unique(queries)) {
        res <- mygene::query(id, species=species, email=email, fields=fields)$hits
        ## damn that is a nested data frame.
        new <- t(as.data.frame(unlist(res), stringsAsFactors=FALSE))
        rownames(new) <- id
        responses[[id]] <- as.data.frame(new, stringsAsFactors=FALSE)
        
        ##res <- query(id, species=species, email=email, fields=fields)$hits
        ##res[["id"]] <- id
        ##input <- list(
    }
    result <- data.table::rbindlist(responses, fill=TRUE, use.names=TRUE)
    ##response <- one_way[["response"]]
    return(result)
}

