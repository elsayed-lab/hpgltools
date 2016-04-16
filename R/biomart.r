#' Extract annotation information from biomart without having to remember the stupid biomart parameters
#'
#' @param species currently only hsapiens
#' @param overwrite overwrite a savefile if it exists
#' @param do_save create a savefile of the annotations
#' @examples
#' \dontrun{
#'  tt = get_biomarg_annotations()
#' }
#' @export
get_biomart_annotations <- function(species="hsapiens", overwrite=FALSE, do_save=FALSE) {
    savefile <- "biomart_annotations.rda"
    if (file.exists(savefile) & overwrite == FALSE) {
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=globalenv())")
        eval(parse(text=load_string))
    }
    dataset <- paste0(species, "_gene_ensembl")
    ##mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset)
    mart <- NULL
    mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
    ensembl <- biomaRt::useDataset(dataset, mart=mart)
    ## The following was stolen from Laura's logs for human annotations.
    ## To see possibilities for attributes, use head(listAttributes(ensembl), n=20L)
    biomart_annotations <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype"), mart=ensembl)
    colnames(biomart_annotations) <- c("ID", "hgnc_symbol", "Description", "Type")
    ## In order for the return from this function to work with other functions in this, the rownames must be set.
    rownames(biomart_annotations) <- make.names(biomart_annotations[, "ID"], unique=TRUE)
    message("Finished downloading ensembl annotations.")
    if (isTRUE(do_save)) {
        message(paste0("Saving annotations to ", savefile, "."))
        save(list=ls(c("biomart_annotations", "go_annotations"), envir=globalenv()), envir=globalenv(), file=savefile)
        message("Finished save().")
    }
    return(biomart_annotations)
}


get_biomart_ontology <- function(species="hsapiens", overwrite=FALSE, do_save=FALSE) {
    savefile <- "biomart_ontology.rda"
    if (file.exists(savefile) & overwrite == FALSE) {
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=globalenv())")
        eval(parse(text=load_string))
    }
    dataset <- paste0(species, "_gene_ensembl")
    ##mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset)
    mart <- NULL
    mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
    ensembl <- biomaRt::useDataset(dataset, mart=mart)
    ## The following was stolen from Laura's logs for human annotations.
    ## To see possibilities for attributes, use head(listAttributes(ensembl), n=20L)
    go_annotations <- biomaRt::getBM(attributes = c("ensembl_gene_id","go_id"), mart=ensembl)
    message(paste0("Finished downloading ensembl go annotations, saving to ", savefile, "."))
    if (isTRUE(do_save)) {
        save(list=ls(c("biomart_annotations", "go_annotations"), envir=globalenv()), envir=globalenv(), file=savefile)
    }
    colnames(go_annotations) <- c("ID","GO")
    return(go_annotations)
}

translate_ids_querymany <- function(queries, from="ensembl", to="entrez", species="human") {
    scopes <- "entrezgene"
    if (from == "ensembl") {
        from_field <- "ensembl.gene"
    } else if (from == "entrez") {
        from_field <- "entrezgene"
    }

    if (to == "entrez") {
        to <- "entrezgene"
    } else if (to == "ensembl") {
        to <- "ensembl.gene"
    }

    one_way <- mygene::queryMany(queries, scopes=from_field, fields=c("uniprot","ensembl.gene","entrezgene", "go"), species=species)
    print(head(one_way))
    queries <- as.data.frame(queries)
    ret <- merge(queries, one_way, by.x="queries", by.y="query", all.x=TRUE)
    return(ret)
}
