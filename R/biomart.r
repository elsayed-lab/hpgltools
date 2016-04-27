## Time-stamp: <Mon Apr 25 17:02:22 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Extract annotation information from biomart without having to remember the stupid biomart parameters
#'
#' @param species currently only hsapiens
#' @param overwrite overwrite a savefile if it exists
#' @param do_save create a savefile of the annotations
#' @param host ensembl hostname to use
#' @param trymart biomart has become a circular dependency, this makes me sad, now to list the marts, you need to have a mart loaded...
#' @examples
#' \dontrun{
#'  tt = get_biomart_annotations()
#' }
#' @export
get_biomart_annotations <- function(species="hsapiens", overwrite=FALSE, do_save=FALSE, host="dec2015.archive.ensembl.org", trymart="ENSEMBL_MART_ENSEMBL") {
    savefile <- "biomart_annotations.rda"
    biomart_annotations <- NULL
    if (file.exists(savefile) & overwrite == FALSE) {
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=globalenv())")
        eval(parse(text=load_string))
        return(biomart_annotations)
    }
    dataset <- paste0(species, "_gene_ensembl")
    ##mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset)
    mart <- NULL
    mart <- try(biomaRt::useMart(biomart=trymart, host=host))
    if (class(mart) == 'try-error') {
        message(paste0("Unable to perform useMart, perhaps the host/mart is incorrect: ", host, " ", trymart, "."))
        marts <- biomaRt::listMarts(host=host)
        mart_names <- as.character(marts[[1]])
        message(paste0("The available marts are: "))
        message(mart_names)
        message("Trying the first one.")
        mart <- biomaRt::useMart(biomart=marts[[1,1]], host=host)
    }
    ensembl <- try(biomaRt::useDataset(dataset, mart=mart))
    if (class(ensembl) == 'try-error') {
        message(paste0("Unable to perform useDataset, perhaps the given dataset is incorrect: ", dataset, "."))
        datasets <- biomaRt::listDatasets(mart=mart)
        print(datasets)
        return(NULL)
    }
    ## The following was stolen from Laura's logs for human annotations.
    ## To see possibilities for attributes, use head(listAttributes(ensembl), n=20L)
    biomart_annotations <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype"), mart=ensembl)
    colnames(biomart_annotations) <- c("ID", "hgnc_symbol", "Description", "Type")
    ## In order for the return from this function to work with other functions in this, the rownames must be set.
    rownames(biomart_annotations) <- make.names(biomart_annotations[, "ID"], unique=TRUE)
    message("Finished downloading ensembl annotations.")
    if (isTRUE(do_save)) {
        message(paste0("Saving annotations to ", savefile, "."))
        save(list=ls(c("biomart_annotations"), envir=globalenv()), envir=globalenv(), file=savefile)
        message("Finished save().")
    }
    return(biomart_annotations)
}

#' Extract gene ontology information from biomart without having to remember the stupid biomart parameters
#'
#' @param species currently only hsapiens
#' @param overwrite overwrite a savefile if it exists
#' @param do_save create a savefile of the annotations
#' @param host ensembl hostname to use
#' @param trymart default mart to try
#' @examples
#' \dontrun{
#'  tt = get_biomart_ontology()
#' }
#' @export
get_biomart_ontology <- function(species="hsapiens", overwrite=FALSE, do_save=FALSE, host="dec2015.archive.ensembl.org", trymart="ENSEMBL_MART_ENSEMBL") {
    savefile <- "biomart_ontology.rda"
    if (file.exists(savefile) & overwrite == FALSE) {
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=globalenv())")
        eval(parse(text=load_string))
        return(go_annotations)
    }
    dataset <- paste0(species, "_gene_ensembl")
    ##mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset)
    mart <- NULL
    mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host=host)
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
