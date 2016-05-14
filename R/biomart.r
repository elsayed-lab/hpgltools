## Time-stamp: <Sat May 14 12:47:22 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Extract annotation information from biomart.
#'
#' Biomart is an amazing resource of information, but using it is a bit annoying.  This function
#' hopes to alleviate some common headaches.
#'
#' @param species Choose a species.
#' @param overwrite Overwite an existing save file?
#' @param do_save Create a savefile of annotations for future runs?
#' @param host Ensembl hostname to use.
#' @param trymart Biomart has become a circular dependency, this makes me sad, now to list the
#'     marts, you need to have a mart loaded...
#' @param include_lengths Also perform a search on structural elements in the genome?
#' @return Df of some (by default) human annotations.
#' @examples
#' \dontrun{
#'  tt = get_biomart_annotations()
#' }
#' @export
get_biomart_annotations <- function(species="hsapiens", overwrite=FALSE, do_save=TRUE,
                                    host="dec2015.archive.ensembl.org",
                                    trymart="ENSEMBL_MART_ENSEMBL", include_lengths=TRUE) {
    savefile <- "biomart_annotations.rda"
    biomart_annotations <- NULL
    if (file.exists(savefile) & overwrite == FALSE) {
        fresh <- new.env()
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=fresh)")
        eval(parse(text=load_string))
        biomart_annotations <- fresh[["biomart_annotations"]]
        return(biomart_annotations)
    }
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
    dataset <- paste0(species, "_gene_ensembl")
    ensembl <- try(biomaRt::useDataset(dataset, mart=mart))
    if (class(ensembl) == 'try-error') {
        message(paste0("Unable to perform useDataset, perhaps the given dataset is incorrect: ", dataset, "."))
        datasets <- biomaRt::listDatasets(mart=mart)
        print(datasets)
        return(NULL)
    }
    ## The following was stolen from Laura's logs for human annotations.
    ## To see possibilities for attributes, use head(listAttributes(ensembl), n=20L)
    biomart_annotation <- biomaRt::getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol",
                                                    "description", "gene_biotype"),
                                       mart=ensembl)
    message("Finished downloading ensembl feature annotations.")
    biomart_annotations <- NULL
    if (isTRUE(include_lengths)) {
        biomart_structure <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "cds_length", "chromosome_name",
                                                         "strand", "start_position", "end_position"),
                                            mart=ensembl)
        message("Finished downloading ensembl structure annotations.")
        biomart_annotations <- merge(biomart_annotation, biomart_structure,
                                     by.x="ensembl_transcript_id", by.y="ensembl_transcript_id",
                                     all.x=TRUE)
        colnames(biomart_annotations) <- c("transcriptID", "geneID", "hgnc_symbol", "Description",
                                           "Type", "length", "chromosome", "strand", "start", "end")
    } else {
        biomart_annotations <- biomart_annotation
        colnames(biomart_annotations) <- c("geneID", "transcriptID", "hgnc_symbol", "Description",
                                           "Type")
    }
    rownames(biomart_annotations) <- make.names(biomart_annotations[, "transcriptID"], unique=TRUE)
    ## In order for the return from this function to work with other functions in this, the rownames must be set.

    if (isTRUE(do_save)) {
        message(paste0("Saving annotations to ", savefile, "."))
        save(list=ls(pattern="biomart_annotations"), file=savefile)
        message("Finished save().")
    }
    return(biomart_annotations)
}

#' Extract gene ontology information from biomart.
#'
#' I perceive that every time I go to acquire annotation data from biomart, they have changed
#' something important and made it more difficult for me to find what I want. I recently found the
#' *.archive.ensembl.org, and so this function uses that to try to keep things predictable, if not
#' consistent.
#'
#' @param species Species to query.
#' @param overwrite Overwrite existing savefile?
#' @param do_save Create a savefile of the annotations?
#' @param host Ensembl hostname to use.
#' @param trymart Default mart to try, newer marts use a different notation.
#' @param secondtry The newer mart name.
#' @return Df of geneIDs and GOIDs.
#' @seealso \link[biomaRt]{getBM}
#' @examples
#' \dontrun{
#'  tt = get_biomart_ontology()
#' }
#' @export
get_biomart_ontology <- function(species="hsapiens", overwrite=FALSE, do_save=FALSE,
                                 host="dec2015.archive.ensembl.org", trymart="ENSEMBL_MART_ENSEMBL",
                                 secondtry="_gene") {
    savefile <- "go_annotations.rda"
    secondtry <- paste0(species, secondtry)
    go_annotations <- NULL
    if (file.exists(savefile) & overwrite == FALSE) {
        fresh <- new.env()
        message("The biomart annotations file already exists, loading from it.")
        load_string <- paste0("load('", savefile, "', envir=fresh)")
        eval(parse(text=load_string))
        go_annotations <- fresh[["go_annotations"]]
        return(go_annotations)
    }
    dataset <- paste0(species, "_gene_ensembl")
    mart <- NULL
    mart <- try(biomaRt::useMart(biomart=trymart, host=host))
    if (class(mart) == "try-error") {
        message(paste0("Unable to perform useMart, perhaps the host/mart is incorrect: ", host, " ", trymart, "."))
        marts <- biomaRt::listMarts(host=host)
        mart_names <- as.character(marts[[1]])
        message(paste0("The available marts are: "))
        message(mart_names)
        message("Trying the first one.")
        mart <- biomaRt::useMart(biomart=marts[[1,1]], host=host)
    }
    ensembl <- biomaRt::useDataset(dataset, mart=mart)
    if (class(ensembl) == 'try-error') {
        message(paste0("Unable to perform useDataset, perhaps the given dataset is incorrect: ", dataset, "."))
        datasets <- biomaRt::listDatasets(mart=mart)
        print(datasets)
        return(NULL)
    }
    go_annotations <- biomaRt::getBM(attributes = c("ensembl_gene_id","go_id"), mart=ensembl)
    message(paste0("Finished downloading ensembl go annotations, saving to ", savefile, "."))
    if (isTRUE(do_save)) {
        message(paste0("Saving ontologies to ", savefile, "."))
        save(list=ls(pattern="biomart_annotations"), file=savefile)
        message("Finished save().")
    }
    colnames(go_annotations) <- c("ID","GO")
    return(go_annotations)
}

#' Use mygene's queryMany to translate gene ID types
#'
#' Juggling between entrez, ensembl, etc can be quite a hassel.  This hopes to make it easier.
#'
#' @param queries Gene IDs to translate.
#' @param from Database to translate IDs from.
#' @param to Database to translate IDs into.
#' @param species Human readable species for translation (Eg. 'human' instead of 'hsapiens'.)
#' @return Df of translated IDs/accessions
#' @seealso \link[mygene]{queryMany}
#' @examples
#' \dontrun{
#'  data <- translate_ids_querymany(genes)
#' }
#' @export
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

#' Use biomart to get orthologs between supported species.
#'
#' Biomart's function getLDS is incredibly powerful, but it makes me think very polite people are
#' going to start knocking on my door, and it fails weirdly pretty much always. This function
#' attempts to alleviate some of that frustration.
#'
#' @param gene_ids List of gene IDs to translate.
#' @param first_species Linnean species name for one species.
#' @param second_species Linnean species name for the second species.
#' @param host Ensembl server to query.
#' @param trymart Assumed mart name to use.
#' @return Df of orthologs.
#' @export
biomart_orthologs <- function(gene_ids, first_species='hsapiens', second_species='mmusculus',
                              host="dec2015.archive.ensembl.org", trymart="ENSEMBL_MART_ENSEMBL") {
    first_mart <- NULL
    first_mart <- try(biomaRt::useMart(biomart=trymart, host=host))
    if (class(first_mart) == 'try-error') {
        message(paste0("Unable to perform useMart, perhaps the host/mart is incorrect: ", host, " ", trymart, "."))
        first_marts <- biomaRt::listMarts(host=host)
        first_mart_names <- as.character(first_marts[[1]])
        message(paste0("The available first_marts are: "))
        message(first_mart_names)
        message("Trying the first one.")
        first_mart <- biomaRt::useMart(biomart=first_marts[[1,1]], host=host)
    }
    first_dataset <- paste0(first_species, "_gene_ensembl")
    first_ensembl <- try(biomaRt::useDataset(first_dataset, mart=first_mart))
    if (class(first_ensembl) == 'try-error') {
        message(paste0("Unable to perform useDataset, perhaps the given dataset is incorrect: ", first_ensembl, "."))
        datasets <- biomaRt::listDatasets(mart=first_mart)
        print(datasets)
        return(NULL)
    }

    second_mart <- NULL
    second_mart <- try(biomaR::useMart(biomart=trymart, host=host))
    if (class(second_mart) == 'try-error') {
        message(paste0("Unable to perform useMart, perhaps the host/mart is incorrect: ", host, " ", trymart, "."))
        second_marts <- biomaRt::listMarts(host=host)
        second_mart_names <- as.character(second_marts[[1]])
        message(paste0("The available second_marts are: "))
        message(second_mart_names)
        message("Trying the first one.")
        second_mart <- biomaRt::useMart(biomart=second_marts[[1,1]], host=host)
    }
    second_dataset <- paste0(second_species, "_gene_ensembl")
    second_ensembl <- try(biomaRt::useDataset(second_dataset, mart=second_mart))
    if (class(second_ensembl) == 'try-error') {
        message(paste0("Unable to perform useDataset, perhaps the given dataset is incorrect: ", second_ensembl, "."))
        datasets <- biomaRt::listDatasets(mart=second_mart)
        print(datasets)
        return(NULL)
    }

    linked_genes <- biomaRt::getLDS(attributes=c('ensembl_gene_id','hgnc_symbol'), values=gene_ids,
                                    mart=first_ensembl, attributesL=c('ensembl_gene_id'), martL=second_ensembl)

    return(linked_genes)
}

## EOF
