#' Extract annotation information from biomart.
#'
#' Biomart is an amazing resource of information, but using it is a bit annoying.  This function
#' hopes to alleviate some common headaches.
#'
#' Tested in test_40ann_biomart.R
#' This goes to some lengths to find the relevant tables in biomart.  But biomart is incredibly
#' complex and one should carefully inspect the output if it fails to see if there are more
#' appropriate marts, datasets, and columns to download.
#'
#' @param species  Choose a species.
#' @param overwrite  Overwite an existing save file?
#' @param do_save  Create a savefile of annotations for future runs?
#' @param host  Ensembl hostname to use.
#' @param drop_haplotypes  Some chromosomes have stupid names because they are from non-standard
#'   haplotypes and they should go away.  Setting this to false stops that.
#' @param trymart  Biomart has become a circular dependency, this makes me sad, now to list the
#'  marts, you need to have a mart loaded.
#' @param trydataset  Choose the biomart dataset from which to query.
#' @param gene_requests  Set of columns to query for description-ish annotations.
#' @param length_requests  Set of columns to query for location-ish annotations.
#' @param include_lengths  Also perform a search on structural elements in the genome?
#' @return Df of some (by default) human annotations.
#' @seealso \pkg{biomaRt}
#'  \code{\link[biomaRt]{listDatasets}} \code{\link[biomaRt]{getBM}}
#' @examples
#' \dontrun{
#'  tt = get_biomart_annotations()
#' }
#' @export
load_biomart_annotations <- function(species="hsapiens", overwrite=FALSE, do_save=TRUE,
                                     host="dec2016.archive.ensembl.org", drop_haplotypes=TRUE,
                                     trymart="ENSEMBL_MART_ENSEMBL",
                                     trydataset=NULL,
                                     gene_requests=c("ensembl_gene_id",
                                                     "version",
                                                     "ensembl_transcript_id",
                                                     "transcript_version", "hgnc_symbol",
                                                     "description", "gene_biotype"),
                                     length_requests=c("ensembl_transcript_id",
                                                       "cds_length", "chromosome_name",
                                                       "strand", "start_position", "end_position"),
                                     include_lengths=TRUE) {
  savefile <- paste0(species, "_biomart_annotations.rda")
  biomart_annotations <- NULL
  if (file.exists(savefile) & overwrite == FALSE) {
    fresh <- new.env()
    message("The biomart annotations file already exists, loading from it.")
    load_string <- paste0("load('", savefile, "', envir=fresh)")
    eval(parse(text=load_string))
    biomart_annotations <- fresh[["biomart_annotations"]]
    retlist <- list(
      "annotation" = biomart_annotations,
      "mart" = "savefile",
      "host" = "savefile",
      "mart_name" = "savefile",
      "rows" = "savefile",
      "dataset" = "savefile"
    )
    return(retlist)
  }
  mart <- NULL
  used_mart <- NULL
  mart <- try(biomaRt::useMart(biomart=trymart, host=host), silent=TRUE)
  if (class(mart) == "try-error") {
    message("Unable to perform useMart, perhaps the host/mart is incorrect: ",
                   host, " ", trymart, ".")
    marts <- biomaRt::listMarts(host=host)
    mart_names <- as.character(marts[[1]])
    message("The available marts are: ")
    message(toString(mart_names))
    message("Trying the first one.")
    mart <- biomaRt::useMart(biomart=marts[[1, 1]], host=host)
    used_mart <- marts[[1, 1]]
  } else {
    used_mart <- trymart
  }
  chosen_dataset <- NULL
  dataset <- NULL
  if (!is.null(trydataset)) {
    dataset <- trydataset
  } else {
    dataset <- paste0(species, "_gene_ensembl")
  }
  second_dataset <- paste0(species, "_eg_gene")
  ensembl <- try(biomaRt::useDataset(dataset, mart=mart), silent=TRUE)
  if (class(ensembl) == "try-error") {
    ensembl <- try(biomaRt::useDataset(second_dataset, mart=mart), silent=TRUE)
    if (class(ensembl) == "try-error") {
      message("Unable to perform useDataset, perhaps the given dataset is incorrect: ", dataset, ".")
      datasets <- biomaRt::listDatasets(mart=mart)
      message(toString(datasets))
      return(NULL)
    } else {
      message("Successfully connected to the ", second_dataset, " database.")
      chosen_dataset <- second_dataset
    }
  } else {
    message("Successfully connected to the ", dataset, " database.")
    chosen_dataset <- dataset
  }
  ## The following was stolen from Laura's logs for human annotations.
  ## To see possibilities for attributes, use head(listAttributes(ensembl), n=20L)
  chosen_annotations <- c()
  available_attribs <- biomaRt::listAttributes(ensembl)[["name"]]
  found_attribs <- gene_requests %in% available_attribs
  if (length(gene_requests) != sum(found_attribs)) {
    message(strwrap(prefix=" ", initial="", "Some attributes in your request list were not in the
 ensembl database. At some point I will show them here."))
    gene_requests <- gene_requests[found_attribs]
  }
  gene_annotations <- biomaRt::getBM(attributes=gene_requests,
                                     mart=ensembl)
  chosen_annotations <- c(gene_requests)
  message("Finished downloading ensembl gene annotations.")
  biomart_annotations <- NULL
  if (isTRUE(include_lengths)) {
    found_attribs <- length_requests %in% available_attribs
    if (length(length_requests) != sum(found_attribs)) {
      message(strwrap(prefix=" ", initial="", "Some attributes in your request list were not in the
 ensembl database. At some point I will show them here."))
      length_requests <- length_requests[found_attribs]
    }
    structure_annotations <- biomaRt::getBM(attributes=length_requests,
                                            mart=ensembl)
    message("Finished downloading ensembl structure annotations.")
    tmp_annot <- data.table::as.data.table(gene_annotations)
    tmp_struct <- data.table::as.data.table(structure_annotations)
    biomart_annotations <- merge(tmp_annot, tmp_struct,
                                 by.x="ensembl_transcript_id", by.y="ensembl_transcript_id",
                                 all.x=TRUE)
    biomart_annotations <- as.data.frame(biomart_annotations)
    chosen_annotations <- c(chosen_annotations, length_requests)
  } else {
    ## Do not include the lengths
    biomart_annotations <- as.data.frame(biomart_annotations)
  }
  ## rownames(biomart_annotations) <- make.names(biomart_annotations[, "transcriptID"], unique=TRUE)
  ## It is not valid to arbitrarily set it to 'transcriptID' because we cannot guarantee that will
  ## be the column name, but I think we can safely assume it will be the 1st column.
  rownames(biomart_annotations) <- make.names(biomart_annotations[, 1], unique=TRUE)

  ## In order for the return from this function to work with other functions in
  ## this, the rownames must be set.

  ## Set strand to +/-
  if (!is.null(biomart_annotations[["strand"]])) {
    biomart_annotations[["strand"]] <- ifelse(biomart_annotations[["strand"]] == "1", "+", "-")
  }

  ## Steve has some excellent questions regarding multiple gene IDs for a single known locus:
  ## Example: hgnc_id:TNF which lies on chromosome 6 has 8 gene IDs, 7 of which appear to be
  ## associated with haplotype chromosomes with names like 'CHR_HSCHR6_MHC_DBB_CTG1' as opposed
  ## to the rather simpler chromosome name '6'.
  ## Thus, if one wishes to get rid of these putatively spurious annotations, we should be
  ## able to grep -v chromosomes with MHC in them.
  if (isTRUE(drop_haplotypes)) {
    message("Dropping haplotype chromosome annotations, set drop_haplotypes=FALSE if this is bad.")
    good_idx <- grepl(x=biomart_annotations[["chromosome_name"]], pattern="^[[:alnum:]]{1,2}$")
    biomart_annotations <- biomart_annotations[good_idx, ]
  }

  if (isTRUE(do_save)) {
    message("Saving annotations to ", savefile, ".")
    save(list=ls(pattern="biomart_annotations"), file=savefile)
    message("Finished save().")
  }
  retlist <- list(
    "annotation" = biomart_annotations,
    "mart" = ensembl,
    "host" = host,
    "mart_name" = used_mart,
    "rows" = chosen_annotations,
    "dataset" = chosen_dataset
  )
  return(retlist)
}

#' Extract gene ontology information from biomart.
#'
#' I perceive that every time I go to acquire annotation data from biomart, they have changed
#' something important and made it more difficult for me to find what I want. I recently found the
#' *.archive.ensembl.org, and so this function uses that to try to keep things predictable, if not
#' consistent.
#'
#' Tested in test_40ann_biomart.R
#' This function makes a couple of attempts to pick up the correct tables from biomart.  It is worth
#' noting that it uses the archive.ensembl host(s) because of changes in table organization after
#' December 2015 as well as an attempt to keep the annotation sets relatively consistent.
#'
#' @param species Species to query.
#' @param overwrite Overwrite existing savefile?
#' @param do_save Create a savefile of the annotations? (if not false, then a filename.)
#' @param host Ensembl hostname to use.
#' @param trymart Default mart to try, newer marts use a different notation.
#' @param secondtry The newer mart name.
#' @param dl_rows  List of rows from the final biomart object to download.
#' @param dl_rowsv2  A second list of potential rows.
#' @return Df of geneIDs and GOIDs.
#' @seealso \pkg{biomaRt}
#'  \code{\link[biomaRt]{listMarts}} \code{\link[biomaRt]{useDataset}} \code{\link[biomaRt]{getBM}}
#' @examples
#' \dontrun{
#'  tt = get_biomart_ontologies()
#' }
#' @export
load_biomart_go <- function(species="hsapiens", overwrite=FALSE, do_save=TRUE,
                            host="dec2015.archive.ensembl.org", trymart="ENSEMBL_MART_ENSEMBL",
                            secondtry="_gene", dl_rows=c("ensembl_gene_id", "go_accession"),
                            dl_rowsv2=c("ensembl_gene_id", "go_id")) {
  secondtry <- paste0(species, secondtry)

  savefile <- paste0(species, "_go_annotations.rda")
  if (!identical(FALSE, do_save)) {
    if (class(do_save) == "character") {
      savefile <- do_save
      do_save <- TRUE
    }
  }

  if (file.exists(savefile) & overwrite == FALSE) {
    fresh <- new.env()
    message("The biomart annotations file already exists, loading from it.")
    load_string <- paste0("load('", savefile, "', envir=fresh)")
    eval(parse(text=load_string))
    biomart_go <- fresh[["biomart_go"]]
    retlist <- list(
      "go" = biomart_go,
      "mart" = "savefile",
      "host" = "savefile",
      "mart_name" = "savefile",
      "rows" = "savefile",
      "dataset" = "savefile"
    )
    return(retlist)
  }

  used_mart <- NULL
  mart <- try(biomaRt::useMart(biomart=trymart, host=host), silent=TRUE)
  if (class(mart) == "try-error") {
    message("Unable to perform useMart, perhaps the host/mart is incorrect: ",
            host, " ", trymart, ".")
    marts <- biomaRt::listMarts(host=host)
    mart_names <- as.character(marts[[1]])
    message("The available marts are: ")
    message(mart_names)
    message("Trying the first one.")
    mart <- biomaRt::useMart(biomart=marts[[1, 1]], host=host)
    used_mart <- marts[[1, 1]]
  } else {
    used_mart <- trymart
  }

  chosen_dataset <- NULL
  dataset <- paste0(species, "_gene_ensembl")
  ensembl <- try(biomaRt::useDataset(dataset, mart=mart), silent=TRUE)
  if (class(ensembl) == "try-error") {
    message("Unable to perform useDataset, perhaps the given dataset is incorrect: ",
            dataset, ".")
    datasets <- biomaRt::listDatasets(mart=mart)
    try_again <- paste0(species, "_eg_gene")
    message("Trying instead to use the dataset: ", try_again)
    ensembl <- biomaRt::useDataset(try_again, mart=mart)
    if (class(ensembl)[[1]] == "Mart") {
      message("That seems to have worked, extracting the resulting annotations.")
      chosen_dataset <- try_again
    } else {
      message("The second attempt failed as well, the following are the available datasets:")
      message(toString(datasets))
      return(NULL)
    }
  } else {
    chosen_dataset <- dataset
  }

  biomart_go <- try(biomaRt::getBM(attributes=dl_rows, mart=ensembl), silent=TRUE)
  if (class(biomart_go) == "try-error") {
    biomart_go <- try(biomaRt::getBM(attributes=dl_rowsv2, mart=ensembl), silent=TRUE)
    dl_rows <- dl_rowsv2
  }
  if (class(biomart_go) == "try-error") {
    message("Unable to download annotation data.")
    return(NULL)
  }
  message("Finished downloading ensembl go annotations, saving to ", savefile, ".")

  if (length(colnames(biomart_go)) == 2) {
    colnames(biomart_go) <- c("ID", "GO")
  }
  if (isTRUE(do_save)) {
    message("Saving ontologies to ", savefile, ".")
    save(list=ls(pattern="biomart_go"), file=savefile)
    message("Finished save().")
  }

  retlist <- list(
    "go" = biomart_go,
    "mart" = ensembl,
    "host" = host,
    "mart_name" = used_mart,
    "rows" = dl_rows,
    "dataset" = chosen_dataset
  )
  return(biomart_go)
}

#' Use biomart to get orthologs between supported species.
#'
#' Biomart's function getLDS is incredibly powerful, but it makes me think very polite people are
#' going to start knocking on my door, and it fails weirdly pretty much always. This function
#' attempts to alleviate some of that frustration.
#'
#' Tested in test_40ann_biomart.R
#' As with my other biomart functions, this one grew out of frustrations when attempting to work
#' with the incredibly unforgiving biomart service.  It does not attempt to guarantee a useful
#' biomart connection, but will hopefully point out potentially correct marts and attributes to use
#' for a successful query.  I can say with confidence that it works well between mice and humans.
#'
#' @param gene_ids List of gene IDs to translate.
#' @param first_species Linnean species name for one species.
#' @param second_species Linnean species name for the second species.
#' @param host Ensembl server to query.
#' @param trymart Assumed mart name to use.
#' @param attributes Key to query
#' @return list of 4 elements:  The first is the set of all ids, as getLDS seems
#'   to always send them all; the second is the subset corresponding to the
#'   actual ids of interest, and the 3rd/4th are other, optional ids from other datasets.
#' @seealso \pkg{biomaRt} \code{\link[biomaRt]{getLDS}}
#'  \code{\link[biomaRt]{useMart}}
#' @examples
#' \dontrun{
#'  mouse_genes <- biomart_orthologs(some_ids)
#'  ## Hopefully the defaults are sufficient to translate from human to mouse.
#'  yeast_genes <- biomart_orthologs(some_ids, first_species='mmusculus', second_species='scerevisiae')
#' }
#' @export
load_biomart_orthologs <- function(gene_ids, first_species="hsapiens",
                                   second_species="mmusculus",
                                   host="dec2016.archive.ensembl.org",
                                   trymart="ENSEMBL_MART_ENSEMBL",
                                   attributes="ensembl_gene_id") {
  first_mart <- NULL
  first_mart <- try(biomaRt::useMart(biomart=trymart, host=host), silent=TRUE)
  if (class(first_mart) == "try-error") {
    message("Unable to perform useMart, perhaps the host/mart is incorrect: ",
            host, " ", trymart, ".")
    first_marts <- biomaRt::listMarts(host=host)
    first_mart_names <- as.character(first_marts[[1]])
    message("The available first_marts are: ")
    message(first_mart_names)
    message("Trying the first one.")
    first_mart <- biomaRt::useMart(biomart=first_marts[[1, 1]], host=host)
  }
  first_dataset <- paste0(first_species, "_gene_ensembl")
  first_ensembl <- try(biomaRt::useDataset(first_dataset, mart=first_mart), silent=TRUE)
  if (class(first_ensembl) == "try-error") {
    message("Unable to perform useDataset, perhaps the given dataset is incorrect: ",
            first_ensembl, ".")
    datasets <- biomaRt::listDatasets(mart=first_mart)
    message(toString(datasets))
    return(NULL)
  }

  second_mart <- NULL
  second_mart <- try(biomaRt::useMart(biomart=trymart, host=host), silent=TRUE)
  if (class(second_mart) == "try-error") {
    message("Unable to perform useMart, perhaps the host/mart is incorrect: ",
            host, " ", trymart, ".")
    second_marts <- biomaRt::listMarts(host=host)
    second_mart_names <- as.character(second_marts[[1]])
    message("The available second_marts are: ")
    message(second_mart_names)
    message("Trying the first one.")
    second_mart <- biomaRt::useMart(biomart=second_marts[[1, 1]], host=host)
  }
  second_dataset <- paste0(second_species, "_gene_ensembl")
  second_ensembl <- try(biomaRt::useDataset(second_dataset, mart=second_mart), silent=TRUE)
  if (class(second_ensembl) == "try-error") {
    message("Unable to perform useDataset, perhaps the given dataset is incorrect: ",
            second_ensembl, ".")
    datasets <- biomaRt::listDatasets(mart=second_mart)
    message(toString(datasets))
    return(NULL)
  }

  possible_first_attributes <- biomaRt::listAttributes(first_ensembl)
  possible_second_attributes <- biomaRt::listAttributes(second_ensembl)

  ## That is right, I had forgotten but it seems to me that no matter
  ## what list of genes I give this stupid thing, it returns all genes.

  ## Note: As of 2018-03 getLDS is more stringent in the queries it allows.  One must choose the same
  ## attributes from the first and second marts, otherwise it throws an error which looks like:
  ## "The query to the BioMart webservice returned an invalid result: the number
  ## of columns in the result table does not equal the number of attributes in
  ## the query. Please report this to the mailing list."
  ## Therefore I am dropping the arguments first_attributes/second_attributes
  ## and just leaving behind 'attributes'.
  linked_genes <- biomaRt::getLDS(attributes=attributes,
                                  values=gene_ids,
                                  mart=first_ensembl,
                                  attributesL=attributes,
                                  martL=second_ensembl)
  kept_idx <- linked_genes[[1]] %in% gene_ids
  kept_genes <- linked_genes[kept_idx, ]
  new_colnames <- colnames(linked_genes)
  new_colnames[[1]] <- first_species
  second_position <- length(attributes) + 1
  new_colnames[[second_position]] <- second_species
  colnames(kept_genes) <- new_colnames
  colnames(linked_genes) <- new_colnames

  linked_genes <- list(
    "all_linked_genes" = linked_genes,
    "subset_linked_genes" = kept_genes,
    "first_attribs" = possible_first_attributes,
    "second_attribs" = possible_second_attributes)
  return(linked_genes)
}

## EOF
