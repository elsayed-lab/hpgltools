#' Load organism annotation data (mouse/human).
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param keytype a, umm keytype? I need to properly read this code.
#' @param chromosome_column  Which column contains the chromosome information?
#' @param strand_column Ibid' but for the gene's strand.
#' @param start_column  Ibid' but to find where each gene starts.
#' @param end_column  Ibid' but to find where each gene ends.
#' @param description_column  Ibid' but to find gene descriptions.
#' @param fields Columns to include in the output.
#' @param biomart_dataset Name of the biomaRt dataset to query for gene type.
#' @return a table of gene information
#' @seealso \pkg{AnnotationDbi} \pkg{dplyr} \pkg{biomaRt}
#'  \code{\link[AnnotationDbi]{select}} \code{\link[AnnotationDbi]{keytypes}}
#' @examples
#' \dontrun{
#'  host <- load_host_annotations(org, c("a","b"))
#' }
#' @export
load_host_annotations <- function(orgdb=NULL, gene_ids=NULL, keytype="ensembl",
                                  chromosome_column="txchrom", strand_column="txstrand",
                                  start_column="txstart", end_column="txend",
                                  description_column="genename", fields=c("geneid"),
                                  biomart_dataset=NULL) {
  message("Consider deprecating this in favor of load_orgdb_annotations().")
  if (is.null(orgdb)) {
    org_pkgstring <- "library(Homo.sapiens); orgdb <- Homo.sapiens"
    eval(parse(text=org_pkgstring))
  } else if (class(orgdb) == "character") {
    org_pkgstring <- paste0("library(", orgdb, "); orgdb <- ", orgdb)
    eval(parse(text=org_pkgstring))
  }
  keytype <- toupper(keytype)
  chromosome_column <- toupper(chromosome_column)
  strand_column <- toupper(strand_column)
  start_column <- toupper(start_column)
  end_column <- toupper(end_column)
  description_column <- toupper(description_column)
  fields <- toupper(fields)

  ## Gene info
  if (is.null(gene_ids)) {
    gene_ids <- try(AnnotationDbi::keys(orgdb, keytype=keytype))
    if (class(gene_ids) == "try-error") {
      if (grepl(x=gene_ids[[1]], pattern="Invalid keytype")) {
        valid_keytypes <- AnnotationDbi::keytypes(orgdb)
        stop(paste0("Try using valid keytypes: ", toString(valid_keytypes)))
      } else {
        stop("There was an error getting the gene ids.")
      }
    } else {
      message("Extracted all gene ids.")
    }
  }
  ## Note querying by "GENEID" will exclude noncoding RNAs
  fields <- c(keytype, fields, chromosome_column,
              strand_column, start_column, end_column,
              description_column)
  gene_info <- AnnotationDbi::select(
                                orgdb,
                                keys=gene_ids,
                                keytype=keytype,
                                columns=fields)
  .data <- NULL
  ## Convert to tbl_df and reorganize
  gene_info <- dplyr::tbl_df(gene_info) %>%
    dplyr::mutate(
             "transcript_length"=abs(.data[[end_column]] - .data[[start_column]]) + 1)
  ## Remove any entries which are either deprecated (e.g. ENSG00000006074) or
  ## missing transcript information (e.g. ENSG00000050327).
  gene_info <- gene_info[!is.na(gene_info[[start_column]]), ]
  ## This process might be easier if we used 'transcripts(), cds(), and friends.'

  if (!is.null(biomart_dataset)) {
    ## Get gene biotype
    ## Main server temporarily unavailable (2015/11/09)
    ## ensembl_mart <- useMart(biomart="ensembl")
    ## biomart_dataset="hsapiens_gene_ensembl"
    ensembl_mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                     host="www.ensembl.org", biomart_dataset)
    biomart <- biomaRt::useDataset(biomart_dataset, mart=ensembl_mart)
    biomart_genes <- biomaRt::getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=biomart)
    gene_info[["type"]] <- biomart_genes$gene_biotype[match(gene_info$gene_id,
                                                            biomart_genes$ensembl_gene_id)]
  }
  return(gene_info)
}

#' Load organism annotation data.
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' Tested in test_45ann_organdb.R
#' This defaults to a few fields which I have found most useful, but the brave can pass it 'all'.
#'
#' @param orgdb OrganismDb instance.
#' @param include_go  Ask the Dbi for gene ontology information?
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param keytype mmm the key type used?
#' @param fields Columns included in the output.
#' @param sum_exons Perform a sum of the exons in the data set?
#' @return Table of geneids, chromosomes, descriptions, strands, types, and lengths.
#' @seealso \pkg{AnnotationDbi} \pkg{GenomicFeatures} \pkg{BiocGenerics}
#'  \code{\link[AnnotationDbi]{columns}} \code{\link[AnnotationDbi]{keytypes}}
#'  \code{\link[AnnotationDbi]{select}} \code{\link[GenomicFeatures]{exonsBy}}
#' @examples
#' \dontrun{
#'  one_gene <- load_orgdb_annotations(org, c("LmJF.01.0010"))
#' }
#' @export
load_orgdb_annotations <- function(orgdb=NULL, gene_ids=NULL, include_go=FALSE, keytype="ensembl",
                                   strand_column="txstrand", start_column="txstart",
                                   end_column="txend",  chromosome_column="chr",
                                   type_column="type", name_column="genename",
                                   fields=NULL, sum_exons=FALSE) {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    org_pkgstring <- "library(Homo.sapiens); orgdb <- Homo.sapiens"
    eval(parse(text=org_pkgstring))
  } else if (class(orgdb) == "character") {
    org_pkgstring <- paste0("library(", orgdb, "); orgdb <- ", orgdb)
    eval(parse(text=org_pkgstring))
  }
  keytype <- toupper(keytype)
  strand_column <- toupper(strand_column)
  start_column <- toupper(start_column)
  end_column <- toupper(end_column)
  chromosome_column <- toupper(chromosome_column)
  type_column <- toupper(type_column)
  name_column <- toupper(name_column)
  fields <- toupper(fields)
  all_fields <- AnnotationDbi::columns(orgdb)

  if (! name_column %in% all_fields) {
    a_name <- grepl(pattern="NAME", x=all_fields)
    new_name_column <- all_fields[a_name][1]
    message(paste0("Unable to find ", name_column, ", setting it to ", new_name_column, "."))
    name_column <- new_name_column
  }
  if (! type_column %in% all_fields) {
    message(paste0("Unable to find ", type_column, " in the db, removing it."))
    type_column <- NULL
  }
  if (! chromosome_column %in% all_fields) {
    message(paste0("Unable to find ", chromosome_column, " in the db, removing it."))
    chromosome_column <- NULL
  }
  if (! strand_column %in% all_fields) {
    message(paste0("Unable to find ", strand_column, " in the db, removing it."))
    strand_column <- NULL
  }
  if (! start_column %in% all_fields) {
    message(paste0("Unable to find ", start_column, " in the db, removing it."))
    start_column <- NULL
  }
  if (! end_column %in% all_fields) {
    message(paste0("Unable to find ", end_column, " in the db, removing it."))
    end_column <- NULL
  }

  if (is.null(fields)) {
    fields <- c(name_column, type_column, chromosome_column, strand_column, start_column, end_column)
  } else {
    fields <- c(name_column, type_column, chromosome_column, strand_column, start_column, end_column, fields)
  }

  if (fields[1] == "all") {
    message(paste0("Selecting the following fields, this might be too many: \n",
                   toString(all_fields)))
    fields <- all_fields
  }

  if (sum(fields %in% all_fields) != length(fields)) {
    missing_idx <- ! fields %in% all_fields
    missing_fields <- fields[missing_idx]
    found_fields <- fields %in% all_fields
    fields <- fields[found_fields]
    message(paste0("Some requested columns are not available: ", toString(missing_fields), "."))
    message(paste0("The following are available: ", toString(all_fields)))
  }

  ## Gene IDs
  if (is.null(gene_ids)) {
    gene_ids <- try(AnnotationDbi::keys(orgdb, keytype=keytype))
    if (class(gene_ids) == "try-error") {
      if (grepl(x=gene_ids[[1]], pattern="Invalid keytype")) {
        valid_keytypes <- AnnotationDbi::keytypes(orgdb)
        stop(paste0("Try using valid keytypes: ", toString(valid_keytypes)))
      } else {
        stop("There was an error getting the gene ids.")
      }
    } else {
      message("Extracted all gene ids.")
    }
  }
  ## Note querying by "GENEID" will exclude noncoding RNAs

  gene_info <- AnnotationDbi::select(
                                x=orgdb,
                                keys=gene_ids,
                                keytype=keytype,
                                columns=fields)

  ## Compute total transcript lengths (for all exons)
  ## https://www.biostars.org/p/83901/
  gene_exons <- try(GenomicFeatures::exonsBy(orgdb, by="gene"), silent=TRUE)
  if (class(gene_exons) == "try-error") {
    gene_exons <- NULL
  }
  transcripts <- try(GenomicFeatures::transcripts(orgdb), silent=TRUE)
  if (class(transcripts) == "try-error") {
    transcripts <- NULL
  }
  colnames(gene_info) <- tolower(colnames(gene_info))
  if (isTRUE(sum_exons)) {
    message("Summing exon lengths, this takes a while.")
    lengths <- lapply(gene_exons, function(x) {
      sum(BiocGenerics::width(GenomicRanges::reduce(x)))
    })
    message("Adding exon lengths to the gene_exons.")
    lengths <- as.data.frame(unlist(lengths), stringsAsFactors=FALSE)
    colnames(lengths) <- "transcript_length"
    gene_info <- merge(gene_info, lengths, by.x=keytype, by.y="row.names")
  }
  rownames(gene_info) <- make.names(gene_info[[1]], unique=TRUE)

  retlist <- list(
    "genes" = gene_info,
    "transcripts" = transcripts)
  return(retlist)
}

#' Retrieve GO terms associated with a set of genes.
#'
#' AnnotationDbi provides a reasonably complete set of GO mappings between gene ID and
#' ontologies.  This will extract that table for a given set of gene IDs.
#'
#' Tested in test_45ann_organdb.R
#' This is a nice way to extract GO data primarily because the Orgdb data sets are extremely fast
#' and flexible, thus by changing the keytype argument, one may use a lot of different ID types
#' and still score some useful ontology data.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Identifiers of the genes to retrieve annotations.
#' @param keytype  The mysterious keytype returns yet again to haunt my dreams.
#' @param columns  The set of columns to request.
#' @return Data frame of gene IDs, go terms, and names.
#' @seealso \pkg{AnnotationDbi} \pkg{GO.db} \pkg{magrittr}
#'  \code{\link[AnnotationDbi]{select}} \code{\link[dplyr]{tbl_df}}
#' @examples
#' \dontrun{
#'  go_terms <- load_go_terms(org, c("a","b"))
#' }
#' @export
load_orgdb_go <- function(orgdb=NULL, gene_ids=NULL, keytype="ensembl",
                          columns=c("go", "goall", "goid")) {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    org_pkgstring <- "library(Homo.sapiens); orgdb <- Homo.sapiens"
    eval(parse(text=org_pkgstring))
  } else if (class(orgdb) == "character") {
    org_pkgstring <- paste0("library(", orgdb, "); orgdb <- ", orgdb)
    eval(parse(text=org_pkgstring))
  }
  tt <- sm(requireNamespace("GO.db"))
  keytype <- toupper(keytype)
  columns <- toupper(columns)
  if (is.null(gene_ids)) {
    gene_ids <- try(AnnotationDbi::keys(orgdb, keytype=keytype), silent=TRUE)
    if (class(gene_ids) == "try-error") {
      avail_types <- AnnotationDbi::keytypes(orgdb)
      if ("GID" %in% avail_types) {
        message(paste0("The chosen keytype was not available.  Using 'GID'."))
        keytype <- "GID"
        gene_ids <- AnnotationDbi::keys(orgdb, keytype=keytype)
      } else {
        keytype <- avail_types[[1]]
        message(paste0("Neither the chosen keytype, nor 'GID' was available.
The available keytypes are: ", toString(avail_types), "choosing ", keytype, "."))
        gene_ids <- AnnotationDbi::keys(orgdb, keytype=keytype)
      }
    }
  }
  if (class(orgdb)[[1]] == "OrganismDb") {
    message("This is an organismdbi, that should be ok.")
  } else if (class(orgdb)[[1]] == "OrgDb" | class(orgdb)[[1]] == "orgdb") {
    message("This is an orgdb, good.")
  } else {
    stop(paste0("This requires either an organismdbi or orgdb instance, not ", class(orgdb)[[1]]))
  }
  available_columns <- AnnotationDbi::columns(orgdb)
  chosen_columns <- c()
  for (col in columns) {
    if (col %in% available_columns) {
      chosen_columns <- c(chosen_columns, col)
    }
  }
  if (is.null(chosen_columns)) {
    stop(paste0("Did not find any of: ", toString(columns),
                " in the set of available columns: ", toString(available_columns)))
  }
  go_terms <- try(sm(AnnotationDbi::select(x=orgdb,
                                           keys=gene_ids,
                                           keytype=keytype,
                                           columns=chosen_columns)))
  if (class(go_terms) == "try-error") {
    if (grep(pattern="Invalid keytype", x=go_terms[[1]])) {
      message("Here are the possible keytypes:")
      message(toString(AnnotationDbi::keytypes(orgdb)))
      stop()
    }
  }
  ## Deduplicate
  go_terms <- go_terms[!duplicated(go_terms), ]
  if ("GO" %in% chosen_columns) {
    go_terms <- go_terms[!is.na(go_terms[["GO"]]), ]
    go_term_names <- sm(AnnotationDbi::select(x=GO.db::GO.db,
                                              keys=unique(go_terms[["GO"]]),
                                              columns=c("TERM", "GOID", "ONTOLOGY")))
    go_terms <- merge(go_terms, go_term_names, by.x="GO", by.y="GOID")
  }

  ## Remove redundant annotations which differ only in source/evidence
  ## and rename ONTOLOGYALL column
  ##unique(go_terms %>% rename(ONTOLOGY=ONTOLOGYALL) %>% na.omit())
  go_terms <- unique(dplyr::tbl_df(go_terms) %>% na.omit())
  return(go_terms)
}

#' I see no reason to have load_host_annotations and load_parasite_annotations.
#'
#' Thus I am making them both into aliases to load_annotations.
#'
#' @param ... Arguments to be passed to load_annotations.
#' @export
load_parasite_annotations <- function(...) {
  load_orgdb_annotations(...)
}

#' Load organism annotation data (mouse/human).
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param mapto Key to map the IDs against.
#' @param keytype  Choose a keytype, this will yell if it doesn't like your choice.
#' @return a table of gene information
#' @seealso \pkg{AnnotationDbi}
#'  \code{\link[AnnotationDbi]{select}} \code{\link[AnnotationDbi]{keytypes}}
#' @examples
#' \dontrun{
#'  host <- map_orgdb_ids(org, c("a","b"))
#' }
#' @export
map_orgdb_ids <- function(orgdb, gene_ids=NULL, mapto=c("ensembl"), keytype="geneid") {
  mapto <- toupper(mapto)
  keytype <- toupper(keytype)
  avail_keytypes <- AnnotationDbi::keytypes(orgdb)
  found_keys <- sum(mapto %in% avail_keytypes)
  if (found_keys < length(mapto)) {
    warning(paste0("The chosen keytype ", mapto, " is not in this orgdb."))
    warning("Try some of the following instead: ", toString(avail_keytypes), ".")
    warning("Going to pull all the availble keytypes, which is probably not what you want.")
    mapto <- avail_keytypes
  }

  test_masterkey <- sum(keytype %in% avail_keytypes)
  if (test_masterkey != 1) {
    warning(paste0("The chosen master key ", keytype, " is not in this orgdb."))
    warning("Try some of the following instead: ", toString(avail_keytypes), ".")
    warning("I am going to choose one arbitrarily, which is probably not what you want.")
    if ("ENTREZID" %in% avail_keytypes) {
      keytype <- "ENTREZID"
      message("Using entrezid as the master key.")
    } else if ("ENSEMBLID" %in% avail_keytypes) {
      keytype <- "ENSEMBLID"
      message("Using ensemblid as the master key.")
    } else
      stop("Could not think of a usable master key.")
  }

  ## If no gene ids were chosen, grab them all.
  if (is.null(gene_ids)) {
    gene_ids <- AnnotationDbi::keys(orgdb, keytype=keytype)
  }
  ## Gene info
  ## Note querying by "GENEID" will exclude noncoding RNAs
  gene_info <- AnnotationDbi::select(x=orgdb, keytype=keytype, keys=gene_ids, columns=mapto)
  colnames(gene_info) <- tolower(colnames(gene_info))
  return(gene_info)
}

#' Generate a set of joins suitable for the creation of an organismdbi package
#'
#' The graph data required in an organismdbi is pretty specific, this function creates it!
#' It does so by iterating through all keytype pairs between the two packages and looking for
#' matching keys, whichever keys have the most matches win.  It is therefore rather slow.
#'
#' @param first_name  Name of the first package to search
#' @param second_name  Name of the second package to search
#' @param starting  What number join to start from
#' @return  A list named join# where the number is the nth join discovered and the elements
#'   are non-zero matches between the sqlite packages described by first_name and second_name.
orgdb_match_keytypes <- function(first_name, second_name, starting=1, exclude=NULL) {
  tt <- sm(requireNamespace(first_name))
  tt <- sm(try(attachNamespace(first_name), silent=TRUE))
  tt <- sm(requireNamespace(second_name))
  tt <- sm(try(attachNamespace(second_name), silent=TRUE))
  org_pkg <- NULL
  tx_pkg <- NULL
  org_pkgstring <- paste0("org_pkg <- ", first_name)
  eval(parse(text=org_pkgstring))
  tx_pkgstring <- paste0("tx_pkg <- ", second_name)
  eval(parse(text=tx_pkgstring))
  org_keytypes <- AnnotationDbi::keytypes(org_pkg)
  if (!is.null(exclude)) {
    keepers <- ! org_keytypes %in% exclude
    org_keytypes <- org_keytypes[keepers]
  }
  tx_keytypes <- AnnotationDbi::keytypes(tx_pkg)
  if (!is.null(exclude)) {
    keepers <- ! tx_keytypes %in% exclude
    tx_keys <- tx_keytypes[keepers]
  }

  key_matches <- list()
  join_number <- starting
  for (orgk in org_keytypes) {
    org_keys <- AnnotationDbi::keys(x=org_pkg, keytype=orgk)
    max_matched <- 0
    matching_keys <- NULL
    for (txk in tx_keytypes) {
      tx_keys <- AnnotationDbi::keys(x=tx_pkg, keytype=txk)
      matches <- sum(org_keys %in% tx_keys)
      ## message(paste0("Found ", matches, " with ", txk, "."))
      if (matches > max_matched) {
        max_matched <- matches
        matching_keys <- c(orgdb=orgk, txdb=txk)
        names(matching_keys) <- c(first_name, second_name)
      }
    }
    if (max_matched > 0) {
      element_name <- paste0("join", join_number)
      join_number <- join_number + 1
      key_matches[[element_name]] <- matching_keys
    }
  }
  return(key_matches)
}

#' Create an orgdb from an taxonID
#'
#' This function is a bit more fragile than I would like.  I am not completely sold
#' on AnnotationHub yet.
#'
#' @param ahid  TaxonID from AnnotationHub
#' @param title  Title for the annotation hub instance
#' @param species  Species to download
#' @param type  Datatype to download
#' @return An Orgdb instance
#' @seealso \pkg{AnnotationHub} \pkg{S4Vectors}
#' @examples
#' \dontrun{
#'  orgdbi <- mytaxIdToOrgDb(taxid)
#' }
#' @export
take_from_ah <- function(ahid=NULL, title=NULL, species=NULL, type="OrgDb") {
  ## Other available types:
  tt <- sm(loadNamespace("AnnotationHub"))
  ah <- sm(AnnotationHub::AnnotationHub())
  message(paste0("Available types: \n", toString(levels(as.factor(ah$rdataclass)))))

  if (!is.null(type)) {
    ah <- AnnotationHub::query(x=ah, pattern=type)
  }
  if (is.null(title) & is.null(species) & is.null(ahid)) {
    ahid <- "AH57973"  ## org.Hs.eg.db.sqlite
  } else if (is.null(ahid) & is.null(title)) {
    ## Then we got a species
    possible <- ah$species
    titles <- ah$title
    hits_idx <- grepl(pattern=species, x=possible)
    first_true <- which.max(hits_idx)
    first_true_name <- titles[first_true]
    hits <- names(ah)[hits_idx]
    message(paste0("The possible hits are: \n",
                   toString(hits), "\nchoosing: ", hits[1],
                   "\nwhich is ", first_true_name))
    ahid <- hits[1]
  } else if (is.null(ahid) & is.null(species)) {
    ## We got a title
    possible <- ah$title
    hits_idx <- grepl(pattern=title, x=possible)
    first_true <- which.max(hits_idx)
    first_true_name <- possible[first_true]
    hits <- names(ah)[hits_idx]
    message(paste0("The possible hits are: \n",
                   toString(hits), "\nchoosing: ", hits[1],
                   "\nwhich is ", first_true_name))
    ahid <- hits[1]
  }

  ah_names <- names(ah)
  ah_titles <- ah$title
  hit_idx <- ah_names == ahid
  hit_num <- which.max(hit_idx)
  hit_title <- ah_titles[hit_num]
  message(paste0("Chose ", ahid, " which is ", hit_title, "."))
  res <- ah[[ahid]]
  return(res)
}

## EOF
