#' Returns metadata for all eupathdb organisms.
#'
#' @param overwrite  Overwrite existing data?
#' @param webservice  Optional alternative webservice for hard-to-find species.
#' @param dir  Where to put the json.
#' @param use_savefile  Make a savefile of the data for future reference.
#' @return  Dataframe with lots of rows for the various species in eupathdb.
#' @author  Keith Hughitt
#' @export
download_eupath_metadata <- function(overwrite=FALSE, webservice="eupathdb",
                                     dir="eupathdb", use_savefile=TRUE) {
  ## Get EuPathDB version (same for all databases)
  savefile <- paste0(webservice, "_metadata-v", format(Sys.time(), "%Y%m"), ".rda")

  if (!file.exists(dir)) {
    dir.create(dir, recursive=TRUE)
  }
  if (isTRUE(use_savefile)) {
    savefile <- file.path(dir, savefile)
    if (isTRUE(overwrite)) {
      file.remove(savefile)
    }
    if (file.exists(savefile)) {
      metadata <- new.env()
      loaded <- load(savefile, envir=metadata)
      metadata <- metadata[["metadata"]]
      return(metadata)
    }
  }

  db_version <- readLines("http://tritrypdb.org/common/downloads/Current_Release/Build_number")
  .data <- NULL  ## To satisfy R CMD CHECK
  shared_tags <- c("Annotation", "EuPathDB", "Eukaryote", "Pathogen", "Parasite")
  tags <- list(
    "AmoebaDB" = c(shared_tags, "Amoeba"),
    "CryptoDB" = c(shared_tags, "Cryptosporidium"),
    "FungiDB" = c(shared_tags, "Fungus", "Fungi"),
    "GiardiaDB" = c(shared_tags, "Giardia"),
    "MicrosporidiaDB" = c(shared_tags, "Microsporidia"),
    "PiroplasmaDB" = c(shared_tags, "Piroplasma"),
    "PlasmoDB" = c(shared_tags, "Plasmodium"),
    "ToxoDB" = c(shared_tags, "Toxoplasmosis"),
    "TrichDB" = c(shared_tags, "Trichomonas"),
    "TriTrypDB" = c(shared_tags, "Trypanosome", "Kinetoplastid", "Leishmania"))
  tag_strings <- lapply(tags, function(x) {
    paste(x, collapse=",") })

  ## construct API request URL
  base_url <- paste0("https://", webservice, ".org/", webservice, "/webservices/")
  query_string <- "OrganismQuestions/GenomeDataTypes.json?o-fields=all"
  request_url <- paste0(base_url, query_string)

  ## retrieve organism metadata from EuPathDB
  metadata_json <- paste0(dir, "/metadata.json")
  file <- download.file(url=request_url, destfile=metadata_json)
  result <- jsonlite::fromJSON(metadata_json)
  records <- result[["response"]][["recordset"]][["records"]]
  message(paste0("Downloaded: ", request_url))

  ## convert to a dataframe
  dat <- data.frame(t(sapply(records[["fields"]], function(x) {
    x[, "value"] })),
    stringsAsFactors=FALSE)
  colnames(dat) <- records[["fields"]][[1]][["name"]]

  ## shared metadata
  ## I wish I were this confident with %>% and transmute, I always get confused by them
  shared_metadata <- dat %>% dplyr::transmute(
                                      "BiocVersion"=as.character(BiocInstaller::biocVersion()),
                                      "Genome"=sub(".gff", "", basename(.data[["URLgff"]])),
                                      "NumGenes"=.data[["genecount"]],
                                      "NumOrthologs"=.data[["orthologcount"]],
                                      "SourceType"="GFF",
                                      "SourceUrl"=.data[["URLgff"]],
                                      "SourceVersion"=db_version,
                                      "Species"=.data[["organism"]],
                                      "TaxonomyId"=.data[["ncbi_tax_id"]],
                                      "Coordinate_1_based"=TRUE,
                                      "DataProvider"=.data[["project_id"]],
                                      "Maintainer"="Keith Hughitt <khughitt@umd.edu>")

  ## Add project-specific tags for each entry
  shared_metadata[["Tags"]] <- sapply(shared_metadata[["DataProvider"]],
                                      function(x) {
                                        tag_strings[[x]] })

  ## replace missing taxonomy ids with NAs
  shared_metadata[["TaxonomyId"]][shared_metadata[["TaxonomyId"]] == ""] <- NA

  ## overide missing taxonomy ids for strains where it can be assigned; ideally
  ## OrgDb and GRanges objects should not depend on taxonomy id information since
  ## this precludes the inclusion of a lot of prokaryotic resources.
  known_taxon_ids <- data.frame(
    species=c("Ordospora colligata OC4",
              "Trypanosoma cruzi CL Brener Esmeraldo-like",
              "Trypanosoma cruzi CL Brener Non-Esmeraldo-like"),
    taxonomy_id=c("1354746", "353153", "353153"),
    stringsAsFactors=FALSE)

  taxon_mask <- shared_metadata[["Species"]] %in% known_taxon_ids[["species"]]
  ind <- match(shared_metadata[taxon_mask, "Species"], known_taxon_ids[["species"]])
  shared_metadata[taxon_mask, ][["TaxonomyId"]] <- as.character(known_taxon_ids[["taxonomy_id"]][ind])

  ## exclude remaining species which are missing taxonomy information from
  ## metadata; cannot construct GRanges/OrgDb instances for them since they are
  ## have no known taxonomy id, and are not in available.species()
  na_ind <- is.na(shared_metadata[["TaxonomyId"]])
  ## I think I will try to hack around this problem.
  ##message(sprintf("- Excluding %d organisms for which no taxonomy id could be assigned (%d remaining)",
  ##                sum(na_ind), sum(!na_ind)))
  ##shared_metadata <- shared_metadata[!na_ind,]

  ## convert remaining taxonomy ids to numeric
  shared_metadata[["TaxonomyId"]] <- as.numeric(shared_metadata[["TaxonomyId"]])

  ## remove any organisms for which no GFF is available
  ## gff_exists <- sapply(shared_metadata[["SourceUrl"]],
  ##                      function(url) { httr::HEAD(url)[["status_code"]] == 200 })
  ## remove any organisms for which no GFF is available
  ## Once again, I will attempt a workaround, probably via bioconductor.
  ## gff_exists <- sapply(shared_metadata$SourceUrl, function(url) { HEAD(url)$status_code == 200 })
  ##message(sprintf("- Excluding %d organisms for which no GFF file is available (%d remaining)",
  ##                sum(!gff_exists), sum(gff_exists)))
  ##shared_metadata <- shared_metadata[gff_exists,]

  ## generate separate metadata table for OrgDB and GRanges targets
  granges_metadata <- shared_metadata %>%
    dplyr::mutate(
             Title=sprintf("Transcript information for %s", .data[["Species"]]),
             Description=sprintf("%s %s transcript information for %s",
                                 .data[["DataProvider"]], .data[["SourceVersion"]],
                                 .data[["Species"]]),
             RDataClass="GRanges",
             DispatchClass="GRanges",
             ResourceName=sprintf("GRanges.%s.%s%s.rda", gsub("[ /.]+", "_", .data[["Species"]]),
                                  tolower(.data[["DataProvider"]]), .data[["SourceVersion"]], "rda")) %>%
    dplyr::mutate(DataPath=file.path("EuPathDB", "GRanges", .data[["BiocVersion"]], .data[["ResourceName"]]))

  metadata <- shared_metadata %>%
    dplyr::mutate(
             "Title"=sprintf("Genome wide annotations for %s", .data[["Species"]]),
             "Description"=sprintf("%s %s annotations for %s",
                                   .data[["DataProvider"]],
                                   .data[["SourceVersion"]],
                                   .data[["Species"]]),
             "RDataClass"="OrgDb",
             "DispatchClass"="SQLiteFile",
             "ResourceName"=sprintf("org.%s.%s.db.sqlite", gsub("[ /.]+", "_", .data[["Species"]]),
                                    tolower(substring(.data[["DataProvider"]], 1, nchar(.data[["DataProvider"]]) - 2)))
           ) %>%
    dplyr::mutate("RDataPath"=file.path("EuPathDB", "OrgDb",
                                        .data[["BiocVersion"]],
                                        .data[["ResourceName"]]))

  if (isTRUE(use_savefile)) {
    if (isTRUE(overwrite) | !file.exists(savefile)) {
      saved <- save(list="metadata", file=savefile)
    }
  }

  return(metadata)
}

#' The new eupath system provides 3 output types for downloading data.  This uses the raw one.
#'
#' For the life of me, I could not figure out how to query the big text tables as the
#' tabular format.  Every query I sent came back telling me I gave it incorrect parameter
#' despite the fact that I was copy/pasting the example given me by the eupathdb maintainers.
#' So, I got mad and asked it for the raw format, and so this function was born.
#'
#' @param entry  Annotation entry for a given species
#' @param question  Which query to try?  Molecular weight is the easiest, as it was their example.
#' @param table_name  Used to make sure all columns are unique by prefixing them with the table name.
#' @param parameters  Query parameters when posting
#' @param columns  Columns for which to ask.
#' @param minutes  How long to wait until giving up and throwing an error.
#' @return  A hopefully huge table of eupath data.
post_eupath_raw <- function(entry, question="GeneQuestions.GenesByMolecularWeight",
                            table_name=NULL, parameters=NULL, columns="primary_key",
                            minutes=40) {
  species <- entry[["Species"]]
  provider <- tolower(entry[["DataProvider"]])
  ## determine appropriate prefix to use
  prefix_mapping <- list(
    "amoebadb" = "amoeba",
    "microbiomedb" = "mbio",
    "microsporidiadb" = "micro",
    "piroplasmadb" = "piro",
    "plasmodb" = "plasmo",
    "schistodb" = "schisto",
    "toxodb" = "toxo"
  )
  uri_prefix <- provider
  if (uri_prefix %in% names(prefix_mapping)) {
    uri_prefix <- prefix_mapping[[uri_prefix]]
  }
  if (is.null(parameters)) {
    parameters <- list("organism" = jsonlite::unbox(species))
  }

  ##query_body <- list(
  ##  ## 3 elements, answerSpec, formatting, format.
  ##  "answerSpec" = list(
  ##    "questionName" = jsonlite::unbox(question),
  ##    "parameters" = parameters,
  ##    "viewFilters" = list(),
  ##    "filters" = list()
  ##  ),
  ##  "formatting" = list(
  ##    "formatConfig" = list(
  ##      "includeHeaders" = jsonlite::unbox("true"),
  ##      "attributes" = columns,
  ##      "attachmentType" = jsonlite::unbox("plain")
  ##    ),
  ##    "format" = jsonlite::unbox("fullRecord")
  ##  ))
  answerlist <- list(
    "questionName" = jsonlite::unbox(question),
    "parameters" = parameters,
    "viewFilters" = list(),
    "filters" = list())
  formattinglist <- list(
    "formatConfig" = list(
      "includeHeaders" = jsonlite::unbox("true"),
      "attributes" = columns,
      "attachmentType" = jsonlite::unbox("plain")),
    "format" = jsonlite::unbox("fullRecord"))
  query_body <- list(
    "answerSpec" = answerlist,
    "formatting" = formattinglist)
  body <- jsonlite::toJSON(query_body)

  api_uri <- sprintf("http://%s.org/%s/service/answer", provider, uri_prefix)
  result <- httr::POST(
                    url=api_uri,
                    body=body,
                    httr::content_type("application/json"),
                    httr::timeout(minutes * 60))
  if (result[["status_code"]] == "422") {
    warning("The provided species does not have a table of weights.")
    return(data.frame())
  } else if (result[["status_code"]] == "400") {
    warning("Status 400 was returned, likely a bad formatConfig.")
    ## Interesting, querying a tritrypdb entry with the above query_body works fine,
    ## but querying microsporidiadb.org with it fails with 'Could not configure reporter 'fullRecord' with passed formatConfig.
    ## I wonder what is different between them?
    ## Weirdly, if I remove the format fields from the formatting list, then I get a nested response,
    ## But this response requires a subset of the original columns and is in yet another format.
    ##columns <- c("primary_key", "sequence_id", "chromosome", "organism", "gene_type", "gene_location_text",
    ##             "gene_name", "gene_exon_count", "is_pseudo", "gene_transcript_count", "gene_ortholog_number",
    ##             "gene_orthomcl_name", "gene_entrez_id", "transcript_length", "exon_count", "strand",
    ##             "cds_length", "tm_count", "molecular_weight", "isoelectric_point", "signalp_scores",
    ##             "signalp_peptide", "annotated_go_function", "annotated_go_process", "annotated_go_component",
    ##             "annotated_go_id_function", "annotated_go_id_process", "annotated_go_id_component",
    ##             "predicted_go_id_function", "predicted_go_id_process", "predicted_go_id_component",
    ##             "ec_numbers", "ec_numbers_derived", "five_prime_utr_length", "three_prime_utr_length",
    ##             "location_text", "gene_previous_ids", "protein_sequence", "cds")
    ##answerlist <- list(
    ##"questionName" = jsonlite::unbox(question),
    ##"parameters" = parameters,
    ##"viewFilters" = list(),
    ##"filters" = list())
    ##formattinglist <- list(
    ##  "formatConfig" = list(
    ##    "includeHeaders" = jsonlite::unbox("true"),
    ##    "attributes" = columns,
    ##    "attachmentType" = jsonlite::unbox("plain")),
    ##  "format" = jsonlite::unbox("tableTabular"))
    ##query_body <- list(
    ##  "answerSpec" = answerlist,
    ##  "formatting" = formattinglist,
    ##  "format" = jsonlite::unbox("fullRecord"))
    ##)
    ##body <- jsonlite::toJSON(query_body)
    ##api_uri <- sprintf("http://%s.org/%s/service/answer", provider, uri_prefix)
    ##result <- httr::POST(
    ##                  url=api_uri,
    ##                  body=body,
    ##                  httr::content_type("application/json"),
    ##                  httr::timeout(minutes * 60))
    ##cont <- httr::content(result)
    return(data.frame())
  } else if (result[["status_code"]] != "200") {
    warning("An error status code was returned.")
    return(data.frame())
  } else if (length(result[["content"]]) < 100) {
    warning("A minimal amount of content was returned.")
  }

  ## Get the content, this will take a while, as the result from eupathdb might be > 50 Mb of stuff.
  cont <- httr::content(result)
  ## Sadly, most of that stuff is completely unwanted.  This is because we are using the
  ## 'fullRecord' format, as it is the only format I have been able to get to work so far.
  ## This format is newline separated fields with entries separated by 4 returns with dashes...
  ## Ergo the following line, which separates the entries by the dashes/returns into individual
  ## strings with the newlines remaining inside them.  So we will need to use some regular
  ## expressions in order to extract the column names and data.
  entries <- strsplit(x=cont, split="\n\n------------------------------------------------------------\n\n")[[1]]
  ## We will read the first entry in order to extract the column names.
  stuff <- read.delim(textConnection(entries[1]), sep="\n", header=FALSE)
  ## My regular expression pattern needs to by greedy in the correct places
  ## because for reasons passing all understanding, some fields have colons inside them...
  mypattern <- "^(.+?)\\: (.+)?$"
  ## If I am going to make column names, I need first to get the first part of stuff: otherstuff
  stupid_column_names <- gsub(pattern=mypattern, replacement="\\1", x=stuff[["V1"]], perl=TRUE)
  column_names <- columns
  ## Then get rid of any punctuation, as there is a column '# TM domains' -- that is bad.
  ##
  ## column_names <- gsub(pattern="[[:punct:]]", replacement="", x=column_names)
  ##
  ## Get rid of any extraneous spaces from removing punctuation, but since I cannot be certain
  ## that there is no punctuation in the middle of words, just look at the beginning of the strings.
  ##
  ## column_names <- gsub(pattern="^ +", replacement="", x=column_names)
  ##
  ## Finally, I do not accept column names with spaces.
  ##
  ## column_names <- gsub(pattern=" ", replacement="_", x=column_names)
  ##
  column_names[1] <- "GID"
  information <- data.frame(row.names=1:length(entries))
  for (col in column_names) {
    new_col <- rep(NA, times=length(entries))
    information <- cbind(information, new_col)
  }
  colnames(information) <- column_names
  bar <- utils::txtProgressBar(style=3)
  ## Now fill in the data using the other side of my regular expression.
  for (c in 1:length(entries)) {
    pct_done <- c / length(entries)
    setTxtProgressBar(bar, pct_done)
    stuff <- read.delim(textConnection(entries[c]), sep="\n", header=FALSE)
    material <- gsub(pattern="^(.+?)\\: (.+)?$", replacement="\\2", x=stuff[["V1"]])
    information[c, ] <- material
  }
  close(bar)
  ## remove duplicated rows
  information <- information[!duplicated(information), ]
  ## In some cases we will want to prefix the columns with the table name...
  if (!is.null(table_name)) {
    for (c in 2:length(colnames(information))) {
      col_name <- colnames(information)[c]
      prefix_string <- paste0(toupper(table_name), "_")
      ## Use if() test this to avoid column names like 'GO_GO_ID'
      foundp <- grepl(pattern=paste0("^", prefix_string), x=toupper(col_name))
      if (!foundp) {
        new_col <- paste0(toupper(table_name), "_", toupper(col_name))
        colnames(information)[c] <- new_col
      }
    }
  }
  return(information)
}

#' Queries one of the EuPathDB APIs using a POST request and returns a
#' dataframe representation of the result.
#' Note: As of 2017/07/13, POST requests are not yet supported on EuPathDB.
#' Note: 2017/07/13 POST queries can only use the new API
#'
#' @param query_body String of additional query arguments
#' @param species Species name if missing an entry
#' @param entry The single metadatum containing the base url of the provider, species, etc.
#' @param metadata  If no entry is provided, then it may be retrieved given a species and this.
#' @param table_name  The name of the table to extract, this is provided to make for prettier labeling.
#' @param minutes  A timeout when querying the eupathdb.
#' @param ...  Extra arguments for stuff like download_metadtata()
#' @return list containing response from API request.
#'
#' More information
#' ----------------
#' 1. https://tritrypdb.org/tritrypdb/serviceList.jsp
#' @author Keith Hughitt
#' @export
post_eupath_table <- function(query_body, species=NULL, entry=NULL, metadata=NULL,
                              table_name=NULL, minutes=30, ...) {

  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- download_eupath_metadata(dir=dir, ...)
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }

  ## determine appropriate prefix to use
  species <- entry[["Species"]]
  provider <- tolower(entry[["DataProvider"]])
  prefix_mapping <- list(
    "amoebadb" = "amoeba",
    "microbiomedb" = "mbio",
    "microsporidiadb" = "micro",
    "piroplasmadb" = "piro",
    "plasmodb" = "plasmo",
    "schistodb" = "schisto",
    "toxodb" = "toxo"
  )
  uri_prefix <- provider
  if (uri_prefix %in% names(prefix_mapping)) {
    uri_prefix <- prefix_mapping[[uri_prefix]]
  }

  ## construct API query
  api_uri <- sprintf("http://%s.org/%s/service/answer", provider, uri_prefix)
  ##result <- httr::POST(api_uri, body=tt,
  body <- jsonlite::toJSON(query_body)
  ##message("The request is:")
  ##message(body)
  result <- httr::POST(
                    url=api_uri,
                    body=body,
                    httr::content_type("application/json"),
                    httr::timeout(minutes * 60))
  if (result[["status_code"]] == "422") {
    warning("The provided species does not have a table of weights.")
    return(data.frame())
  } else if (result[["status_code"]] != "200") {
    warning("An error status code was returned.")
    return(data.frame())
  } else if (length(result[["content"]]) < 100) {
    warning("A minimal amount of content was returned.")
  }

  result <- httr::content(result)
  result <- read.delim(textConnection(result), sep="\t")
  ## If nothing was received, return nothing.
  if (nrow(result) == 0) {
    return(data.frame())
  }

  ## If a column is just 'X', then I think it can go away.
  non_stupid_columns <- colnames(result) != "X"
  result <- result[, non_stupid_columns]

  ## simplify column names, the are downloaded with
  ## annoyingly stupid names like:
  ## > colnames(dat)
  ## [1] "X.Gene.ID."                        "X.pathway_source_id."
  ## [3] "X.Pathway."                        "X.Pathway.Source."
  ## [5] "X.EC.Number.Matched.in.Pathway."   "X.expasy_url."
  ## [7] "X...Reactions.Matching.EC.Number."
  new_colnames <- toupper(colnames(result))
  ## Get rid of dumb X. prefix
  new_colnames <- gsub("^X\\.+", replacement="", x=new_colnames)
  ## Get rid of spurious end .
  new_colnames <- gsub("\\.$", replacement="", x=new_colnames)
  ## Get rid of internal .'s
  new_colnames <- gsub("\\.", replacement="_", x=new_colnames)
  colnames(result) <- new_colnames
  colnames(result)[1] <- "GID"
  ## remove duplicated rows
  result <- result[!duplicated(result), ]
  if (!is.null(table_name)) {
    for (c in 2:length(colnames(result))) {
      col_name <- colnames(result)[c]
      new_col <- paste0(toupper(table_name), "_", toupper(col_name))
      colnames(result)[c] <- new_col
    }
  }
  return(result)
}

#' Gather all available annotation data for a given eupathdb species.
#'
#' This function fills in the parameters to post_eupath_raw() so that one can download
#' all the available data for a given parasite into one massive table.  It should also
#' provide some constraints to the data rather than leaving it all as characters.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @param ...  Used for downloading metadata.
#' @return  A big honking table.
post_eupath_annotations <- function(species="Leishmania major", entry=NULL,
                                    metadata=NULL, dir="eupathdb", ...) {
  savefile <- file.path(dir, "annotations.rda")
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
      ## metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
    species <- entry[["Species"]]
  }

  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    result <- new.env()
    load(savefile, envir=result)
    result <- result[["result"]]
    return(result)
  }

  ## query body as a structured list
  ## This list was generated by going to:
  ## view-source:http://tritrypdb.org/webservices/GeneQuestions/GenesByMolecularWeight.wadl
  ## scrolling down to the 'o-fields' section, and writing down the most likely useful column names.
  ## These are written one per line in an attempt to make looking for new/changed columns from one
  ## eupathdb release to the next easier.
  field_list <- list(
    "primary_key" = "Gene ID",
    "transcript_link" = "Transcript ID",
    "sequence_id" = "Genomic Sequence ID",
    "chromosome" = "Chromosome",
    "gene_product" = "Product Description",
    "organism" = "Organism",
    "gene_type" = "Gene Type",
    "gene_location_text" = "Genomic Location (Gene)",
    "gene_name" = "Gene Name or Symbol",
    "gene_exon_count" = "# Exons in Gene",
    "is_pseudo" = "Is Pseudo",
    "gene_transcript_count" = "# Transcripts",
    "gene_ortholog_number" = "Ortholog count",
    "gene_paralog_number" = "Paralog Count",
    "gene_orthomcl_name" = "Ortholog Group",
    "gene_total_hts_snps" = "Total SNPs All Strains",
    "gene_hts_nonsynonymous_snps" = "NonSynonymous SNPs All Strains",
    "gene_hts_synonymous_snps" = "Synonymous SNPs All Strains",
    "gene_hts_noncoding_snps" = "Non-Coding SNPs All Strains",
    "gene_hts_stop_codon_snps" = "SNPs with Stop Codons All Strains",
    "gene_hts_nonsyn_syn_ratio" = "NonSyn/Syn SNP Ratio All Strains",
    "uniprot_id" = "UniProt ID",
    "gene_entrez_id" = "Entrez Gene ID",
    "transcript_product" = "Transcript Product Description",
    "transcript_length" = "Transcript Length",
    "exon_count" = "# Exons in Transcript",
    "strand" = "Gene Strand",
    "cds_length" = "CDS Length",
    "tm_count" = "# TM Domains",
    "molecular_weight" = "Molecular Weight",
    "isoelectric_point" = "Isoelectric Point",
    "signalp_scores" = "Signalp Scores",
    "signalp_peptide" = "Signalp Peptide",
    "annotated_go_function" = "Curated GO Functions",
    "annotated_go_process" = "Curated GO Processes",
    "annotated_go_component" = "Curated GO Components",
    "annotated_go_id_function" = "Curated GO Function IDs",
    "annotated_go_id_process" = "Curated GO Process IDs",
    "annotated_go_id_component" = "Curated GO Component IDs",
    "predicted_go_id_function" = "Computed GO Function IDs",
    "predicted_go_id_process" = "Computed GO Process IDs",
    "predicted_go_id_component" = "Computed GO Component IDs",
    "ec_numbers" = "EC numbers",
    "ec_numbers_derived" = "EC numbers from OrthoMCL",
    "five_prime_utr_length" = "Annotated 5' UTR length",
    "three_prime_utr_length" = "Annotated 3' UTR length",
    "location_text" = "Genomic Location (Transcript)",
    "gene_previous_ids" = "Previous ID(s)",
    "transcripts_found_per_gene" = "# Transcripts that Met Search Critera",
    "transcript_index_per_gene" = "Transcript local index",
    "transcript_sequence" = "Predicted RNA/mRNA Sequence (introns spliced out)",
    "protein_sequence" = "Predicted Protein Sequence",
    "cds" = "Coding Sequence",
    "uri" = "Image",
    "gene_source_id" = "gene_source_id",
    "source_id" = "source_id")
  parameters <- list(
    "organism" = jsonlite::unbox(species),
    "min_molecular_weight" = jsonlite::unbox("1"),
    "max_molecular_weight" = jsonlite::unbox("10000000000000000")
  )
  result <- post_eupath_raw(entry,
                            question="GeneQuestions.GenesByMolecularWeight",
                            parameters=parameters,
                            table="annot",
                            columns=names(field_list))
  colnames(result) <- tolower(colnames(result))
  numeric_columns <- c(
    "annot_gene_exon_count",
    "annot_gene_transcript_count",
    "annot_gene_ortholog_number",
    "annot_gene_paralog_number",
    "annot_gene_total_hts_snps",
    "annot_gene_hts_nonsynonymous_snps",
    "annot_gene_hts_synonymous_snps",
    "annot_gene_hts_noncoding_snps",
    "annot_gene_hts_stop_codon_snps",
    "annot_gene_hts_nonsyn_syn_ratio",
    "annot_transcript_length",
    "annot_exon_count",
    "annot_cds_length",
    "annot_tm_count",
    "annot_molecular_weight",
    "annot_isoelectric_point",
    "annot_five_prime_utr_length",
    "annot_three_prime_utr_length")
  for (col in numeric_columns) {
    result[[col]] <- as.numeric(result[[col]])
  }
  factor_columns <- c(
    "annot_chromosome",
    "annot_gene_type",
    "annot_is_pseudo",
    "annot_strand",
    "annot_tm_count",
    "annot_exon_count")
  for (col in factor_columns) {
    result[[col]] <- as.factor(result[[col]])
  }
  colnames(result) <- toupper(colnames(result))

  message(paste0("Saving annotations to ", savefile))
  save(result, file=savefile)
  return(result)
}

#'  Use the post interface to get GO data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @param ... Extra options when downloading metadata.
#' @return  A big honking table.
post_eupath_go_table <- function(species="Leishmania major", entry=NULL,
                                 metadata=NULL, dir="eupathdb", ...) {
  savefile <- file.path(dir, "go_table.rda")
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
      ## metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  species <- entry[["Species"]]
  ## query body as a structured list
  ## Parameters taken from the pdf "Exporting Data - Web Services.pdf" received
  ## from Cristina

  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    result <- new.env()
    load(savefile, envir=result)
    result <- result[["result"]]
    return(result)
  }

  query_body <- list(
    ## 3 elements, answerSpec, formatting, format.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox("GeneQuestions.GenesByTaxonGene"),
      "parameters" = list("organism" = jsonlite::unbox(species)),
      "viewFilters" = list(),
      "filters" = list()
    ),
    "formatting" = list(
      "formatConfig" = list(
        "tables" = "GOTerms",
        "includeEmptyTables" = jsonlite::unbox("true"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="go")
  message(paste0("Saving annotations to ", savefile))
  save(result, file=savefile)
  return(result)
}

#'  Use the post interface to get ortholog data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @param ... Extra options for downloading metadata.
#' @return  A big honking table.
post_eupath_ortholog_table <- function(species="Leishmania major", entry=NULL,
                                       metadata=NULL, dir="eupathdb", ...) {
  savefile <- file.path(dir, "ortholog_table.rda")
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      ## metadata <- sm(download_eupath_metadata(dir=dir, ...))
      metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  species <- entry[["Species"]]

  ## query body as a structured list
  ## Parameters taken from the pdf "Exporting Data - Web Services.pdf" received
  ## from Cristina
  query_body <- list(
    ## 3 elements, answerSpec, formatting, format.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox("GeneQuestions.GenesByTaxonGene"),
      "parameters" = list("organism" = jsonlite::unbox(species)),
      "viewFilters" = list(),
      "filters" = list()
    ),
    "formatting" = list(
      "formatConfig" = list(
        "tables" = "Orthologs",
        ## "includeEmptyTables" = jsonlite::unbox("true"),
        "includeEmptyTables" = jsonlite::unbox("false"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    result <- new.env()
    load(savefile, envir=result)
    result <- result[["result"]]
    return(result)
  }

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="orthologs")

  message(paste0("Saving annotations to ", savefile))
  save(result, file=savefile)
  return(result)
}

#'  Use the post interface to get interpro data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @param ... Extra options when downloading metadata.
#' @return  A big honking table.
post_eupath_interpro_table <- function(species="Leishmania major strain Friedlin", entry=NULL,
                                       metadata=NULL, dir="eupathdb", ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      ## metadata <- sm(download_eupath_metadata(dir=dir, ...))
      metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }

  species <- entry[["Species"]]
  ## query body as a structured list
  ## Parameters taken from the pdf "Exporting Data - Web Services.pdf" received
  ## from Cristina
  query_body <- list(
    ## 3 elements, answerSpec, formatting, format.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox("GeneQuestions.GenesByTaxonGene"),
      "parameters" = list("organism" = jsonlite::unbox(species)),
      "viewFilters" = list(),
      "filters" = list()
    ),
    "formatting" = list(
      "formatConfig" = list(
        "tables" = "InterPro",
        "includeEmptyTables" = jsonlite::unbox("true"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  savefile <- file.path(dir, "interpro_table")
  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    result <- new.env()
    load(savefile, envir=result)
    result <- result[["result"]]
    return(result)
  }


  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="interpro")

  message(paste0("Saving annotations to ", savefile))
  save(result, file=savefile)
  return(result)
}

#'  Use the post interface to get pathway data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @param ... Extra options when downloading metadata
#' @return  A big honking table.
post_eupath_pathway_table <- function(species="Leishmania major", entry=NULL,
                                      metadata=NULL, dir="eupathdb", ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
      ## metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  species <- entry[["Species"]]
  ## query body as a structured list
  ## Parameters taken from the pdf "Exporting Data - Web Services.pdf" received
  ## from Cristina
  query_body <- list(
    ## 2 elements, answerSpec, formatting.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox("GeneQuestions.GenesByTaxonGene"),
      "parameters" = list("organism" = jsonlite::unbox(species)),
      "viewFilters" = list(),
      "filters" = list()
    ),
    "formatting" = list(
      "formatConfig" = list(
        "tables" = "MetabolicPathways",
        "includeEmptyTables" = jsonlite::unbox("false"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  savefile <- file.path(dir, "pathway_table.rda")
  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    result <- new.env()
    load(savefile, envir=result)
    result <- result[["result"]]
    return(result)
  }

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="pathway")

  message(paste0("Saving annotations to ", savefile))
  save(result, file=savefile)
  return(result)
}

#' Query ortholog tables from the eupathdb one gene at a time.
#'
#' Querying the full ortholog table at eupathdb.org fails mysteriously.
#' This is a horrible brute-force approach to get around this.
#'
#' @param species  What species to query
#' @param dir  Directory to which to save intermediate data (currently unused)
#' @param entry  An entry from the eupathdb metadata to use for other parameters.
#' @param metadata  The set of eupathdb metadata from which to query.
#' @param ...  Extra parameters for downloading eupathdb metadata.
#' @export
get_orthologs_all_genes <- function(species="Leishmania major", dir="eupathdb",
                                    entry=NULL, metadata=NULL, ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
      ## metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
    species <- entry[["Species"]]
  }
  ## query body as a structured list
  field_list <- c(
    "primary_key")
  parameters <- list(
    "organism" = jsonlite::unbox(species),
    "min_molecular_weight" = jsonlite::unbox("1"),
    "max_molecular_weight" = jsonlite::unbox("10000000000000000")
  )
  message("Getting the set of possible genes.")
  result <- post_eupath_raw(entry,
                            question="GeneQuestions.GenesByMolecularWeight",
                            parameters=parameters,
                            columns=field_list)

  savefile <- file.path(dir, "ortholog_table.rda")
  if (file.exists(savefile)) {
    message("We can save some time by reading the savefile.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    all_orthologs <- new.env()
    load(savefile, envir=all_orthologs)
    all_orthologs <- all_orthologs[["all_orthologs"]]
    return(all_orthologs)
  }

  all_orthologs <- data.frame()
  message("Downloading orthologs, one gene at a time.")
  bar <- utils::txtProgressBar(style=3)
  for (c in 1:length(result)) {
    pct_done <- c / length(result)
    setTxtProgressBar(bar, pct_done)
    id <- result[c]
    orthos <- sm(get_orthologs_one_gene(species=species, gene=id, entry=entry))
    all_orthologs <- rbind(all_orthologs, orthos)
  }
  close(bar)

  message(paste0("Saving annotations to ", savefile))
  save(all_orthologs, file=savefile)
  return(all_orthologs)
}

get_orthologs_one_gene <- function(species="Leishmania major", gene="LmjF.01.0010",
                               dir="eupathdb", entry=NULL, metadata=NULL, ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
      ## metadata <- sm(download_eupath_metadata(dir=dir))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  species <- entry[["Species"]]
  provider <- tolower(entry[["DataProvider"]])
  prefix_mapping <- list(
    "amoebadb" = "amoeba",
    "microbiomedb" = "mbio",
    "microsporidiadb" = "micro",
    "piroplasmadb" = "piro",
    "plasmodb" = "plasmo",
    "schistodb" = "schisto",
    "toxodb" = "toxo"
  )
  uri_prefix <- provider
  if (uri_prefix %in% names(prefix_mapping)) {
    uri_prefix <- prefix_mapping[[uri_prefix]]
  }

  question <- "GenesOrthologousToAGivenGene"
  ## I am not sure what I was doing with the following 5 lines of code.
  ## It looks like I was intending to gather the set of required parameters programmatically.
  ## However I did not finish the logic.
  params_uri <- sprintf("http://%s.org/%s/webservices/GeneQuestions/%s.wadl",
                        provider, uri_prefix, question)
  result <- xml2::read_html(params_uri)
  test <- rvest::html_nodes(result, "param")
  param_string <- rvest::html_attr(x=test, name="default")[[1]]

  parameters <- list(
    "organism" = jsonlite::unbox(species),
    "single_gene_id" = jsonlite::unbox(gene))
  columns <- c("primary_key", "organism", "orthomcl_link",
               "gene_ortholog_number", "gene_paralog_number")
  query_body <- list(
    ## 3 elements, answerSpec, formatting, format.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox(paste0("GeneQuestions.", question)),
      "parameters" = parameters,
      "viewFilters" = list(),
      "filters" = list()
    ),
    "formatting" = list(
      "formatConfig" = list(
        "includeHeaders" = jsonlite::unbox("true"),
        "attributes" = columns,
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("fullRecord")
    ))
  api_uri <- sprintf("http://%s.org/%s/service/answer", provider, uri_prefix)
  body <- jsonlite::toJSON(query_body)
  result <- httr::POST(
                    url=api_uri,
                    body=body,
                    httr::content_type("application/json"))

  if (result[["status_code"]] == "422") {
    stop("There is a missing parameter.")
  } else if (result[["status_code"]] == "400") {
    stop("An invalid format configuration was provided.")
  } else if (result[["status_code"]] != "200") {
    warning("An error status code was returned.")
    return(data.frame())
  } else if (length(result[["content"]]) < 100) {
    warning("A minimal amount of content was returned.")
  }

  cont <- httr::content(result)
  entries <- strsplit(x=cont, split="\n\n------------------------------------------------------------\n\n")[[1]]
  stuff <- read.delim(textConnection(entries[1]), sep="\n", header=FALSE)
  mypattern <- "^(.+?)\\: (.+)?$"
  ## If I am going to make column names, I need first to get the first part of stuff: otherstuff
  column_names <- gsub(pattern=mypattern, replacement="\\1", x=stuff[["V1"]], perl=TRUE)
  ## Then get rid of any punctuation, as there is a column '# TM domains' -- that is bad.
  column_names <- gsub(pattern="[[:punct:]]", replacement="", x=column_names)
  ## Get rid of any extraneous spaces from removing punctuation, but since I cannot be certain
  ## that there is no punctuation in the middle of words, just look at the beginning of the strings.
  column_names <- gsub(pattern="^ +", replacement="", x=column_names)
  ## Finally, I do not accept column names with spaces.
  column_names <- gsub(pattern=" ", replacement="_", x=column_names)
  column_names[1] <- "Ortholog_ID"
  information <- data.frame(row.names=1:length(entries))
  for (col in column_names) {
    new_col <- rep(NA, times=length(entries))
    information <- cbind(information, new_col)
  }
  colnames(information) <- column_names
  ## Now fill in the data using the other side of my regular expression.
  for (c in 1:length(entries)) {
    stuff <- read.delim(textConnection(entries[c]), sep="\n", header=FALSE)
    material <- gsub(pattern="^(.+?)\\: (.+)?$", replacement="\\2", x=stuff[["V1"]])
    information[c, ] <- material
  }
  ## remove duplicated rows
  information <- information[!duplicated(information), ]
  ## Now fill in the original ID
  information[["GID"]] <- gene
  ## The dplyr way of moving a column to the front.
  ## information <- information %>% dplyr::select(GID, everything())
  ## The base way of moving a column to the front
  new_order <- c(which(colnames(information) == "GID"), which(colnames(information) != "GID"))
  information <- information[, new_order]
  colnames(information) <- toupper(colnames(information))

  return(information)
}

## EOF
