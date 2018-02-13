#' Search the eupathdb metadata for a given species substring.
#'
#' If the specific species is not found, look for a reasonably approximation.
#' stop() if nothing is found.
#'
#' @param species  Guess
#' @param metadata  Eupathdb metadata.
#' @param ...  Extra arguments passed to download_eupath_metadata()
#' @return  A single entry from the eupathdb metadata.
check_eupath_species <- function(species="Leishmania major strain Friedlin", metadata=NULL, ...) {
  if (is.null(metadata)) {
    message("Starting metadata download.")
    metadata <- download_eupath_metadata(...)
    ## metadata <- download_eupath_metadata(dir=dir)
    message("Finished metadata download.")
  }
  all_species <- metadata[["Species"]]
  entry <- NULL
  grep_hits <- grepl(species, all_species)
  grepped_hits <- all_species[grep_hits]
  if (species %in% all_species) {
    entry <- metadata[metadata[["Species"]] == species, ]
    message(paste0("Found: ", entry[["Species"]]))
  } else if (sum(grep_hits > 0)) {
    species <- grepped_hits[[1]]
    entry <- metadata[metadata[["Species"]] == species, ]
    message("Found the following hits: ", toString(grepped_hits), ", choosing the first.")
  } else {
    message(paste0("Here are the possible species: ", toString(all_species)))
    stop("Did not find your species.")
  }
  return(entry)
}

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
                                     dir="eupathdb", use_savefile=TRUE, ...) {
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
  ##base_url <- "http://eupathdb.org/eupathdb/webservices/"
  base_url <- paste0("http://", webservice, ".org/", webservice, "/webservices/")
  query_string <- "OrganismQuestions/GenomeDataTypes.json?o-fields=all"
  request_url <- paste0(base_url, query_string)

  ## retrieve organism metadata from EuPathDB
  metadata_json <- paste0(dir, "/metadata.json")
  file <- download.file(url=request_url, destfile=metadata_json, method="curl", quiet=FALSE)
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

#' Given 2 species names from the eupathdb, make orthology tables betwixt them.
#'
#' The eupathdb provides such a tremendous wealth of information.  For me
#' though, it is difficult sometimes to boil it down into just the bits of
#' comparison I want for 1 species or between 2 species.  A singularly common
#' question I am asked is: "What are the most similar genes between species x
#' and y among these two arbitrary parasites?"  There are lots of ways to poke
#' at this question: run BLAST/fasta36, use biomart, query the ortholog tables
#' from the eupathdb, etc...  However, in all these cases, it is not trivial to
#' ask the next question:  What about: a:b and b:a?
#' This function attempts to address that for the case of two eupath species
#' from the same domain. (tritrypdb/fungidb/etc.)  It does however assume that
#' the sqlite package has been installed locally, if not it suggests you run the
#' make_organismdbi function in order to do that.
#'
#' @param first  Species name from one eupath database.
#' @param second Another species name from the same eupath database.
#' @param abbrev  It is likely you do not care about the reciprocals and such,
#'   in that case the second species is irrelvant and just a species
#'   abbreviation is sufficient.
#' @param xref  Cross reference the stuff in the sqlite vs. the abbreviation?
#'   If no, then this will just return the entire ortholog table.
#'
#' The xref parameter reminds me, this function does make a peculiar attempt to
#'   handle a lot of corner cases without failure, but instead to try to always
#'   provide some useful information, even if that information is not precisely
#'   what you want, but the only thing I can find.
#' With that in mind, one could easily argue that this function is too long and
#'   it should be split into a couple pieces, but hopefully it is sufficiently
#'   clear that this is not necessary -- or I will do it later.
#' @return List containing the information available given the species provided
#'   and/or the abbreviation/xref.
#' @export
extract_eupath_orthologs <- function(first, second, abbrev=NULL, xref=TRUE, ...) {
  arglist <- list(...)

  ortholog_column <- "ORTHOLOGS_ORTHOLOG"
  if (!is.null(arglist[["ortholog_column"]])) {
    ortholog_column <- arglist[["ortholog_column"]]
  }

  ##metadata <- download_eupath_metadata(...)
  metadata <- download_eupath_metadata()
  first_name <- check_eupath_species(first)[["Species"]]
  second_name <- check_eupath_species(second)[["Species"]]

  first_pkg <- make_eupath_pkgnames(species=first_name, metadata=metadata)[["orgdb"]]
  second_pkg <- make_eupath_pkgnames(species=second_name, metadata=metadata)[["orgdb"]]

  num_found <- 0
  pkgs <- list()
  tt <- try(do.call("library", as.list(first_pkg)), silent=TRUE)
  if (class(tt) == "try-error") {
    message(paste0("Did not find the package: ",
                   first_pkg,
                   ". Will not be able to do reciprocal hits."))
    message(paste0("Perhaps try invoking the following: hpgltools::make_eupath_organismdbi('",
                   first_name, "')"))
  } else {
    num_found <- num_found + 1
    message(paste0("Loaded: ", first_pkg))
    pkgs[[first_name]] <- get(first_pkg)
  }

  tt <- try(do.call("library", as.list(second_pkg)), silent=TRUE)
  if (class(tt) == "try-error") {
    message(paste0("Did not find the package: ",
                   second_pkg,
                   ". Will not be able to do reciprocal hits."))
    message(paste0("Perhaps try invoking the following: hpgltools::make_eupath_organismdbi('",
                   second_name, "')"))
  } else {
    num_found <- num_found + 1
    message(paste0("Loaded: ", second_pkg))
    pkgs[[second_name]] <- get(second_pkg)
  }
  if (num_found == 0) {
    stop("Found neither package.  Cannot continue.")
  }

  all_orthos <- list()
  gene_lists <- list()
  for (i in 1:length(pkgs)) {
    name <- names(pkgs)[[i]]
    pkg <- pkgs[[i]]
    gene_lists[[name]] <- AnnotationDbi::keys(pkg)
    all_orthos[[name]] <- AnnotationDbi::select(x=pkg,
                                                keytype="GID",
                                                keys=gene_lists[[name]],
                                                columns=ortholog_column)
  }

  ## If a single pkg is found, then we cannot do a cross reference of two species,
  ## but we _can_ extract all genes which match a given species prefix pattern.
  retlist <- list()
  if (num_found == 1) {
    retlist[["all_genes"]] <- gene_lists[[1]]
    if (is.null(abbrev)) {
      message(paste0("Extracting all orthologs for ", names(pkgs)[[1]], "."))
      retlist[["orthologs"]] <- all_orthos[[1]]
      return(retlist)
    } else {
      ortho_idx <- grepl(pattern=paste0("^", abbrev), x=all_orthos[[1]][[ortholog_column]])
      if (sum(ortho_idx) == 0) {
        warning(paste0("The abbreviation: ", abbrev, " did not match any of the ortholog table."))
        warning("Returning the entire ortholog table so that you can figure out where things went wrong.")
        retlist[["orthologs"]] <- all_orthos[[1]]
      } else {
        tmp_table <- all_orthos[[1]]
        colnames(tmp_table) <- c(names(pkgs)[1], abbrev)
        retlist[["orthologs"]] <- tmp_table[ortho_idx, ]
      }
      return(retlist)
    }
  }

  ## We have finished handling the case where we have 1 table to search in
  ## Now what to do when we have two?
  ## Well, the way I envision this:
  ## 1. Use table 1's keys to set the prefix for table 2
  ## 2. Use table 2's ... for table 1
  ## 3. Index each of the two tables using the other key's prefix

  ## The following code should only therefore get hit if the above did _not_.
  ## In order to get the appropriate eupathdb-specific prefix, we want to gsub out the
  ## beginning of the genes from each package until either a '.' or '_'.
  first_string <- gene_lists[[1]][[1]]  ## Grab the first gene in each set of genes.
  second_string <- gene_lists[[2]][[1]]

  ## My pattern needs to greedily(use the ?) match at the beginning of the string(^)
  ## A set of characters(()) which matches anything(.) at least for one character(+)
  ## Then it needs to match a set(()) of either(|) a dot(\\.) or underscore(_)
  ## not greedily(no ?) followed by anything(.*) until the end of the string($).
  idx_pattern <- "^(.+?)(\\.|_).*$"
  ## When performing the replacement, pull out the first set of the 2 sets above(\\1)
  ## And it needs to acquire this information from the first gene in both sets of genes.
  first_idx_pattern <- gsub(pattern=idx_pattern,
                            replacement="\\1",
                            x=first_string,
                            perl=TRUE)
  second_idx_pattern  <- gsub(pattern=idx_pattern,
                              replacement="\\1",
                              x=second_string,
                              perl=TRUE)
  ## When doing the grep, make sure it searches from the beginning of each string(^).
  ## Also make sure to pull it from the second column of each of the two dataframes of orthologs.
  colnames(all_orthos[[1]]) <- c("first_gene", "second_ortholog") ## Explicitly rename the columns
  ## so that it is easier to see where I am pulling the index information.
  colnames(all_orthos[[2]]) <- c("second_gene", "first_ortholog")
  second_species_idx_in_first_db <- grepl(pattern=paste0("^", second_idx_pattern),
                                          x=all_orthos[[1]][, "second_ortholog"])
  first_species_idx_in_second_db <- grepl(pattern=paste0("^", first_idx_pattern),
                                       x=all_orthos[[2]][, "first_ortholog"])
  first_species_orthologs <- all_orthos[[1]][second_species_idx_in_first_db, ]
  second_species_orthologs <- all_orthos[[2]][first_species_idx_in_second_db, ]

  ## Now test for reciprocal best hittedness (yeah, that is not a word).
  ## To do this I will merge the two tables of orthologs and see where the ortholog makes a
  ## bridge from geneW_speciesA -> orthologX_speciesB -> geneW_speciesA and/or
  ##             geneY_speciesB -> orthologZ_speciesA -> geneY_speciesB
  shared_genes_first_second <- merge(first_species_orthologs,
                                     second_species_orthologs,
                                     by.x="first_gene",
                                     by.y="first_ortholog")
  shared_genes_second_first <- merge(second_species_orthologs,
                                     first_species_orthologs,
                                     by.x="second_gene",
                                     by.y="second_ortholog")
  reciprocal_idx_first <- shared_genes_first_second[["first_gene"]] ==
    shared_genes_first_second[["second_gene"]]
  reciprocals_first <- shared_genes_first_second[reciprocal_idx_first, ]
  reciprocal_idx_second <- shared_genes_second_first[["first_gene"]] ==
    shared_genes_second_first[["second_gene"]]
  reciprocals_second <- shared_genes_second_first[reciprocal_idx_second, ]

  ## While we are at it, lets get the set of genes in speciesA which do not show up in the
  ## orthologous list for speciesB and vice-versa
  first_non_orthologous_idx <- !(gene_lists[[1]] %in% all_orthos[[2]][["first_ortholog"]])
  first_non_orthologous <- gene_lists[[1]][first_non_orthologous_idx]
  second_non_orthologous_idx <- !(gene_lists[[2]] %in% all_orthos[[1]][["second_ortholog"]])
  second_non_orthologous <- gene_lists[[2]][first_non_orthologous_idx]
  retlist <- list(
    "non_orthologs_first_species" = first_non_orthologous,
    "non_orthologs_second_species" = second_non_orthologous,
    "reciprocal_match_first_second" = reciprocals_first,
    "reciprocal_match_second_first" = reciprocals_second,
    "first_genes_second_orthologs" = first_species_orthologs,
    "second_genes_first_orthologs" = second_species_orthologs,
    "first_orthologs_merged" = shared_genes_first_second,
    "second_orthologs_merged" = shared_genes_second_first)
  return(retlist)
}

#' Generate a BSgenome package from the eupathdb.
#'
#' Since we go to the trouble to try and generate nice orgdb/txdb/organismdbi packages, it
#' seems to me that we ought to also be able to make a readable genome package.  I should
#' probably use some of the logic from this to make the organismdbi generator smarter.
#'
#' @param species  Species to create.
#' @param entry  Single eupathdb metadata entry.
#' @param dir  Working directory.
#' @param reinstall  Rewrite an existing package directory.
#' @param metadata  Eupathdb metadata dataframe.
#' @param ... Extra arguments for downloading metadata when not provided.
#' @return  List of package names generated (only 1).
#' @export
make_eupath_bsgenome <- function(species="Leishmania major strain Friedlin", entry=NULL,
                                 dir="eupathdb", reinstall=FALSE, metadata=NULL, ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- download_eupath_metadata(dir=dir, ...)
      ## metadata <- download_eupath_metadata(dir=dir)
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  taxa <- make_taxon_names(entry)
  pkgnames <- make_eupath_pkgnames(species=species, metadata=metadata)
  pkgname <- pkgnames[["bsgenome"]]
  if (pkgname %in% installed.packages() & !isTRUE(reinstall)) {
    message(paste0(pkgname, " is already installed, set reinstall=TRUE if you wish to reinstall."))
    retlist <- list(
      "bsgenome_name" = pkgname
    )
    return(retlist)
  }

  ## Figure out the version numbers and download urls.
  db_version <- entry[["SourceVersion"]]
  fasta_start <- entry[["SourceUrl"]]
  fasta_starturl <- sub(pattern="gff",
                        replacement="fasta",
                        x=fasta_start,
                        perl=TRUE)
  fasta_url <- sub(pattern="\\.gff", replacement="_Genome\\.fasta",
                   x=fasta_starturl)
  fasta_hostname <- sub(pattern="http://(.*)\\.org.*$",
                        replacement="\\1",
                        x=fasta_start)
  genome_filename <- file.path(dir, paste0(pkgname, ".fasta"))

  ## Find a spot to dump the fasta files
  bsgenome_base <- file.path(dir)
  bsgenome_dir <- file.path(bsgenome_base, pkgname)
  if (!file.exists(bsgenome_dir)) {
    created <- dir.create(bsgenome_dir, recursive=TRUE)
  }
  ## Download them to this directory.
  downloaded <- download.file(url=fasta_url, destfile=genome_filename, quiet=FALSE)
  ## And extract all the individual chromosomes into this directory.
  input <- Biostrings::readDNAStringSet(genome_filename)
  output_list <- list()
  sequence_names <- "c("
  message("Writing chromosome files, this may take some time for especially fragmented scaffolds.")
  for (index in 1:length(input)) {
    chr <- names(input)[index]
    chr_name <- strsplit(chr, split=" ")[[1]][1]
    chr_file <- file.path(bsgenome_dir, paste0(chr_name, ".fa"))
    output <- Biostrings::writeXStringSet(input[index], chr_file, append=FALSE,
                                          compress=FALSE, format="fasta")
    output_list[[chr_name]] <- chr_file
    sequence_names <- paste0(sequence_names, '"', chr_name, '", ')
  }
  sequence_names <- gsub(pattern=", $", replacement=")", x=sequence_names)

  ## Now start creating the DESCRIPTION file
  desc_file <- file.path(bsgenome_dir, "DESCRIPTION")
  descript <- desc::description$new("!new")
  descript$set(Package=pkgname)
  author <- "Ashton Trey Belew <abelew@umd.edu>"
  title <- paste0(taxa[["genus"]], " ", taxa[["species"]], " strain ", taxa[["strain"]],
                  " version ", db_version)
  descript$set(Title=title)
  descript$set(Author=author)
  version_string <- format(Sys.time(), "%Y.%m")
  descript$set(Version=version_string)
  descript$set(Maintainer=author)
  descript$set(Description=paste0("A full genome from the eupathdb/tritrypdb for ", title, "."))
  descript$set(License="Artistic-2.0")
  descript$set(URL="http://eupathdb.org")
  descript$set(BugReports="http://github.com/elsayed-lab")
  descript$set(seqs_srcdir=bsgenome_dir)
  descript$set(seqnames=sequence_names)
  descript$set(organism=taxa[["taxon"]])
  descript$set(common_name=taxa[["genus_species"]])
  descript$set(provider=fasta_hostname)
  descript$set(provider_version=paste0(fasta_hostname, " ", db_version))
  descript$set(release_date=format(Sys.time(), "%Y%m%d"))
  descript$set(BSgenomeObjname=paste0(taxa[["genus_species"]], "_", taxa[["strain"]]))
  descript$set(release_name=db_version)
  descript$set(organism_biocview=paste0(taxa[["genus_species"]], "_", taxa[["strain"]]))
  descript$del("LazyData")
  descript$del("Authors@R")
  descript$del("URL")
  descript$del("BugReports")
  descript$del("Encoding")
  description_file <- file.path(bsgenome_dir, "DESCRIPTION")
  descript$write(description_file)

  ## Generate the package, this puts it into the cwd.
  message("Starting forgeBSgenomeDataPkg().")
  annoying <- sm(BSgenome::forgeBSgenomeDataPkg(description_file))

  if (class(annoying) != "try-error") {
    inst <- sm(try(devtools::install(pkgname, quiet=TRUE)))
  }

  retlist <- list()
  if (class(inst) != "try-error") {
    retlist[["bsgenome_name"]] <- pkgname
    ## Clean up a little.
    deleted <- unlink(x=bsgenome_dir, recursive=TRUE, force=TRUE)
    built <- sm(try(devtools::build(pkgname, quiet=TRUE)))
    if (class(built) != "try-error") {
      moved <- file.rename(paste0(pkgname, "_", version_string, ".tar.gz"),
                           paste0(bsgenome_dir, "_", version_string, ".tar.gz"))
      final_deleted <- unlink(x=pkgname, recursive=TRUE, force=TRUE)
    }
  } else {
    retlist <- inst
  }
  return(retlist)
}

#' Generate standardized package names for the various eupathdb species.
#'
#' In my test directory, I have a little for loop which randomly chooses a couple
#' of eupathdb species for which to try and generate genome/annotation packages.
#' I am on my 4th or 5th iteration of passing that loop and in all of them
#' I have found some new and exciting exception to how a strain should be named.
#' The default argument for this function shows the funniest one so far.
#' With that in mind, this function should provide consistent, valid package names.
#'
#' @param species  Guess.
#' @param metadata  Eupathdb metadata.
#' @param ...  Further arguments to pass to download_eupath_metadata()
#' @return  List of package names and some booleans to see if they have already been installed.
make_eupath_pkgnames <- function(species="Coprinosis.cinerea.okayama7#130",
                                 metadata=NULL, ...) {
  if (is.null(metadata)) {
    message("Starting metadata download.")
    metadata <- sm(download_eupath_metadata(...))
    ## metadata <- download_eupath_metadata(dir=dir)
    message("Finished metadata download.")
  }
  all_species <- metadata[["Species"]]
  entry <- NULL
  grep_hits <- grepl(species, all_species)
  grepped_hits <- all_species[grep_hits]
  if (species %in% all_species) {
    entry <- metadata[metadata[["Species"]] == species, ]
  } else if (sum(grep_hits > 0)) {
    species <- grepped_hits[[1]]
    entry <- metadata[metadata[["Species"]] == species, ]
    message("Found the following hits: ", toString(grepped_hits), ", choosing the first.")
  } else {
    message(paste0("Here are the possible species: ", toString(all_species)))
    stop("Did not find your species.")
  }

  taxa <- make_taxon_names(entry)
  first_char <- strsplit(taxa[["genus"]], split="")[[1]][[1]]
  pkg_list <- list(
    "bsgenome" = paste0("BSGenome.", taxa[["taxon"]], ".v", entry[["SourceVersion"]]),
    "bsgenome_installed" = FALSE,
    "organismdbi" = paste0(taxa[["taxon"]], ".v", entry[["SourceVersion"]]),
    "organismdbi_installed" = FALSE,
    "orgdb" = paste0("org.", first_char, taxa[["species_strain"]], ".v",
                     entry[["SourceVersion"]], ".eg.db"),
    "orgdb_installed" = FALSE,
    "txdb" = paste0("TxDb.", taxa[["genus"]], ".", taxa[["species_strain"]],
                    ".", entry[["DataProvider"]], ".v", entry[["SourceVersion"]]),
    "txdb_installed" = FALSE
  )

  inst <- as.data.frame(installed.packages())
  if (pkg_list[["bsgenome"]] %in% inst[["Package"]]) {
    pkg_list[["bsgenome_installed"]] <- TRUE
  }
  if (pkg_list[["organismdbi"]] %in% inst[["Package"]]) {
    pkg_list[["organismdbi_installed"]] <- TRUE
  }
  if (pkg_list[["orgdb"]] %in% inst[["Package"]]) {
    pkg_list[["orgdb_installed"]] <- TRUE
  }
  if (pkg_list[["txdb"]] %in% inst[["Package"]]) {
    pkg_list[["txdb_installed"]] <- TRUE
  }

  return(pkg_list)
}

#' Create an organismDbi instance for an eupathdb organism.
#'
#' @param species  A species in the eupathDb metadata.
#' @param entry  A row from the eupathdb metadataframe.
#' @param dir  Directory in which to build the packages.
#' @param reinstall  Overwrite existing data files?
#' @param metadata  Use a pre-existing metadata table, or download a new one.
#' @param kegg_abbreviation  For when we cannot automagically find the kegg species id.
#' @param ...  Extra arguments when downloading metadata.
#' @return  The result of attempting to install the organismDbi package.
#' @author  Keith Hughitt
#' @export
make_eupath_organismdbi <- function(species="Leishmania major strain Friedlin", entry=NULL,
                                    dir="eupathdb", reinstall=FALSE, metadata=NULL,
                                    kegg_abbreviation=NULL, exclude_join="ENTREZID",
                                    ...) {
  arglist <- list(...)
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- download_eupath_metadata(dir=dir, ...)
      ## metadata <- download_eupath_metadata(dir=dir)
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  taxa <- make_taxon_names(entry)
  pkgnames <- make_eupath_pkgnames(species=species, metadata=metadata)
  pkgname <- pkgnames[["organismdbi"]]
  if (isTRUE(pkgnames[["organismdbi_installed"]]) & !isTRUE(reinstall)) {
    message(paste0(pkgname, " is already installed, set reinstall=TRUE if you wish to reinstall."))
    retlist <- list(
      "organdb_name" = pkgname
    )
    return(retlist)
  }
  orgdb_name <- pkgnames[["orgdb"]]
  txdb_name <- pkgnames[["txdb"]]
  orgdb_ret <- make_eupath_orgdb(
    species=species,
    entry=entry,
    metadata=metadata,
    dir=dir,
    kegg_abbreviation=kegg_abbreviation,
    reinstall=reinstall)
  if (is.null(orgdb_ret)) {
    return(NULL)
  }
  txdb_ret <- make_eupath_txdb(
    species=species,
    entry=entry,
    metadata=metadata,
    dir=dir,
    reinstall=reinstall)
  if (is.null(txdb_ret)) {
    return(NULL)
  }

  tt <- sm(requireNamespace(orgdb_name))
  tt <- sm(requireNamespace(txdb_name))
  ##libstring <- paste0("library(", orgdb_name, ")")
  ##library_ret <- try(eval(parse(text=libstring)))
  ##if (class(library_ret) == "try-error") {
  ##  message("The orgdb failed to load.")
  ##  return(NULL)
  ##}
  test <- sm(do.call("library", as.list(orgdb_name)))
  ##libstring <- paste0("library(", txdb_name, ")")
  ##library_ret <- try(eval(parse(text=libstring)))
  ##if (class(library_ret) == "try-error") {
  ##  message("The orgdb failed to load.")
  ##  return(NULL)
  ##}
  test <- sm(do.call("library", as.list(txdb_name)))
  organism <- taxa[["taxon"]]
  required <- sm(requireNamespace("OrganismDbi"))

  message("Joining the txdb and orgdb objects.")
  version <- paste0(entry[["SourceVersion"]], ".0")
  count <- 0
  graph_data <- list()
  geneids_found <- "GID" %in% AnnotationDbi::keytypes(get(orgdb_name)) &&
    "GENEID" %in% AnnotationDbi::keytypes(get(txdb_name))
  if (isTRUE(geneids_found)) {
    count <- count + 1
    name <- paste0("join", count)
    graph_data[[name]] <- c(orgdb="GID",  txdb="GENEID")
    names(graph_data[[name]]) <- c(orgdb_name, txdb_name)
  }
  required <- sm(requireNamespace("GO.db"))
  required <- sm(try(attachNamespace("GO.db"), silent=TRUE))
  ## FIXME Theoretically we should no longer have columns with names like
  ## GO_GO
  orgdb_go_col <- "GO_GO_ID"
  goids_found <- "GOID" %in% AnnotationDbi::keytypes(get("GO.db")) &&
    orgdb_go_col %in% AnnotationDbi::keytypes(get(orgdb_name))
  if (isTRUE(goids_found)) {
    count <- count + 1
    name <- paste0("join", count)
    graph_data[[name]] <- c(GO.db="GOID", orgdb=orgdb_go_col)
    names(graph_data[[name]]) <- c("GO.db", orgdb_name)
  }
  ## FIXME Theoretically we should no longer have columns with names like
  ## PATHWAY_PATHWAY
  required <- sm(requireNamespace("reactome.db"))
  required <- sm(try(attachNamespace("reactome.db"), silent=TRUE))
  orgdb_path_col <- "PATHWAY_PATHWAY"
  reactomeids_found <- "REACTOMEID" %in% AnnotationDbi::keytypes(get("reactome.db")) &&
    orgdb_path_col %in% AnnotationDbi::keytypes(get(orgdb_name))
  if (isTRUE(reactomeids_found)) {
    count <- count + 1
    name <- paste0("join", count)
    graph_data[[name]] <- c(reactome.db="REACTOMEID", orgdb=orgdb_path_col)
    names(graph_data[[name]]) <- c("reactome.db", orgdb_name)
  }

  author <- as.character(entry[["Maintainer"]])
  maintainer <- as.character(entry[["Maintainer"]])
  final_dir <- file.path(dir, pkgname)
  if (file.exists(final_dir)) {
    if (isTRUE(reinstall)) {
      unlinkret <- unlink(x=final_dir, recursive=TRUE)
    } else {
      if (file.exists(paste0(final_dir, ".bak"))) {
        unlinkret <- unlink(x=paste0(final_dir, ".bak"),
                            recursive=TRUE)
      }
      renamed <- file.rename(from=final_dir, to=paste0(final_dir, ".bak"))
    }
  }
  version <- format(as.numeric(version), nsmall=1)

  tmp_pkg_dir <- file.path(dir)
  if (!file.exists(tmp_pkg_dir)) {
    dir.create(tmp_pkg_dir, recursive=TRUE)
  }
  version_string <- format(Sys.time(), "%Y.%m")
  organdb <- OrganismDbi::makeOrganismPackage(
                            pkgname=pkgname,
                            graphData=graph_data,
                            organism=organism,
                            version=version_string,
                            maintainer=maintainer,
                            author=author,
                            destDir=tmp_pkg_dir,
                            license="Artistic-2.0"
                          )
  organdb_path <- sm(clean_pkg(final_dir))
  organdb_path <- sm(clean_pkg(organdb_path, removal="_", replace="", sqlite=FALSE))
  organdb_path <- sm(clean_pkg(organdb_path, removal="_like", replace="like", sqlite=FALSE))
  if (class(organdb) == "list") {
    inst <- sm(try(devtools::install(organdb_path)))
    if (class(inst) != "try-error") {
      built <- sm(try(devtools::build(organdb_path, quiet=TRUE)))
      if (class(built) != "try-error") {
        final_deleted <- unlink(x=organdb_path, recursive=TRUE, force=TRUE)
      }
    }
  }
  final_organdb_name <- basename(organdb_path)
  retlist <- list(
    "orgdb_name" = orgdb_name,
    "txdb_name" = txdb_name,
    "organdb_name" = final_organdb_name
  )
  return(retlist)
}

#' Create an orgdb SQLite database from the tables in eupathdb.
#'
#' This now uses the new POST version of the eupathdb.  Theoretically it is better,
#' I am not yet convinced, but the QUERY version of the eupathdb apparently will not
#' be supported over time.
#'
#' @param species  A specific species ID to query
#' @param entry  If not provided, then species will get this, it contains all the information.
#' @param dir  Where to put all the various temporary files.
#' @param kegg_abbreviation  If known, provide the kegg abbreviation.
#' @param reinstall  Re-install an already existing orgdb?
#' @param metadata  Use an existing metadata table to get the entry?
#' @param ...  Extra parameters when searching for metadata
#' @return  Currently only the name of the installed package.  This should probably change.
make_eupath_orgdb <- function(species=NULL, entry=NULL, dir="eupathdb",
                              kegg_abbreviation=NULL, reinstall=FALSE, metadata=NULL, ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- download_eupath_metadata(dir=dir, ...)
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  taxa <- make_taxon_names(entry)
  pkgnames <- make_eupath_pkgnames(species=species, metadata=metadata)
  pkgname <- pkgnames[["orgdb"]]
  if (isTRUE(pkgnames[["orgdb_installed"]]) & !isTRUE(reinstall)) {
    message(paste0(pkgname, " is already installed, set reinstall=TRUE if you wish to reinstall."))
    retlist <- list(
      "orgdb_name" = pkgname
    )
    return(retlist)
  }

  if (!file.exists(dir)) {
    created <- dir.create(dir, recursive=TRUE)
  }

  do_kegg <- TRUE
  if (is.null(kegg_abbreviation)) {
    kegg_abbreviation <- get_kegg_orgn(paste0(taxa[["genus"]], " ", taxa[["species"]]))
    if (length(kegg_abbreviation) == 0) {
      do_kegg <- FALSE
    }
  }

  chosen_provider <- entry[["DataProvider"]]
  chosen_species <- entry[["Species"]]

  message("Downloading full annotation table.")
  gene_table <- try(post_eupath_annotations(species=chosen_species, entry=entry, dir=dir))
  if (class(gene_table) == "try-error") {
    gene_table <- data.frame()
  }
  message("Downloading GO table.")
  go_table <- try(post_eupath_go_table(species=chosen_species, entry=entry, dir=dir))
  if (class(go_table) == "try-error") {
    go_table <- data.frame()
  }
  message("Downloading ortholog table.")
  ortholog_table <- try(post_eupath_ortholog_table(species=chosen_species, entry=entry, dir=dir))
  if (class(ortholog_table) == "try-error") {
    ortholog_table <- data.frame()
  }
  message("Downloading interpro table.")
  interpro_table <- try(post_eupath_interpro_table(species=chosen_species, entry=entry, dir=dir))
  if (class(interpro_table) == "try-error") {
    interpro_table <- data.frame()
  }
  message("Downloading pathway table.")
  pathway_table <- try(post_eupath_pathway_table(species=chosen_species, entry=entry, dir=dir))
  if (class(pathway_table) == "try-error") {
    pathway_table <- data.frame()
  }
  message("Downloading kegg table.")
  kegg_table <- try(load_kegg_annotations(species=taxa[["genus_species"]],
                                          flatten=FALSE, abbreviation=kegg_abbreviation))
  if (class(kegg_table) == "try-error") {
    kegg_table <- data.frame()
  }
  if (nrow(gene_table) == 0) {
    warning("Unable to create an orgdb for this species.")
    return(NULL)
  }
  chromosome_table <- gene_table[, c("GID", "GENOMIC_SEQUENCE_ID")]
  colnames(chromosome_table) <- c("GID", "CHR_ID")
  type_table <- gene_table[, c("GID", "GENE_TYPE")]
  colnames(type_table) <- c("GID", "GENE_TYPE")

  ## Compile list of arguments for makeOrgPackage call
  version_string <- format(Sys.time(), "%Y.%m")
  orgdb_args <- list(
    "gene_info"  = gene_table,
    "chromosome" = chromosome_table,
    "type" = type_table,
    "version" = version_string,
    "author" = entry[["Maintainer"]],
    "maintainer" = entry[["Maintainer"]],
    "tax_id" = as.character(entry[["TaxonomyId"]]),
    "genus" = taxa[["genus"]],
    "species" = paste0(taxa[["species_strain"]], ".v", entry[["SourceVersion"]]),
    "outputDir" = dir)

  ## add non-empty tables
  if (is.null(go_table)) {
    message("This should not be possible, but the go table is still null.")
  } else if (nrow(go_table) > 0) {
    orgdb_args[["go"]] <- go_table
  }
  if (is.null(ortholog_table)) {
    message("This should not be possible, but the ortholog table is still null.")
  } else if (nrow(ortholog_table) > 0) {
    orgdb_args[["orthologs"]] <- ortholog_table
  }
  if (is.null(interpro_table)) {
    message("This should not be possible, but the interpro table is still null.")
  } else if (nrow(interpro_table) > 0) {
    orgdb_args[["interpro"]] <- interpro_table
  }
  if (is.null(pathway_table)) {
    message("This should not be possible, but the pathway table is still null.")
  } else if (nrow(pathway_table) > 0) {
    orgdb_args[["pathway"]] <- pathway_table
  }
  if (is.null(kegg_table)) {
    message("This should not be possible, but the kegg table is still null.")
  } else if (nrow(kegg_table) > 0) {
    orgdb_args[["kegg"]] <- kegg_table
  }

  ## Make sure no duplicated stuff snuck through, or makeOrgPackage throws an error.
  ## Make sure that every GID field is character, too
  ## -- otherwise you get 'The type of data in the 'GID' columns must be the same for all data.frames.'
  used_columns <- c()
  for (i in 1:length(orgdb_args)) {
    argname <- names(orgdb_args)[i]
    if (class(orgdb_args[[i]])[1] == "data.frame") {
      ## Make sure that the column names in this data frame are unique.
      ## This starts at 2 because the first column should _always_ by 'GID'
      for (cn in 2:length(colnames(orgdb_args[[i]]))) {
        colname <- colnames(orgdb_args[[i]])[cn]
        if (colname %in% used_columns) {
          new_colname <- paste0(toupper(argname), "_", colname)
          colnames(orgdb_args[[i]])[cn] <- new_colname
          used_columns <- c(used_columns, new_colname)
        } else {
          used_columns <- c(used_columns, colname)
        }
      }
      ## First swap out NA to ""
      na_tmp <- orgdb_args[[i]]
      na_set <- is.na(na_tmp)
      na_tmp[na_set] <- ""
      orgdb_args[[i]] <- na_tmp
      ## Then remove duplicated elements.
      orgdb_dups <- duplicated(orgdb_args[[i]])
      if (sum(orgdb_dups) > 0) {
        tmp <- orgdb_args[[i]]
        tmp <- tmp[!orgdb_dups, ]
        orgdb_args[[i]] <- tmp
      }
      ## Finally, make sure all GID columns are characters
      orgdb_args[[i]][["GID"]] <- as.character(orgdb_args[[i]][["GID"]])
    } ## End checking for data.frames
  }

  ## The following lines are because makeOrgPackage fails stupidly if the directory exists.
  backup_path <- file.path(dir, paste0(pkgname, ".bak"))
  first_path <- file.path(dir, pkgname)
  if (file.exists(backup_path)) {
    message(paste0(backup_path, " already exists, deleting it."))
    ret <- unlink(backup_path, recursive=TRUE)
  }
  if (file.exists(first_path)) {
    message(paste0(first_path, " already exists,
 backing it up."))
    ret <- file.rename(first_path, backup_path)
  }

  tt <- sm(require("AnnotationForge"))
  message(sprintf("- Calling makeOrgPackage for %s", entry[["Species"]]))
  orgdb_path <- sm(try(do.call("makeOrgPackage", orgdb_args)))
  if (class(orgdb_path) == "try-error") {
    return(NULL)
  }

  ## Fix name in sqlite metadata table
  dbpath <- file.path(
    orgdb_path, "inst/extdata", sub(".db", ".sqlite", basename(orgdb_path)))
  message(sprintf("- Fixing sqlite Orgdb sqlite database %s", dbpath))

  ## make sqlite database editable
  Sys.chmod(dbpath, mode="0644")
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=dbpath)
  ## update SPECIES field
  query <- sprintf('UPDATE metadata SET value="%s" WHERE name="SPECIES";', entry[["Species"]])
  sq_result <- RSQLite::dbSendQuery(conn=db, query)
  ## update ORGANISM field
  query <- sprintf('UPDATE metadata SET value="%s" WHERE name="ORGANISM";', entry[["Species"]])
  sq_result <- RSQLite::dbSendQuery(conn=db, query)
  ## lock it back down
  Sys.chmod(dbpath, mode="0444")

  ## Clean up any strangeness in the DESCRIPTION file
  orgdb_path <- sm(clean_pkg(orgdb_path))
  orgdb_path <- sm(clean_pkg(orgdb_path, removal="_", replace=""))
  orgdb_path <- sm(clean_pkg(orgdb_path, removal="_like", replace="like"))
  testthat::expect_equal(first_path, orgdb_path)
  ## And install the resulting package.
  inst <- sm(try(devtools::install(orgdb_path, quiet=TRUE)))
  if (class(inst) != "try-error") {
    built <- sm(try(devtools::build(orgdb_path, quiet=TRUE)))
    if (class(built) != "try-error") {
      final_deleted <- unlink(x=orgdb_path, recursive=TRUE, force=TRUE)
    }
  }

  ## return the path to the sqlite database
  retlist <- list(
    "orgdb_name" = pkgname
  )
  return(retlist)

}

#' Generate TxDb for EuPathDB organism
#'
#' @param species  guess
#' @param entry  One dimensional dataframe with organism metadata.
#' @param dir  Base directory for building the package.
#' @param reinstall  Overwrite an existing installed package?
#' @param metadata  dataframe of eupathdb metadata.
#' @param ...  Extra arguments for getting metadata.
#' @return TxDb instance name.
#' @author atb
#' @export
make_eupath_txdb <- function(species=NULL, entry=NULL, dir="eupathdb",
                             reinstall=FALSE, metadata=NULL, ...) {
  if (is.null(entry) & is.null(species)) {
    stop("Need either an entry or species.")
  } else if (is.null(entry)) {
    if (is.null(metadata)) {
      metadata <- sm(download_eupath_metadata(dir=dir, ...))
    }
    entry <- check_eupath_species(species=species, metadata=metadata)
  }
  taxa <- make_taxon_names(entry)
  pkgnames <- make_eupath_pkgnames(species=species, metadata=metadata)
  pkgname <- pkgnames[["txdb"]]

  if (isTRUE(pkgnames[["txdb_installed"]]) & !isTRUE(reinstall)) {
    message(paste0(pkgname, " is already installed, set reinstall=TRUE if you wish to reinstall."))
    retlist <- list(
      "txdb_name" = pkgname
    )
    return(retlist)
  }

  ## save gff as tempfile
  input_gff <- file.path(dir, paste0(pkgname, ".gff"))
  tt <- sm(download.file(url=entry[["SourceUrl"]], destfile=input_gff,
                         method="internal", quiet=FALSE))

  chr_entries <- read.delim(file=input_gff, header=FALSE, sep="")
  chromosomes <- chr_entries[["V1"]] == "##sequence-region"
  chromosomes <- chr_entries[chromosomes, c("V2", "V3", "V4")]
  colnames(chromosomes) <- c("ID", "start", "end")
  chromosome_info <- data.frame(
    "chrom" = chromosomes[["ID"]],
    "length" = as.numeric(chromosomes[["end"]]),
    "is_circular" = NA,
    stringsAsFactors=FALSE)

  txdb_metadata <- as.data.frame(t(entry))
  txdb_metadata[["name"]] <- rownames(txdb_metadata)
  colnames(txdb_metadata) <- c("value", "name")
  txdb_metadata <- txdb_metadata[, c("name", "value")]
  txdb <- sm(try(GenomicFeatures::makeTxDbFromGFF(
                                    file=input_gff,
                                    format="gff",
                                    chrominfo=chromosome_info,
                                    dataSource=entry[["SourceUrl"]],
                                    organism=paste0(taxa[["genus"]], " ", taxa[["species"]]),
                                    ## metadata=t(entry))
                                    )))
  if (class(txdb) == "try-error") {
    message("The txdb creation failed.")
    return(NULL)
  }

  ## This is the section I yanked
  provider <- GenomicFeatures:::.getMetaDataValue(txdb, "Data source")
  providerVersion <- GenomicFeatures:::.getTxDbVersion(txdb)
  dbType <- GenomicFeatures:::.getMetaDataValue(txdb, "Db type")
  authors <- GenomicFeatures:::.normAuthor(entry[["Maintainer"]], entry[["Maintainer"]])
  template_path <- system.file("txdb-template", package = "GenomicFeatures")
  version_string <- format(Sys.time(), "%Y.%m")
  symvals <- list(
    "PKGTITLE" = paste("Annotation package for", dbType, "object(s)"),
    "PKGDESCRIPTION" = paste("Exposes an annotation databases generated from",
                             GenomicFeatures:::.getMetaDataValue(txdb, "Data source"),
                             "by exposing these as", dbType, "objects"),
    "PKGVERSION" = version_string,
    "AUTHOR" = paste(authors, collapse = ", "),
    "MAINTAINER" = as.character(GenomicFeatures:::.getMaintainer(authors)),
    "GFVERSION" = GenomicFeatures:::.getMetaDataValue(txdb, "GenomicFeatures version at creation time"),
    "LIC" = "Artistic-2.0",
    "DBTYPE" = dbType,
    "ORGANISM" = GenomicFeatures:::.getMetaDataValue(txdb, "Organism"),
    "SPECIES" = GenomicFeatures:::.getMetaDataValue(txdb, "Organism"),
    "PROVIDER" = provider,
    "PROVIDERVERSION" = providerVersion,
    "RELEASEDATE" = GenomicFeatures:::.getMetaDataValue(txdb, "Creation time"),
    ## SOURCEURL = GenomicFeatures:::.getMetaDataValue(txdb, "Resource URL"),
    "SOURCEURL" = entry[["SourceUrl"]],
    "ORGANISMBIOCVIEW" = gsub(" ", "_", GenomicFeatures:::.getMetaDataValue(txdb, "Organism")),
    "TXDBOBJNAME" = pkgname)
  if (any(duplicated(names(symvals)))) {
    str(symvals)
    stop("'symvals' contains duplicated symbols")
  }
  is_OK <- sapply(symvals, S4Vectors::isSingleString)
  if (!all(is_OK)) {
    bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
    stop("values for symbols ", bad_syms, " are not single strings")
  }
  if (!file.exists(dir)) {
    tt <- dir.create(dir, recursive=TRUE)
  }

  pkg_list <- sm(Biobase::createPackage(
                            pkgname=pkgname,
                            destinationDir=dir,
                            originDir=template_path,
                            symbolValues=symvals,
                            unlink=TRUE))
  db_dir <- file.path(dir, pkgname, "inst", "extdata")
  if (!file.exists(db_dir)) {
    tt <- dir.create(db_dir, recursive=TRUE)
  }
  db_path <- file.path(db_dir, paste(pkgname, "sqlite", sep="."))

  obj <- try(AnnotationDbi::saveDb(txdb, file=db_path))
  if (class(obj) == "try-error") {
    warning("Failed to save the txdb object.")
  }

  install_dir <- file.path(dir, pkgname)
  install_dir <- sm(clean_pkg(install_dir))
  install_dir <- sm(clean_pkg(install_dir, removal="_", replace=""))
  install_dir <- sm(clean_pkg(install_dir, removal="_like", replace="like"))

  inst <- sm(try(devtools::install(install_dir, quiet=TRUE)))
  if (class(inst) != "try-error") {
    built <- sm(try(devtools::build(install_dir, quiet=TRUE)))
    if (class(built) != "try-error") {
      final_deleted <- unlink(x=install_dir, recursive=TRUE, force=TRUE)
    }
  }

  retlist <- list(
    "txdb_name" = pkgname)
  return(retlist)
}

#' Iterate through the various ways of representing taxon names
#'
#' Spend some time making sure they are valid, too.  Thus we want to get rid of weird characters
#' like hash marks, pipes, etc.
#'
#' @param entry  An entry of the eupathdb metadata.
#' @return  A list of hopefully valid nomenclature names to be used elsewhere in this family.
make_taxon_names <- function(entry) {
  taxon <- entry[["Species"]]
  species_parts <- unlist(strsplit(taxon, " "))
  genus <- species_parts[1]
  first <- toupper(substring(genus, 1, 1))
  species <- species_parts[2]
  strain <- ""
  if (length(species_parts) >= 3) {
    for (part in 3:length(species_parts)) {
      strain <- paste0(strain, ".", species_parts[part])
    }
  }
  strain <- gsub(pattern="^\\.", replacement="", x=strain)

  ## Start by getting rid of annoying text
  strain_pattern <- "(_| |\\.)strain(_| |\\.)"
  taxon <- gsub(pattern=strain_pattern, replacement="\\1", x=taxon)
  genus <- gsub(pattern=strain_pattern, replacement="\\1", x=genus)
  species <- gsub(pattern=strain_pattern, replacement="\\1", x=species)
  strain <- gsub(pattern=strain_pattern, replacement="\\1", x=strain)

  ## Replace some annoying characters with .
  annoying_pattern <- "(_|-|#| )"
  taxon <- gsub(pattern=annoying_pattern, replacement="\\.", x=taxon)
  genus <- gsub(pattern=annoying_pattern, replacement="\\.", x=genus)
  species <- gsub(pattern=annoying_pattern, replacement="\\.", x=species)
  strain <- gsub(pattern=annoying_pattern, replacement="\\.", x=strain)

  ## Completely remove the truly stupid characters.
  stupid_pattern <- "(\\,|\\/|\\?|\\|\\[|\\])"
  taxon <- gsub(pattern=stupid_pattern, replacement="", x=taxon)
  genus <- gsub(pattern=stupid_pattern, replacement="", x=genus)
  species <- gsub(pattern=stupid_pattern, replacement="", x=species)
  strain <- gsub(pattern=stupid_pattern, replacement="", x=strain)

  ## There are a few extra-weirdos with double-.s
  silly_pattern <- "\\.\\."
  taxon <- gsub(pattern=silly_pattern, replacement="\\.", x=taxon)
  genus <- gsub(pattern=silly_pattern, replacement="\\.", x=genus)
  species <- gsub(pattern=silly_pattern, replacement="\\.", x=species)
  strain <- gsub(pattern=silly_pattern, replacement="\\.", x=strain)

  species_strain <- paste0(species, ".", strain)
  genus_species <- paste0(genus, ".", species)

  species_strain <- paste(unlist(strsplit(taxon, split="\\."))[-1], collapse=".")
  genus_species <- paste0(genus, ".", species)

  gspecies <- paste0(first, species)
  gsstrain <- paste0(gspecies, strain)

  taxa <- list(
    "taxon" = taxon,
    "genus" = genus,
    "species" = species,
    "strain" = strain,
    "species_strain" = species_strain,
    "genus_species" = genus_species,
    "gspecies" = gspecies,
    "gsstrain" = gsstrain)
  return(taxa)
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

  query_body <- list(
    ## 3 elements, answerSpec, formatting, format.
    "answerSpec" = list(
      "questionName" = jsonlite::unbox(question),
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
  body <- jsonlite::toJSON(query_body)

  api_uri <- sprintf("http://%s.org/%s/service/answer", provider, uri_prefix)
  result <- httr::POST(
                    api_uri,
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
  column_names <- gsub(pattern=mypattern, replacement="\\1", x=stuff[["V1"]], perl=TRUE)
  ## Then get rid of any punctuation, as there is a column '# TM domains' -- that is bad.
  column_names <- gsub(pattern="[[:punct:]]", replacement="", x=column_names)
  ## Get rid of any extraneous spaces from removing punctuation, but since I cannot be certain
  ## that there is no punctuation in the middle of words, just look at the beginning of the strings.
  column_names <- gsub(pattern="^ +", replacement="", x=column_names)
  ## Finally, I do not accept column names with spaces.
  column_names <- gsub(pattern=" ", replacement="_", x=column_names)
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
#' @param provider Name of data provider to query (e.g. 'TriTrypDB')
#' @param query_body String of additional query arguments
#' @return list containing response from API request.
#'
#' More information
#' ----------------
#' 1. http://tritrypdb.org/tritrypdb/serviceList.jsp
#' @author Keith Hughitt
#' @export
post_eupath_table <- function(query_body, species=NULL, entry=NULL, metadata=NULL,
                              table_name=NULL, minutes=20) {

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
  result <- httr::POST(
                    api_uri,
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
#' @return  A big honking table.
post_eupath_annotations <- function(species="Leishmania major", entry=NULL,
                                    metadata=NULL, dir="eupathdb") {
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
    "primary_key", "sequence_id", "chromosome", "organism", "gene_type", "gene_location_text",
    "gene_name", "gene_exon_count", "is_pseudo", "gene_transcript_count", "gene_ortholog_number",
    "gene_paralog_number", "gene_orthomcl_name", "gene_total_hts_snps",
    "gene_hts_nonsynonymous_snps", "gene_hts_noncoding_snps", "gene_hts_stop_codon_snps",
    "gene_hts_nonsyn_syn_ratio", "uniprot_id", "gene_entrez_id", "transcript_length",
    "exon_count", "strand", "cds_length", "tm_count", "molecular_weight", "isoelectric_point",
    "signalp_scores", "signalp_peptide", "annotated_go_function", "annotated_go_process",
    "annotated_go_component", "annotated_go_id_function", "annotated_go_id_process",
    "annotated_go_id_component", "predicted_go_id_function", "predicted_go_id_process",
    "predicted_go_id_component", "ec_numbers", "ec_numbers_derived", "five_prime_utr_length",
    "three_prime_utr_length", "location_text", "gene_previous_ids", "protein_sequence", "cds")
  parameters <- list(
    "organism" = jsonlite::unbox(species),
    "min_molecular_weight" = jsonlite::unbox("1"),
    "max_molecular_weight" = jsonlite::unbox("10000000000000000")
  )
  result <- post_eupath_raw(entry,
                            question="GeneQuestions.GenesByMolecularWeight",
                            parameters=parameters,
                            columns=field_list)
  colnames(result) <- tolower(colnames(result))
  numeric_columns <- c(
    "exons_in_gene", "transcripts", "ortholog_count", "paralog_count", "total_snps_all_strains",
    "nonsynonymous_snps_all_strains", "noncoding_snps_all_strains", "snps_with_stop_codons_all_strains",
    "nonsynsyn_snp_ratio_all_strains", "transcript_length", "cds_length",
    "molecular_weight", "isoelectric_point", "annotated_5_utr_length", "annotated_3_utr_length")
  for (col in numeric_columns) {
    result[[col]] <- as.numeric(result[[col]])
  }
  factor_columns <- c(
    "chromosome", "gene_type", "is_pseudo", "gene_strand", "tm_domains", "exons_in_transcript")
  for (col in factor_columns) {
    result[[col]] <- as.factor(result[[col]])
  }
  colnames(result) <- toupper(colnames(result))
  return(result)
}

#'  Use the post interface to get GO data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @return  A big honking table.
post_eupath_go_table <- function(species="Leishmania major", entry=NULL,
                                 metadata=NULL, dir="eupathdb") {
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
  return(result)
}

#'  Use the post interface to get ortholog data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @return  A big honking table.
post_eupath_ortholog_table <- function(species="Leishmania major", entry=NULL,
                                       metadata=NULL, dir="eupathdb") {
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
        "includeEmptyTables" = jsonlite::unbox("true"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="orthologs")
  return(result)
}

#'  Use the post interface to get interpro data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @return  A big honking table.
post_eupath_interpro_table <- function(species="Leishmania major strain Friedlin", entry=NULL,
                                       metadata=NULL, dir="eupathdb") {
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

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="interpro")
  return(result)
}

#'  Use the post interface to get pathway data.
#'
#' @param species  guess.
#' @param entry  The full annotation entry.
#' @param metadata  A metadata table from which to get some annotation data.
#' @param dir  FIXME: I want to write some intermediate data to dir in case of transient error.
#' @return  A big honking table.
post_eupath_pathway_table <- function(species="Leishmania major", entry=NULL,
                                      metadata=NULL, dir="eupathdb") {
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
        "tables" = "MetabolicPathways",
        "includeEmptyTables" = jsonlite::unbox("true"),
        "attachmentType" = jsonlite::unbox("plain")
      ),
      "format" = jsonlite::unbox("tableTabular")
    ))

  result <- post_eupath_table(query_body, species=species, entry=entry, table_name="pathway")
  return(result)
}

## EOF
