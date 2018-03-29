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

#' Given 2 species names from the eupathdb, make orthology tables betwixt them.
#'
#' The eupathdb provides such a tremendous wealth of information.  For me
#' though, it is difficult sometimes to boil it down into just the bits of
#' comparison I want for 1 species or between 2 species.  A singularly common
#' question I am asked is: "What are the most similar genes between species x
#' and y among these two arbitrary parasites?"  There are lots of ways to poke
#' at this question: run BLAST/fasta36, use biomart, query the ortholog tables
#' from the eupathdb, etc.  However, in all these cases, it is not trivial to
#' ask the next question:  What about: a:b and b:a?
#' This function attempts to address that for the case of two eupath species
#' from the same domain. (tritrypdb/fungidb/etc.)  It does however assume that
#' the sqlite package has been installed locally, if not it suggests you run the
#' make_organismdbi function in order to do that.
#'
#' @param db  Species name (subset) from one eupath database.
#' @param master  Primary keytype to use for indexing the various tables.
#' @param query_species A list of exact species names to search for.  If uncertain
#'   about them, add print_speciesnames=TRUE and be ready for a big blob of
#'   text.  If left null, then it will pull all species.
#' @param id_column What column in the database provides the set of ortholog IDs?
#' @param org_column What column provides the species name?
#' @param url_column What column provides the orthomcl group ID?
#' @param count_column Name of the column with the count of species represented.
#' @param print_speciesnames Dump the species names for diagnostics?
#' @return A big table of orthoMCL families, the columns are:
#'   1.  GID: The gene ID
#'   2.  ORTHOLOG_ID: The gene ID of the associated ortholog.
#'   3.  ORTHOLOG_SPECIES: The species of the associated ortholog.
#'   4.  ORTHOLOG_URL: The OrthoMCL group ID's URL.
#'   5.  ORTHOLOG_COUNT: The number of all genes from all species represented in
#'   this group.
#'   6.  ORTHOLOG_GROUP: The family ID
#'   7.  QUERIES_IN_GROUP: How many of the query species are represented in this
#'   group?
#'   8.  GROUP_REPRESENTATION: ORTHOLOG_COUNT / the number of possible species.
#' @export
extract_eupath_orthologs <- function(db, master="GID",
                                     query_species=NULL,
                                     id_column="ORTHOLOG_ID",
                                     org_column="ORTHOLOGS_ORGANISM",
                                     url_column="ORTHOLOGS_ORTHOLOG_GROUP",
                                     count_column="ORTHOLOGS_ORTHOLOG_COUNT",
                                     print_speciesnames=FALSE) {

  load_pkg <- function(name, ...) {
    metadata <- download_eupath_metadata(...)
    ## metadata <- download_eupath_metadata()
    first_name <- check_eupath_species(name)[["Species"]]
    first_pkg <- make_eupath_pkgnames(species=first_name, metadata=metadata)[["orgdb"]]
    tt <- try(do.call("library", as.list(first_pkg)), silent=TRUE)
    if (class(tt) == "try-error") {
      message(paste0("Did not find the package: ",
                     first_pkg,
                     ". Will not be able to do reciprocal hits."))
      message(paste0("Perhaps try invoking the following: hpgltools::make_eupath_organismdbi('",
                     first_name, "')"))
      pkg <- NULL
    } else {
      message(paste0("Loaded: ", first_pkg))
      pkg <- get(first_pkg)
    }
    return(pkg)
  }

  pkg <- NULL
  if (class(db)[1] == "OrgDb") {
    pkg <- db
  } else if (class(first)[1] == "character") {
    pkg <- load_pkg(db)
  } else {
    stop("I only understand orgdbs or the name of a species.")
  }

  columns <- c(id_column, org_column, url_column, count_column)
  gene_set <- AnnotationDbi::keys(pkg, keytype=master)
  all_orthos <- AnnotationDbi::select(x=pkg,
                                      keytype=master,
                                      keys=gene_set,
                                      columns=columns)
  all_orthos[["ORTHOLOG_GROUP_ID"]] <- gsub(pattern="^.*>(.*)<\\/a>$",
                                            replacement="\\1",
                                            x=all_orthos[[url_column]])

  num_possible <- 1
  species_names <- unique(all_orthos[[org_column]])
  if (is.null(query_species)) {
    query_species <- species_names
  }
  num_possible <- length(species_names)
  message(paste0("There are ", num_possible, " possible species in this group."))

  if (isTRUE(print_speciesnames)) {
    print(toString(species_names))
    return(invisible())
  }

  ## Now pull out the species of interest
  found_species <- 0
  for (sp in query_species) {
    if (sp %in% all_orthos[[org_column]]) {
      message(paste0("Found species: ", sp))
    } else {
      message(paste0("Did not find species: ", sp))
    }
  }
  kept_orthos_idx <- all_orthos[[org_column]] %in% query_species
  kept_orthos <- all_orthos[kept_orthos_idx, ]
  colnames(kept_orthos) <- c(master, "ORTHOLOG_ID", "ORTHOLOG_SPECIES",
                             "ORTHOLOG_URL", "ORTHOLOG_COUNT", "ORTHOLOG_GROUP")
  kept_orthos[["ORTHOLOG_COUNT"]] <- as.integer(kept_orthos[["ORTHOLOG_COUNT"]])

  kept_orthos_dt <- data.table::as.data.table(kept_orthos) %>%
    dplyr::group_by_(master) %>%
    dplyr::add_count_(master)
  colnames(kept_orthos_dt) <- c(master, "ORTHOLOG_ID", "ORTHOLOG_SPECIES",
                             "ORTHOLOG_URL", "ORTHOLOG_COUNT", "ORTHOLOG_GROUP",
                             "QUERIES_IN_GROUP")

  kept_orthos_dt[["GROUP_REPRESENTATION"]] <- kept_orthos_dt[["ORTHOLOG_COUNT"]] / num_possible
  num_queries <- length(query_species)

  return(kept_orthos_dt)
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
  fasta_hostname <- sub(pattern="https://(.*)\\.org.*$",
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
  descript$set(URL="https://eupathdb.org")
  descript$set(BugReports="https://github.com/elsayed-lab")
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
  ## ortholog_table <- try(post_eupath_ortholog_table(species=chosen_species, entry=entry, dir=dir))
  ortholog_table <- try(get_orthologs_all_genes(species=chosen_species, entry=entry, dir=dir))
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
  ## message("Downloading kegg table.")
  ## kegg_table <- try(load_kegg_annotations(species=taxa[["genus_species"]],
  ##                                         flatten=FALSE, abbreviation=kegg_abbreviation))
  kegg_table <- data.frame()
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

  tt <- sm(require("AnnotationForge"), wrap=FALSE)
  message(sprintf("- Calling makeOrgPackage for %s", entry[["Species"]]))
  orgdb_path <- sm(try(do.call("makeOrgPackage", orgdb_args)), wrap=FALSE)
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

## EOF
