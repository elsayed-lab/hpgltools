#' TriTrypDB gene information table parser
#'
#' An example input file is the T. brucei Lister427 gene information table
#' available at:
#' http://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427/txt/data/TriTrypDB-5.0_TbruceiLister427Gene.txt
#'
#' @author Keith Hughitt
#' @param file Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe of gene info.
#' @export
parse_gene_info_table <- function(file, verbose=FALSE) {
  ## Create empty vector to store dataframe rows
  N <- 1e5
  go_rows <- data.frame("GO" = rep("", N),
                        "ONTOLOGY" = rep("", N),
                        "TERM" = rep("", N),
                        "SOURCE" = rep("", N),
                        "EVIDENCE" = rep("", N),
                        stringsAsFactors = FALSE)

  genedb <- data.frame("chromosome" = rep(NA, N),
                       "start" = rep(NA, N),
                       "stop" = rep(NA, N),
                       "strand" = rep("", N),
                       "type" = rep("", N),
                       "transcript_length" = rep(NA, N),
                       "cds_length" = rep(NA, N),
                       "pseudogene" = rep(NA, N),
                       "description" = rep("", N),
                       "mweight" = rep("", N),
                       "isopoint" = rep("", N),
                       stringsAsFactors = FALSE)
  ## Regular expression to extract location info
  location_regex <- "([0-9,]*) - ([0-9,]*) \\(([-+])\\)"
  gene_ids <- c()
  go_gene_ids <- c()
  gene_num <- 0
  go_num <- 0
  message("Using readr to read the txt file.")
  read_vec <- readr::read_lines(file, progress=TRUE)
  message("Starting to iterate over the txt file, this takes a long time.")
  chromosome <- ""
  seqid <- ""
  for (i in 1:length(read_vec)) {
    line <- read_vec[i]
    ## Example: Gene ID: TcCLB.397923.10
    if (grepl("^Gene ID", line)) {
      gene_num <- gene_num + 1
      gene_id <- local_get_value(line)
      gene_ids[gene_num] <- gene_id
      if (verbose) {
        if ((gene_num %% 100) == 0) {
          message(sprintf('Processing gene %d: %s', gene_num, gene_id))
        }
      }
    } else if(grepl("^Molecular Weight", line)) {
      ## Example: Molecular Weight: 37091
      if (grepl("^Molecular Weight: Not Assigned", line)) {
        mweight <- NA
      } else {
        mweight <- as.numeric(local_get_value(line))
      }
    } else if(grepl("^Isoelectric Point", line)) {
      ## Example: Isoelectric Point: 7.79
      if (grepl("^Isoelectric Point: Not Assigned", line)) {
        ipoint <- NA
      } else {
        ipoint <- as.numeric(local_get_value(line))
      }
    } else if (grepl("^Genomic Sequence ID", line)) {
      ## Example: Genomic Sequence ID: Tcruzi_56
      ## Genomic Sequence ID serves as chromosome in some genomes.
      ## And, since it comes before the Chromosome line, I will need a check in the chromosome assignment below.
      seqid <- as.character(local_get_value(line))
      ## message(paste0("Set seqid to ", seqid))
    } else if (grepl("^Chromosome", line)) { ## Example: Chromosome: Not Assigned
      if (grepl("^Chromosome: Not Assigned", line)) {
        chromosome <- seqid
      } else {
        chromosome <- as.character(local_get_value(line))
      }
    } else if (grepl("^Genomic Location:", line)) {
      ## Example: Genomic location: Tcruzi_56: 2 - 586 (+)
      result <- unlist(regmatches(line, regexec(location_regex, line)))
      gene_start <- as.numeric(gsub(",", "", result[2], fixed=TRUE))
      gene_stop  <- as.numeric(gsub(",", "", result[3], fixed=TRUE))
      strand <- result[4]
    } else if (grepl("^Gene Type", line)) {
      ## Example: Gene type: protein coding
      gene_type <- local_get_value(line)
    } else if (grepl("^Product Description", line)) {
      ## Example: Product Description: mucin-associated surface protein (MASP), putative
      description <- local_get_value(line)
    } else if (grepl("^Transcript Length", line)) {
      ## Example: Transcript length: 585
      transcript_length <- as.numeric(local_get_value(line))
    } else if (grepl("^CDS Length", line)) { ## Example: CDS length: 585
      val <- local_get_value(line)
      if (val == "null") {
        cds_length <- NA
      } else {
        cds_length <- as.numeric(val)
      }
    } else if (grepl("^Is Pseudo:", line)) {
      ## Pseudogene
      is_pseudo <- ifelse((local_get_value(line) == "Yes"), TRUE, FALSE)
    } else if (grepl("^GO:", line)) {
      ## Gene ontology terms
      go_num <- go_num + 1
      go_gene_ids[go_num] <- gene_id
      go_rows[go_num, ] <- c(head(unlist(strsplit(line, '\t')), 5))
      ##} else if (grepl("^PFAM", line)) { ## PFAM IDs
      ##    pfam_id <- unlist(strsplit(line, "\t"))[2]
      ##    pfam_ids <- paste0(pfam_id, " ", pfam_ids)
    } else if (grepl("^---", line)) {
      ## End of a gene's description
      ## message(paste0("Got to end of entry ", chromosome))
      ## Skip gene if it is not assigned to a chromosome
      if (is.na(chromosome)) {
        next
      }
      ## Otherwise add row to dataframe
      ## message("Adding elements:")
      ##element <- c(chromosome, gene_start, gene_stop, strand,
      ##             gene_type, transcript_length, cds_length, is_pseudo,
      ##             description, mweight, ipoint, pfam_ids)
      element <- c(chromosome, gene_start, gene_stop, strand,
                   gene_type, transcript_length, cds_length, is_pseudo,
                   description, mweight, ipoint)
      ## print(element)
      genedb[gene_num, ] <- element
    }  ## END ELSE!
  } ## End parsing the file
  rm(read_vec)

  go_rows <- go_rows[1:go_num, ]
  go_rows <- cbind("GID" = go_gene_ids,
                   go_rows)

  ## get rid of unallocated rows
  genedb <- genedb[1:gene_num, ]
  ## use gene id as row name
  rownames(genedb) <- gene_ids
  ## fix numeric types
  for (colname in c("start", "stop", "transcript_length", "cds_length")) {
    genedb[, colname] <- as.numeric(genedb[, colname])
  }
  for (colname in c("chromosome", "description")) {
    genedb[, colname] <- as.character(genedb[, colname])
  }
  ## sort data frame
  genedb <- genedb[with(genedb, order(chromosome, start)),]
  ret <- list(
    "genes" = genedb,
    "go" = go_rows)
  return(ret)
}

#' TriTrypDB gene information table GO term parser
#'
#' @author Keith Hughitt
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe where each line includes a gene/GO terms pair
#'         along with some addition information about the GO term. Note that
#'         because each gene may have multiple GO terms, a single gene ID may
#'         appear on multiple lines.
#' @export
parse_gene_go_terms <- function (filepath, verbose=FALSE) {
  ##require(tools)
  if (tools::file_ext(filepath) == "gz") {
    fp <- gzfile(filepath, open="rb")
  } else {
    fp <- file(filepath, open="r")
  }

  ## Create empty vector to store dataframe rows
  N <- 1e5
  gene_ids <- c()
  go_rows <- data.frame("GO"=rep("", N),
                        "ONTOLOGY"=rep("", N), TERM=rep("", N),
                        "SOURCE"=rep("", N), EVIDENCE=rep("", N),
                        stringsAsFactors=FALSE)
  ## Counter to keep track of row number
  i <- j <- 1
  ## Iterate through lines in file
  while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
    ## Gene ID
    if (grepl("^Gene ID", x)) {
      gene_id <- local_get_value(x)
      if (verbose) {
        message(sprintf("Processing gene %d: %s", i, gene_id))
      }
      i <- i + 1
    }

    ## Gene Ontology terms
    else if (grepl("^GO:", x)) {
      gene_ids[j] <- gene_id
      go_rows[j, ] <- c(head(unlist(strsplit(x, "\t")), 5))
      j <- j + 1
    }
  }
  ## get rid of unallocated rows
  go_rows <- go_rows[1:j-1, ]
  ## add gene id column
  go_rows <- cbind(GID=gene_ids, go_rows)
  ## close file pointer
  close(fp)
  return(go_rows)
}

#' Download the various data files from http://tritrypdb.org/
#'
#' The tritrypdb nicely makes their downloads standardized!
#'
#' @param version  What version of the tritrypdb to use?
#' @param species  Human readable species to use.
#' @param strain  Strain of the given species to download.
#' @param dl_dir  Directory into which to download the various files.
#' @param quiet  Print download progress?
#' @return  List of downloaded files.
#' @examples
#' \dontrun{
#'  filenames <- tritryp_downloads(species="lmajor", strain="friedlin", version="28")
#' }
#' @export
download_tritrypdb_text <- function(version="27", species="lmajor",
                                    strain="friedlin", dl_dir="organdb/tritryp", quiet=TRUE) {
  files_downloaded <- 0
  files_found <- 0

  if (!file.exists(dl_dir)) {
    dir.create(dl_dir, recursive=TRUE)
  }

  uc_species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
  uc_strain <- paste0(toupper(substr(strain, 1, 1)), substr(strain, 2, nchar(strain)))
  gff_filename <- paste0("TriTrypDB-", version, "_", uc_species, uc_strain, ".gff")

  gff_path <- paste0(dl_dir, "/", gff_filename)
  gff_url <- paste0("http://tritrypdb.org/common/downloads/release-", version,
                    "/", uc_species, uc_strain, "/gff/data/", gff_filename)
  if (file.exists(gff_path)) {
    files_found <- files_found + 1
  } else {
    message(paste0("Downloading ", gff_url))
    path <- try(download.file(url=gff_url, destfile=gff_path, quiet=quiet))
    if (class(path) != "try-error") {
      message("Finished downloading gff.")
      files_downloaded <- files_downloaded + 1
    } else {
      message("Failed downloading gff.")
    }
  }

  txt_filename <- paste0("TriTrypDB-", version, "_", uc_species, uc_strain, "Gene.txt")
  txt_url <- paste0("http://tritrypdb.org/common/downloads/release-", version,
                    "/", uc_species, uc_strain, "/txt/", txt_filename)
  txt_path <- paste0(dl_dir, "/", txt_filename)
  if (file.exists(txt_path)) {
    files_found <- files_found + 1
  } else {
    message(paste0("Downloading ", txt_url))
    path <- try(download.file(url=txt_url, destfile=txt_path, quiet=quiet))
    if (class(path) != "try-error") {
      message("Finished downloading txt file.")
      files_downloaded <- files_downloaded + 1
    } else {
      message("Failed downloading txt file.")
    }
  }

  fasta_filename <- paste0("TriTrypDB-", version, "_", uc_species, uc_strain, "_Genome.fasta")
  fasta_url <- paste0("http://tritrypdb.org/common/downloads/release-", version,
                      "/", uc_species, uc_strain, "/fasta/data/", fasta_filename)
  fasta_path <- paste0(dl_dir, "/", fasta_filename)
  if (file.exists(fasta_path)) {
    files_found <- files_found + 1
  } else {
    message(paste0("Downloading ", fasta_url))
    path <- try(download.file(url=fasta_url, destfile=fasta_path, quiet=quiet))
    if (class(path) != "try-error") {
      message("Finished downloading fasta file.")
      files_downloaded <- files_downloaded + 1
    } else {
      message("Failed downloading fasta file.")
    }
  }

  alias_filename <- paste0("TriTrypDB-", version, "_", uc_species, uc_strain, "_GeneAliases.txt")
  alias_url <- paste0("http://tritrypdb.org/common/downloads/release-", version,
                      "/", uc_species, uc_strain, "/txt/", alias_filename)
  alias_path <- paste0(dl_dir, alias_filename)
  if (file.exists(alias_path)) {
    files_found <- files_found + 1
  } else {
    message(paste0("Downloading ", alias_url))
    path <- try(download.file(url=alias_url, destfile=alias_path, quiet=TRUE))
    if (class(path) != "try-error") {
      message("Finished downloading alias file.")
      files_downloaded <- files_downloaded + 1
    } else {
      message("Failed downloading alias file.")
    }
  }

  ret <- list(
    "gff" = gff_path,
    "txt" = txt_path,
    "genome" = fasta_path,
    "aliases" = alias_path,
    "found" = files_found,
    "downloaded" = files_downloaded)
  return(ret)
}

#' Use taxize to get ncbi taxon IDs
#'
#' taxize looks like it might be awesome, but it is also pretty annoying
#'
#' @param species  Human readable species name
#' @return potential NCBI taxon IDs
#' @seealso \pkg{taxize}
#' @examples
#' \dontrun{
#'  taxonid <- get_ncbi_taxonid(species="Trypanosoma cruzi")
#' }
#' @export
get_ncbi_taxonid <- function(species="Leishmania major") {
  taxid <- taxize::get_eolid(sciname=species)[[1]]
  return(taxid)
}

#' Create an OrganismDbi for a species at the TriTrypDb
#'
#' OrganismDbi instances are pretty neat, they pull together OrgDb and TxDb.  With any luck, this
#' function provides the ability to pull together all the data from the TriTrypDb, GO.db, and
#' KEGGREST in order to accomplish these peculiar tasks.
#'
#' @param id  Unique tritrypdb identifier.
#' @param cfg  A configuration dataframe, when null it will be replaced by reading a csv file in
#'  inst/extdata.
#' @param output_dir  The directory into which to put the various intermediate files, including
#'  downloads from the TriTrypdb, the created OrgDb and TxDb instances, and the final
#'  OrganismDbi.
#' @param ...  Extra arguments including a boolean for whether to include kegg.
#' @return A path, some data files, and a kitty!
#' @seealso \pkg{AnnotationForge} \pkg{OrganismDbi}
#' @examples
#' \dontrun{
#'  crazytown <- make_organismdbi()  ## wait a loong time
#' }
#' @export
make_tritrypdb_organismdbi <- function(id="lmajor_friedlin",
                                       cfg=NULL, output_dir="organdb/tritryp", ...) {
  warning("make_tritrypdb_organismdbi is deprecated in favor of make_eupath_organismdbi.")
  warning("it works for now, but will die a sad, lonely death.")
  arglist <- list(...)
  reinstall <- TRUE
  if (!is.null(arglist[["reinstall"]])) {
    reinstall <- arglist[["reinstall"]]
  }
  kegg <- arglist[["kegg"]]
  cfg <- get_eupath_config(cfg)
  cfg_line <- cfg[id, ]
  version <- as.character(cfg_line[["db_version"]])
  strain <- as.character(cfg_line[["strain"]])
  shortname <- as.character(cfg_line[["shortname"]])
  files <- download_tritrypdb_text(version=version, species=shortname, strain=strain)
  ##files <- tritryp_downloads(version=version, species=shortname, strain=strain, ...)
  ##savefile <- paste0(output_dir, "/", cfg_line[["id"]], ".rda")
  orgdb_info <- make_orgdb_info(files[["gff"]], txt=files[["txt"]])
  message("Starting make_orgdb.")
  orgdb_result <- make_tritrypdb_orgdb(orgdb_info, id=id, cfg=cfg_line,
                                       output_dir=output_dir, kegg=kegg)
  orgdb_package <- orgdb_result[["package_name"]]
  message("Starting make_txdb")
  txdb_result <- make_tritrypdb_txdb(orgdb_info, cfg_line=cfg_line,
                                     gff=files[["gff"]], output_dir=output_dir)
  txdb_package <- txdb_result[["package_name"]]

  graph_data <- list(
    "join1" = c(GO.db="GOID", orgdb="GO"),
    "join2" = c(orgdb="GID",  txdb="GENEID")
  )

  names(graph_data[["join1"]]) = c("GO.db", orgdb_package)
  names(graph_data[["join2"]]) = c(orgdb_package, txdb_package)

  requireNamespace(orgdb_package)
  requireNamespace(txdb_package)
  libstring <- paste0("library(", orgdb_package, ")")
  eval(parse(text=libstring))
  libstring <- paste0("library(", txdb_package, ")")
  eval(parse(text=libstring))
  organism <- paste0(cfg_line[["genus"]], " ",
                     gsub(pattern="[[:punct:]]", replacement="", x=cfg_line[["species"]]),
                     " ",
                     gsub(pattern="[[:punct:]]", replacement="", x=cfg_line[["strain"]]))
  organism <- gsub(pattern="like", replacement="", x=organism)
  requireNamespace("OrganismDbi")
  pkgname <- paste0(as.character(cfg_line[["organismdb_name"]]), ".",
                    cfg_line[["db_name"]], cfg_line[["db_version"]])
  author <- as.character(cfg_line[["author"]])
  maintainer <- as.character(cfg_line[["maintainer"]])

  destination <- paste0(output_dir, "/organismdbi")
  if (file.exists(destination)) {
    unlink(x=destination, recursive=TRUE)
  }
  dir.create(destination, recursive=TRUE)
  version <- format(as.numeric(version), nsmall=1)
  organdb <- OrganismDbi::makeOrganismPackage(
                            pkgname = pkgname,
                            graphData = graph_data,
                            organism = organism,
                            version = version,
                            maintainer = maintainer,
                            author = author,
                            destDir = destination,
                            license = "Artistic-2.0"
                          )
  organdb_path <- paste0(destination, "/", pkgname)
  organdb_path <- clean_pkg(organdb_path)
  organdb_path <- clean_pkg(organdb_path, removal="_like", replace="like")
  if (class(organdb) == "list") {
    inst <- devtools::install(organdb_path)
  }
  return(inst)
}

#' Cleans up illegal characters in packages generated by make_organismdbi(), make_orgdb(),
#' and make_txdb(). This attempts to fix some of the common problems therein.
#'
#' OrganismDbi instances are pretty neat, they pull together OrgDb and TxDb.  With any luck, this
#' function provides the ability to pull together all the data from the TriTrypDb, GO.db, and
#' KEGGREST in order to accomplish these peculiar tasks.
#'
#' @param path  Location for the original Db/Dbi instance.
#' @param removal  String to remove from the instance.
#' @param replace  What to replace removal with, when necessary.
#' @return  A new OrgDb/TxDb/OrganismDbi
#' @examples
#' \dontrun{
#'  crazytown <- make_organismdbi()  ## wait a loong time
#' }
#' @export
clean_pkg <- function(path, removal="-like", replace="", sqlite=TRUE) {
  ## This is because TxDb creation fails if you have an author like 'abelew <abelew@gmail.com>'
  ##at_cmd <- paste0("sed -i 's/ at /\\@/g' ", path, "/DESCRIPTION")
  basedir <- dirname(path)
  dir <- basename(path)
  full_path <- file.path(basedir, dir)
  at_cmd <- paste0("perl -p -i -e 's/ at /\\@/g' ", full_path, "/DESCRIPTION")
  message(paste0("Rewriting DESCRIPTION: ", at_cmd))
  system(command=at_cmd)
  ## Since I changed @ to at I figured . could be dot too
  ##dot_cmd <- paste0("sed -i 's/ dot /\\./g' ", path, "/DESCRIPTION")
  dot_cmd <- paste0("perl -p -i -e 's/ dot /\\./g' ", full_path, "/DESCRIPTION")
  message(paste0("Rewriting DESCRIPTION to remove dot: ", dot_cmd))
  system(dot_cmd)
  new_dir <- dir
  new_path <- file.path(basedir, new_dir)
  if (grepl(pattern=removal, x=dir)) {
    ## Get rid of the -like in the path name
    new_dir <- gsub(pattern=removal, replacement=replace, x=dir)
    new_path <- file.path(basedir, new_dir)
    ## And rename the directory
    mv_cmd <- paste0("mv ", path, " ", new_path)
    message(paste0("moving orgdb: ", mv_cmd))
    system(mv_cmd)
    ## Collect the text files in the new package and remove all -like instances in them
    ##find_cmd <- paste0("sed -i 's/", removal, "/", replace,
    find_cmd <- paste0("perl -p -i -e 's/",
                       removal, "/", replace,
                       "/g' $(find ", new_path,
                       " -type f | grep -v 'sqlite' | grep -v 'zzz' | grep -v 'rda')")
    message(paste0("rewriting orgdb files: ", find_cmd))
    system(find_cmd)

    if (isTRUE(sqlite)) {
      ## Move the sqlite file, now the directory has been renamed.
      ## So when we go to move it we need to take that into account.
      old_sqlite_base <- gsub(pattern=".db", replacement="", x=dir)
      sqlite_basename <- basename(dir)
      sqlite_basename <- gsub(pattern=".sqlite", replacement="", x=sqlite_basename)
      ## sqlite_name <- gsub(pattern=".db", replacement=".sqlite", x=sqlite_basename)
      old_sqlite_file <- file.path(new_dir, "inst", "extdata", paste0(old_sqlite_base, ".sqlite"))
      old_sqlite <- file.path(basedir, old_sqlite_file)
      new_sqlite_file <- gsub(pattern=removal, replacement=replace, x=old_sqlite_file)
      new_sqlite <- file.path(basedir, new_sqlite_file)
      sqlite_mv_cmd <- paste0("mv ", old_sqlite, " ", new_sqlite)
      message(paste0("moving sqlite file: ", sqlite_mv_cmd))
      system(sqlite_mv_cmd)
      ## orgdb_dir <- new_dir
      new_pkg_name <- gsub(pattern=removal, replacement=replace, x=sqlite_basename)
      ## Update the orgdb sqlite file to reflect the new name
      final_sqlite_cmd <- paste0("chmod +w ", new_sqlite, " ; sqlite3 ", new_sqlite,
                                 " \"UPDATE metadata SET value='", new_pkg_name,
                                 "' WHERE name='SPECIES';\" ; chmod -w ", new_sqlite)
      message(paste0("rewriting sqlite db:", final_sqlite_cmd))
      system(final_sqlite_cmd)
    }
  }
  message("Finished clean_pkg.")
  return(new_path)
}

#' Make an orgDb object from some information provided by make_orgdb_info()
#'
#' An orgDb object should provide some useful annotation data including fun stuff like gene
#' ontology, kegg, etc.  In the case of the species at the TriTrypDb, much of this information
#' is available in the species .txt file.  This function takes that data and collates it into the
#' final orgDb objects using AnnotationForge.  It then makes some attempts to ensure that the
#' resulting material created in the filesystem conforms to specifications which allow one to have
#' multiple strains, etc.  Finally, if everything goes according to plan, it calls
#' devtools::install() and installs the resulting package.
#'
#' @param orgdb_info  List of data frames generated by make_orgdb_info()
#' @param id  Human readable species identifier, keys off the cfg data frame.
#' @param cfg  Configuration data extracted either from inst/eupath_configuration.csv or provided by
#'     the user.
#' @param kegg  Attempt adding kegg data?
#' @param output_dir  Base output directory for the resulting packages.
#' @param ... Args to pass through.
#' @return List of the resulting package name(s) and whether they installed.
#' @seealso \pkg{AnnotationForge} \pkg{devtools}
#'  \code{\link[AnnotationForge]{makeOrgPackage}}
#' @examples
#' \dontrun{
#'  orgdb_installedp <- make_orgdb(id="tcruzi_clbrener")
#' }
#' @export
make_tritrypdb_orgdb <- function(orgdb_info, id="lmajor_friedlin", cfg=NULL,
                                 kegg=TRUE, output_dir="organismdbi", ...) {
  arglist <- list(...)
  orgdb_pre <- paste0(output_dir, "/orgdb")
  if (!file.exists(orgdb_pre)) {
    dir.create(orgdb_pre, recursive=TRUE)
  }

  gene_info <- orgdb_info[["gene_info"]]
  gene_info[["GID"]] <- as.character(gene_info[["GID"]])
  rownames(gene_info) <- gene_info[["GID"]]
  chr_info <- orgdb_info[["chr_info"]]  ## Maps gene IDs to chromosomes
  chr_info[["GID"]] <- as.character(chr_info[["GID"]])
  go_info <- orgdb_info[["go_info"]]
  go_info[["GID"]] <- as.character(go_info[["GID"]])
  go_info <- go_info[, c("GID", "GO", "EVIDENCE")]
  gene_types <- orgdb_info[["gene_types"]]
  gene_types[["GID"]] <- as.character(gene_types[["GID"]])
  chromosome_info <- orgdb_info[["chromosome_info"]]  ## The lengths of the chromosomes

  shortname <- stringi::stri_trans_general(str=cfg[["shortname"]], id="Title")
  orgdb_base_name <- paste0("org.",
                            gsub(x=shortname, pattern="[[:punct:]]", replacement=""),
                            ".",
                            gsub(x=cfg[["strain"]], pattern="[[:punct:]]", replacement=""),
                            ".eg")
  
  orgdb_pkg_name <- paste0(orgdb_base_name, ".db")
  orgdb_sqlite_name <- paste0(orgdb_base_name, ".sqlite")
  assumed_dir <- paste0(orgdb_pre, "/", orgdb_pkg_name)
  message(paste0("The directory will be: ", assumed_dir))
  if (file.exists(assumed_dir)) {
    try(unlink(x=orgdb_pre, recursive=TRUE), silent=TRUE)
    try(unlink(x=orgdb_sqlite_name), silent=TRUE)
    try(dir.create(orgdb_pre, recursive=TRUE), silent=TRUE)
  }

  ## We need to ensure that none of the inputs for makeOrgPackage have duplicated rows.
  ## It looks like go/kegg are the most likely candidates for this particular problem.
  ## Also, make all the NAs to ""
  test_gene_info <- duplicated(gene_info)
  gene_info <- gene_info[!test_gene_info, ]
  test_gene_info <- is.na(gene_info)
  gene_info[test_gene_info] <- ""
  test_chr_info <- duplicated(chr_info)
  chr_info <- chr_info[!test_chr_info, ]
  test_go_info <- duplicated(go_info)
  go_info <- go_info[!test_go_info, ]
  test_types <- duplicated(gene_types)
  gene_types <- gene_types[!test_types, ]
  orgdb_dir <- NULL

  clean_genus <- as.character(gsub(pattern="[[:punct:]]", replacement="", x=cfg[["genus"]]))
  clean_species <- paste0(as.character(gsub(pattern="[[:punct:]]", replacement="", x=cfg[["species"]])),
                          ".",
                          as.character(gsub(pattern="[[:punct:]]", replacement="", x=cfg[["strain"]]))
                          )
  if (isTRUE(kegg)) {
    kegg_species <- paste0(cfg[["genus"]], " ", cfg[["species"]])
    kegg_info <- get_kegg_genepaths(species=kegg_species)
    test_kegg_info <- duplicated(kegg_info)
    kegg_info <- kegg_info[!test_kegg_info, ]
    kegg_info[["GID"]] <- as.character(kegg_info[["GID"]])
    orgdb_dir <- AnnotationForge::makeOrgPackage(
                                    gene_info = gene_info,
                                    chromosome = chr_info,
                                    go = go_info,
                                    kegg_info = kegg_info,
                                    type = gene_types,
                                    version = format(as.numeric(cfg[["db_version"]]), nsmall=1),
                                    author = as.character(cfg[["author"]]),
                                    maintainer = as.character(cfg[["maintainer"]]),
                                    ## Maybe use taxize for this and remove from the csv?
                                    tax_id = as.character(cfg[["tax_id"]]),
                                    genus = clean_genus,
                                    species = clean_species,
                                    outputDir = orgdb_pre,
                                    goTable = "go")
  } else {
    orgdb_dir <- AnnotationForge::makeOrgPackage(
                                    gene_info = gene_info,
                                    chromosome = chr_info,
                                    go = go_info,
                                    type = gene_types,
                                    version = format(as.numeric(cfg[["db_version"]]), nsmall=1),
                                    author = as.character(cfg[["author"]]),
                                    maintainer = as.character(cfg[["maintainer"]]),
                                    tax_id = as.character(cfg[["tax_id"]]),
                                    genus = clean_genus,
                                    species = clean_species,
                                    outputDir = orgdb_pre,
                                    goTable = "go")
  }

  orgdb_dir <- clean_pkg(orgdb_dir)
  if (cfg[["strain"]] == "CLBrenerNon-Esmeraldo-like") {
    orgdb_dir <- sm(clean_pkg(orgdb_dir, removal="Non-Esmeraldo", replace="NonEsmeraldo"))
  }

  inst <- FALSE
  if (!is.null(orgdb_dir)) {
    inst <- sm(try(devtools::install(orgdb_dir)))
  }
  orgdb_pkg_name <- basename(orgdb_dir)
  ## The result is the pathname of the created orgdb directory
  ret <- list(
    package_name = orgdb_pkg_name,
    result = inst)
  return(ret)
}

#' Create a TxDb object given data provided by make_orgdb_info()
#'
#' Much like make_orgdb() above, this uses the same data to generate a TxDb object.
#'
#' @param orgdb_info  List of data frames generated by make_orgdb_info().
#' @param cfg_line  Configuration data frame as per make_orgdb.
#' @param gff  File to read
#' @param from_gff  Use a gff file?
#' @param output_dir  Place to put rda intermediates.
#' @param ...   Extra arguments to pass through.
#' @return List of the resulting txDb package and whether it installed.
#' @seealso \pkg{GenomicFeatures} \pkg{Biobase} \pkg{devtools}
#'  \code{\link[Biobase]{createPackage}}
#' @examples
#' \dontrun{
#'  txdb <- make_txdb(orgdb_output)
#' }
#' @export
make_tritrypdb_txdb <- function(orgdb_info, cfg_line, gff=NULL,
                                from_gff=FALSE, output_dir="organismdbi", ...) {
  ## Sections of this were stolen from GenomicFeatures
  ## because it hates me.
  arglist <- list(...)

  destination <- output_dir
  chromosome_info <- orgdb_info[["chromosome_info"]]
  requireNamespace("GenomicFeatures")
  destination <- paste0(destination, "/txdb")
  db_version <- format(as.numeric(cfg_line[["db_version"]]), nsmall=1)
  maintainer <- as.character(cfg_line[["maintainer"]])
  maintainer <- gsub(pattern=" at ", replacement="@", x=maintainer)
  maintainer <- gsub(pattern=" dot ", replacement=".", x=maintainer)
  author <- as.character(cfg_line[["author"]])
  author <- gsub(pattern=" at ", replacement="@", x=author)
  author <- gsub(pattern=" dot ", replacement=".", x=author)
  db_url <- as.character(cfg_line[["db_url"]])
  clean_species <- gsub(pattern="[[:punct:]]", replacement="", x=cfg_line[["species"]])
  clean_strain <- gsub(pattern="[[:punct:]]", replacement="", x=cfg_line[["strain"]])
  package_name <<- paste0("TxDb.",
                          cfg_line[["shortname"]],
                          ".",
                          clean_strain,
                          ".", cfg_line[["db_name"]],
                          cfg_line[["db_version"]])

  txdb <- NULL
  if (!is.null(gff)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(
                               file=gff,
                               format="gff3",
                               chrominfo=chromosome_info,
                               ## exonRankAttributeName=NA,
                               dataSource=paste0(cfg_line[["db_name"]], "_", cfg_line[["db_version"]]),
                               organism=paste0(cfg_line[["genus"]], " ", clean_species))
    package_name <<- paste0("TxDb.", cfg_line[["shortname"]],
                            ".", clean_strain,
                            ".", cfg_line[["db_name"]],
                            cfg_line[["db_version"]])
  }

  ## This is the section I yanked
  provider <- GenomicFeatures:::.getMetaDataValue(txdb, "Data source")
  providerVersion <- GenomicFeatures:::.getTxDbVersion(txdb)
  dbType <- GenomicFeatures:::.getMetaDataValue(txdb, "Db type")
  authors <- GenomicFeatures:::.normAuthor(author, maintainer)
  template_path <- system.file("txdb-template", package = "GenomicFeatures")
  symvals <- list(
    "PKGTITLE" = paste("Annotation package for", dbType, "object(s)"),
    "PKGDESCRIPTION" = paste("Exposes an annotation databases generated from",
                             GenomicFeatures:::.getMetaDataValue(txdb, "Data source"), "by exposing these as", dbType, "objects"),
    "PKGVERSION" = db_version,
    ## "AUTHOR" = as.character(paste(authors, collapse = ", ")),
    "AUTHOR" = "Keith Hughitt, Ashton Trey Belew",
    "AUTHORS@R" =
      'c(person("Ashton", "Belew", email="abelew@gmail.com", role=c("aut","cre")),
                      person("Keith", "Hughitt", email="khughitt@umd.edu", role=c("aut")))',
    "MAINTAINER" = "Ashton Trey Belew <abelew@gmail.com>",
    ## "MAINTAINER" = as.character(GenomicFeatures:::.getMaintainer(authors)),
    "GFVERSION" = GenomicFeatures:::.getMetaDataValue(txdb, "GenomicFeatures version at creation time"),
    "LIC" = "Artistic-2.0",
    "DBTYPE" = dbType,
    "ORGANISM" = GenomicFeatures:::.getMetaDataValue(txdb, "Organism"),
    "SPECIES" = GenomicFeatures:::.getMetaDataValue(txdb, "Organism"),
    "PROVIDER" = provider,
    "PROVIDERVERSION" = providerVersion,
    "RELEASEDATE" = GenomicFeatures:::.getMetaDataValue(txdb, "Creation time"),
    ## SOURCEURL = GenomicFeatures:::.getMetaDataValue(txdb, "Resource URL"),
    "SOURCEURL" = db_url,
    "ORGANISMBIOCVIEW" = gsub(" ", "_", GenomicFeatures:::.getMetaDataValue(txdb, "Organism")),
    "TXDBOBJNAME" = package_name)
  if (any(duplicated(names(symvals)))) {
    str(symvals)
    stop("'symvals' contains duplicated symbols")
  }
  is_OK <- sapply(symvals, S4Vectors::isSingleString)
  if (!all(is_OK)) {
    bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
    stop("values for symbols ", bad_syms, " are not single strings")
  }
  pkg_list <- Biobase::createPackage(
                         pkgname=package_name,
                         destinationDir=destination,
                         originDir=template_path,
                         symbolValues=symvals,
                         unlink=TRUE)
  db_path <- file.path(destination, package_name, "inst", "extdata",
                       paste(package_name, "sqlite", sep="."))
  if (!file.exists(dirname(db_path))) {
    dir.create(dirname(db_path), recursive=TRUE)
  }
  obj <- try(AnnotationDbi::saveDb(txdb, file=db_path))
  if (class(obj) == "try-error") {
    warning("Failed to save the txdb object.")
  }

  install_dir <- paste0(destination, "/", package_name)
  install_dir <- clean_pkg(install_dir)
  ## Things with '-' are problematic.
  if (cfg_line[["strain"]] == "CLBrenerNon-Esmeraldo-like") {
    install_dir <- sm(clean_pkg(install_dir, removal="Non-Esmeraldo", replace="NonEsmeraldo"))
  }
  result <- sm(try(devtools::install(install_dir)))
  package_name <- basename(install_dir)
  ret <- list(
    "package_name" = package_name,
    "result" = result)
  return(ret)
}

#' Grab some configuration data collated and used to make OrganismDbi/OrgDb/TxDb objects.
#'
#' This function uses some data copied into inst/ to decide some parameters used for generating the
#' various packages generated here.
#'
#' Tested in test_46ann_tritrypdb.R
#' This function is sort of stupid and perhaps will be removed.  I keep a small csv file of some
#' TriTrypDB specific metadata, things like data base version number, URL schemes, etc.
#' This reads that and extracts the relevant information.
#'
#' @param cfg  Optional data frame
#' @return Dataframe of configuration data, a few columns are required, run it with no args to see
#'     which ones.
#' @export
get_eupath_config <- function(cfg=NULL) {
  cfg_data <- data.frame()
  if (is.null(cfg)) {
    cfg <- system.file("eupathdb_organdb.csv", package="hpgltools")
    cfg_data <- read.csv(cfg, stringsAsFactors=FALSE)
  } else if (class(cfg) == "data.frame") {
    cfg_data <- cfg
  } else if (class(cfg) == "character") {
    cfg_data <- read.csv(cfg, stringsAsFactors=FALSE)
  }
  rownames(cfg_data) <- cfg_data[["id"]]
  return(cfg_data)
}

#' Generate the (large) set of data frames required to make functional OrgDb/TxDb/OrganismDbi
#' objects.
#'
#' This function should probably be split into a few more pieces as it is pretty unwieldy at the
#' moment.
#' @param gff  File to read gff annotations from.
#' @param txt  File to read txt annotations from.
#' @param kegg  Boolean deciding whether to try for KEGG data.
#' @return  List containing gene information (likely from the txt file), chromosome information
#'   (gff file), gene types (gff file), gene ontology information, and potentially kegg information.
#' @seealso \pkg{rtracklayer} \pkg{GenomicRanges}
#' @examples
#' \dontrun{
#'  orgdb_data <- make_orgdb_info(gff="lmajor.gff", txt="lmajor.txt")
#' }
make_orgdb_info <- function(gff, txt=NULL, kegg=TRUE) {
  savefile <- paste0(txt, ".rda")
  gff_entries <- GenomicRanges::as.data.frame(rtracklayer::import.gff3(gff))
  gene_types <- gff_entries[["type"]] == "gene"
  genes <- gff_entries[gene_types, ]
  gene_info <- genes
  ##gene_info <- as.data.frame(mcols(genes))
  gene_info[["description"]] <- gsub("\\+", " ", gene_info[["description"]])  ## Get rid of stupid characters
  colnames(gene_info) <- toupper(colnames(gene_info))
  colnames(gene_info)[colnames(gene_info) == "ID"] <- "GID"

  chromosome_types <- c("apicoplast_chromosome", "chromosome", "contig",
                        "geneontig", "genecontig", "random_sequence", "supercontig")
  available_types <- gff_entries[["type"]]
  chromosome_entries <- available_types %in% chromosome_types
  chromosomes <- gff_entries[chromosome_entries, ]
  chromosome_info <- data.frame(
    "chrom" = chromosomes[["ID"]],
    "length" = as.numeric(chromosomes[["size"]]),
    "is_circular" = NA)

  gid_index <- grep("GID", colnames(gene_info))
  ## Move gid to the front of the line.
  gene_info <- gene_info[, c(gid_index, (1:ncol(gene_info))[-gid_index])]
  colnames(gene_info) <- paste0("GENE", colnames(gene_info))
  colnames(gene_info)[1] <- "GID"
  num_rows <- nrow(gene_info)
  ## Get rid of character(0) crap and NAs
  is.empty <- function(stuff) {
    (length(stuff) == 0) && (typeof(stuff) == "character")
  }
  empty_colp <- function(x) {
    y <- if (length(x)) {
           do.call("cbind", lapply(x, "is.empty"))
         } else {
           matrix(FALSE, length(row.names(x)), 0)
         }
    if (.row_names_info(x) > 0L) {
      rownames(y) <- row.names(x)
    }
    y
  }
  for (col in colnames(gene_info)) {
    tmp_col <- gene_info[[col]]
    empty_index <- empty_colp(tmp_col)
    tmp_col[empty_index] <- NA
    gene_info[[col]] <- tmp_col
    if (sum(!is.na(gene_info[[col]])) == num_rows) {
      gene_info[, !(colnames(gene_info) %in% col)]
    }
    gene_info[[col]] <- as.character(gene_info[[col]])
  }
  gene_info[is.na(gene_info)] <- "null"
  gene_info[gene_info == "NA"] <- "null"
  gene_info[["GENEALIAS"]] <- gsub(pattern="^c\\(", replacement="", x=gene_info[["GENEALIAS"]])
  gene_info[["GENEALIAS"]] <- gsub(pattern="\\)$", replacement="", x=gene_info[["GENEALIAS"]])
  gene_info[["GENEALIAS"]] <- gsub(pattern='"', replacement="", x=gene_info[["GENEALIAS"]])

  ## This function takes a long time.
  txt_information <- NULL
  if (file.exists(savefile)) {
    message("Reading the txt file takes a long time in R.  Happily we have a savefile for that.")
    message(paste0("Delete the file ", savefile, " to regenerate."))
    txt_information <- new.env()
    load(savefile, envir=txt_information)
    txt_information <- txt_information[["txt_information"]]
  } else {
    message("Reading the txt takes a long time.")
    txt_information <- parse_gene_info_table(file=txt, verbose=TRUE)
    save(txt_information, file=savefile)
  }
  gene_set <- txt_information[["genes"]]
  go_info <- txt_information[["go"]]

  chr_info <- data.frame(
    "GID" = rownames(gene_set),
    "CHR" = gene_set[["chromosome"]])

  gene_types <- data.frame(
    "GID" = rownames(gene_set),
    "TYPE" = gene_set[["type"]])

  ## The information in the following list really should be coming from parse_gene_info_table()
  ret <- list(
    "gene_info" = gene_info,
    "chr_info" = chr_info,
    "gene_types" = gene_types,
    "go_info" = go_info,
    "chromosome_info" = chromosome_info
  )
  return(ret)
}

#' EuPathDB gene information table GO term parser
#'
#' Note: EuPathDB currently includes some GO annotations corresponding to
#' obsolete terms. For example, the L. major gene LmjF.19.1390
#' (http://tritrypdb.org/tritrypdb/showRecord.do?name=GeneRecordClasses.GeneRecordClass&source_id=LmjF.19.1390&project_id=TriTrypDB)
#' includes the term "GO:0003702" on the website and txt data file. The
#' term has been deprecated, and does not have a category associated with it
#' on the website. These will not be included in the final database.
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @return Returns a dataframe where each line includes a gene/GO terms pair
#'         along with some addition information about the GO term. Note that
#'         because each gene may have multiple GO terms, a single gene ID may
#'         appear on multiple lines.
#' @export
parse_go_terms <- function(filepath) {
  fp <- NULL
  if (tools::file_ext(filepath) == "gz") {
    fp <- gzfile(filepath, open="rb")
  } else {
    fp <- file(filepath, open="r")
  }

  ## Create empty vector to store dataframe rows
  N <- 1e5
  gene_ids <- c()
  go_rows <- data.frame(GO=rep("", N),
                        ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
                        SOURCE=rep("", N), EVIDENCE=rep("", N),
                        stringsAsFactors=FALSE)

  ## Counter to keep track of row number
  i <- j <- 1

  ## Iterate through lines in file
  while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
    ## Gene ID
    if (grepl("^Gene ID", x)) {
      gene_id <- local_get_value(x)
      i <- i + 1
    }

    ## Gene Ontology terms
    else if (grepl("^GO:", x)) {
      gene_ids[j] <- gene_id
      go_rows[j, ] <- c(head(unlist(strsplit(x, "\t")), 5))
      j <- j + 1
    }
  }

  ## get rid of unallocated rows
  go_rows <- go_rows[1:j-1, ]
  ## drop unneeded columns
  go_rows <- go_rows[, c("GO", "EVIDENCE")]
  ## add gene id column
  go_rows <- cbind(GID=gene_ids, go_rows)
  ## close file pointer
  close(fp)
  ## TODO: Determine source of non-unique rows in the dataframe
  ## (May have to do with multiple types of evidence?)
  return(unique(go_rows))
}

#' EuPathDB gene information table InterPro domain parser
#'
#' @author Keith Hughitt
#' @param filepath Location of TriTrypDB gene information table.
#' @return Returns a dataframe where each line includes a gene/domain pairs.
#' @export
parse_interpro_domains <- function(filepath) {
  fp <- NULL
  if (tools::file_ext(filepath) == "gz") {
    fp <- gzfile(filepath, open="rb")
  } else {
    fp <- file(filepath, open="r")
  }

  ## Create empty vector to store dataframe rows
  ## N <- 1e5
  gene_ids <- c()
  ##interpro_rows = data.frame(GO=rep("", N),
  ##                     ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
  ##                     SOURCE=rep("", N), EVIDENCE=rep("", N),
  ##                     stringsAsFactors=FALSE)

  ## InterPro table columns
  cols <- c("name", "interpro_id", "primary_id", "secondary_id", "description",
            "start_min", "end_min", "evalue")

  ## Counter to keep track of row number
  i <- 1
  ## Iterate through lines in file
  while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
    ## Gene ID
    if (grepl("^Gene ID", x)) {
      gene_ids[i] <- local_get_value(x)
      i <- i + 1
    }

    ## Parse InterPro table
    else if (grepl("TABLE: InterPro Domains", x)) {
      ## Skip column header row
      trash <- readLines(fp, n=1)
      ## Continue reading until end of table
      raw_table <- ""
      entry <- readLines(fp, n=1)
      while (length(entry) != 0) {
        if (raw_table == "") {
          raw_table <- entry
        } else {
          raw_table <- paste(raw_table, entry, sep="\n")
        }
        entry <- readLines(fp, n=1)
      }

      ## If table length is greater than 0, read ino
      buffer <- textConnection(raw_table)
      interpro_table <- read.delim(buffer, header=FALSE, col.names=cols)
    }
  }

  ## add gene id column
  go_rows <- cbind(GID=gene_ids, go_rows)
  ## close file pointer
  close(fp)
  ## TODO: Determine source of non-unique rows in the dataframe
  return(unique(go_rows))
}

## EOF
