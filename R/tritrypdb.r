#' TriTrypDB gene information table parser
#'
#' An example input file is the T. brucei Lister427 gene information table
#' available at:
#' http://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427/txt/data/TriTrypDB-5.0_TbruceiLister427Gene.txt
#'
#' @author Keith Hughitt
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe of gene info.
parse_gene_info_table <- function(file, verbose=FALSE) {
    ## Create empty vector to store dataframe rows

    .get_value <- function(x) {
        return(gsub(" ","", tail(unlist(strsplit(x, ": ")), n=1), fixed=TRUE))
    }

    N <- 1e5
    genes <- data.frame("chromosome" = rep(NA, N),
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
    go_rows <- data.frame("GO" = rep("", N),
                          "ONTOLOGY" = rep("", N),
                          "TERM" = rep("", N),
                          "SOURCE" = rep("", N),
                          "EVIDENCE" = rep("", N),
                          stringsAsFactors = FALSE)
    ## Regular expression to extract location info
    location_regex <- '([0-9,]*) - ([0-9,]*) \\(([-+])\\)'
    gene_ids <- c()
    go_gene_ids <- c()
    gene_num <- 0
    go_num <- 0
    message("Using readr to read the txt file.")
    read_vec <- readr::read_lines(file, progress=TRUE)
    message("Starting to iterate over the txt file, this takes a long time.")
    for (i in 1:length(read_vec)) {
        line <- read_vec[i]

        ## Gene ID
        if (grepl("^Gene ID", line)) {
            gene_num <- gene_num + 1
            gene_id <- .get_value(line)
            gene_ids[gene_num] <- gene_id
            if (verbose) {
                if ((gene_num %% 1000) == 0) {
                    message(sprintf('Processing gene %d: %s', gene_num, gene_id))
                }
            }
        }

        ##Molecular Weight: 37091
        else if(grepl("^Molecular Weight", line)) {
            if (grepl("^Molecular Weight: Not Assigned", line)) {
                mweight <- NA
            } else {
                mweight <- as.numeric(.get_value(line))
            }
        }

        ##Isoelectric Point: 7.79
        else if(grepl("^Isoelectric Point", line)) {
            if (grepl("^Isoelectric Point: Not Assigned", line)) {
                ipoint <- NA
            } else {
                ipoint <- as.numeric(.get_value(line))
            }
        }

        ## Chromosome number
        else if(grepl("^Chromosome", line)) {
            if (grepl("^Chromosome: Not Assigned", line)) {
                chromosome <- NA
            } else {
                chromosome <- as.numeric(.get_value(line))
            }
        }

        ## Genomic location
        else if (grepl("^Genomic Location:", line)) {
            result <- unlist(regmatches(line, regexec(location_regex, line)))
            gene_start <- as.numeric(gsub(",", "", result[2], fixed=TRUE))
            gene_stop  <- as.numeric(gsub(",", "", result[3], fixed=TRUE))
            strand <- result[4]
        }

        ## Gene type
        else if (grepl("^Gene Type", line)) {
            gene_type <- .get_value(line)
        }

        ## Product Description
        else if (grepl("^Product Description", line)) {
            description <- .get_value(line)
        }

        ## Transcript length
        else if (grepl("^Transcript Length", line)) {
            transcript_length <- as.numeric(.get_value(line))
        }

        ## CDS length
        else if (grepl("^CDS Length", line)) {
            val <- .get_value(line)
            if (val == 'null') {
                cds_length <- NA
            } else {
                cds_length <- as.numeric(val)
            }
        }

        ## Pseudogene
        else if (grepl("^Is Pseudo:", line)) {
            is_pseudo <- ifelse((.get_value(line) == "Yes"), TRUE, FALSE)
        }

        ## Gene ontology terms
        else if (grepl("^GO:", line)) {
            go_num <- go_num + 1
            go_gene_ids[go_num] <- gene_id
            go_rows[go_num, ] <- c(head(unlist(strsplit(line, '\t')), 5))
        }

        ## End of gene description
        else if (grepl("^---", line)) {
            ## Skip gene if it is not assigned to a chromosome
            if (is.na(chromosome)) {
                next
            }

            ## Otherwise add row to dataframe
            genes[gene_num, ] <- c(chromosome, gene_start, gene_stop, strand,
                            gene_type, transcript_length, cds_length, is_pseudo,
                            description, mweight, ipoint)
        }


    } ## End parsing the file
    rm(read_vec)

    go_rows <- go_rows[1:go_num, ]
    go_rows <- cbind("GID" = go_gene_ids,
                     go_rows)

    ## get rid of unallocated rows
    genes <- genes[1:gene_num, ]
    ## use gene id as row name
    rownames(genes) <- gene_ids
    ## fix numeric types
    for (colname in c('chromosome', 'start', 'stop', 'transcript_length', 'cds_length')) {
        genes[, colname] <- as.numeric(genes[, colname])
    }
    # sort data frame
    genes <- genes[with(genes, order(chromosome, start)),]
    ret <- list(
        "genes" = genes,
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
    if (tools::file_ext(filepath) == 'gz') {
        fp <- gzfile(filepath, open='rb')
    } else {
        fp <- file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N <- 1e5
    gene_ids <- c()
    go_rows <- data.frame("GO"=rep("", N),
                          "ONTOLOGY"=rep("", N), TERM=rep("", N),
                          "SOURCE"=rep("", N), EVIDENCE=rep("", N),
                          stringsAsFactors=FALSE)
    # Counter to keep track of row number
    i <- j <- 1
    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id <- .get_value(x)
            if (verbose) {
                print(sprintf('Processing gene %d: %s', i, gene_id))
            }
            i <- i + 1
        }

        # Gene Ontology terms
        else if (grepl("^GO:", x)) {
            gene_ids[j] <- gene_id
            go_rows[j,] <- c(head(unlist(strsplit(x, '\t')), 5))
            j <- j + 1
        }
    }
    # get rid of unallocated rows
    go_rows <- go_rows[1:j-1,]
    # add gene id column
    go_rows <- cbind(GID=gene_ids, go_rows)
    # close file pointer
    close(fp)
    return(go_rows)
}

#
# Parses a key: value string and returns the value
#


#' @export
tritryp_downloads <- function(version="24", species="lmajor", strain="friedlin", dl_dir="organdb/tritryp", quiet=TRUE) {
    files_downloaded <- 0
    files_found <- 0

    if (!file.exists(dl_dir)) {
        dir.create(dl_dir, recursive=TRUE)
    }

    uc_species <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))
    uc_strain <- paste0(toupper(substr(strain, 1, 1)), substr(strain, 2, nchar(strain)))
    gff_filename <- paste0("TriTrypDB-", version, "_", uc_species, uc_strain, ".gff")

    gff_path <- paste0(dl_dir, "/", gff_filename)
    gff_url <- paste0("http://tritrypdb.org/common/downloads/release-", version, "/", uc_species, uc_strain, "/gff/data/", gff_filename)
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
    txt_url <- paste0("http://tritrypdb.org/common/downloads/release-", version, "/", uc_species, uc_strain, "/txt/", txt_filename)
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
    fasta_url <- paste0("http://tritrypdb.org/common/downloads/release-", version, "/", uc_species, uc_strain, "/fasta/data/", fasta_filename)
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
    alias_url <- paste0("http://tritrypdb.org/common/downloads/release-", version, "/", uc_species, uc_strain, "/txt/", alias_filename)
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
