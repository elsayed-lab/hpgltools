#' @export
make_organismdbi <- function(id="lmajor_friedlin", cfg=NULL, output_dir="organdb", ...) {
    arglist <- list(...)
    kegg <- arglist[["kegg"]]
    cfg <- get_eupath_config(cfg)[id, ]
    version <- as.character(cfg[["db_version"]])
    strain <- as.character(cfg[["strain"]])
    shortname <- as.character(cfg[["shortname"]])
    files <- tritryp_downloads(version=version, species=shortname, strain=strain)
    ##files <- tritryp_downloads(version=version, species=shortname, strain=strain, ...)
    savefile <- paste0(output_dir, "/", cfg[["id"]], ".rda")
    orgdb_info <- NULL
    if (file.exists(savefile)) {
        message("Found a previous savefile for this species, loading from it.")
        orgdb_info <- new.env()
        load(savefile, envir=orgdb_info)
        orgdb_info <- orgdb_info[["orgdb_info"]]
        message(paste0("Loaded.  Delete <", savefile, "> to regenerate the rdata file."))
    } else {
        message("Reading from a freshly downloaded txt/gff file.")
        orgdb_info <- make_orgdb_info(files[["gff"]], files[["txt"]])
        save(orgdb_info, file=savefile)
    }

    orgdb_result <- make_orgdb(orgdb_info, id=id, cfg=cfg, output_dir=output_dir, kegg=kegg)
    orgdb_package <- orgdb_result[["package_name"]]
    txdb_result <- make_txdb(orgdb_info, id=id, cfg=cfg, gff=files[["gff"]], output_dir=output_dir)
    txdb_package <- txdb_result[["package_name"]]

    graph_data <- list(
        join1=c(GO.db='GOID', orgdb='GO'),
        join2=c(orgdb='GID',  txdb='GENEID')
    )

    names(graph_data[["join1"]]) = c('GO.db', orgdb_package)
    names(graph_data[["join2"]]) = c(orgdb_package, txdb_package)

    requireNamespace(orgdb_package)
    requireNamespace(txdb_package)
    libstring <- paste0("library(", orgdb_package, ")")
    eval(parse(text=libstring))
    libstring <- paste0("library(", txdb_package, ")")
    eval(parse(text=libstring))
    organism <- paste0(cfg[["genus"]], " ", cfg[["species"]], " ", cfg[["strain"]])
    organism <- gsub(pattern="-like", replacement="", x=organism)
    requireNamespace("OrganismDbi")
    pkgname <- as.character(cfg[["organismdb_name"]])
    author <- as.character(cfg[["author"]])
    maintainer <- as.character(cfg[["maintainer"]])

    destination <- paste0(output_dir, "/organismdbi")
    if (file.exists(destination)) {
        unlink(x=destination, recursive=TRUE)
    }
    dir.create(destination, recursive=TRUE)
    version = format(as.numeric(version), nsmall=1)
    organdb <- OrganismDbi::makeOrganismPackage(
        pkgname = pkgname,
        graphData = graph_data,
        organism = organism,
        version = version,
        maintainer = maintainer,
        author = author,
        destDir = destination,
        license='Artistic-2.0'
    )
    organdb_path <- paste0(destination, "/", pkgname)
    organdb_path <- pkg_cleaner(organdb_path)
    if (class(organdb) == "list") {
        inst <- devtools::install(organdb_path)
    }
    return(inst)
}

#' @export
pkg_cleaner <- function(path) {
    ## This is because TxDb creation fails if you have an author like 'abelew <abelew@gmail.com>'
    at_cmd <- paste0("sed -i 's/ at /\\@/g' ", path, "/DESCRIPTION")
    message(paste0("Rewriting DESCRIPTION: ", at_cmd))
    system(command=at_cmd)
    ## Since I changed @ to at I figured . could be dot too
    dot_cmd <- paste0("sed -i 's/ dot /\\./g' ", path, "/DESCRIPTION")
    message(paste0("Rewriting DESCRIPTION to remove dot: ", dot_cmd))
    system(dot_cmd)
    new_dir <- path
    if (grepl(pattern="-like", x=path)) {
        ## Get rid of the -like in the path name
        new_dir <- gsub(pattern="-like", replacement="", x=path)
        ## And rename the directory
        mv_cmd <- paste0("mv ", path, " ", new_dir)
        message(paste0("moving orgdb: ", mv_cmd))
        system(mv_cmd)
        ## Collect the text files in the new package and remove all -like instances in them
        find_cmd <- paste0("sed -i 's/-like//g' $(find ", new_dir, " -type f | grep -v sqlite)")
        message(paste0("rewriting orgdb files: ", find_cmd))
        system(find_cmd)

        ## Move the sqlite file
        old_sqlite_basename <- basename(path)
        old_sqlite_base <- gsub(pattern=".db", replacement="", x=old_sqlite_basename)
        sqlite_basename <- basename(new_dir)
        sqlite_base <- gsub(pattern=".db", replacement="", x=sqlite_basename)
        sqlite_name <- paste0(sqlite_base, ".sqlite")
        old_sqlite <- paste0(new_dir, "/inst/extdata/", old_sqlite_base, ".sqlite")
        new_sqlite <- gsub(pattern="-like", replacement="", x=old_sqlite)
        sqlite_mv_cmd <- paste0("mv ", old_sqlite, " ", new_sqlite)
        message(paste0("moving sqlite file: ", sqlite_mv_cmd))
        system(sqlite_mv_cmd)
        orgdb_dir <- new_dir
        new_pkg_name <- gsub(pattern="-like", replacement="", x=sqlite_basename)
        ## Update the orgdb sqlite file to reflect the new name
        final_sqlite_cmd <- paste0("chmod +w ", new_sqlite, " ; sqlite3 ", new_sqlite, " \"UPDATE metadata SET value='", new_pkg_name, "' WHERE name='SPECIES';\" ; chmod -w ", new_sqlite)
        message(paste0("rewriting sqlite db:", final_sqlite_cmd))
        system(final_sqlite_cmd)
    }
    message("Finished pkg_cleaner.")
    return(new_dir)
}

#' @export
make_orgdb <- function(orgdb_info, id="lmajor_friedlin", cfg=NULL, kegg=NULL, output_dir="organismdbi", ...) {
    arglist=list(...)
    orgdb_pre <- paste0(output_dir, "/orgdb")
    if (!file.exists(orgdb_pre)) {
        dir.create(orgdb_pre, recursive=TRUE)
    }

    gene_info <- orgdb_info[["gene_info"]]
    gene_info[["GID"]] <- as.character(gene_info[["GID"]])
    rownames(gene_info) <- gene_info[["GID"]]
    chr_info <- orgdb_info[["chr_info"]]
    chr_info[["GID"]] <- as.character(chr_info[["GID"]])
    go_info <- orgdb_info[["go_info"]]
    go_info[["GID"]] <- as.character(go_info[["GID"]])
    go_info <- go_info[, c("GID", "GO", "EVIDENCE")]
    gene_types <- orgdb_info[["gene_types"]]
    gene_types[["GID"]] <- as.character(gene_types[["GID"]])

    orgdb_base_name <- paste0("org.", cfg[["shortname"]], ".", cfg[["strain"]], ".eg")
    orgdb_pkg_name <- paste0(orgdb_base_name, ".db")
    orgdb_sqlite_name <- paste0(orgdb_base_name, ".sqlite")
    assumed_dir <- paste0(output_dir, "/orgdb/", orgdb_pkg_name)

    if (file.exists(assumed_dir)) {
        unlink(x=assumed_dir, recursive=TRUE)
    }
    orgdb_dir <- NULL
    if (is.null(kegg)) {
        orgdb_dir <- AnnotationForge::makeOrgPackage(
            gene_info = gene_info,
            chromosome = chr_info,
            go = go_info,
            type = gene_types,
            version = format(as.numeric(cfg[["db_version"]]), nsmall=1),
            author = as.character(cfg[["author"]]),
            maintainer = as.character(cfg[["maintainer"]]),
            tax_id = as.character(cfg[["tax_id"]]),
            genus = as.character(cfg[["genus"]]),
            species = paste0(as.character(cfg[["species"]]), ".", as.character(cfg[["strain"]])),
            outputDir = orgdb_pre,
            goTable = "go")
    } else {
        kegg_info <- orgdb_kegg(cfg, gene_info)
        orgdb_dir <- makeOrgPackage(
            gene_info = gene_info,
            chromosome = chr_info,
            go = go_info,
            kegg_info = kegg_table,
            type = gene_types,
            version = as.character(cfg[["db_version"]]),
            author = as.character(cfg[["author"]]),
            maintainer = as.character(cfg[["maintainer"]]),
            tax_id = as.character(cfg[["tax_id"]]),
            genus = as.character(cfg[["genus"]]),
            species = as.character(cfg[["species"]]),
            outputDir = orgdb_pre,
            goTable = "go")
    }

    orgdb_dir <- pkg_cleaner(orgdb_dir)

    inst <- FALSE
    if (!is.null(orgdb_dir)) {
        inst <- sm(devtools::install(orgdb_dir))
    }
    orgdb_pkg_name <- gsub(pattern="-like", replacement="", x=orgdb_pkg_name)
    ## The result is the pathname of the created orgdb directory
    ret <- list(
        package_name = orgdb_pkg_name,
        result = inst)
    return(ret)
}

make_txdb <- function(orgdb_info, cfg, gff=NULL, output_dir="organismdbi", ...) {
    arglist <- list(...)

    chromosome_info <- orgdb_info[["chromosome_info"]]
    txdb <- GenomicFeatures::makeTxDbFromGFF(
        file=gff,
        format='gff3',
        chrominfo=chromosome_info,
        ## exonRankAttributeName=NA,
        dataSource=paste0(cfg[["db_name"]], "_", cfg[["db_version"]]),
        organism=paste0(cfg[["genus"]], " ", cfg[["species"]])
    )

    requireNamespace("GenomicFeatures")
    destination <- paste0(output_dir, "/txdb")
    db_version <- format(as.numeric(cfg[["db_version"]]), nsmall=1)
    maintainer <- as.character(cfg[["maintainer"]])
    author <- as.character(cfg[["author"]])
    package_name <- paste0("TxDb.", cfg[["shortname"]], ".", cfg[["strain"]], ".", cfg[["db_name"]], cfg[["db_version"]])

    if (file.exists(destination)) {
        unlink(x=destination, recursive=TRUE)
    }
    dir.create(destination, recursive=TRUE)
    result <- GenomicFeatures::makeTxDbPackage(
        txdb=txdb,
        version=db_version,
        maintainer=maintainer,
        author=author,
        destDir=destination,
        license="Artistic-2.0",
        pkgname=package_name
    )
    ## What in the flying hell means 'Error in cpSubsCon(src[k], destname) : UNRESOLVED SYMBOLS:
    ## Line 5 : @umd.edu>, atb <atb@u'
    ## The entire person/author/maintainer system in R is utterly stupid.
    install_dir <- paste0(destination, "/", package_name)
    install_dir <- pkg_cleaner(install_dir)
    result <- devtools::install(install_dir)
    package_name <- gsub(pattern="-like", replacement="", x=package_name)
    ret <- list(
        package_name = package_name,
        result = result)
    return(ret)
}

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

#' @export
make_orgdb_info <- function(gff, txt) {
    gff_entries <- rtracklayer::import.gff3(gff)

    genes <- gff_entries[gff_entries$type == "gene"]  ## WTF? why does this work?
    gene_info <- as.data.frame(GenomicRanges::mcols(genes))
    gene_info[["description"]] <- gsub("\\+", " ", gene_info[["description"]])  ## Get rid of stupid characters
    colnames(gene_info) <- toupper(colnames(gene_info))
    colnames(gene_info)[colnames(gene_info) == "ID"] <- "GID"

    chromosome_types <- c("apicoplast_chromosome", "chromosome", "contig",
                          "geneontig", "genecontig", "random_sequence", "supercontig")
    chromosomes <- gff_entries[gff_entries$type %in% chromosome_types]
    chromosome_info <- data.frame(
        "chrom" = chromosomes$ID,
        "length" = as.numeric(chromosomes$size),
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
    is.empty.col <- function(x) {
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
        empty_index <- is.empty.col(tmp_col)
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
    txt_information <- parse_gene_info_table(file=txt, verbose=TRUE)
    genes <- txt_information[["genes"]]
    go_info <- txt_information[["go"]]

    chr_info <- data.frame(
        "GID" = rownames(genes),
        "CHR" = genes[["chromosome"]])

    gene_types <- data.frame(
        "GID" = rownames(genes),
        "TYPE" = genes[["type"]])

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

orgdb_kegg <- function(cfg, gene_info) {
    org_abbreviation <- paste0(tolower(substring(cfg[["genus"]], 1, 1)),
                              substring(cfg[["species"]], 1, 2))

    # Overides for cases where KEGG abbreviation differes from the above
    # pattern.

    # L. braziliensis
    if (org_abbreviation == 'lbr') {
        org_abbreviation <- 'lbz'
    }

    # For some species, it is necessary to map gene ids from KEGG to what is
    # currently used on TriTrypDB.
    #
    # TODO: Generalize if possible
    #
    ##if (org_abbreviation == 'tbr') {
        # Load GeneAlias file and convert entry in KEGG results
        # to newer GeneDB/TriTrypDB identifiers.
        ## kegg_id_mapping <- list()
        # example alias file entries
        #Tb927.10.2410  TRYP_x-70a06.p2kb545_720  Tb10.70.5290
        #Tb927.9.15520  Tb09.244.2520  Tb09.244.2520:mRNA
        #Tb927.8.5760   Tb08.26E13.490
        #Tb10.v4.0258   Tb10.1120
        #Tb927.11.7240  Tb11.02.5150  Tb11.02.5150:mRNA  Tb11.02.5150:pep
        ##for (alias in strsplit(alias, split=", ")) {
        ##    # get first and third columns in the alias file
        ##    old_ids <- row[2:length(row)]
        ##
        ##    for (old_id in old_ids[grepl('Tb\\d+\\.\\w+\\.\\d+', old_ids)]) {
        ##        kegg_id_mapping[old_id] <- row[1]
        ##    }
        ##}
    ##
    ##} else if (org_abbreviation == 'lma') {
        # L. major
        #
        # Convert KEGG identifiers to TriTrypDB identifiers
        #
        # Note that this currently skips a few entries with a different
        # format, e.g. "md:lma_M00359", and "bsid:85066"
        #
    ##} else if (org_abbreviation == 'tcr') {
        ##kegg_id_mapping <- list()
        ##
        ##for (row in rows) {
        ##    # get first and third columns in the alias file
        ##    kegg_id_mapping[row[3]] <- row[1]
       ## }

        # Example: "tcr:509463.30" -> ""
        ##convert_kegg_gene_ids <- function(kegg_ids) {
        ##    kegg_to_genedb(kegg_ids, kegg_id_mapping)
        ##}
    ##}

    # data frame to store kegg gene mapping and pathway information
    kegg_mapping <- data.frame()
    kegg_pathways <- data.frame()

    pathways <- unique(keggLink("pathway", org_abbreviation))

    # Iterate over pathways and query genes for each one
    for (pathway in pathways) {
        message(sprintf("Processing genes for KEGG pathway %s", pathway))

        # Get pathway info
        meta <- keggGet(pathway)[[1]]
        pathway_desc  <- ifelse(is.null(meta$DESCRIPTION), '', meta$DESCRIPTION)
        pathway_class <- ifelse(is.null(meta$CLASS), '', meta$CLASS)
        kegg_pathways <- rbind(kegg_pathways,
                               data.frame(
                                   "pathway" = pathway,
                                   "name" = meta$PATHWAY_MAP,
                                   "class" = pathway_class,
                                   "description" = pathway_desc))

        # Get genes in pathway
        kegg_ids <- as.character(keggLink(org_abbreviation, pathway))
        gene_ids <- convert_kegg_gene_ids(kegg_ids)

        # Map old T. brucei gene names
        if (org_abbreviation == 'tbr') {
            old_gene_ids <- gene_ids
            gene_ids <- c()

            for (x in old_gene_ids) {
                if (x %in% names(kegg_id_mapping)) {
                    gene_ids <- append(gene_ids, kegg_id_mapping[[x]])
                } else {
                    gene_ids <- append(gene_ids, x)
                }
            }
        }

        if (!is.null(gene_ids)) {
            kegg_mapping <- unique(rbind(kegg_mapping,
                data.frame(GID=gene_ids, pathway=pathway)))
        }
    }

    # Drop columns with unrecognized identifiers
kegg_mapping <- kegg_mapping[complete.cases(kegg_mapping),]

# Combined KEGG table
kegg_table <- merge(kegg_mapping, kegg_pathways, by='pathway')
colnames(kegg_table) <- c("KEGG_PATH", "GID", "KEGG_NAME", "KEGG_CLASS",
                         "KEGG_DESCRIPTION")

# reorder so GID comes first
kegg_table <- kegg_table[,c(2, 1, 3, 4, 5)]
}

convert_kegg_gene_ids <- function(kegg_ids, kegg_id_mapping) {
    result <- c()
    for (kegg_id in kegg_ids) {
        if (substring(kegg_id, 1, 4) == 'tbr:') {
            # T. brucei
            result <- append(result,
                gsub('tbr:', '', kegg_id))
        } else if (substring(kegg_id, 1, 4) == 'tcr:') {
            # T. cruzi
            result <- append(result,
                gsub('tcr:', 'TcCLB.', kegg_id))
        } else if (substring(kegg_id, 1, 4) == 'tgo:') {
            # T. gondii
            result <- append(result, gsub('tgo:', '', gsub('_', '.', kegg_id)))
        } else if (substring(kegg_id, 1, 9) == 'lbz:LBRM_') {
            # L. braziliensis (lbz:LBRM_01_0080)
            result <- append(result, gsub('LBRM', 'LbrM',
                     gsub("_", "\\.", substring(kegg_id, 5))))
        } else if (substring(kegg_id, 1, 9) == 'lma:LMJF_') {
            # L. major (lma:LMJF_11_0100)
            result <- append(result,
                gsub('LMJF', 'LmjF',
                     gsub("_", "\\.", substring(kegg_id, 5))))
        } else if (substring(kegg_id, 1, 8) == 'lma:LMJF') {
            # L. major (lma:LMJF10_TRNALYS_01)
            parts <- unlist(strsplit(kegg_id, "_"))
            result <- append(result,
                sprintf("LmjF.%s.%s.%s",
                        substring(kegg_id, 9, 10),
                        parts[2], parts[3]))
        } else {
            print(sprintf("Skipping KEGG id: %s", kegg_id))
            result <- append(result, NA)
        }
    }
    return(result)
} ## End convert_kegg_gene_ids

#
# kegg_to_genedb
#
# Takes a list of KEGG gene identifiers and returns a list of GeneDB
# ids corresponding to those genes.
#
kegg_to_genedb = function(kegg_ids, gene_mapping) {
    # query gene ids 10 at a time (max allowed)
    result = c()

    for (x in split(kegg_ids, ceiling(seq_along(kegg_ids) / 10))) {
        query = keggGet(x)
        for (item in query) {
            dblinks = item$DBLINKS
            genedb_id = dblinks[grepl('GeneDB', dblinks)]
            if (length(genedb_id) > 0) {
                # get old-style t. cruzi identifier
                old_id = substring(genedb_id, 9)

                # if possible, map to new id and add to results
                if (!is.null(gene_mapping[[old_id]])) {
                    result = append(result, gene_mapping[[old_id]])
                }
            }
        }
    }
    return(result)
}

#'
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
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe where each line includes a gene/GO terms pair
#'         along with some addition information about the GO term. Note that
#'         because each gene may have multiple GO terms, a single gene ID may
#'         appear on multiple lines.
#'
parse_go_terms = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N = 1e5
    gene_ids = c()
    go_rows = data.frame(GO=rep("", N),
                         ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
                         SOURCE=rep("", N), EVIDENCE=rep("", N),
                         stringsAsFactors=FALSE)

    # Counter to keep track of row number
    i = j = 1

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
            i = i + 1
        }

        # Gene Ontology terms
        else if (grepl("^GO:", x)) {
            gene_ids[j] = gene_id
            go_rows[j,] = c(head(unlist(strsplit(x, '\t')), 5))
            j = j + 1
        }
    }

    # get rid of unallocated rows
    go_rows = go_rows[1:j-1,]

    # drop unneeded columns
    go_rows = go_rows[,c('GO', 'EVIDENCE')]

    # add gene id column
    go_rows = cbind(GID=gene_ids, go_rows)

    # close file pointer
    close(fp)

    # TODO: Determine source of non-unique rows in the dataframe
    # (May have to do with multiple types of evidence?)
    return(unique(go_rows))
}

#'
#' EuPathDB gene information table InterPro domain parser
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#'
#' @return Returns a dataframe where each line includes a gene/domain pairs.
#'
parse_interpro_domains = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    #N = 1e5
    #gene_ids = c()
    #interpro_rows = data.frame(GO=rep("", N),
    #                     ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
    #                     SOURCE=rep("", N), EVIDENCE=rep("", N),
    #                     stringsAsFactors=FALSE)

    # InterPro table columns
    cols = c('name', 'interpro_id', 'primary_id', 'secondary_id', 'description',
             'start_min', 'end_min', 'evalue')

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
        }

        # Parse InterPro table
        else if (grepl("TABLE: InterPro Domains", x)) {
            # Skip column header row
            trash = readLines(fp, n=1)

            # Continue reading until end of table
            raw_table = ""

            entry = readLines(fp, n=1)

            while(length(entry) != 0) {
                if (raw_table == "") {
                    raw_table = entry
                } else {
                    raw_table = paste(raw_table, entry, sep='\n')
                }
                entry = readLines(fp, n=1)
            }

            # If table length is greater than 0, read ino
            buffer = textConnection(raw_table)

            interpro_table = read.delim(buffer, header=FALSE, col.names=cols)

        }
    }

    # add gene id column
    go_rows = cbind(GID=gene_ids, go_rows)

    # close file pointer
    close(fp)

    # TODO: Determine source of non-unique rows in the dataframe
    return(unique(go_rows))
}
