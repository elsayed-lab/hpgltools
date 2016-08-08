
#' @export
orgdb_genes <- function(gff, txt) {
    gff <- "/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener_nonesmeraldo-like/annotation/TriTrypDB-27_TcruziCLBrenerNon-Esmeraldo-like.gff"
    gff_entries <- rtracklayer::import.gff3(gff)
    genes <- gff_entries[gff_entries$type == "gene"]  ## WTF? why does this work?
    gene_info <- as.data.frame(elementMetadata(genes))
    gene_info[["description"]] <- gsub("\\+", " ", gene_info[["description"]])  ## Get rid of stupid characters
    colnames(gene_info) <- toupper(colnames(gene_info))
    colnames(gene_info)[colnames(gene_info) == "ID"] <- "GID"
    gid_index <- grep("GID", colnames(gene_info))
    ## Move gid to the front of the line.
    gene_info <- gene_info[, c(gid_index, (1:ncol(gene_info))[-gid_index])]
    colnames(gene_info) <- paste0("GENE", colnames(gene_info))
    colnames(gene_info)[1] <- "GID"
    num_rows <- nrow(gene_info)
    gene_info[["GENEALIAS"]] <- as.character(gene_info[["GENEALIAS"]])
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

    chr_info <- data.frame(
        "GID" = genes$ID,
        "CHR" = as.character(seqnames(genes)))

    gene_types <- orgdb_gene_types(txt)

    ret <- list(
        "gene_info" = gene_info,
        "chr_info" = chr_info,
        "gene_types" = gene_types,
    )
    return(ret)
}

orgdb_go <- function(go_file) {
    go_table <- parse_go_terms(settings$txt)

    ## Map from non-primary IDs to primary GO ids;
    ## non-primary IDs are filtered out by makeOrgPackage
    problem_rows <- go_table[!go_table$GO %in% keys(GO.db),]
    synonyms <- problem_rows$GO

    ##Create a mapping data frame
    synonym_mapping <- data.frame()

    for (syn in synonyms) {
        if (!is.null(GOSYNONYM[[syn]])) {
            synonym_mapping <- rbind(synonym_mapping, c(syn, GOSYNONYM[[syn]]@GOID))
        }
    }

    ## replace alternative GO term identifiers
    if (nrow(synonym_mapping) > 0) {
        colnames(synonym_mapping) <- c('synonym', 'primary')
        synonym_mapping <- unique(synonym_mapping)

        go_table$GO[!go_table$GO %in% keys(GO.db)] <- synonym_mapping$primary[match(synonyms, synonym_mapping$synonym)]
        go_table <- unique(go_table[complete.cases(go_table),])
    }
        return(go_table)
}

orgdb_kegg <- function(kegg) {
    org_abbreviation <- paste0(tolower(substring(settings$genus, 1, 1)),
                              substring(settings$species, 1, 2))

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
    if (org_abbreviation == 'tbr') {
        # Load GeneAlias file and convert entry in KEGG results
        # to newer GeneDB/TriTrypDB identifiers.
        fp <- file(settings$aliases)
        rows <- strsplit(readLines(fp), "\t")
        close(fp)

        kegg_id_mapping <- list()

        # example alias file entries
        #Tb927.10.2410  TRYP_x-70a06.p2kb545_720  Tb10.70.5290
        #Tb927.9.15520  Tb09.244.2520  Tb09.244.2520:mRNA
        #Tb927.8.5760   Tb08.26E13.490
        #Tb10.v4.0258   Tb10.1120
        #Tb927.11.7240  Tb11.02.5150  Tb11.02.5150:mRNA  Tb11.02.5150:pep
        for (row in rows) {
            # get first and third columns in the alias file
            old_ids <- row[2:length(row)]

            for (old_id in old_ids[grepl('Tb\\d+\\.\\w+\\.\\d+', old_ids)]) {
                kegg_id_mapping[old_id] <- row[1]
            }
        }

    } else if (org_abbreviation == 'lma') {
        # L. major
        #
        # Convert KEGG identifiers to TriTrypDB identifiers
        #
        # Note that this currently skips a few entries with a different
        # format, e.g. "md:lma_M00359", and "bsid:85066"
        #
    } else if (org_abbreviation == 'tcr') {
        fp <- file(settings$aliases)
        rows <- strsplit(readLines(fp), "\t")
        close(fp)

        kegg_id_mapping <- list()

        for (row in rows) {
            # get first and third columns in the alias file
            kegg_id_mapping[row[3]] <- row[1]
        }

        # Example: "tcr:509463.30" -> ""
        ##convert_kegg_gene_ids <- function(kegg_ids) {
        ##    kegg_to_genedb(kegg_ids, kegg_id_mapping)
        ##}
    }

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
                               data.frame(pathway=pathway,
                                         name=meta$PATHWAY_MAP,
                                         class=pathway_class,
                                         description=pathway_desc))

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
    # Save KEGG mapping
    write.csv(kegg_mapping, file=kegg_mapping_file, quote=FALSE,
              row.names=FALSE)
    write.table(kegg_pathways, file=kegg_pathways_file, quote=FALSE,
                row.names=FALSE, sep='\t')

    # Drop columns with unrecognized identifiers
kegg_mapping <- kegg_mapping[complete.cases(kegg_mapping),]

# Combined KEGG table
kegg_table <- merge(kegg_mapping, kegg_pathways, by='pathway')
colnames(kegg_table) <- c("KEGG_PATH", "GID", "KEGG_NAME", "KEGG_CLASS",
                         "KEGG_DESCRIPTION")

# reorder so GID comes first
kegg_table <- kegg_table[,c(2, 1, 3, 4, 5)]
}

make_orgdb_pkg <- function(gene, chr, go, kegg) {
    if (is.null(kegg)) {
    org_result <- makeOrgPackage(
        gene_info  = gene_info,
        chromosome = chr_info,
        go         = go_table,
        type       = gene_types,
        version    = db_version,
        author     = settings$author,
        maintainer = settings$maintainer,
        outputDir  = settings$output_dir,
        tax_id     = settings$tax_id,
        genus      = settings$genus,
        species    = settings$species,
        goTable    = "go")
    } else {
        org_result <- makeOrgPackage(
        gene_info  = gene_info,
        chromosome = chr_info,
        go         = go_table,
        kegg = kegg_table,
        type       = gene_types,
        version    = db_version,
        author     = settings$author,
        maintainer = settings$maintainer,
        outputDir  = settings$output_dir,
        tax_id     = settings$tax_id,
        genus      = settings$genus,
        species    = settings$species,
        goTable    = "go")
    }
    return(org_result)
}

gather_txdb <- function() {
    gff = import.gff3(settings$gff)
    ch = gff[gff$type %in% c('apicoplast_chromosome', 'chromosome', 'contig',
                             'geneontig', 'random_sequence', 'supercontig')]

    ##genes = gff[gff$type == 'gene']
    ##gene_ch = unique(as.character(chrom(genes)))

    chrom_info = data.frame(
        chrom=ch$ID,
        length=as.numeric(ch$size),
        is_circular=NA
    )

    ## 2015/06/16 Switching backt o mRNA entries to construct TxDb -- database is
    ## intended for mRNAs so ncRNAs will be handled separately.

    ##txdb = makeTranscriptDbFromGFF(
    txdb = makeTxDbFromGFF(
        file=settings$gff,
        format='gff3',
        chrominfo=chrom_info,
        ## exonRankAttributeName=NA,
        dataSource=sprintf('%s %s', settings$db_name, settings$db_version),
        organism=paste(settings$genus, settings$species)
    )

    ## Save transcript database
    short_name = paste0(substring(tolower(settings$genus), 1, 1), settings$species)
    saveDb(txdb, file=file.path(build_dir, sprintf("%s.sqlite", short_name)))

    ## R package versions must be of the form "x.y"
    db_version = paste(settings$db_version, '0', sep='.')

    ## Build TxDB package
    makeTxDbPackage(
        txdb,
        destDir=settings$output_dir,
        version=db_version,
        maintainer=settings$maintainer,
        author=settings$author,
        license='Artistic-2.0'
    )
}

make_organdb <- function() {

    names(graph_data$join1) = c('GO.db', settings$orgdb_name)
    names(graph_data$join2) = c(settings$orgdb_name, settings$txdb_name)

    ## R package versions must be of the form "x.y"
    db_version = paste(settings$db_version, '0', sep='.')

    makeOrganismPackage(
        pkgname=settings$organismdb_name,
        graphData=graph_data,
        organism=paste(settings$genus, settings$species),
        version=db_version,
        maintainer=settings$maintainer,
        author=settings$author,
        destDir=settings$output_dir,
        license='Artistic-2.0'
    )
}

install_orgdb <- function() {

}
prepare_dbs <- function() {
    ## This part really does not make sense to me.
    ##printf -v orgdb_name_old 'org.%s%s.eg.db' ${genus:0:1} ${species}
    ##echo "Processing $orgdb_name_old..."
    ## name without the .db suffix
    ##orgdb_name_short_old=${orgdb_name_old/.db/}
    ##orgdb_name_short=${orgdb_name/.db/}
    ## rename and enter directory
    ##mv $output_dir/$orgdb_name_old $output_dir/$orgdb_name
    ##cd $output_dir/$orgdb_name
    ## Fix DESCRIPTION
    ##sed -i "s/$orgdb_name_old/$orgdb_name/g" DESCRIPTION
    ##sed -i "s/species:.*/species: $description/g" DESCRIPTION
    ##sed -i "s/Entrez/$db_name/g" DESCRIPTION
    ## Fix NAMESPACE
    ##sed -i "s/$orgdb_name_short_old/$orgdb_name_short/g" NAMESPACE
    ## Fix sqlite database
    ##dbpath=inst/extdata/${orgdb_name_short}.sqlite
    ##mv inst/extdata/${orgdb_name_short_old}.sqlite $dbpath
    ##chmod +w $dbpath
    ##sqlite3 $dbpath "UPDATE metadata SET value=\"$description\" WHERE name='SPECIES';"
    ##chmod -w $dbpath
    ## Fix manual pages
    ##for suffix in "BASE.Rd" "ORGANISM.Rd" "_dbconn.Rd"; do
    ##mv man/${orgdb_name_short_old}${suffix} man/${orgdb_name_short}${suffix}
    ##sed -i "s/$orgdb_name_short_old/$orgdb_name_short/g" man/${orgdb_name_short}${suffix} 
    ##done
    ## Fix zzz.R
    ##sed -i "s/$orgdb_name_short_old/$orgdb_name_short/g" R/zzz.R
    ##
    ## Generate OrgDb README.md
    ##
    ##cat << EOF > README.md
    ## $orgdb_name
##Genome-wide annotation package for *$description*, based on
##annotations from [$db_name $db_version]($db_url).
##This package was generated using the tools from
##[https://github.com/elsayed-lab/eupathdb-organismdb](github.com/eupathdb-organismdb).
##
##Installation
##------------
##
##You can install the latest version from Github using:
##
##\`\`\` r
##library('devtools')
##install_github('elsayed-lab/$orgdb_name')
##\`\`\`
##
##Usage
##-----
##
##This package is based on the Bioconductor
##[AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
##interface. As such, the methods for interacting with this package are similar
##to the ways one can interact with other commonly-used annotation packages such as
##[org.Hs.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).
##
##Example usage:
##
##\`\`\`r
##library($orgdb_name)
##
## list available fields to query
##columns($orgdb_name)
##
## get first 10 genes
##gene_ids = head(keys($orgdb_name), 10)
##
## gene names and descriptions
##annotations = AnnotationDbi::select($orgdb_name, 
##                                    keys=gene_ids, 
##                                    keytype='GID', 
##                                    columns=c('CHROMOSOME', 'GENENAME'))
##head(annotations)
##
## GO terms
##go_terms = AnnotationDbi::select($orgdb_name, 
##                                 keys=gene_ids, 
##                                 keytype='GID', 
##                                 columns=c('GO', 'ONTOLOGYALL'))
##head(go_terms)
##
## KEGG pathways
##kegg_paths = AnnotationDbi::select($orgdb_name,
##                                   keys=gene_ids, 
##                                   keytype='GID', 
##                                   columns=c('KEGG_NAME', 'KEGG_PATH'))
##head(kegg_paths)
##\`\`\`
##
##For more information, check out the [AnnotationDbi - Introduction to Annotation
##packages vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).
##
##Additional resources that may be helpful:
##
##1. http://www.bioconductor.org/help/workflows/annotation-data/
##2. http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
##3. http://training.bioinformatics.ucdavis.edu/docs/2012/05/DAV/lectures/annotation/annotation.html
##EOF
##
##cd $cwd
##
##
## TranscriptDB
##
##printf -v txdb_name_old 'TxDb.%s%s.%s.%s' ${genus:0:1} ${species} ${db_name} ${db_version}
##echo "Processing $txdb_name_old..."
##
##mv $output_dir/$txdb_name_old $output_dir/$txdb_name
##cd $output_dir/$txdb_name
##
### Fix DESCRIPTION
##sed -i "s/$txdb_name_old/$txdb_name/g" DESCRIPTION
##sed -i "s/species:.*/species: $description/g" DESCRIPTION
##
## Fix NAMESPACE
##sed -i "s/$txdb_name_old/$txdb_name/g" NAMESPACE
##
## Fix sqlite database
##dbpath=inst/extdata/${txdb_name}.sqlite
##mv inst/extdata/${txdb_name_old}.sqlite $dbpath
##
## Fix Manual pages
##sed -i "s/$txdb_name_old/$txdb_name/g" man/package.Rd
##
##
## Generate TxDb README.md
##
##cat << EOF > README.md
## $txdb_name
##
##Transcript annotation package for *$description*, based on
##annotated genes from [$db_name $db_version]($db_url).
##
##This package was generated using the tools from
##[https://github.com/elsayed-lab/eupathdb-organismdb](github.com/eupathdb-organismdb).
##
##Installation
##------------
##
##You can install the latest version from Github using:
##
##\`\`\` r
##library('devtools')
##install_github('elsayed-lab/$txdb_name')
##\`\`\`
##
##Usage
##-----
##
##This package is based on the Bioconductor
##[AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
##interface. As such, the methods for interacting with this package are similar
##to the ways one can interact with other commonly-used annotation packages such as
##[TxDb.Hsapiens.UCSC.hg19.knownGene](http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html).
##
##Example usage:
##
##\`\`\`r
##library($txdb_name)
##
## list available fields to query
##columns($txdb_name)
##
## get first 10 genes
##gene_ids = head(keys($txdb_name), 10)
##
## gene coordinates and strand
##genes = AnnotationDbi::select($txdb_name, 
##                              keys=gene_ids, 
##                              keytype='GENEID', 
##                              columns=c('TXSTART', 'TXEND', 'TXSTRAND'))
##
##head(genes)
##\`\`\`
##
##For more information, check out the [AnnotationDbi - Introduction to Annotation
##packages vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).
##
##Additional resources that may be helpful:
##
##1. http://www.bioconductor.org/help/workflows/annotation-data/
##2. http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
##3. http://training.bioinformatics.ucdavis.edu/docs/2012/05/DAV/lectures/annotation/annotation.html
##EOF
##
## Install OrgDB and TxDb
##echo "Installing databases"
##cd $cwd
##
##R CMD INSTALL $output_dir/$orgdb_name
##R CMD INSTALL $output_dir/$txdb_name
##
##echo "Done!"
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












#'
#' EuPathDB gene information table gene type parser
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @return Returns a dataframe mapping gene ids to gene types
#' @export
orgdb_gene_types <- function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp <- gzfile(filepath, open='rb')
    } else {
        fp <- file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N <- 1e5
    gene_ids <- c()
    gene_types <- c()

    # Counter to keep track of row number
    i = 1

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
        }
        # Gene type
        else if (grepl("^Gene Type:", x)) {
            gene_ids[i]   = gene_id
            gene_types[i] = .get_value(x)
            i = i + 1
        }
    }
    # close file pointer
    close(fp)

    return(data.frame(GID=gene_ids, TYPE=gene_types))
}





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

#
# Parses a key: value string and returns the value
#
.get_value = function(x) {
    return(gsub(" ","", tail(unlist(strsplit(x, ': ')), n=1), fixed=TRUE))
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
