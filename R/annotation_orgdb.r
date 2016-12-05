## Everything in this was written by Keith, I stole it with his permission and incorporated it here.
## I might change a few function declarations in this: tbl_df
## for example is not included in any import declarations and so I will likely
## re-call it with dplyr::tbl_df() -- or add it to the imports.

#' I see no reason to have load_host_annotations and load_parasite_annotations.
#'
#' Thus I am making them both into aliases to load_annotations.
#'
#' @param ... Arguments to be passed to load_annotations.
#' @export
load_parasite_annotations <- function(...) {
    load_annotations(...)
}

#' Load organism annotation data (parasite).
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb OrganismDb instance.
#' @param include_go  Ask the Dbi for gene ontology information?
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param keytype mmm the key type used?
#' @param fields Columns included in the output.
#' @param sum_exons Perform a sum of the exons in the data set?
#' @return Table of geneids, chromosomes, descriptions, strands, types, and lengths.
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' one_gene <- load_annotations(org, c("LmJF.01.0010"))
#' }
#' @export
load_annotations <- function(orgdb, gene_ids=NULL, include_go=FALSE, keytype="ENSEMBL",
                                      ## fields=c("CHR", "GENENAME", "TXSTRAND",
                                      fields=NULL, sum_exons=FALSE) {
    ## "TXSTART", "TXEND", "TYPE")) {
    keytype <- toupper(keytype)
    all_fields <- AnnotationDbi::columns(orgdb)
    if (is.null(fields)) {
        fields <- c("CHR", "GENENAME", "TXSTRAND", "TXSTART", "TXEND", "TYPE")
    } else if (fields[[1]] == "all") {
        fields <- all_fields
    } else {
        fields <- toupper(fields)
    }

    if (sum(fields %in% all_fields) != length(fields)) {
        message(toString(fields %in% all_fields))
        stop(paste0("Some of the fields requested are not in the available, we do have the following ", toString(all_fields)))
    }
    ## fields <- c("CHR", "GENENAME", "TXSTRAND", "TXSTART", "TXEND", "TYPE")

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

    transcript_length <- NULL
    gene_info <- AnnotationDbi::select(orgdb,
                                       keys=gene_ids,
                                       keytype=keytype,
                                       columns=fields)

    ## Compute total transcript lengths (for all exons)
    ## https://www.biostars.org/p/83901/
    gene_exons <- GenomicFeatures::exonsBy(orgdb, by="gene")
    transcripts <- GenomicFeatures::transcripts(orgdb)
    colnames(gene_info) <- tolower(colnames(gene_info))
    gene_info_ret <- NULL
    if (isTRUE(sum_exons)) {
        message("Summing exon lengths, this takes a while.")
        lengths <- lapply(gene_exons, function(x) {
            sum(BiocGenerics::width(GenomicRanges::reduce(x)))
        })
        message("Adding exon lengths to the gene_exons.")
        lengths <- as.data.frame(unlist(lengths))
        colnames(lengths) <- "transcript_length"
        gene_info <- merge(gene_info, lengths, by.x=keytype, by.y="row.names")
    }
    rownames(gene_info) <- make.names(gene_info[[1]], unique=TRUE)

    retlist <- list(
        "genes" = gene_info,
        "transcripts" = transcripts)
    return(retlist)
}

#' Load organism annotation data (mouse/human).
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param keytype a, umm keytype? I need to properly read this code.
#' @param fields Columns to include in the output.
#' @param biomart_dataset Name of the biomaRt dataset to query for gene type.
#' @return a table of gene information
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' host <- load_host_annotations(org, c("a","b"))
#' }
#' @export
load_host_annotations <- function(orgdb, gene_ids=NULL, keytype="ENSEMBL",
                                 fields=c("TXCHROM", "GENENAME", "TXSTRAND",
                                          "TXSTART", "TXEND"),
                                 biomart_dataset="hsapiens_gene_ensembl") {
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
    gene_info <- dplyr::select_(orgdb,
                        keys=gene_ids,
                        keytype=keytype,
                        columns=fields)
    ## Remove any entries which are either deprecated (e.g. ENSG00000006074) or
    ## missing transcript information (e.g. ENSG00000050327).
    ## gene_info <- gene_info[!is.na(gene_info$TXSTART),]
    transcript_length <- NULL
    ## Convert to tbl_df and reorganize
    gene_info <- dplyr::tbl_df(gene_info) %>%
        dplyr::mutate_(transcript_length=abs("TXEND" - "TXSTART") + 1) %>%

    ## filter(keytype %in% gene_ids) %>%
    ## Are TXSTRAND and friends quotable?
    dplyr::select_(
        gene_id=get(keytype),
        chromosome="TXCHROM",
        description="GENENAME",
        strand="TXSTRAND",
        transcript_length=transcript_length)

    ## Get gene biotype
    ## Main server temporarily unavailable (2015/11/09)
    ## ensembl_mart <- useMart(biomart="ensembl")
    ensembl_mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                     host="www.ensembl.org", biomart_dataset)
    biomart <- biomaRt::useDataset(biomart_dataset, mart=ensembl_mart)
    biomart_genes <- biomaRt::getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=biomart)
    gene_info$type <- biomart_genes$gene_biotype[match(gene_info$gene_id,
                                                       biomart_genes$ensembl_gene_id)]
    return(gene_info)
}

#' Retrieve GO terms associated with a set of genes.
#'
#' AnnotationDbi provides a reasonably complete set of GO mappings between gene ID and
#' ontologies.  This will extract that table for a given set of gene IDs.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Identifiers of the genes to retrieve annotations.
#' @param keytype the mysterious keytype returns yet again to haunt my dreams
#' @return Data frame of gene IDs, go terms, and names.
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' go_terms <- load_go_terms(org, c("a","b"))
#' }
#' @export
load_go_terms <- function(orgdb, gene_ids, keytype="ENSEMBL") {
    requireNamespace("GO.db")
    requireNamespace("magrittr")
    if (is.null(gene_ids)) {
        stop("gene_ids may not be null.")
    }
    go_terms <- suppressWarnings(
        AnnotationDbi::select(orgdb,
                              "keys" = gene_ids,
                              "keytype" = keytype,
                              "columns" = c("GO"))[,c(1,2)]
    )
    ## Deduplicate
    go_terms <- go_terms[!duplicated(go_terms), ]
    go_terms <- go_terms[!is.na(go_terms[["GO"]]), ]
    go_term_names <- suppressWarnings(
        AnnotationDbi::select(GO.db::GO.db,
                               keys=unique(go_terms$GO),
                               columns=c("TERM", "GOID", "ONTOLOGY")))
    go_terms <- merge(go_terms, go_term_names, by.x="GO", by.y="GOID")

    # Remove redundant annotations which differ only in source/evidence
    # and rename ONTOLOGYALL column
    #unique(go_terms %>% rename(ONTOLOGY=ONTOLOGYALL) %>% na.omit())
    go_terms <- unique(dplyr::tbl_df(go_terms) %>% na.omit())
    return(go_terms)
}

#' Creates a gene/KEGG mapping dataframe.
#'
#' In much the same way AnnotationDbi provides GO data, it also provides KEGG data.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Identifiers of the genes to retrieve annotations.
#' @param keytype the keytype, damn I really need to read this code
#' @return Df of kegg mappings
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' kegg_data <- load_kegg_mapping(org, c("a","b"))
#' }
#' @export
load_kegg_mapping <- function(orgdb, gene_ids, keytype='ENSEMBL') {
    kegg_mapping <- suppressWarnings(
        dplyr::tbl_df(AnnotationDbi::select(orgdb, keys=gene_ids,
                            keytype=keytype,
                            columns=c('KEGG_PATH'))) %>% na.omit()
    )
    colnames(kegg_mapping) <- c('gene', 'category')
    ## goseq does not support tbl_df instances
    return(as.data.frame(kegg_mapping))
}

#' Creates a KEGG pathway/description mapping dataframe.
#'
#' Use AnnotationDbi to map descriptions of KEGG pathways to gene IDs.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Identifiers of the genes to retrieve annotations.
#' @param keytype as per the previous functions, I don't know what this does yet
#' @return Character list of pathways.
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' pathnames <- load_kegg_pathways(org, c("a","b","c")
#' }
#' @export
load_kegg_pathways <- function(orgdb, gene_ids, keytype='ENSEMBL') {
    kegg_pathways <- suppressWarnings(
        dplyr::tbl_df(AnnotationDbi::select(orgdb, keys=gene_ids,
                            keytype=keytype,
                            columns=c('KEGG_PATH', 'KEGG_NAME', 'KEGG_CLASS',
                                        'KEGG_DESCRIPTION')))
    )
    kegg_pathways <- kegg_pathways %>%
        na.omit() %>%
        ## AnnotationDbi::select(KEGG_PATH, KEGG_NAME, KEGG_CLASS, KEGG_DESCRIPTION)
        ## I think these should be quoted
        dplyr::select_("KEGG_PATH", "KEGG_NAME", "KEGG_CLASS", "KEGG_DESCRIPTION")
        #select(-get(keytype))
    colnames(kegg_pathways) <- c('pathway', 'name', 'class', 'description')
    return(kegg_pathways)
}

#' Maps KEGG identifiers to ENSEMBL gene ids.
#'
#' Takes a list of KEGG gene identifiers and returns a list of ENSEMBL
#' ids corresponding to those genes.
#'
#' @param kegg_ids List of KEGG identifiers to be mapped.
#' @return Ensembl IDs as a character list.
#' @seealso \link[KEGGREST]{keggGet}
#' @examples
#' \dontrun{
#' ensembl_list <- kegg_to_ensembl("a")
#' }
#' @export
kegg_to_ensembl <- function(kegg_ids) {
    ## query gene ids 10 at a time (max allowed)
    result <- c()
    for (x in split(kegg_ids, ceiling(seq_along(kegg_ids) / 3))) {
        print(x)
        query <- KEGGREST::keggGet(x)
        for (item in query) {
            dblinks <- item$DBLINKS
            ensembl_id <- dblinks[grepl('Ensembl', dblinks)]
            if (length(ensembl_id) > 0) {
                result <- append(result, substring(ensembl_id, 10))
                ## TESTING
                if (length(ensembl_id) > 1) {
                    warning(sprintf("One to many KEGG mapping for gene %s", x))
                }
            }
        }
    }
    return(result)
}

#' Generate GENE/KEGG mapping.
#'
#' This uses KEGGREST and related function kegg_to_ensembl() to associate genes to KEGG pathways.
#'
#' @param pathways Vector of KEGG pathway IDs returned from call to keggLink()
#'        e.g. "path:mmu05134".
#' @param org_abbreviation KEGG identifier for the species of interest (e.g.
#'        "hsa" for Homo sapiens).
#' @param verbose talky talky?
#' @return Df mapping kegg and gene IDs.
#' @seealso \link[KEGGREST]{keggLink}
#' @examples
#' \dontrun{
#' kegg_df <- generate_gene_kegg_mapping(path, org)
#' }
#' @export
generate_gene_kegg_mapping <- function(pathways, org_abbreviation, verbose=FALSE) {
    ## data frame to store kegg gene mapping
    kegg_mapping <- data.frame()
    ## Iterate over pathways and query genes for each one
    for (pathway in pathways) {
        if (verbose) {
            message(sprintf("Processing for pathway %s", pathway))
        }
        ## Get genes in pathway
        kegg_ids <- as.character(KEGGREST::keggLink(org_abbreviation, pathway))
        gene_ids <- kegg_to_ensembl(kegg_ids)
        kegg_mapping <- unique(rbind(kegg_mapping,
                                     data.frame(
                                         "gene" = gene_ids,
                                         "pathway" = pathway)))
    }
    return(kegg_mapping)
}

#' Generate a KEGG PATHWAY / description mapping.
#'
#' Make an easier to use df of KEGG -> descriptions using keggGet.
#'
#' @param pathways Vector of KEGG pathway identifiers.
#' @param verbose talk talk?
#' @return Data frame describing some kegg pathways
#' @seealso \link[KEGGREST]{keggGet}
#' @examples
#'  mapping <- generate_kegg_pathway_mapping(c("hsa00040", "hsa00100"))
#' @export
generate_kegg_pathway_mapping <- function(pathways, verbose=FALSE) {
    ## result data frame
    kegg_pathways <- data.frame()
    for (pathway in pathways) {
        if (verbose) {
            message(sprintf("Processing for pathway %s", pathway))
        }
        ## Get pathway info

        ##        meta <- try(KEGGREST::keggGet(pathway, "kgml")[[1]])
        meta <- try(KEGGREST::keggGet(pathway)[[1]])
        if (class(meta) != "try-error") {
            pathway_desc  <- ifelse(is.null(meta[["DESCRIPTION"]]), '', meta[["DESCRIPTION"]])
            pathway_class <- ifelse(is.null(meta[["CLASS"]]), '', meta[["CLASS"]])
            kegg_pathways <- rbind(kegg_pathways,
                                   data.frame("pathway" = pathway,
                                              "name" = meta[["PATHWAY_MAP"]],
                                              "class" = pathway_class,
                                              "description" = pathway_desc))
        }
    }
    return(kegg_pathways)
}

#' Load the appropriate orgDb environment for a given species.
#'
#' Ok, so these are a bit more complex than I realized.  The heirarchy as I now understand it
#' (probably wrong) is that orgdb objects provide ID mappings among the various DBs.  txdb objects
#' provide the actual annotation information, and organismdbs acquire both (but only exist for a few
#' species). Let's face it, I will never remember that the yeast orgdb is 'org.Sc.sgd.something'.  This
#' function is intended to make that process easier.  Feed it a species name which makes sense:
#' 'homo_sapiens' and it will assume you mean orgdb.whatever and load that into your
#' environment. This should also make a reasonable attempt at installing the appropriate orgdb if it
#' is not already in your R library tree.
#'
#' @param species Human readable species name
#' @return orgdb object for the relevant species, or an error if I don't have a mapping for it.
#' @examples
#'  object <- choose_orgdb("homo_sapiens")
#' @export
choose_orgdb <- function(species="saccharomyces_cerevisiae") {
    used_species <- list(
        "homo_sapiens" = c("org.Hs.eg.db", "bioconductor"),
        "mus_musculus" = c("org.Mm.eg.db", "bioconductor"),
        "leishmania_major" = c("org.LmjF.tritryp.db", "elsayed-lab"),
        "trypanosoma_cruzi_clb" = c("org.TcCLB.clb.tritryp.db", "elsayed-lab"),
        "trypanosoma_cruzi_esmer" = c("org.TcCLB.esmer.tritryp.db", "elsayed-lab"),
        "trypanosoma_cruzi_nonesmer" = c("org.TcCLB.nonesmer.tritryp.db", "elsayed-lab"),
        "drosophila_melanogaster" = c("org.Dm.eg.db", "bioconductor"),
        "saccharomyces_cerevisiae" = c("org.Sc.sgd.db", "bioconductor")
    )
    found <- 0
    org <- NULL
    avail_keys <- NULL
    if (species %in% names(used_species)) {
        message("I have a species name mapped to orgdb for ", species, ".")
    } else {
        stop("I have not yet assigned an orgdb to this species.")
    }
    try_list <- used_species[[species]]
    try_orgdb <- try_list[[1]]
    try_source <- try_list[[2]]
    org <- try(loadNamespace(try_orgdb))  ## try loading it into the environment named 'org'
    if (class(org) == "try-error") {  ## if that was not found, install it.
        if (try_source == "bioconductor") {
            installedp <- require.auto(try_orgdb)  ## install from bioconductor
        } else {
            installedp <- require.auto(paste0(try_source, "/", try_orgdb))  ## install from github
        }
        org <- loadNamespace(try_orgdb)

    }
    org <- org[[try_orgdb]]
    if (is.null(org)) {
        stop("Did not extract the relevant orgDb.")
    } else {
        avail_keys <- AnnotationDbi::keytypes(org)
        message("Attached the orgDb for ", species, " with keys: ", toString(avail_keys), ".")
    }
    return(org)
}

#' Load the appropriate TxDb environment for a given species.
#'
#' Ok, so these are a bit more complex than I realized.  The heirarchy as I now understand it
#' (probably wrong) is that orgdb objects provide ID mappings among the various DBs.  txdb objects
#' provide the actual annotation information, and organismdbs acquire both (but only exist for a few
#' species). Let's face it, I will never remember that the yeast orgdb is 'org.Sc.sgd.something'.  This
#' function is intended to make that process easier.  Feed it a species name which makes sense:
#' 'homo_sapiens' and it will assume you mean orgdb.whatever and load that into your
#' environment. This should also make a reasonable attempt at installing the appropriate orgdb if it
#' is not already in your R library tree.
#'
#' @param species Human readable species name
#' @return orgdb object for the relevant species, or an error if I don't have a mapping for it.
#' @examples
#'  object <- choose_txdb("homo_sapiens")
#' @export
choose_txdb <- function(species="saccharomyces_cerevisiae") {
    used_species <- list(
        "homo_sapiens" = c("TxDb.Hsapiens.UCSC.hg38.knownGene", "bioconductor"),
        "mus_musculus" = c("TxDb.Mmusculus.UCSC.mm10.knownGene", "bioconductor"),
        "leishmania_major" = c("TxDb.LmajorFriedlin.tritryp27.genes", "elsayed-lab"),
        "trypanosoma_cruzi_clb" = c("TxDb.TcruziCLBrener.tritryp27.genes", "elsayed-lab"),
        "trypanosoma_cruzi_sylvio" = c("TxDb.TcruziSylvio.tritryp27.genes", "elsayed-lab"),
        "trypanosoma_cruzi_esmer" = c("TxDb.TcruziCLBrenerEsmer.tritryp27.genes", "elsayed-lab"),
        "trypanosoma_cruzi_nonesmer" = c("TxDb.TcruziCLBrenerNonEsmer.tritryp27.genes", "elsayed-lab"),
        "drosophila_melanogaster" = c("TxDb.Dmelanogaster.UCSC.dm6.ensGene", "bioconductor"),
        "saccharomyces_cerevisiae" = c("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "bioconductor")
    )
    found <- 0
    tx <- NULL
    avail_keys <- NULL
    if (species %in% names(used_species)) {
        message("I have a species name mapped to txdb for ", species, ".")
    } else {
        stop("I have not yet assigned an txdb to this species.")
    }
    try_list <- used_species[[species]]
    try_txdb <- try_list[[1]]
    try_source <- try_list[[2]]
    tx <- try(loadNamespace(try_txdb))  ## try loading it into the environment named 'tx'
    if (class(tx) == "try-error") {  ## if that was not found, install it.
        if (try_source == "bioconductor") {
            installedp <- require.auto(try_txdb)  ## install from bioconductor
        } else {
            installedp <- require.auto(paste0(try_source, "/", try_txdb))  ## install from github
        }
        tx <- loadNamespace(try_txdb)

    }
    tx <- tx[[try_txdb]]
    ## avail_namespaces <- ls(paste0("package:", try_txdb))
    if (is.null(tx)) {
        stop("Did not extract the relevant txDb.")
    } else {
        avail_keys <- AnnotationDbi::keytypes(tx)
        message("Attached the txDb for ", species, " with keys: ", toString(avail_keys), ".")
    }
    return(tx)
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
#' @seealso \link[AnnotationDbi]{select}
#' @examples
#' \dontrun{
#' host <- load_host_annotations(org, c("a","b"))
#' }
#' @export
orgdb_idmap <- function(orgdb, gene_ids=NULL, mapto=c('ensembl'), keytype="geneid") {
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
    gene_info <- AnnotationDbi::select(orgdb, keytype=keytype, keys=gene_ids, columns=mapto)
    colnames(gene_info) <- tolower(colnames(gene_info))
    return(gene_info)
}

## EOF
