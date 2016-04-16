## Time-stamp: <Sat Apr 16 00:36:57 2016 Ashton Trey Belew (abelew@gmail.com)>
## Everything in this was written by Keith, I stole it with his permission and incorporated it here.
## I might change a few function declarations in this: tbl_df
## for example is not included in any import declarations and so I will likely
## re-call it with dplyr::tbl_df() -- or add it to the imports.

#' Load organism annotation data (parasite)
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype mmm the key type used
#' @param fields Columns to include in the output
#' @return a table of geneids, chromosomes, descriptions, strands, types, and lengths
#' @export
load_parasite_annotations <- function(orgdb, gene_ids, keytype='ENSEMBL',
                            fields=c('CHR', 'GENENAME', 'TXSTRAND',
                                     'TXSTART', 'TXEND', 'TYPE')) {
    ## Gene info
    ## Note querying by "GENEID" will exclude noncoding RNAs
    transcript_length <- NULL
    gene_info <- AnnotationDbi::select(orgdb,
                                       keys=gene_ids,
                                       keytype=keytype,
                                       columns=fields)

    ## Compute total transcript lengths (for all exons)
    ## https://www.biostars.org/p/83901/
    gene_exons <- GenomicFeatures::exonsBy(orgdb, by='gene')
    lengths <- lapply(gene_exons, function(x) {
        sum(BiocGenerics::width(GenomicRanges::reduce(x)))
    })
    lengths <- as.data.frame(unlist(lengths))
    colnames(lengths) <- 'transcript_length'

    gene_info <- merge(gene_info, lengths, by.x=keytype, by.y='row.names')

    ## Convert to tbl_df and reorganize
    ret <- dplyr::tbl_df(gene_info) %>%
        AnnotationDbi::select(
            gene_id     = get(keytype),
            chromosome  = CHR,
            description = GENENAME,
            strand      = TXSTRAND,
            type        = TYPE,
            transcript_length = transcript_length
        )
    return(ret)
}

#' Load organism annotation data (mouse/human)
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the OrganismDbi interface.
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype a, umm keytype? I need to properly read this code.
#' @param fields Columns to include in the output
#' @param biomart_dataset Name of the biomaRt dataset to query for gene type
#' @return a table of gene information
#' @export
load_host_annotations <- function(orgdb, gene_ids, keytype='ENSEMBL',
                                 fields=c('TXCHROM', 'GENENAME', 'TXSTRAND',
                                          'TXSTART', 'TXEND'),
                                 biomart_dataset='hsapiens_gene_ensembl') {
    ## Gene info
    ## Note querying by "GENEID" will exclude noncoding RNAs
    gene_info <- AnnotationDbi::select(orgdb,
                        keys=gene_ids,
                        keytype=keytype,
                        columns=fields)
    ## Remove any entries which are either deprecated (e.g. ENSG00000006074) or
    ## missing transcript information (e.g. ENSG00000050327).
    ## gene_info <- gene_info[!is.na(gene_info$TXSTART),]
    transcript_length <- NULL
    ## Convert to tbl_df and reorganize
    gene_info <- dplyr::tbl_df(gene_info) %>%
        dplyr::mutate(transcript_length=abs(TXEND - TXSTART) + 1) %>%

    ## filter(keytype %in% gene_ids) %>%
    ## Are TXSTRAND and friends quotable?
    AnnotationDbi::select(
        gene_id     = get(keytype),
        chromosome  = TXCHROM,
        description = GENENAME,
        strand      = TXSTRAND,
        transcript_length = transcript_length
    )

    ## Get gene biotype
    ## Main server temporarily unavailable (2015/11/09)
    ## ensembl_mart <- useMart(biomart="ensembl")
    ensembl_mart <- bioMart::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                     host="www.ensembl.org", biomart_dataset)
    biomart <- bioMart::useDataset(biomart_dataset, mart=ensembl_mart)
    biomart_genes <- bioMart::getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=biomart)
    gene_info$type <- biomart_genes$gene_biotype[match(gene_info$gene_id,
                                                       biomart_genes$ensembl_gene_id)]
    return(gene_info)
}

#' Retrieve GO terms associated with a set of genes
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype the mysterious keytype returns yet again to haunt my dreams
#' @return a data frame of gene IDs, go terms, and names.
#' @export
load_go_terms <- function(orgdb, gene_ids, keytype='ENSEMBL') {
    go_terms <- suppressWarnings(
        #tbl_df(AnnotationDbi::select(orgdb, keys=gene_ids,
        #                    keytype=keytype,
        #                    columns=c('GO', 'TERM', 'ONTOLOGYALL')))
        AnnotationDbi::select(orgdb, keys=gene_ids,
                                     keytype=keytype,
                                     columns=c("GO"))[,c(1,2)]
    )
    # Deduplicate
    go_terms <- go_terms[!duplicated(go_terms),]
    requireNamespace("GO.db")
    requireNamespace("magrittr")
    go_term_names <- suppressWarnings(
        AnnotationDbi::select(GO.db::GO.db, keys=unique(go_terms$GO),
                              columns=c("TERM", "GOID", "ONTOLOGY"))
    )
    requireNamespace("magrittr")
    #go_terms$TERM <- go_term_names$TERM[match(go_terms$GO, go_term_names$GOID)]
    go_terms <- merge(go_terms, go_term_names, by.x='GO', by.y='GOID')

    # Remove redundant annotations which differ only in source/evidence
    # and rename ONTOLOGYALL column
    #unique(go_terms %>% rename(ONTOLOGY=ONTOLOGYALL) %>% na.omit())
    unique(dplyr::tbl_df(go_terms) %>% na.omit())
    return(go_terms)
}

#' Creates a gene/KEGG mapping dataframe
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype the keytype, damn I really need to read this code
#' @return a df of kegg mappings
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

#' Creates a KEGG pathway/description mapping dataframe
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype as per the previous functions, I don't know what this does yet
#' @return a character list of pathways
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
        AnnotationDbi::select("KEGG_PATH", "KEGG_NAME", "KEGG_CLASS", "KEGG_DESCRIPTION")
        #select(-get(keytype))
    colnames(kegg_pathways) <- c('pathway', 'name', 'class', 'description')
    return(kegg_pathways)
}

#' Maps KEGG identifiers to ENSEMBL gene ids
#'
#' Takes a list of KEGG gene identifiers and returns a list of ENSEMBL
#' ids corresponding to those genes.
#'
#' @param kegg_ids List of KEGG identifiers to be mapped
#' @return some ensembl IDs
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

#' Generate GENE/KEGG mapping
#'
#' @param pathways Vector of KEGG pathway IDs returned from call to keggLink()
#'        e.g. "path:mmu05134"
#' @param org_abbreviation KEGG identifier for the species of interest (e.g.
#'        "hsa" for Homo sapiens)
#' @param verbose talky talky?
#' @return a df mapping kegg and gene IDs
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
                                     data.frame(gene=gene_ids, pathway=pathway)))
    }
    return(kegg_mapping)
}

#' Generate a KEGG PATHWAY / description mapping
#'
#' @param pathways Vector of KEGG pathway identifiers
#' @param verbose talk talk?
#' @return a data frame describing some kegg pathways
#' @export
generate_kegg_pathway_mapping <- function(pathways, verbose=FALSE) {
    ## result data frame
    kegg_pathways <- data.frame()
    for (pathway in pathways) {
        if (verbose) {
            message(sprintf("Processing for pathway %s", pathway))
        }
        ## Get pathway info
        meta <- KEGGREST::keggGet(pathway)[[1]]
        pathway_desc  <- ifelse(is.null(meta$DESCRIPTION), '', meta$DESCRIPTION)
        pathway_class <- ifelse(is.null(meta$CLASS), '', meta$CLASS)
        kegg_pathways <- rbind(kegg_pathways,
                               data.frame("pathway" = pathway,
                                          "name" = meta$PATHWAY_MAP,
                                          "class" = pathway_class,
                                          "description" = pathway_desc))
    }
    return(kegg_pathways)
}

## EOF
