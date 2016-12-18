#' Perform the array of analyses in the 2016-04 version of clusterProfiler
#'
#' The new version of clusterProfiler has a bunch of new toys.  However, it is more stringent in
#' terms of input in that it now explicitly expects to receive annotation data in terms of a orgdb
#' object.  This is mostly advantageous, but will probably cause some changes in the other ontology
#' functions in the near future.  This function is an initial pass at making something similar to my
#' previous 'simple_clusterprofiler()' but using these new toys.
#'
#' @param sig_genes  Dataframe of genes deemed 'significant.'
#' @param all_genes  Dataframe of all genes in the analysis, primarily for gse analyses.
#' @param orgdb  Name of the orgDb used for gathering annotation data.
#' @param orgdb_from  Name of a key in the orgdb used to cross reference to entrez IDs.
#' @param orgdb_to  List of keys to grab from the orgdb for cross referencing ontologies.
#' @param internal  Used by the 'use_internal_data' flag.
#' @param go_level  How deep into the ontology tree should this dive for over expressed categories.
#' @param pcutoff  P-value cutoff for 'significant' analyses.
#' @param qcutoff  Q-value cutoff for 'significant' analyses.
#' @param fc_column  When extracting vectors of all genes, what column should be used?
#' @param updown  Include the less than expected ontologies?
#' @param permutations  How many permutations for GSEA-ish analyses?
#' @param min_groupsize  Minimum size of an ontology before it is included.
#' @param kegg_prefix  Many KEGG ids need a prefix before they will cross reference.
#' @param mings  What is the minimum ontology group's size?
#' @param kegg_organism  Choose the 3 letter KEGG organism name here.
#' @param categories  How many categories should be plotted in bar/dot plots?
#' @param parallel  Perform slow operations in parallel?
#' @return a list
#' @export
simple_clusterprofiler <- function(sig_genes, all_genes, orgdb="org.Dm.eg.db",
                                   orgdb_from="FLYBASE", orgdb_to="ENTREZID", internal=TRUE,
                                   go_level=3, pcutoff=0.05, qcutoff=0.1, fc_column="logFC", updown="up",
                                   permutations=100, min_groupsize=5, kegg_prefix="Dmel_", mings=5,
                                   kegg_organism="dme", categories=12, parallel=TRUE) {
    requireNamespace(package="clusterProfiler", quietly=TRUE)
    requireNamespace(orgdb)
    org <- loadNamespace(orgdb) ## put the orgDb instance into an environment
    org <- org[[orgdb]] ## Then extract it
    mapper_keys <- AnnotationDbi::keytypes(org)
    all_genenames <- rownames(all_genes)
    orgdb_from <- toupper(orgdb_from)
    orgdb_to <- toupper(orgdb_to)
    ## Interestingly, these bitr calls fail on travis but work fine on my system.
    ## It looks like the version on travis does not require the orgdb.
    all_genes_df <- try(clusterProfiler::bitr(all_genenames, fromType=orgdb_from,
                                              toType=orgdb_to, OrgDb=org))
    if (class(all_genes_df) == "try-error") {
        all_genes_df <- try(clusterProfiler::bitr(all_genenames, fromType=orgdb_from,
                                                  toType=orgdb_to, annoDb=orgdb))
    }
    sig_genenames <- rownames(sig_genes)
    sig_genes_df <- try(clusterProfiler::bitr(sig_genenames, fromType=orgdb_from,
                                              toType=orgdb_to, OrgDb=org))
    if (class(sig_genes_df) == "try-error") {
        sig_genes_df <- try(clusterProfiler::bitr(sig_genenames, fromType=orgdb_from,
                                                  toType=orgdb_to, annoDb=orgdb))
    }
    universe <- AnnotationDbi::keys(org, keytype=orgdb_to)
    all_genes_df <- merge(all_genes, all_genes_df, by.x="row.names", by.y=orgdb_from)
    ## Rename the first column
    colnames(all_genes_df)[1] <- orgdb_from
    if (is.null(all_genes_df[[fc_column]])) {
        stop("The fold change column appears to provide no genes, try another column in the data set.")
    }

    ## Why did I do this? ## Ahh for the GSE analyses, they want ordered gene IDs
    if (updown == "up") {
        all_genes_df <- all_genes_df[ order(all_genes_df[[fc_column]], decreasing=TRUE), ]
    } else {
        all_genes_df <- all_genes_df[ order(all_genes_df[[fc_column]], decreasing=FALSE), ]
    }

    message("Calculating GO groups.")
    ggo_mf <- clusterProfiler::groupGO(gene=sig_genes_df[[orgdb_to]], OrgDb=org,
                                       ont="MF", level=go_level, readable=TRUE)
    ggo_bp <- clusterProfiler::groupGO(gene=sig_genes_df[[orgdb_to]], OrgDb=org,
                                       ont="BP", level=go_level, readable=TRUE)
    ggo_cc <- clusterProfiler::groupGO(gene=sig_genes_df[[orgdb_to]], OrgDb=org,
                                       ont="CC", level=go_level, readable=TRUE)
    group_go <- list(
        "MF" = as.data.frame(DOSE::summary(ggo_mf)),
        "BP" = as.data.frame(DOSE::summary(ggo_bp)),
        "CC" = as.data.frame(DOSE::summary(ggo_cc)))

    message("Calculating enriched GO groups.")
    enrich_results <- list(
        "all_mf" = NULL,
        "sig_mf" = NULL,
        "all_bp" = NULL,
        "sig_bp" = NULL,
        "all_cc" = NULL,
        "sig_cc" = NULL)
    ego_all_mf <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="MF",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=1.0, readable=TRUE)
    ego_sig_mf <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="MF",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=pcutoff, readable=TRUE)
    ego_all_bp <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="BP",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=1.0, readable=TRUE)
    ego_sig_bp <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="BP",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=pcutoff, readable=TRUE)
    ego_all_cc <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="CC",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=1.0, readable=TRUE)
    ego_sig_cc <- clusterProfiler::enrichGO(gene=sig_genes_df[[orgdb_to]], universe=universe,
                                            OrgDb=org, ont="CC",
                                            minGSSize=mings, pAdjustMethod="BH",
                                            pvalueCutoff=pcutoff, readable=TRUE)
    ##} ## End else if we should do a parallel search

    enrich_go <- list(
        "MF_all" = as.data.frame(DOSE::summary(ego_all_mf)),
        "MF_sig" = as.data.frame(DOSE::summary(ego_sig_mf)),
        "BP_all" = as.data.frame(DOSE::summary(ego_all_bp)),
        "BP_sig" = as.data.frame(DOSE::summary(ego_sig_bp)),
        "CC_all" = as.data.frame(DOSE::summary(ego_all_cc)),
        "CC_sig" = as.data.frame(DOSE::summary(ego_sig_cc)))

    message("Performing GSE analyses of gene lists (this is slow).")
    genelist <- as.vector(all_genes_df[[fc_column]])
    names(genelist) <- all_genes_df[[orgdb_to]]
    gse_all_mf <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="MF",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=1.0)
    gse_sig_mf <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="MF",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=pcutoff)
    gse_all_bp <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="BP",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=1.0)
    gse_sig_bp <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="BP",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=pcutoff)
    gse_all_cc <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="CC",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=1.0)
    gse_sig_cc <- clusterProfiler::gseGO(geneList=genelist, OrgDb=org, ont="CC",
                                         nPerm=permutations, minGSSize=min_groupsize,
                                         pvalueCutoff=pcutoff)
    gse_go <- list(
        "MF_all" = as.data.frame(DOSE::summary(gse_all_mf)),
        "MF_sig" = as.data.frame(DOSE::summary(gse_sig_mf)),
        "BP_all" = as.data.frame(DOSE::summary(gse_all_bp)),
        "BP_sig" = as.data.frame(DOSE::summary(gse_sig_bp)),
        "CC_all" = as.data.frame(DOSE::summary(gse_all_cc)),
        "CC_sig" = as.data.frame(DOSE::summary(gse_sig_cc)))

    message("Performing KEGG analyses.")
    ## The help documentation _says_ this uses entrez gene ids, but looking at the environment
    ## created by clusterProfiler belies this.
    kegg_sig_ids <- sig_genes_df[[orgdb_from]]
    all_kegg <- clusterProfiler::enrichKEGG(kegg_sig_ids, organism=kegg_organism,
                                            pvalueCutoff=1.0, use_internal_data=internal)
    enrich_kegg <- clusterProfiler::enrichKEGG(gene=kegg_sig_ids, organism=kegg_organism,
                                               pvalueCutoff=pcutoff, use_internal_data=internal)

    kegg_genelist <- as.vector(all_genes_df[[fc_column]])
    names(kegg_genelist) <- all_genes_df[[orgdb_from]]
    gse_all_kegg <- clusterProfiler::gseKEGG(geneList=kegg_genelist, organism=kegg_organism,
                                             nPerm=permutations, minGSSize=min_groupsize,
                                             pvalueCutoff=1.0, use_internal_data=internal)
    gse_sig_kegg <- clusterProfiler::gseKEGG(geneList=kegg_genelist, organism=kegg_organism,
                                             nPerm=permutations, minGSSize=min_groupsize,
                                             pvalueCutoff=pcutoff, use_internal_data=internal)
    gse_all_mkegg <- NULL
    gse_sig_mkegg <- NULL
    ## The following does not work anylonger with weird errors that I am not inclined to chase down.
    ##Sys.sleep(3)
    ##gse_all_mkegg <- clusterProfiler::gseMKEGG(geneList=kegg_genelist, organism=kegg_organism,
    ##                                           pvalueCutoff=1.0)
    ##Sys.sleep(3)
    ##gse_sig_mkegg <- clusterProfiler::gseMKEGG(geneList=kegg_genelist, organism=kegg_organism,
    ##pvalueCutoff=pcutoff)
    kegg_data <- list(
        "kegg_all" = as.data.frame(DOSE::summary(all_kegg)),
        "kegg_sig" = as.data.frame(DOSE::summary(enrich_kegg)),
        "kegg_gse_all" = as.data.frame(DOSE::summary(gse_all_kegg)),
        "kegg_gse_sig" = as.data.frame(DOSE::summary(gse_sig_kegg)))
##        "kegg_gsem_all" <- as.data.frame(summary(gse_all_mkegg)),
##        "kegg_gsem_sig" <- as.data.frame(summary(gse_sig_mkegg)))

    message("Plotting results.")
    map_sig_mf <- map_sig_bp <- map_sig_cc <- NULL
    tt <- try(DOSE::enrichMap(ego_sig_mf))
    if (class(tt) != "try-error") {
        map_sig_mf <- recordPlot()
    }
    tt <- try(DOSE::enrichMap(ego_sig_bp))
    if (class(tt) != "try-error") {
        map_sig_bp <- recordPlot()
    }
    tt <- try(DOSE::enrichMap(ego_sig_cc))
    if (class(tt) != "try-error") {
        map_sig_cc <- recordPlot()
    }
    net_sig_mf <- net_sig_bp <- net_sig_cc <- NULL
    tt <- try(DOSE::cnetplot(ego_sig_mf, categorySize="pvalue", foldChange=genelist))
    if (class(tt) != "try-error") {
        net_sig_mf <- recordPlot()
    }
    tt <- try(DOSE::cnetplot(ego_sig_bp, categorySize="pvalue", foldChange=genelist))
    if (class(tt) != "try-error") {
        net_sig_bp <- recordPlot()
    }
    tt <- try(DOSE::cnetplot(ego_sig_cc, categorySize="pvalue", foldChange=genelist))
    if (class(tt) != "try-error") {
        net_sig_cc <- recordPlot()
    }
    tree_sig_mf <- tree_sig_bp <- tree_sig_cc <- NULL
    tree_mf <- try(clusterProfiler::plotGOgraph(ego_sig_mf), silent=TRUE)
    if (class(tree_mf) != "try-error") {
        tree_sig_mf <- recordPlot()
    }
    tree_bp <- try(clusterProfiler::plotGOgraph(ego_sig_bp), silent=TRUE)
    if (class(tree_bp) != "try-error") {
        tree_sig_bp <- recordPlot()
    }
    tree_cc <- try(clusterProfiler::plotGOgraph(ego_sig_cc), silent=TRUE)
    if (class(tree_cc) != "try-error") {
        tree_sig_cc <- recordPlot()
    }

    plotlist <- list(
        "ggo_mf_bar" = barplot(ggo_mf, drop=TRUE, showCategory=categories),
        "ggo_bp_bar" = barplot(ggo_bp, drop=TRUE, showCategory=categories),
        "ggo_cc_bar" = barplot(ggo_cc, drop=TRUE, showCategory=categories),
        "ego_all_mf" = barplot(ego_all_mf, showCategory=categories, drop=TRUE),
        "ego_all_bp" = barplot(ego_all_bp, showCategory=categories, drop=TRUE),
        "ego_all_cc" = barplot(ego_all_cc, showCategory=categories, drop=TRUE),
        "ego_sig_mf" = barplot(ego_sig_mf, showCategory=categories, drop=TRUE),
        "ego_sig_bp" = barplot(ego_sig_bp, showCategory=categories, drop=TRUE),
        "ego_sig_cc" = barplot(ego_sig_cc, showCategory=categories, drop=TRUE),
        "dot_all_mf" = clusterProfiler::dotplot(ego_all_mf),
        "dot_all_bp" = clusterProfiler::dotplot(ego_all_bp),
        "dot_all_cc" = clusterProfiler::dotplot(ego_all_cc),
        "dot_sig_mf" = clusterProfiler::dotplot(ego_sig_mf),
        "dot_sig_bp" = clusterProfiler::dotplot(ego_sig_bp),
        "dot_sig_cc" = clusterProfiler::dotplot(ego_sig_cc),
        "map_sig_mf" = map_sig_mf,
        "map_sig_bp" = map_sig_bp,
        "map_sig_cc" = map_sig_cc,
        "net_sig_mf" = net_sig_mf,
        "net_sig_bp" = net_sig_bp,
        "net_sig_cc" = net_sig_cc,
        "tree_sig_mf" = tree_sig_mf,
        "tree_sig_bp" = tree_sig_bp,
        "tree_sig_cc" = tree_sig_cc)

    retlist <- list(
        "all_mappings" = all_genes_df,
        "sig_mappings" = sig_genes_df,
        "group_go" = group_go,
        "enrich_go" = enrich_go,
        "gse_go" = gse_go,
        "kegg_data" = kegg_data,
        "plots" = plotlist)
    return(retlist)
}

#' Set up appropriate option sets for clusterProfiler
#'
#' This hard-sets some defaults for orgdb/kegg databases when using clusterProfiler.
#'
#' @param species  Currently it only works for humans and fruit flies.
cp_options <- function(species) {
    if (species == "dmelanogaster") {
        options <- list(
            orgdb = "org.Dm.eg.db",
            orgdb_from = "FLYBASE",
            orgdb_to = c("ENSEMBL","SYMBOL","ENTREZID"),
            kegg_prefix = "Dmel_",
            kegg_organism = "dme",
            kegg_id_column = "FLYBASECG")
    } else if (species == "hsapiens") {
        options <- list(
            orgdb = "org.Hs.eg.db",
            orgdb_from = "ENSEMBL",
            orgdb_to = c("ENSEMBL","SYMBOL","ENTREZID"),
            kegg_prefix = "Hsa_",
            kegg_organism = "hsa",
            kegg_id_column = "")
    }
    return(options)
}

#' Generic enrichment using clusterProfiler.
#'
#' culsterProfiler::enricher provides a quick and easy enrichment analysis given a set of
#' siginficant' genes and a data frame which connects each gene to a category.
#'
#' @param sig_genes Set of 'significant' genes as a table.
#' @param de_table All genes from the original analysis.
#' @param goids_df Dataframe of GO->ID matching the gene names of sig_genes to GO categories.
#' @return Table of 'enriched' categories.
simple_cp_enricher <- function(sig_genes, de_table, goids_df=NULL) {
    all_genenames <- rownames(de_table)
    sig_genenames <- rownames(sig_genes)
    enriched <- clusterProfiler::enricher(sig_genenames, TERM2GENE=goids_df)
    retlist <- list(
        "enriched" = as.data.frame(DOSE::summary(enriched)))
    return(retlist)
}

## EOF
