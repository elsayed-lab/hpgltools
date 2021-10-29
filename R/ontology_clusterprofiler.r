## ontology_clusterprofiler.r: Methods to simplify using clusterProfiler.  I
## think clusterprofiler is probably the most complete and comprehensive GSEA
## toolkit.  It is not necessarily as easy to use as I might prefer.  This seeks
## to fill in some corner case.

#' Perform the array of analyses in the 2016-04 version of clusterProfiler
#'
#' The new version of clusterProfiler has a bunch of new toys.  However, it is
#' more stringent in terms of input in that it now explicitly expects to receive
#' annotation data in terms of a orgdb object.  This is mostly advantageous, but
#' will probably cause some changes in the other ontology functions in the near
#' future.  This function is an initial pass at making something similar to my
#' previous 'simple_clusterprofiler()' but using these new toys.
#'
#' @param sig_genes Dataframe of genes deemed 'significant.'
#' @param de_table Dataframe of all genes in the analysis, primarily for gse
#'  analyses.
#' @param orgdb Name of the orgDb used for gathering annotation data.
#' @param orgdb_from Name of a key in the orgdb used to cross reference to entrez IDs.
#' @param orgdb_to List of keys to grab from the orgdb for cross referencing
#'  ontologies.
#' @param go_level How deep into the ontology tree should this dive for over
#'  expressed categories.
#' @param pcutoff P-value cutoff for 'significant' analyses.
#' @param qcutoff Q-value cutoff for 'significant' analyses.
#' @param fc_column When extracting vectors of all genes, what column should be used?
#' @param second_fc_column When extracting vectors of all genes, what column
#'  should be tried the second time around?
#' @param updown Include the less than expected ontologies?
#' @param permutations How many permutations for GSEA-ish analyses?
#' @param min_groupsize Minimum size of an ontology before it is included.
#' @param kegg_prefix Many KEGG ids need a prefix before they will cross reference.
#' @param kegg_organism Choose the 3 letter KEGG organism name here.
#' @param do_gsea Perform gsea searches?
#' @param categories How many categories should be plotted in bar/dot plots?
#' @param excel Print the results to an excel file?
#' @param do_david Attempt to use the DAVID database for a search?
#' @param david_id Which column to use for cross-referencing to DAVID?
#' @param david_user Default registered username to use.
#' @return a list
#' @seealso [clusterProfiler] [AnnotationDbi] [KEGGREST]
#' @examples
#' \dontrun{
#'  holyasscrackers <- simple_clusterprofiler(gene_list, all_genes, "org.Dm.eg.db")
#' }
#' @export
simple_clusterprofiler <- function(sig_genes, de_table = NULL, orgdb = "org.Dm.eg.db",
                                   orgdb_from = NULL, orgdb_to = "ENTREZID",
                                   go_level = 3, pcutoff = 0.05,
                                   qcutoff = 0.1, fc_column = "logFC",
                                   second_fc_column = "limma_logfc",
                                   updown = "up", permutations = 1000, min_groupsize = 5,
                                   kegg_prefix = NULL, kegg_organism = NULL, do_gsea = TRUE,
                                   categories = 12, excel = NULL, do_david = FALSE,
                                   david_id = "ENTREZ_GENE_ID",
                                   david_user = "unknown@unknown.org") {
  tt <- sm(requireNamespace(package = "clusterProfiler", quietly = TRUE))
  tt <- sm(requireNamespace(package = "DOSE", quietly = TRUE))
  org <- NULL

  ## Start off by figuring out what was given, an OrgDb or the name of one.
  if (class(orgdb)[[1]] == "OrgDb") {
    org <- orgdb
  } else if ("character" %in% class(orgdb)) {
    tt <- sm(requireNamespace(orgdb))
    org <- loadNamespace(orgdb) ## put the orgDb instance into an environment
    org <- org[[orgdb]] ## Then extract it
  } else {
    stop("Need either the name of an orgdb package or the orgdb itself.")
  }

  ## It is likely that we will need to query multiple different keys from the OrgDb
  ## So gather them now for later reference.
  mapper_keys <- AnnotationDbi::keytypes(org)
  ## If we must, we can extract the set of all genes from the orgdb
  ## However, this means we may not do a GSEA analysis
  universe_genes <- AnnotationDbi::keys(org)
  if (is.null(de_table)) {
    do_gsea <- FALSE
    all_genenames <- universe_genes
  } else {
    all_genenames <- rownames(de_table)
  }
  sig_genenames <- rownames(sig_genes)
  orgdb_to <- toupper(orgdb_to)

  de_table_namedf <- NULL
  sig_genes_namedf <- NULL
  test_genes_df <- NULL
  test_sig_df <- NULL
  num_sig <- 0
  num_hits <- 0
  orgdb_sig_from <- orgdb_from
  if (is.null(orgdb_from)) {
    message("Testing available OrgDb keytypes for the best mapping to entrez.")
    for (k in mapper_keys) {
      test_genes_df <- sm(try(clusterProfiler::bitr(all_genenames, fromType = k,
                                                    toType = orgdb_to, OrgDb = org), silent = TRUE))
      test_sig_df <- sm(try(clusterProfiler::bitr(sig_genenames, fromType = k,
                                                  toType = orgdb_to, OrgDb = org), silent = TRUE))
      if (class(test_genes_df) == "try-error") {
        test_genes_df <- data.frame()
      }
      if (class(test_sig_df) == "try-error") {
        test_sig_df <- data.frame()
      }
      test_num_hits <- nrow(test_genes_df)
      if (test_num_hits > num_hits) {
        orgdb_from <- k
        num_hits <- test_num_hits
        de_table_namedf <- test_genes_df
      }
      test_sig_hits <- nrow(test_sig_df)
      if (test_sig_hits > num_sig) {
        orgdb_sig_from <- k
        num_sig <- test_sig_hits
        sig_genes_namedf <- test_sig_df
      }
    }
    message("Chose keytype: ", orgdb_from, " for all genes because it had ", num_hits,
            " out of ", length(all_genenames), " genes.")
    message("Chose keytype: ", orgdb_sig_from, " for sig genes because it had ", num_sig,
            " out of ", length(sig_genenames), " genes.")
  } else { ## If we do have a column for the OrgDB
    de_table_namedf <- sm(try(clusterProfiler::bitr(all_genenames, fromType = orgdb_from,
                                                    toType = orgdb_to, OrgDb = org), silent = TRUE))
    sig_genes_namedf <- sm(try(clusterProfiler::bitr(sig_genenames, fromType = orgdb_from,
                                                     toType = orgdb_to, OrgDb = org), silent = TRUE))
  }

  if (is.null(sig_genes[[fc_column]]) & is.null(sig_genes[[second_fc_column]])) {
    stop("The fold change column provided no genes, try another column in the data set.")
  } else if (is.null(sig_genes[[fc_column]])) {
    fc_column <- second_fc_column
  }

  gsea_fc_column <- fc_column
  if (is.null(de_table[[gsea_fc_column]]) & is.null(de_table[[second_fc_column]])) {
    message("Unable to find the fold-change column in the de table, not doing gsea.")
    do_gsea <- FALSE
  } else if (is.null(de_table[[gsea_fc_column]])) {
    gsea_fc_column <- second_fc_column
  }

  ## Acquire the set of IDs against which all queries need to be made
  universe_to <- AnnotationDbi::keys(org, keytype = orgdb_to)
  ## And the set of similar IDs mapped against the significance table.
  all_gene_list <- de_table_namedf[[orgdb_to]]
  all_gene_drop <- !is.na(all_gene_list)
  sig_gene_list <- sig_genes_namedf[[orgdb_to]]
  sig_gene_drop <- !is.na(sig_gene_list)
  sig_gene_list <- sig_gene_list[sig_gene_drop]
  if (is.null(sig_gene_list)) {
    stop("No genes were found between the significant genes and the universe.")
  }

  ## Now we have a universe of geneIDs and significant IDs
  ## Let us perform some analyses...
  message("Calculating GO groups.")
  ggo_mf <- ggo_bp <- ggo_cc <- NULL
  ggo_mf <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "MF", level = go_level))
  ggo_bp <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "BP", level = go_level))
  ggo_cc <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "CC", level = go_level))

  group_go <- list(
    "MF" = as.data.frame(ggo_mf, stringsAsFactors = FALSE),
    "BP" = as.data.frame(ggo_bp, stringsAsFactors = FALSE),
    "CC" = as.data.frame(ggo_cc, stringsAsFactors = FALSE))
  message("Found ", nrow(group_go[["MF"]]),
          " MF, ", nrow(group_go[["BP"]]),
          " BP, and ", nrow(group_go[["CC"]]), " CC hits.")

  message("Calculating enriched GO groups.")
  enrich_results <- list(
    "all_mf" = NULL,
    "sig_mf" = NULL,
    "all_bp" = NULL,
    "sig_bp" = NULL,
    "all_cc" = NULL,
    "sig_cc" = NULL)
  ego_all_mf <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "MF", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = 1.0)
  ego_sig_mf <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "MF", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = pcutoff)
  ego_all_bp <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "BP", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = 1.0)
  ego_sig_bp <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "BP", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = pcutoff)
  ego_all_cc <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "CC", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = 1.0)
  ego_sig_cc <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "CC", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = "BH",
                                          pvalueCutoff = pcutoff)
  enrich_go <- list(
    "MF_all" = as.data.frame(ego_all_mf, stringsAsFactors = FALSE),
    "MF_sig" = as.data.frame(ego_sig_mf, stringsAsFactors = FALSE),
    "BP_all" = as.data.frame(ego_all_bp, stringsAsFactors = FALSE),
    "BP_sig" = as.data.frame(ego_sig_bp, stringsAsFactors = FALSE),
    "CC_all" = as.data.frame(ego_all_cc, stringsAsFactors = FALSE),
    "CC_sig" = as.data.frame(ego_sig_cc, stringsAsFactors = FALSE))
  message("Found ", nrow(enrich_go[["MF_sig"]]),
          " MF, ", nrow(enrich_go[["BP_sig"]]),
          " BP, and ", nrow(enrich_go[["CC_sig"]]), " CC enriched hits.")

  gse_go <- list()
  de_table_merged <- NULL
  if (isTRUE(do_gsea)) {
    ## Why did I do this? ## Ahh for the GSE analyses, they want ordered gene IDs
    ## Add the entrezIDs to the end
    de_table_merged <- merge(de_table, de_table_namedf, by.x = "row.names", by.y = 1)
    if (updown == "up") {
      new_order <- order(de_table_merged[[gsea_fc_column]], decreasing = TRUE)
      de_table_merged <- de_table_merged[new_order, ]
    } else {
      new_order <- order(de_table_merged[[gsea_fc_column]], decreasing = FALSE)
      de_table_merged <- de_table_merged[new_order, ]
    }

    ## Hmm this is odd, in the previous calls, I used orgdb_to, but in this set
    ## I am using orgdb_from...
    message("Performing GSE analyses of gene lists (this is slow).")
    genelist <- as.vector(de_table_merged[[gsea_fc_column]])
    names(genelist) <- de_table_merged[[orgdb_to]]
    ## 2020 04: Adding a pvalue cutoff argument causes an error, I do not know why.
    gse <- sm(clusterProfiler::gseGO(geneList = genelist, OrgDb = org,
                                     ont = "ALL", nPerm = permutations,
                                     minGSSize = min_groupsize))
    gse_go <- as.data.frame(gse)
    message("Found ", nrow(gse_go), " enriched hits.")
  } else {
    genelist <- as.vector(sig_genes[[fc_column]])
    names(genelist) <- rownames(sig_genes)
  }

  ## Now extract the kegg organism/gene IDs.
  if (is.null(kegg_organism)) {
    org_meta <- AnnotationDbi::metadata(org)
    org_row <- org_meta[["name"]] == "ORGANISM"
    organism <- org_meta[org_row, "value"]

    ## Only grab the first of potentially multiple outputs.
    kegg_organism <- get_kegg_orgn(species = organism)
    do_kegg <- TRUE
    if (length(kegg_organism) > 0) {
      kegg_organism <- kegg_organism[[1]]
    } else {
      do_kegg <- FALSE
      kegg_organism <- NULL
    }
  }

  all_kegg <- enrich_kegg <- NULL
  if (isTRUE(do_kegg)) {
    kegg_universe <- KEGGREST::keggConv(kegg_organism, "ncbi-geneid")
    kegg_sig_names <- glue("ncbi-geneid:{sig_gene_list}")
    kegg_sig_intersect <- kegg_sig_names %in% names(kegg_universe)
    message("Found ", sum(kegg_sig_intersect),
            " matches between the significant gene list and kegg universe.")

    if (sum(kegg_sig_intersect) > 0) {
      all_names <- names(kegg_universe)
      small_universe <- kegg_universe[intersect(kegg_sig_names, names(kegg_universe))]
      kegg_sig_ids <- unique(as.character(small_universe))
      ##kegg_sig_ids <- unique(as.character(kegg_universe[kegg_sig_intersect]))
      kegg_sig_ids <- gsub(pattern = glue("{kegg_organism}:"),
                           replacement = "", x = kegg_sig_ids)

      message("Performing KEGG analyses.")
      all_kegg <- clusterProfiler::enrichKEGG(kegg_sig_ids, organism = kegg_organism,
                                              keyType = "kegg",
                                              pvalueCutoff = 1.0)
      enrich_kegg <- sm(clusterProfiler::enrichKEGG(kegg_sig_ids, organism = kegg_organism,
                                                    keyType = "kegg",
                                                    pvalueCutoff = pcutoff))
    }
  }

  if (is.null(all_kegg)) {
    do_gsea <- FALSE
  }

  gse_all_kegg <- NULL
  gse_sig_kegg <- NULL
  gse_all_mkegg <- NULL
  gse_sig_mkegg <- NULL
  if (isTRUE(do_gsea)) {
    lastcol <- ncol(de_table_merged)
    kegg_genelist <- as.vector(de_table_merged[[fc_column]])
    names(kegg_genelist) <- de_table_merged[[lastcol]]

    kegg_all_names <- glue("ncbi-geneid:{names(kegg_genelist)}")
    kegg_all_intersect <- kegg_all_names %in% names(kegg_universe)
    message("Found ", sum(kegg_all_intersect),
            " matches between the gene list and kegg universe.")
    all_names <- names(kegg_universe)
    large_universe <- kegg_universe[intersect(kegg_all_names, names(kegg_universe))]
    kegg_all_ids <- unique(as.character(large_universe))
    kegg_all_ids <- gsub(pattern = glue("{kegg_organism}:"), replacement = "", x = kegg_all_ids)
    names(kegg_genelist) <- kegg_all_ids

    internal <- FALSE
    gse_all_kegg <- sm(
      clusterProfiler::gseKEGG(geneList = kegg_genelist, organism = kegg_organism,
                               nPerm = permutations, minGSSize = min_groupsize,
                               pvalueCutoff = 1.0, use_internal_data = internal))
    gse_sig_kegg <- sm(
      clusterProfiler::gseKEGG(geneList = kegg_genelist, organism = kegg_organism,
                               nPerm = permutations, minGSSize = min_groupsize,
                               pvalueCutoff = pcutoff, use_internal_data = internal))
  }

  kegg_data <- list(
    "kegg_all" = as.data.frame(all_kegg, stringsAsFactors = FALSE),
    "kegg_sig" = as.data.frame(enrich_kegg, stringsAsFactors = FALSE),
    "kegg_gse_all" = as.data.frame(gse_all_kegg, stringsAsFactors = FALSE),
    "kegg_gse_sig" = as.data.frame(gse_sig_kegg, stringsAsFactors = FALSE))
  message("Found ", nrow(kegg_data[["kegg_sig"]]), " KEGG enriched hits.")

  david_data <- NULL
  tt <- sm(please_install("RDAVIDWebService"))
  if (isTRUE(do_david)) {
    message("Attempting DAVID search.")
    david_search <- try(clusterProfiler::enrichDAVID(
                                           gene = sig_gene_list,
                                           minGSSize = min_groupsize,
                                           idType = david_id,
                                           david.user = david_user), silent = TRUE)
    if (class(david_search)[[1]] == "try-error") {
      david_data <- NULL
    } else {
      david_data <- as.data.frame(david_search, stringsAsFactors = FALSE)
    }
    message("Found ", nrow(david_data), " DAVID hits.")
  }

  ##testing <- clusterProfiler::enrichPathway(sig_gene_list)

  message("Plotting results.")
  map_sig_mf <- try(clusterProfiler::emapplot(ego_sig_mf), silent = TRUE)
  map_sig_bp <- try(clusterProfiler::emapplot(ego_sig_bp), silent = TRUE)
  map_sig_cc <- try(clusterProfiler::emapplot(ego_sig_cc), silent = TRUE)
  net_sig_mf <- try(
    clusterProfiler::cnetplot(ego_sig_mf, categorySize = "pvalue",
                              foldChange = genelist), silent = TRUE)
  net_sig_bp <- try(
    clusterProfiler::cnetplot(ego_sig_bp, categorySize = "pvalue",
                              foldChange = genelist), silent = TRUE)
  net_sig_cc <- try(
    clusterProfiler::cnetplot(ego_sig_cc, categorySize = "pvalue",
                              foldChange = genelist), silent = TRUE)

  tree_sig_mf <- tree_sig_bp <- tree_sig_cc <- NULL
  tree_mf <- sm(try(clusterProfiler::plotGOgraph(ego_sig_mf), silent = TRUE))
  if (class(tree_mf)[[1]] != "try-error") {
    tree_sig_mf <- recordPlot()
  }
  tree_bp <- sm(try(clusterProfiler::plotGOgraph(ego_sig_bp), silent = TRUE))
  if (class(tree_bp)[[1]] != "try-error") {
    tree_sig_bp <- recordPlot()
  }
  tree_cc <- sm(try(clusterProfiler::plotGOgraph(ego_sig_cc), silent = TRUE))
  if (class(tree_cc)[[1]] != "try-error") {
    tree_sig_cc <- recordPlot()
  }

  pvalue_plotlist <- list(
    ## I want to split the following list, but I am not sure which belong here.
  )
  ggo_mf_bar <- try(barplot(ggo_mf, drop = TRUE,
                            showCategory = categories), silent = TRUE)
  if (class(ggo_mf_bar)[[1]] == "try-error") {
    ggo_mf_bar <- NULL
  }
  ggo_bp_bar <- try(barplot(ggo_bp, drop = TRUE,
                            showCategory = categories), silent = TRUE)
  if (class(ggo_bp_bar)[[1]] == "try-error") {
    ggo_bp_bar <- NULL
  }
  ggo_cc_bar <- try(barplot(ggo_cc, drop = TRUE,
                            showCategory = categories), silent = TRUE)
  if (class(ggo_cc_bar)[[1]] == "try-error") {
    ggo_cc_bar <- NULL
  }
  ego_all_mf_bar <- try(barplot(ego_all_mf,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_all_mf)[[1]] == "try-error") {
    ego_all_mf <- NULL
  }
  ego_all_bp_bar <- try(barplot(ego_all_bp,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_all_bp)[[1]] == "try-error") {
    ego_all_bp <- NULL
  }
  ego_all_cc_bar <- try(barplot(ego_all_cc,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_all_cc)[[1]] == "try-error") {
    ego_all_cc <- NULL
  }
  ego_sig_mf_bar <- try(barplot(ego_sig_mf,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_sig_mf)[[1]] == "try-error") {
    ego_sig_mf <- NULL
  }
  ego_sig_bp_bar <- try(barplot(ego_sig_bp,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_sig_bp)[[1]] == "try-error") {
    ego_sig_bp <- NULL
  }
  ego_sig_cc_bar <- try(barplot(ego_sig_cc,
                                showCategory = categories, drop = TRUE), silent = TRUE)
  if (class(ego_sig_cc)[[1]] == "try-error") {
    ego_sig_cc <- NULL
  }
  dot_all_mf <- sm(try(clusterProfiler::dotplot(ego_all_mf), silent = TRUE))
  if (class(dot_all_mf)[[1]] == "try-error") {
    dot_all_mf <- NULL
  }
  dot_all_bp <- try(clusterProfiler::dotplot(ego_all_bp), silent = TRUE)
  if (class(dot_all_bp)[[1]] == "try-error") {
    dot_all_bp <- NULL
  }
  dot_all_cc <- try(clusterProfiler::dotplot(ego_all_cc), silent = TRUE)
  if (class(dot_all_cc)[[1]] == "try-error") {
    dot_all_cc <- NULL
  }
  dot_sig_mf <- try(clusterProfiler::dotplot(ego_sig_mf), silent = TRUE)
  if (class(dot_sig_mf)[[1]] == "try-error") {
    dot_sig_mf <- NULL
  }
  dot_sig_bp <- try(clusterProfiler::dotplot(ego_sig_bp), silent = TRUE)
  if (class(dot_sig_bp)[[1]] == "try-error") {
    dot_sig_bp <- NULL
  }
  dot_sig_cc <- try(clusterProfiler::dotplot(ego_sig_cc), silent = TRUE)
  if (class(dot_sig_cc)[[1]] == "try-error") {
    dot_sig_cc <- NULL
  }

  plotlist <- list(
    "ggo_mf_bar" = ggo_mf_bar,
    "ggo_bp_bar" = ggo_bp_bar,
    "ggo_cc_bar" = ggo_cc_bar,
    "ego_all_mf" = ego_all_mf_bar,
    "ego_all_bp" = ego_all_bp_bar,
    "ego_all_cc" = ego_all_cc_bar,
    "ego_sig_mf" = ego_sig_mf_bar,
    "ego_sig_bp" = ego_sig_bp_bar,
    "ego_sig_cc" = ego_sig_cc_bar,
    "dot_all_mf" = dot_all_mf,
    "dot_all_bp" = dot_all_bp,
    "dot_all_cc" = dot_all_cc,
    "dot_sig_mf" = dot_sig_mf,
    "dot_sig_bp" = dot_sig_bp,
    "dot_sig_cc" = dot_sig_cc,
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
    "all_mappings" = de_table_namedf,
    "sig_mappings" = sig_genes_namedf,
    "group_go" = group_go,
    "enrich_go" = enrich_go,
    "gse_go" = gse_go,
    "kegg_data" = kegg_data,
    "david_data" = david_data,
    "plots" = plotlist,
    "pvalue_plots" = plotlist)
  class(retlist) <- c("clusterprofiler_result", "list")
  if (!is.null(excel)) {
    message("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_cp_data(retlist, excel = excel)))
    if (class(excel_ret) == "try-error") {
      message("Writing the data to excel failed.")
    }
  }
  return(retlist)
}

#' Set up appropriate option sets for clusterProfiler
#'
#' This hard-sets some defaults for orgdb/kegg databases when using
#' clusterProfiler.
#'
#' @param species  Currently it only works for humans and fruit flies.
cp_options <- function(species) {
  if (species == "dmelanogaster") {
    options <- list(
      orgdb = "org.Dm.eg.db",
      orgdb_from = "FLYBASE",
      orgdb_to = c("ENSEMBL", "SYMBOL", "ENTREZID"),
      kegg_prefix = "Dmel_",
      kegg_organism = "dme",
      kegg_id_column = "FLYBASECG")
  } else if (species == "hsapiens") {
    options <- list(
      orgdb = "org.Hs.eg.db",
      orgdb_from = "ENSEMBL",
      orgdb_to = c("ENSEMBL", "SYMBOL", "ENTREZID"),
      kegg_prefix = "Hsa_",
      kegg_organism = "hsa",
      kegg_id_column = "")
  }
  return(options)
}

#' Generic enrichment using clusterProfiler.
#'
#' culsterProfiler::enricher provides a quick and easy enrichment analysis given
#' a set of siginficant' genes and a data frame which connects each gene to a
#' category.
#'
#' @param sig_genes Set of 'significant' genes as a table.
#' @param de_table All genes from the original analysis.
#' @param go_db Dataframe of GO->ID matching the gene names of sig_genes to GO
#'  categories.
#' @return Table of 'enriched' categories.
simple_cp_enricher <- function(sig_genes, de_table, go_db = NULL) {
  all_genenames <- rownames(de_table)
  sig_genenames <- rownames(sig_genes)
  enriched <- clusterProfiler::enricher(sig_genenames, TERM2GENE = go_db)
  retlist <- list(
    "enriched" = as.data.frame(enriched, stringsAsFactors = FALSE))
  return(retlist)
}

## EOF
