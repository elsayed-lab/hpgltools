## Time-stamp: <Mon Jan  4 12:47:10 2016 Ashton Trey Belew (abelew@gmail.com)>

#' A simplification function for gostats, in the same vein as those written for clusterProfiler, goseq, and topGO.
#'
#' GOstats has a couple interesting peculiarities:  Chief among them: the gene IDs must be integers.
#' As a result, I am going to have this function take a gff file in order to get the go ids and
#' gene ids on the same page.
#'
#' @param gff The annotation information for this genome
#' @param de_genes The set of differentially expressed genes in the limma format as before
#' @param goids The set of GOids, as before in the format ID/GO
#'
#' @return dunno yet
#' @seealso \code{\link{GOstats}}
#' @export
simple_gostats = function(de_genes, gff, goids, universe_merge="ID", second_merge_try="locus_tag", organism="fun", pcutoff=0.10, direction="over", conditional=FALSE, categorysize=NULL, gff_type="CDS") {
    ## The import(gff) is being used for this primarily because it uses integers for the rownames and because it (should)
    ## contain every gene in the 'universe' used by GOstats, as much it ought to be pretty much perfect.
    annotation = hpgltools:::gff2df(gff)
    ## This is similar to logic in ontology_goseq and is similarly problematic.
    annotation = annotation[annotation$type == gff_type, ]
    annotation[, universe_merge] = make.names(annotation[, universe_merge], unique=TRUE)
    if (universe_merge %in% names(annotation)) {
        universe = annotation[,c(universe_merge, "width")]
    } else if (second_merge_try %in% names(annotation)) {
        universe = annotation[,c(second_merge_try, "width")]
    } else if ("transcript_name" %in% names(annotation)) {
        universe = annotation[,c("transcript_name", "width")]
    } else {
        stop("Unable to cross reference the annotations into a universe for background checks.")
    }
    ## This section is a little odd
    ## The goal is to collect a consistent set of numeric gene IDs
    ## In addition, one must cross reference those IDs consistently with the universe of all genes.
    ## Thus in a few linues I will be doing a merge of all genes against the de_genes and another merge
    ## of the gene<->go mappings, finally extracting the portions of the resulting dataframe into a format suitable for
    ## casting as a GOFrame/GOAllFrame
    colnames(universe) = c("geneid","width")
    universe$id = rownames(universe)
    universe = universe[complete.cases(universe),]

    if (is.null(de_genes$ID)) {
        de_genes$ID = rownames(de_genes)
    }
    universe_cross_de = merge(universe, de_genes, by.x="geneid", by.y="ID")
    degenes_ids = universe_cross_de$id
    universe_ids = universe$id
    gostats_go = merge(universe, goids, by.x="geneid", by.y="ID")
    if (nrow(gostats_go) == 0) {
        stop("The merging of the universe vs. goids failed.")
    }
    gostats_go$frame.Evidence = "TAS"
    colnames(gostats_go) = c("sysName","width", "frame.gene_id", "frame.go_id", "frame.Evidence")
    gostats_go = gostats_go[,c("frame.go_id", "frame.Evidence", "frame.gene_id")]
    gostats_frame = GOFrame(gostats_go, organism=organism)
    gostats_all = GOAllFrame(gostats_frame)
    require.auto("GSEABase", verbose=FALSE)
    message("Creating the gene set collection.")
    gsc = GeneSetCollection(gostats_all, setType=GOCollection())

    mf_over = bp_over = cc_over = NULL
    mf_under = bp_under = cc_under = NULL
    message("Performing gene set enrichment of molecular function over representation.")
    mf_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="MF", pvalueCutoff=pcutoff,
                                   conditional=conditional, testDirection="over")
    mf_over = hyperGTest(mf_params)
    message(paste0("Found ", nrow(GOstats::summary(mf_over)), " over represented molecular function categories."))
    message("Performing gene set enrichment of biological process over representation.")
    bp_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="BP", pvalueCutoff=pcutoff,
                                   conditional=FALSE, testDirection="over")
    bp_over = hyperGTest(bp_params)
    message(paste0("Found ", nrow(GOstats::summary(bp_over)), " over represented biological process categories."))
    message("Performing gene set enrichment of cellular component over representation.")
    cc_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="CC", pvalueCutoff=pcutoff,
                                   conditional=FALSE, testDirection="over")
    cc_over = hyperGTest(cc_params)
    message(paste0("Found ", nrow(GOstats::summary(cc_over)), " over represented cellular component categories."))
    message("Performing gene set enrichment of molecular function under representation.")
    mf_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="MF", pvalueCutoff=pcutoff,
                                   conditional=conditional, testDirection="under")
    mf_under = hyperGTest(mf_params)
    message(paste0("Found ", nrow(GOstats::summary(mf_under)), " under represented molecular function categories."))
    message("Performing gene set enrichment of biological process under representation.")
    bp_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="BP", pvalueCutoff=pcutoff,
                                   conditional=FALSE, testDirection="under")
    bp_under = hyperGTest(bp_params)
    message(paste0("Found ", nrow(GOstats::summary(bp_under)), " under represented biological process categories."))
    message("Performing gene set enrichment of cellular component under representation.")
    cc_params = GSEAGOHyperGParams(name=paste("GSEA of ", organism, sep=""), geneSetCollection=gsc,
                                   geneIds=degenes_ids, universeGeneIds=universe_ids,
                                   ontology="CC", pvalueCutoff=pcutoff,
                                   conditional=FALSE, testDirection="under")
    cc_under = hyperGTest(cc_params)
    message(paste0("Found ", nrow(GOstats::summary(cc_under)), " under represented cellular component categories."))
    mf_over_table = bp_over_table = cc_over_table = NULL
    mf_under_table = bp_under_table = cc_under_table = NULL
    mf_over_table = GOstats::summary(mf_over, pvalue=1.0, htmlLinks=TRUE)
    bp_over_table = GOstats::summary(bp_over, pvalue=1.0, htmlLinks=TRUE)
    cc_over_table = GOstats::summary(cc_over, pvalue=1.0, htmlLinks=TRUE)
    mf_under_table = GOstats::summary(mf_under, pvalue=1.0, htmlLinks=TRUE)
    bp_under_table = GOstats::summary(bp_under, pvalue=1.0, htmlLinks=TRUE)
    cc_under_table = GOstats::summary(cc_under, pvalue=1.0, htmlLinks=TRUE)
    if (!is.null(dim(mf_over_table))) {
        mf_over_table$qvalue = qvalue(mf_over_table$Pvalue)$qvalues
    }
    if (!is.null(dim(bp_over_table))) {
        bp_over_table$qvalue = qvalue(bp_over_table$Pvalue)$qvalues
    }
    if (!is.null(dim(cc_over_table))) {
        cc_over_table$qvalue = qvalue(cc_over_table$Pvalue)$qvalues
    }
    if (!is.null(dim(mf_under_table))) {
        mf_under_table$qvalue = qvalue(mf_under_table$Pvalue)$qvalues
    }
    if (!is.null(dim(bp_under_table))) {
        bp_under_table$qvalue = qvalue(bp_under_table$Pvalue)$qvalues
    }
    if (!is.null(dim(cc_under_table))) {
        cc_under_table$qvalue = qvalue(cc_under_table$Pvalue)$qvalues
    }

    if (is.null(categorysize)) {
        mf_over_sig = GOstats::summary(mf_over)
        bp_over_sig = GOstats::summary(bp_over)
        cc_over_sig = GOstats::summary(cc_over)
        mf_under_sig = GOstats::summary(mf_under)
        bp_under_sig = GOstats::summary(bp_under)
        cc_under_sig = GOstats::summary(cc_under)
    } else {
        mf_over_sig = GOstats::summary(mf_over, categorySize=categorysize)
        bp_over_sig = GOstats::summary(bp_over, categorySize=categorysize)
        cc_over_sig = GOstats::summary(cc_over, categorySize=categorysize)
        mf_under_sig = GOstats::summary(mf_under, categorySize=categorysize)
        bp_under_sig = GOstats::summary(bp_under, categorySize=categorysize)
        cc_under_sig = GOstats::summary(cc_under, categorySize=categorysize)
    }
    if (!is.null(dim(mf_over_sig))) {
        mf_over_sig$definition = try(godef(mf_over_sig$GOMFID), silent=TRUE)
    } else {
        mf_over_sig = NULL
    }
    if (!is.null(dim(bp_over_sig))) {
        bp_over_sig$definition = try(godef(bp_over_sig$GOBPID), silent=TRUE)
    } else {
        bp_over_sig = NULL
    }
    if (!is.null(dim(cc_over_sig))) {
        cc_over_sig$definition = try(godef(cc_over_sig$GOCCID), silent=TRUE)
    } else {
        bp_over_sig = NULL
    }
    if (!is.null(dim(mf_under_sig))) {
        mf_under_sig$definition = try(godef(mf_under_sig$GOMFID), silent=TRUE)
    } else {
        mf_under_sig = NULL
    }
    if (!is.null(dim(bp_under_sig))) {
        bp_under_sig$definition = try(godef(bp_under_sig$GOBPID), silent=TRUE)
    } else {
        bp_under_sig = NULL
    }
    if (!is.null(dim(cc_under_sig))) {
        cc_under_sig$definition = try(godef(cc_under_sig$GOCCID), silent=TRUE)
    } else {
        bp_under_sig = NULL
    }

    gostats_p_mf_over = try(hpgl_histogram(mf_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_mf_under = try(hpgl_histogram(mf_under_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_bp_over = try(hpgl_histogram(bp_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_bp_under = try(hpgl_histogram(bp_under_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_cc_over = try(hpgl_histogram(cc_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_cc_under = try(hpgl_histogram(cc_under_table$Pvalue, bins=20), silent=TRUE)

##    mf_over_table = as.data.frame(mf_over_table)
##    bp_over_table = as.data.frame(bp_over_table)
#    cc_over_table = as.data.frame(cc_over_table)
#    mf_under_table = as.data.frame(mf_under_table)
#    bp_under_table = as.data.frame(bp_under_table)
#    cc_under_table = as.data.frame(cc_under_table)

    ret_list = list(mf_over_all=mf_over_table, bp_over_all=bp_over_table, cc_over_all=cc_over_table,
        mf_under_all=mf_under_table, bp_under_all=bp_under_table, cc_under_all=cc_under_table,
        mf_over_enriched=mf_over_sig, bp_over_enriched=bp_over_sig, cc_over_enriched=cc_over_sig,
        mf_under_enriched=mf_under_sig, bp_under_enriched=bp_under_sig, cc_under_enriched=cc_under_sig,
        gostats_mfp_over=gostats_p_mf_over, gostats_bpp_over=gostats_p_bp_over, gostats_ccp_over=gostats_p_cc_over,
        gostats_mfp_under=gostats_p_mf_under, gostats_bpp_under=gostats_p_bp_under, gostats_ccp_under=gostats_p_cc_under)

    pvalue_plots = try(hpgltools:::gostats_pval_plots(ret_list))
    ret_list$pvalue_plots = pvalue_plots

    return(ret_list)
}

## Take gostats data and print it on a tree as topGO does
#' Make fun trees a la topgo from goseq data.
#'
#' @param de_genes some differentially expressed genes
#' @param mf_over/bp_over/cc_over/mf_under/bp_under/cc_under over/under expression data
#' @param goid_map a mapping of IDs to GO in the Ramigo expected format
#' @param score_limit maximum score to include as 'significant'
#' @param goids_df a dataframe of available goids (used to generate goid_map)
#' @param overwrite overwrite the goid_map?
#' @param selector a function to choose differentially expressed genes in the data
#' @param pval_column a column in the data to be used to extract pvalue scores
#'
#' @return plots! Trees! oh my!
#' @seealso \code{\link{topGO}}
#' @export
gostats_trees = function(de_genes, mf_over, bp_over, cc_over, mf_under, bp_under, cc_under, goid_map="reference/go/id2go.map", score_limit=0.01, goids_df=NULL, overwrite=FALSE, selector="topDiffGenes", pval_column="adj.P.Val") {
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO = topGO::readMappings(file=goid_map)
    annotated_genes = names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID = make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes = factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) = annotated_genes
    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata = new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=annFun.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=annFun.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=annFun.gene2GO, gene2GO=geneID2GO)
    } else {
        pvals = as.vector(de_genes[[pval_column]])
        names(pvals) = rownames(de_genes)
        mf_GOdata = new("topGOdata", description="MF", ontology="MF", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", description="BP", ontology="BP", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", description="CC", ontology="CC", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    }
    mf_over_enriched_ids = mf_over$GOMFID
    bp_over_enriched_ids = bp_over$GOBPID
    cc_over_enriched_ids = cc_over$GOCCID
    mf_under_enriched_ids = mf_under$GOMFID
    bp_under_enriched_ids = bp_under$GOBPID
    cc_under_enriched_ids = cc_under$GOCCID
    mf_over_enriched_scores = mf_over$Pvalue
    names(mf_over_enriched_scores) = mf_over_enriched_ids
    bp_over_enriched_scores = bp_over$Pvalue
    names(bp_over_enriched_scores) = bp_over_enriched_ids
    cc_over_enriched_scores = cc_over$Pvalue
    names(cc_over_enriched_scores) = cc_over_enriched_ids
    mf_under_enriched_scores = mf_under$Pvalue
    names(mf_under_enriched_scores) = mf_under_enriched_ids
    bp_under_enriched_scores = bp_under$Pvalue
    names(bp_under_enriched_scores) = bp_under_enriched_ids
    cc_under_enriched_scores = cc_under$Pvalue
    names(cc_under_enriched_scores) = cc_under_enriched_ids

    mf_avail_nodes = as.list(mf_GOdata@graph@nodes)
    names(mf_avail_nodes) = mf_GOdata@graph@nodes
    mf_over_nodes = mf_over_enriched_scores[names(mf_over_enriched_scores) %in% names(mf_avail_nodes)]
    mf_under_nodes = mf_under_enriched_scores[names(mf_under_enriched_scores) %in% names(mf_avail_nodes)]
    mf_over_included = length(which(mf_over_nodes <= score_limit))
    mf_under_included = length(which(mf_under_nodes <= score_limit))
    mf_over_tree_data = try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_over_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=mf_over_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    mf_under_tree_data = try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_under_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=mf_under_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(mf_over_tree_data) == 'try-error') {
        message("There was an error generating the over MF tree.")
        mf_over_tree = NULL
    } else {
        mf_over_tree = recordPlot()
    }
    if (class(mf_under_tree_data) == 'try-error') {
        message("There was an error generating the under MF tree.")
        mf_under_tree = NULL
    } else {
        mf_under_tree = recordPlot()
    }

    bp_avail_nodes = as.list(bp_GOdata@graph@nodes)
    names(bp_avail_nodes) = bp_GOdata@graph@nodes
    bp_over_nodes = bp_over_enriched_scores[names(bp_over_enriched_scores) %in% names(bp_avail_nodes)]
    bp_under_nodes = bp_under_enriched_scores[names(bp_under_enriched_scores) %in% names(bp_avail_nodes)]
    bp_over_included = length(which(bp_over_nodes <= score_limit))
    bp_under_included = length(which(bp_under_nodes <= score_limit))
    bp_over_tree_data = try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_over_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=bp_over_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    bp_under_tree_data = try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_under_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=bp_under_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(bp_over_tree_data) == 'try-error') {
        message("There was an error generating the over BP tree.")
        bp_over_tree = NULL
    } else {
        bp_over_tree = recordPlot()
    }
    if (class(bp_under_tree_data) == 'try-error') {
        message("There was an error generating the under BP tree.")
        bp_under_tree = NULL
    } else {
        bp_under_tree = recordPlot()
    }

    cc_avail_nodes = as.list(cc_GOdata@graph@nodes)
    names(cc_avail_nodes) = cc_GOdata@graph@nodes
    cc_over_nodes = cc_over_enriched_scores[names(cc_over_enriched_scores) %in% names(cc_avail_nodes)]
    cc_under_nodes = cc_under_enriched_scores[names(cc_under_enriched_scores) %in% names(cc_avail_nodes)]
    cc_over_included = length(which(cc_over_nodes <= score_limit))
    cc_under_included = length(which(cc_under_nodes <= score_limit))
    cc_over_tree_data = try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_over_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=cc_over_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    cc_under_tree_data = try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_under_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=cc_under_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(cc_over_tree_data) == 'try-error') {
        message("There was an error generating the over CC tree.")
        cc_over_tree = NULL
    } else {
        cc_over_tree = recordPlot()
    }
    if (class(cc_under_tree_data) == 'try-error') {
        message("There was an error generating the under CC tree.")
        cc_under_tree = NULL
    } else {
        cc_under_tree = recordPlot()
    }

    trees = list(
        MF_over=mf_over_tree, BP_over=bp_over_tree, CC_over=cc_over_tree,
        MF_overdata=mf_over_tree_data, BP_overdata=bp_over_tree_data, CC_overdata=cc_over_tree_data,
        MF_under=mf_under_tree, BP_under=bp_under_tree, CC_under=cc_under_tree,
        MF_underdata=mf_under_tree_data, BP_underdata=bp_under_tree_data, CC_underdata=cc_under_tree_data
    )
    return(trees)
}

#' Make a pvalue plot similar to that from clusterprofiler from gostats data
#'
#' clusterprofiler provides beautiful plots describing significantly overrepresented categories.
#' This function attempts to expand the repetoire of data available to them to include data from gostats.
#'
#' @param mf_over molecular function data overrepresented
#' @param bp_over biological process data overrepresented
#' @param cc_over cellular component data overrepresented
#' @param mf_under molecular function data underrepresented
#' @param bp_under biological process data underrepresented
#' @param cc_under cellular component data underrepresented
#' @param cutoff 0.1 what is the maximum pvalue allowed
#' @param wrapped_width how big to make the text so that it is legible
#' @param n 10 how many groups to include in the plot
#' @section
#' warning hey, the pval_plot function upon which this is based now has a bunch of new helpers now that I understand how the ontology trees work better, this should take advantage of that, but currently does not.
#'
#' @return plots!
#' @seealso \code{\link{clusterProfiler}} \code{\link{pval_plot}}
#' @export
gostats_pval_plots = function(gs_result, wrapped_width=20, cutoff=0.1, n=12, group_minsize=5) {
    ## TODO: replace the subset calls
    mf_over = gs_result$mf_over_enriched
    mf_under = gs_result$mf_under_enriched
    bp_over = gs_result$bp_over_enriched
    bp_under = gs_result$bp_under_enriched
    cc_over = gs_result$cc_over_enriched
    cc_under = gs_result$cc_under_enriched

    plotting_mf_over = mf_over
    mf_pval_plot_over = NULL
    if (is.null(mf_over)) {
        plotting_mf_over = NULL
    } else {
        plotting_mf_over$score = plotting_mf_over$ExpCount
        plotting_mf_over = subset(plotting_mf_over, Term != "NULL")
        plotting_mf_over = subset(plotting_mf_over, Pvalue <= cutoff)
        plotting_mf_over = subset(plotting_mf_over, Size >= group_minsize)
        plotting_mf_over = plotting_mf_over[order(plotting_mf_over$Pvalue),]
        plotting_mf_over = head(plotting_mf_over, n=n)
        plotting_mf_over = plotting_mf_over[,c("Term","Pvalue","score")]
        colnames(plotting_mf_over) = c("term","pvalue","score")
    }
    if (nrow(plotting_mf_over) > 0) {
        mf_pval_plot_over = pval_plot(plotting_mf_over, ontology="MF")
    }
    plotting_mf_under = mf_under
    mf_pval_plot_under = NULL
    if (is.null(mf_under)) {
        plotting_mf_under = NULL
    } else {
        plotting_mf_under$score = plotting_mf_under$ExpCount
        plotting_mf_under = subset(plotting_mf_under, Term != "NULL")
        plotting_mf_under = subset(plotting_mf_under, Pvalue <= cutoff)
        plotting_mf_under = subset(plotting_mf_under, Size >= group_minsize)
        plotting_mf_under = plotting_mf_under[order(plotting_mf_under$Pvalue),]
        plotting_mf_under = head(plotting_mf_under, n=n)
        plotting_mf_under = plotting_mf_under[,c("Term","Pvalue","score")]
        colnames(plotting_mf_under) = c("term","pvalue","score")
    }
    if (nrow(plotting_mf_under) > 0) {
        mf_pval_plot_under = pval_plot(plotting_mf_under, ontology="MF")
    }
    plotting_bp_over = bp_over
    bp_pval_plot_over = NULL
    if (is.null(bp_over)) {
        plotting_bp_over = NULL
    } else {
        plotting_bp_over$score = plotting_bp_over$ExpCount
        plotting_bp_over = subset(plotting_bp_over, Term != "NULL")
        plotting_bp_over = subset(plotting_bp_over, Pvalue <= 0.1)
        plotting_bp_over = subset(plotting_bp_over, Size > 10)
        plotting_bp_over = plotting_bp_over[order(plotting_bp_over$Pvalue),]
        plotting_bp_over = head(plotting_bp_over, n=n)
        plotting_bp_over = plotting_bp_over[,c("Term","Pvalue","score")]
        colnames(plotting_bp_over) = c("term","pvalue","score")
    }
    if (nrow(plotting_bp_over) > 0) {
        bp_pval_plot_over = pval_plot(plotting_bp_over, ontology="BP")
    }
    plotting_bp_under = bp_under
    bp_pval_plot_under = NULL
    if (is.null(bp_under)) {
        plotting_bp_under = NULL
    } else {
        plotting_bp_under$score = plotting_bp_under$ExpCount
        plotting_bp_under = subset(plotting_bp_under, Term != "NULL")
        plotting_bp_under = subset(plotting_bp_under, Pvalue <= 0.1)
        plotting_bp_under = subset(plotting_bp_under, Size > 10)
        plotting_bp_under = plotting_bp_under[order(plotting_bp_under$Pvalue),]
        plotting_bp_under = head(plotting_bp_under, n=n)
        plotting_bp_under = plotting_bp_under[,c("Term","Pvalue","score")]
        colnames(plotting_bp_under) = c("term","pvalue","score")
    }
    if (nrow(plotting_bp_under) > 0) {
        bp_pval_plot_under = pval_plot(plotting_bp_under, ontology="BP")
    }
    plotting_cc_over = cc_over
    cc_pval_plot_over = NULL
    if (is.null(cc_over)) {
        plotting_cc_over = NULL
    } else {
        plotting_cc_over$score = plotting_cc_over$ExpCount
        plotting_cc_over = subset(plotting_cc_over, Term != "NULL")
        plotting_cc_over = subset(plotting_cc_over, Pvalue <= 0.1)
        plotting_cc_over = subset(plotting_cc_over, Size > 10)
        plotting_cc_over = plotting_cc_over[order(plotting_cc_over$Pvalue),]
        plotting_cc_over = head(plotting_cc_over, n=n)
        plotting_cc_over = plotting_cc_over[,c("Term","Pvalue","score")]
        colnames(plotting_cc_over) = c("term","pvalue","score")
    }
    if (nrow(plotting_cc_over) > 0) {
        cc_pval_plot_over = pval_plot(plotting_cc_over, ontology="CC")
    }
    plotting_cc_under = cc_under
    cc_pval_plot_under = NULL
    if (is.null(cc_under)) {
        plotting_cc_under = NULL
    } else {
        plotting_cc_under$score = plotting_cc_under$ExpCount
        plotting_cc_under = subset(plotting_cc_under, Term != "NULL")
        plotting_cc_under = subset(plotting_cc_under, Pvalue <= 0.1)
        plotting_cc_under = subset(plotting_cc_under, Size > 10)
        plotting_cc_under = plotting_cc_under[order(plotting_cc_under$Pvalue),]
        plotting_cc_under = head(plotting_cc_under, n=n)
        plotting_cc_under = plotting_cc_under[,c("Term","Pvalue","score")]
        colnames(plotting_cc_under) = c("term","pvalue","score")
    }
    if (nrow(plotting_cc_under) > 0) {
        cc_pval_plot_under = pval_plot(plotting_cc_under, ontology="CC")
    }

    pval_plots = list(mfp_plot_over=mf_pval_plot_over, bpp_plot_over=bp_pval_plot_over, ccp_plot_over=cc_pval_plot_over,
        mf_subset_over=plotting_mf_over, bp_subset_over=plotting_bp_over, cc_subset_over=plotting_cc_over,
        mfp_plot_under=mf_pval_plot_under, bpp_plot_under=bp_pval_plot_under, ccp_plot_under=cc_pval_plot_under,
        mf_subset_under=plotting_mf_under, bp_subset_under=plotting_bp_under, cc_subset_under=plotting_cc_under)

    return(pval_plots)
}

## EOF
