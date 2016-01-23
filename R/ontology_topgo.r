## Time-stamp: <Sat Jan 23 15:30:27 2016 Ashton Trey Belew (abelew@gmail.com)>

#' simple_topgo()  Perform a simplified topgo analysis
#'
#' This will attempt to make it easier to run topgo on a set of genes.
#'
#' @param de_genes  a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param goid_map default='reference/go/id2go.map'  a file containing mappings of genes to goids in the format expected by topgo
#' @param goids_df default=NULL  a data frame of the goids which may be used to make the goid_map
#' @param pvals default=NULL  a set of pvalues in the DE data which may be used to improve the topgo results
#' @param limitby default='fisher'  test to index the results by
#' @param limit default=0.1  ontology pvalue to use as the lower limit
#' @param signodes default=100  I don't remember right now
#' @param sigforall default=TRUE  provide the significance for all nodes?
#' @param numchar default=300  character limit for the table of results
#' @param selector default='topDiffGenes'  a function name for choosing genes to include
#' @param overwrite default=FALSE  yeah I do not remember this one either
#' @param densities default=FALSE  the densities, yeah, the densities
#' @param pval_plots default=TRUE  include pvalue plots of the results a la clusterprofiler
#'
#' @return a big list including the various outputs from topgo
#' @export
simple_topgo <- function(de_genes, goid_map="reference/go/id2go.map", goids_df=NULL,
                         pvals=NULL, limitby="fisher", limit=0.1, signodes=100,
                         sigforall=TRUE, numchar=300, selector="topDiffGenes",
                         pval_column="adj.P.Val", overwrite=FALSE, densities=FALSE,
                         pval_plots=TRUE, ...) {
### Some neat ideas from the topGO documentation:
### geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")
### A variant of these operations make it possible to give topGO scores so that
### a larger array of tests may be performed
### x <- topDiffGenes(geneList)
### sum(x) ## the number of selected genes
### If we do something like above to give scores to all the 'DEgenes', then we set up the GOdata object like this:
### mf_GOdata = new("topGOdata", description="something", ontology="BP", allGenes = entire_geneList, geneSel=topDiffGenes, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=2)
    ## The following library invocation is in case it was unloaded for pathview
    require.auto("topGO")
    require.auto("Hmisc")
    gomap_info <- make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    message(paste0("simple_topgo(): Found ID->GO map: ", gomap_info))
    geneID2GO <- topGO::readMappings(file=goid_map)
    annotated_genes <- names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID <- make.names(rownames(de_genes), unique=TRUE)
    }
    ##    interesting_genes = factor(as.integer(annotated_genes %in% de_genes$ID))
    fisher_interesting_genes <- as.factor(as.integer(annotated_genes %in% de_genes$ID))
    names(fisher_interesting_genes) <- annotated_genes
    ks_interesting_genes <- as.integer(annotated_genes %nin% de_genes$ID)
    if (!is.null(de_genes$P.Value)) {
        ## I think this needs to include the entire gene universe, not only the set of x differentially expressed genes
        ## Making this an explicit as.vector(as.numeric()) because it turns out the values from DESeq are characters.
        pvals <- as.vector(as.numeric(de_genes[[pval_column]]))
        names(pvals) <- rownames(de_genes)
        for (p in 1:length(pvals)) {
            name <- names(pvals)[p]
            ks_interesting_genes[[name]] <- pvals[p]
        }
    }
    ks_interesting_genes <- as.vector(ks_interesting_genes)
    names(ks_interesting_genes) <- annotated_genes

    fisher_mf_GOdata <- new("topGOdata", ontology="MF", allGenes=fisher_interesting_genes,
                            annot=annFUN.gene2GO, gene2GO=geneID2GO)
    fisher_bp_GOdata <- new("topGOdata", ontology="BP", allGenes=fisher_interesting_genes,
                            annot=annFUN.gene2GO, gene2GO=geneID2GO)
    fisher_cc_GOdata <- new("topGOdata", ontology="CC", allGenes=fisher_interesting_genes,
                            annot=annFUN.gene2GO, gene2GO=geneID2GO)
    ks_mf_GOdata <- new("topGOdata", description="MF", ontology="MF", allGenes=ks_interesting_genes,
                        geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    ks_bp_GOdata <- new("topGOdata", description="BP", ontology="BP", allGenes=ks_interesting_genes,
                        geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    ks_cc_GOdata <- new("topGOdata", description="CC", ontology="CC", allGenes=ks_interesting_genes,
                        geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)

    test_stat <- new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
    mf_fisher_result <- topGO::getSigGroups(fisher_mf_GOdata, test_stat)
    bp_fisher_result <- topGO::getSigGroups(fisher_bp_GOdata, test_stat)
    cc_fisher_result <- topGO::getSigGroups(fisher_cc_GOdata, test_stat)
    test_stat <- new("classicScore", testStatistic=GOKSTest, name="KS tests")
    mf_ks_result <- topGO::getSigGroups(ks_mf_GOdata, test_stat)
    bp_ks_result <- topGO::getSigGroups(ks_bp_GOdata, test_stat)
    cc_ks_result <- topGO::getSigGroups(ks_cc_GOdata, test_stat)
    test_stat <- new("elimScore", testStatistic=GOKSTest, name="Fisher test", cutOff=0.01)
    mf_el_result <- topGO::getSigGroups(fisher_mf_GOdata, test_stat)
    bp_el_result <- topGO::getSigGroups(fisher_bp_GOdata, test_stat)
    cc_el_result <- topGO::getSigGroups(fisher_cc_GOdata, test_stat)
    test_stat <- new("weightCount", testStatistic=GOFisherTest, name="Fisher test", sigRatio="ratio")
    mf_weight_result <- topGO::getSigGroups(fisher_mf_GOdata, test_stat)
    bp_weight_result <- topGO::getSigGroups(fisher_bp_GOdata, test_stat)
    cc_weight_result <- topGO::getSigGroups(fisher_cc_GOdata, test_stat)

    mf_fisher_pdist <- try(hpgltools::hpgl_histogram(mf_fisher_result@score, bins=20))
    mf_ks_pdist <- try(hpgltools::hpgl_histogram(mf_ks_result@score, bins=20))
    mf_el_pdist <- try(hpgltools::hpgl_histogram(mf_el_result@score, bins=20))
    mf_weight_pdist <- try(hpgltools::hpgl_histogram(mf_weight_result@score, bins=20))
    bp_fisher_pdist <- try(hpgltools::hpgl_histogram(bp_fisher_result@score, bins=20))
    bp_ks_pdist <- try(hpgltools::hpgl_histogram(bp_ks_result@score, bins=20))
    bp_el_pdist <- try(hpgltools::hpgl_histogram(bp_el_result@score, bins=20))
    bp_weight_pdist <- try(hpgltools::hpgl_histogram(bp_weight_result@score, bins=20))
    cc_fisher_pdist <- try(hpgltools::hpgl_histogram(cc_fisher_result@score, bins=20))
    cc_ks_pdist <- try(hpgltools::hpgl_histogram(cc_ks_result@score, bins=20))
    cc_el_pdist <- try(hpgltools::hpgl_histogram(cc_el_result@score, bins=20))
    cc_weight_pdist <- try(hpgltools::hpgl_histogram(cc_weight_result@score, bins=20))
    p_dists <- list(mf_fisher=mf_fisher_pdist, bp_fisher=bp_fisher_pdist,
                    cc_fisher=cc_fisher_pdist, mf_ks=mf_ks_pdist, bp_ks=bp_ks_pdist,
                    cc_ks=cc_ks_pdist, mf_el=mf_el_pdist, bp_el=bp_el_pdist,
                    cc_el=cc_el_pdist, mf_weight=mf_weight_pdist,
                    bp_weight=bp_weight_pdist, cc_weight=cc_weight_pdist)

    results <- list(fmf_godata=fisher_mf_GOdata, fbp_godata=fisher_bp_GOdata,
                    fcc_godata=fisher_cc_GOdata, kmf_godata=ks_mf_GOdata,
                    kbp_godata=ks_bp_GOdata, kcc_godata=ks_cc_GOdata,
                    mf_fisher=mf_fisher_result, bp_fisher=bp_fisher_result,
                    cc_fisher=cc_fisher_result, mf_ks=mf_ks_result,
                    bp_ks=bp_ks_result, cc_ks=cc_ks_result,
                    mf_el=mf_el_result, bp_el=bp_el_result,
                    cc_el=cc_el_result, mf_weight=mf_weight_result,
                    bp_weight=bp_weight_result, cc_weight=cc_weight_result)

    tables <- try(topgo_tables(results, limitby=limitby, limit=limit))
    if (class(tables)[1] == 'try-error') {
        tables <- NULL
    }

    mf_densities <- bp_densities <- cc_densities <- list()
    if (isTRUE(densities)) {
        mf_densities <- hpgltools:::plot_topgo_densities(fisher_mf_GOdata, tables$mf)
        bp_densities <- hpgltools:::plot_topgo_densities(fisher_bp_GOdata, tables$bp)
        cc_densities <- hpgltools:::plot_topgo_densities(fisher_cc_GOdata, tables$cc)
    } else {
        message("simple_topgo(): Set densities=TRUE for ontology densitya plots.")
    }

    information <- list(
        mf_godata=fisher_mf_GOdata, bp_godata=fisher_bp_GOdata, cc_godata=fisher_cc_GOdata,
        kmf_godata=ks_mf_GOdata, kbp_godata=ks_bp_GOdata, kcc_godata=ks_cc_GOdata,
        results=results, tables=tables, mf_densities=mf_densities,
        bp_densities=bp_densities, cc_densities=cc_densities, pdists=p_dists)

    if (isTRUE(pval_plots)) {
        information$pvalue_plots <- topgo_pval_plot(information)
    }

    return(information)
}

#' plot_topgo_densities()  Plot the density of categories vs. the possibilities
#'
#' This can make a large number of plots
#'
#' @param godata  the result from topgo
plot_topgo_densities <- function(godata, table) {
    ret <- list()
    for (id in table$GO.ID) {
        message(id)
        print(hpgltools:::hpgl_GroupDensity(godata, id, ranks=TRUE))
        added_plot <- recordPlot()
        ret[[id]] <- added_plot
    }
    return(ret)
}

#' topgo_tables() Make pretty tables out of topGO data
#'
#' The topgo function GenTable is neat, but it needs some simplification to not be obnoxious
#'
#' @param result a topgo result
#' @param limit a pvalue limit defining 'significant'
#' @param limitby fisher - what type of test to perform
#' @param numchar 300 how many characters to allow in the description
#' @param orderby classic which of the available columns to order the table by?
#' @param ranksof classic which of the available columns are used to rank the data?
#' @export
topgo_tables <- function(result, limit=0.1, limitby="fisher", numchar=300, orderby="classic", ranksof="classic") {
    ## The following if statement could be replaced by get(limitby)
    ## But I am leaving it as a way to ensure that no shenanigans ensue
    mf_allRes <- bp_allRes <- cc_allRes <- mf_interesting <- bp_interesting <- cc_interesting <- NULL
    if (limitby == "fisher") {
        mf_siglist <- names(which(result$mf_fisher@score <= limit))
        bp_siglist <- names(which(result$bp_fisher@score <= limit))
        cc_siglist <- names(which(result$bp_fisher@score <= limit))
    } else if (limitby == "KS") {
        mf_siglist <- names(which(result$mf_ks@score <= limit))
        bp_siglist <- names(which(result$bp_ks@score <= limit))
        cc_siglist <- names(which(result$bp_ks@score <= limit))
    } else if (limitby == "EL") {
        mf_siglist <- names(which(result$mf_el@score <= limit))
        bp_siglist <- names(which(result$bp_el@score <= limit))
        cc_siglist <- names(which(result$bp_el@score <= limit))
    } else if (limitby == "weight") {
        mf_siglist <- names(which(result$mf_weight@score <= limit))
        bp_siglist <- names(which(result$bp_weight@score <= limit))
        cc_siglist <- names(which(result$bp_weight@score <= limit))
    } else {
        stop("I can only limit by: fisher, KS, EL, or weight.")
    }
    mf_topnodes <- length(mf_siglist)
    if (mf_topnodes > 0) {
        mf_allRes <- try(topGO::GenTable(result$fmf_godata, classic=result$mf_fisher, KS=result$mf_ks,
                                         EL=result$mf_el, weight=result$mf_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=mf_topnodes, numChar=numchar))
        if (class(mf_allRes) != 'try-error') {
            mf_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$mf_fisher))$qvalues)
            mf_allRes <- merge(mf_allRes, mf_qvalues, by.x="GO.ID", by.y="row.names")
            mf_allRes$classic <- as.numeric(mf_allRes$classic)
            mf_allRes <- mf_allRes[with(mf_allRes, order(classic)), ]
            colnames(mf_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            mf_interesting <- subset(mf_allRes, get(limitby) <= limit)
            rownames(mf_interesting) <- NULL
            mf_interesting$ont <- "MF"
            mf_interesting <- mf_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                                "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    bp_topnodes <- length(bp_siglist)
    if (bp_topnodes > 0) {
        bp_allRes <- try(topGO::GenTable(result$fbp_godata, classic=result$bp_fisher, KS=result$bp_ks,
                                         EL=result$bp_el, weight=result$bp_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=bp_topnodes, numChar=numchar))
        if (class(bp_allRes) != 'try-error') {
            bp_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$bp_fisher))$qvalues)
            bp_allRes <- merge(bp_allRes, bp_qvalues, by.x="GO.ID", by.y="row.names", all.x=TRUE)
            bp_allRes$classic <- as.numeric(bp_allRes$classic)
            bp_allRes <- bp_allRes[with(bp_allRes, order(classic)), ]
            colnames(bp_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            bp_interesting <- subset(bp_allRes, get(limitby) <= limit)
            rownames(bp_interesting) <- NULL
            bp_interesting$ont <- "BP"
            bp_interesting = bp_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                               "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    cc_topnodes <- length(cc_siglist)
    if (cc_topnodes > 0) {
        cc_allRes <- try(topGO::GenTable(result$fcc_godata, classic=result$cc_fisher, KS=result$cc_ks,
                                         EL=result$cc_el, weight=result$cc_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=cc_topnodes, numChar=numchar))
        if (class(cc_allRes) != 'try-error') {
            cc_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$cc_fisher))$qvalues)
            cc_allRes <- merge(cc_allRes, cc_qvalues, by.x="GO.ID", by.y="row.names")
            cc_allRes$classic <- as.numeric(cc_allRes$classic)
            cc_allRes <- cc_allRes[with(cc_allRes, order(classic)), ]
            colnames(cc_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            cc_interesting <- subset(cc_allRes, get(limitby) <= limit)
            rownames(cc_interesting) <- NULL
            cc_interesting$ont <- "CC"
            cc_interesting <- cc_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                                "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    tables <- list(mf=mf_allRes, bp=bp_allRes, cc=cc_allRes,
                   mf_interesting=mf_interesting, bp_interesting=bp_interesting,
                   cc_interesting=cc_interesting)
    return(tables)
}

#' Print trees from topGO
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param goid_map a file containing mappings of genes to goids in the format expected by topgo
#'
#' @return a big list including the various outputs from topgo
#' @export
topgo_trees <- function(tg, score_limit=0.01, sigforall=TRUE, do_mf_fisher_tree=TRUE,
                        do_bp_fisher_tree=TRUE, do_cc_fisher_tree=TRUE, do_mf_ks_tree=FALSE,
                        do_bp_ks_tree=FALSE, do_cc_ks_tree=FALSE, do_mf_el_tree=FALSE,
                        do_bp_el_tree=FALSE, do_cc_el_tree=FALSE, do_mf_weight_tree=FALSE,
                        do_bp_weight_tree=FALSE, do_cc_weight_tree=FALSE) {
    mf_fisher_nodes <- mf_fisher_tree <- NULL
    if (isTRUE(do_mf_fisher_tree)) {
        included <- length(which(topGO::score(tg$results$mf_fisher) <= score_limit))
        mf_fisher_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_fisher),
                                                                      useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                      useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_fisher_nodes)[1] != 'try-error') {
            mf_fisher_tree <- try(recordPlot())
        }
    }
    bp_fisher_nodes <- bp_fisher_tree <- NULL
    if (isTRUE(do_bp_fisher_tree)) {
        included <- length(which(topGO::score(tg$results$bp_fisher) <= score_limit))
        bp_fisher_nodes <- try(suppressWarnings(showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_fisher),
                                                               useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_fisher_nodes)[1] != 'try-error') {
            bp_fisher_tree <- try(recordPlot())
        }
    }
    cc_fisher_nodes <- cc_fisher_tree <- NULL
    if (isTRUE(do_cc_fisher_tree)) {
        included <- length(which(topGO::score(tg$results$cc_fisher) <= score_limit))
        cc_fisher_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_fisher),
                                                                      useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                      useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_fisher_nodes)[1] != 'try-error') {
            cc_fisher_tree <- try(recordPlot())
        }
    }
    mf_ks_nodes <- mf_ks_tree <- NULL
    if (isTRUE(do_mf_ks_tree)) {
        included <- length(which(topGO::score(tg$results$mf_ks) <= score_limit))
        mf_ks_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_ks),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_ks_nodes)[1] != 'try-error') {
            mf_ks_tree <- try(recordPlot())
        }
    }
    bp_ks_nodes <- bp_ks_tree <- NULL
    if (isTRUE(do_bp_ks_tree)) {
        included <- length(which(topGO::score(tg$results$bp_ks) <= score_limit))
        bp_ks_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_ks),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_ks_nodes)[1] != 'try-error') {
            bp_ks_tree <- try(recordPlot())
        }
    }
    cc_ks_nodes <- cc_ks_tree <- NULL
    if (isTRUE(do_cc_ks_tree)) {
        included <- length(which(topGO::score(tg$results$cc_ks) <= score_limit))
        cc_ks_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_ks),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_ks_nodes)[1] != 'try-error') {
            cc_ks_tree <- try(recordPlot())
        }
    }
    mf_el_nodes <- mf_el_tree <- NULL
    if (isTRUE(do_mf_el_tree)) {
        included <- length(which(topGO::score(tg$results$mf_el) <= score_limit))
        mf_el_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_el),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_el_nodes)[1] != 'try-error') {
            mf_el_tree <- try(recordPlot())
        }
    }
    bp_el_nodes <- bp_el_tree <- NULL
    if (isTRUE(do_bp_el_tree)) {
        included <- length(which(topGO::score(tg$results$bp_el) <= score_limit))
        bp_el_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_el),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_el_nodes)[1] != 'try-error') {
            bp_el_tree <- try(recordPlot())
        }
    }
    cc_el_nodes <- cc_el_tree <- NULL
    if (isTRUE(do_cc_el_tree)) {
        included <- length(which(topGO::score(tg$results$cc_el) <= score_limit))
        cc_el_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_el),
                                                                  useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                  useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_el_nodes)[1] != 'try-error') {
            cc_el_tree <- try(recordPlot())
        }
    }
    mf_weight_nodes <- mf_weight_tree <- NULL
    if (isTRUE(do_mf_weight_tree)) {
        included <- length(which(topGO::score(tg$results$mf_weight) <= score_limit))
        mf_weight_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_weight),
                                                                      useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                      useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_weight_nodes)[1] != 'try-error') {
            mf_weight_tree <- try(recordPlot())
        }
    }
    bp_weight_nodes <- bp_weight_tree <- NULL
    if (isTRUE(do_bp_weight_tree)) {
        included <- length(which(topGO::score(tg$results$bp_weight) <= score_limit))
        bp_weight_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_weight),
                                                                      useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                      useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_weight_nodes)[1] != 'try-error') {
            bp_weight_tree <- try(recordPlot())
        }
    }
    cc_weight_nodes <- cc_weight_tree <- NULL
    if (isTRUE(do_cc_weight_tree)) {
        included <- length(which(topGO::score(tg$results$cc_weight) <= score_limit))
        cc_weight_nodes <- try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_weight),
                                                                      useInfo="all", sigForAll=sigforall, firstSigNodes=included,
                                                                      useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_weight_nodes)[1] != 'try-error') {
            cc_weight_tree <- try(recordPlot())
        }
    }

    trees <- list(
        mf_fisher_nodes=mf_fisher_nodes, bp_fisher_nodes=bp_fisher_nodes, cc_fisher_nodes=cc_fisher_nodes,
        mf_ks_nodes=mf_ks_nodes, bp_ks_nodes=bp_ks_nodes, cc_ks_nodes=cc_ks_nodes,
        mf_el_nodes=mf_el_nodes, bp_el_nodes=bp_el_nodes, cc_el_nodes=cc_el_nodes,
        mf_weight_nodes=mf_weight_nodes, bp_weight_nodes=bp_weight_nodes, cc_weight_nodes=cc_weight_nodes,
        mf_fisher_tree=mf_fisher_tree, bp_fisher_tree=bp_fisher_tree, cc_fisher_tree=cc_fisher_tree,
        mf_ks_tree=mf_ks_tree, bp_ks_tree=bp_ks_tree, cc_ks_tree=cc_ks_tree,
        mf_el_tree=mf_el_tree, bp_el_tree=bp_el_tree, cc_el_tree=cc_el_tree,
        mf_weight_tree=mf_weight_tree, bp_weight_tree=bp_weight_tree, cc_weight_tree=cc_weight_tree)
    return(trees)
}

#' Make a go mapping from IDs in a format suitable for topGO
#'
#' @param goid_map A topGO mapping file
#' @param goids_df If there is no goid_map, create it with this
#' @param overwrite A boolean, if it already exists, rewrite the mapping file?
#'
#' @return a summary of the new goid table
#'
#' @export
make_id2gomap <- function(goid_map="reference/go/id2go.map", goids_df=NULL, overwrite=FALSE) {
    id2go_test <- file.info(goid_map)
    goids_dir <- dirname(goid_map)
    new_go <- NULL
    if (!file.exists(goids_dir)) {
        dir.create(goids_dir, recursive=TRUE)
    }
    if (isTRUE(overwrite)) {
        if (is.null(goids_df)) {
            stop("There is neither a id2go file nor a data frame of goids.")
        } else {
            message("Attempting to generate a id2go file in the format expected by topGO.")
            new_go <- plyr::ddply(goids_df, .(ID), summarise, GO=paste(unique(GO), collapse=','))
            write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
            rm(id2go_test)
        }
    } else { ## overwrite is not true
        if (is.na(id2go_test$size)) {
            if (is.null(goids_df)) {
                stop("There is neither a id2go file nor a data frame of goids.")
            } else {
                message("Attempting to generate a id2go file in the format expected by topGO.")
                new_go <- plyr::ddply(goids_df, .(ID), summarise, GO=paste(unique(GO), collapse=','))
                write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
                id2go_test <- file.info(goid_map)
            }
        } else { ## There already exists a file, so return its stats
            ## new_go = id2go_test
            new_go <- goid_map
        }
    }
    return(new_go)
}

hpgl_topdiffgenes <- function(scores, df=de_genes, direction="up") {
    ## Testing parameters
    ##scores = pvals
    ##df = epi_cl14clbr_high
    ## Here is the original topDiffGenes
    ## topDiffGenes <- function(allScore) {
    ##   return(allScore < 0.01)
    ##}
    ## my version of this will expect a limma result table from which I will extract the entries with low p-values
    ## and logFCs which are high or low
    quartiles <- summary(df)
}

#' A very simple selector of strong scoring genes (by p-value)
#'
#' This function was provided in the topGO documentation, but not defined.
#' It was copied/pasted here.  I have ideas for including up/down expression
#' but have so far deemed them not needed because I am feeding topGO
#' already explicit lists of genes which are up/down/whatever.
#' But it still is likely to be useful to be able to further subset the data.
#'
#' @param allScore The scores of the genes
#' @export
topDiffGenes <- function(allScore) { return(allScore < 0.01) }

#' Make a pvalue plot from topgo data
#'
#' @param topgo_data some data from topgo!
#'
#' @return a plot!
#' @seealso \code{\link{goseq}}
#' @export
topgo_pval_plot <- function(topgo, wrapped_width=20, cutoff=0.1, n=12, type="fisher") {
    mf_newdf <- topgo$tables$mf[,c("GO.ID", "Term", "Annotated","Significant", type)]
    mf_newdf$term <- as.character(lapply(strwrap(mf_newdf$Term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    mf_newdf$pvalue <- as.numeric(mf_newdf[[type]])
    mf_newdf <- subset(mf_newdf, get(type) < cutoff)
    mf_newdf <- mf_newdf[order(mf_newdf$pvalue, mf_newdf[[type]]),]
    mf_newdf <- head(mf_newdf, n=n)
    mf_newdf$score <- mf_newdf$Significant / mf_newdf$Annotated
    mf_pval_plot <- pval_plot(mf_newdf, ontology="MF")

    bp_newdf <- topgo$tables$bp[,c("GO.ID", "Term", "Annotated","Significant",type)]
    bp_newdf$term <- as.character(lapply(strwrap(bp_newdf$Term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    bp_newdf$pvalue <- as.numeric(bp_newdf[[type]])
    bp_newdf <- subset(bp_newdf, get(type) < cutoff)
    bp_newdf <- bp_newdf[order(bp_newdf$pvalue, bp_newdf[[type]]),]
    bp_newdf <- head(bp_newdf, n=n)
    bp_newdf$score <- bp_newdf$Significant / bp_newdf$Annotated
    bp_pval_plot <- pval_plot(bp_newdf, ontology="MF")

    cc_newdf <- topgo$tables$cc[,c("GO.ID", "Term", "Annotated","Significant",type)]
    cc_newdf$term <- as.character(lapply(strwrap(cc_newdf$Term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    cc_newdf$pvalue <- as.numeric(cc_newdf[[type]])
    cc_newdf <- subset(cc_newdf, get(type) < cutoff)
    cc_newdf <- cc_newdf[order(cc_newdf$pvalue, cc_newdf[[type]]),]
    cc_newdf <- head(cc_newdf, n=n)
    cc_newdf$score <- cc_newdf$Significant / cc_newdf$Annotated
    cc_pval_plot <- pval_plot(cc_newdf, ontology="CC")

    pval_plots <- list(MF=mf_pval_plot, BP=bp_pval_plot, CC=cc_pval_plot)
    return(pval_plots)
}


## this function will plot the GO DAG or parts of it
## sigNodes:     a named vector of terms p-values, the names are the GO terms
## wantedNodes:  the nodes that we want to find, we will plot this nodes with
##               a different color. The vector contains the names pf the nodes
## oldSigNodes:  used to plot the (new) sigNodes in the same collor range
##               as the old ones
## export.to.dot.file: is a global variable given the name of the output .dot file

#' @export
getEdgeWeights <- function(graph) {
    weightsList <- graph::edgeWeights(graph)
    to <- lapply(weightsList, names)
    from <- nodes(graph)

    if (any(is.na(unlist(to))) || any(is.na(from))) {
        stop("Edge names do not match node names.")
    }

    edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
    edge.weights <- unlist(weightsList)
    names(edge.weights) <- edge.names
    return(edge.weights)
}

#' A minor hack in the topGO GOplot function
#'
#' @export
hpgl_GOplot <- function(dag, sigNodes, dag.name='GO terms', edgeTypes=TRUE,
                        nodeShape.type=c('box','circle','ellipse','plaintext')[3],
                        genNodes=NULL, wantedNodes=NULL, showEdges=TRUE, useFullNames=TRUE,
                        oldSigNodes=NULL, nodeInfo=nodeInfo, maxchars=30) {
    if(!missing(sigNodes)) {
        sigNodeInd = TRUE
    } else {
        sigNodeInd = FALSE
    }

    ## we set the global Graphviz attributes
    ##  graphAttrs <- getDefaultAttrs(layoutType = 'dot')
    graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
    graphAttrs$cluster <- NULL
    graphAttrs$edge$arrowsize = "0.4"
    graphAttrs$edge$weight = "0.01"

    ##graphAttrs$graph$splines <- FALSE
    graphAttrs$graph$size = "12.0,12.0"
    graphAttrs$graph$margin = "0.0,0.0"
    ##  graphAttrs$graph$ranksep = "0.02"
    ##  graphAttrs$graph$nodesep = "0.30"

    ## set the node shape
    graphAttrs$node$shape <- nodeShape.type
    ##graphAttrs$node$fixedsize <- FALSE
    ## set the fontsize for the nodes labels
    graphAttrs$node$fontsize <- '20.0'
    graphAttrs$node$height <- '2.0'
    graphAttrs$node$width <- '3.0'
    graphAttrs$graph$size = "12,12"
    graphAttrs$node$color = "lightblue"
    graphAttrs$node$fontname = "arial"
    graphAttrs$node$style = "invis"

    ## set the local attributes lists
    nodeAttrs <- list()
    edgeAttrs <- list()

    ## try to use adaptive node size
    ##nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))

    if(is.null(nodeInfo)) {
        nodeInfo <- character(numNodes(dag))
        names(nodeInfo) <- nodes(dag)
    } else {
##        print(class(nodeInfo))
##        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
        nodeInfo = gsub("(\\w.{18}).*(\\\\\\n)","\\1\\2", nodeInfo, perl=TRUE)
        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
    }
    ##teststring = paste("test:", nodeInfo)
    ##print(teststring)

  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
    node.names <- nodes(dag)
    if(!useFullNames) {
        nodeAttrs$label <- sapply(node.names,
                                  function(x) {
                                      return(paste(substr(x, 4, nchar(node.names[1])),
                                                   nodeInfo[x], sep = ''))
                                  })
    } else {
        nodeAttrs$label <- paste(node.names, nodeInfo, sep = '')
        names(nodeAttrs$label) <- node.names
    }

  ## we will change the shape and the color of the nodes that generated the dag
    if(!is.null(wantedNodes)) {
        diffNodes <- setdiff(wantedNodes, genNodes)
        if(length(diffNodes) > 0) {
            nodeAttrs$color[diffNodes] <- rep('lightblue', .ln <- length(diffNodes))
            nodeAttrs$shape[diffNodes] <- rep('circle', .ln)
            nodeAttrs$height[diffNodes] <- rep('0.45', .ln)
            ## nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
            ## nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
        }
    }

    ## we will change the shape and the color of the nodes we want back
    if(!is.null(genNodes)) {
        nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
        nodeAttrs$shape[genNodes] <- rep('box', .ln)
        ## nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)
    }

    ## we will use different fillcolors for the nodes
    if(sigNodeInd) {
        if(!is.null(oldSigNodes)) {
            old.logSigNodes <- log10(sort(oldSigNodes[nodes(dag)]))
            old.range <- range(old.logSigNodes)
            logSigNodes <- log10(sort(sigNodes))
            logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
            logSigNodes[logSigNodes > old.range[2]] <- old.range[2]
            ## debug:  old.range == range(logSigNodes)
            ## if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)) {
            ## print(old.range)
            ## print(range(logSigNodes))
            ## stop('some stupid error here :)')
            ## }
        } else {
            old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
        }

        sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
        old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)

        mm <- max(sigColor, old.sigColor)
        sigColor <- sigColor + (mm - max(sigColor))

        colorMap <- heat.colors(mm)
        nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
    }

    if(!showEdges) {
        graphAttrs$edge$color <- 'white'
    } else {
        ## if we want to differentiate between 'part-of' and 'is-a' edges
        if (edgeTypes) {
            ## 0 for a is_a relation,  1 for a part_of relation
            ## edgeAttrs$color <- ifelse(getEdgeWeights(dag) == 0, 'black', 'red')
            edgeAttrs$color <- ifelse(hpgltools::getEdgeWeights(dag) == 0, 'black', 'black')
        }
        ## plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)
    }

    return(agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                  nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}


GOplot.orig <- function(dag, sigNodes, dag.name='GO terms', edgeTypes=TRUE,
                   nodeShape.type=c('box', 'circle', 'ellipse', 'plaintext')[3],
                   genNodes=NULL, wantedNodes=NULL, showEdges=TRUE, useFullNames=FALSE,
                   oldSigNodes=NULL, nodeInfo=NULL) {

    if(!missing(sigNodes))
        sigNodeInd = TRUE
    else
        sigNodeInd = FALSE

    ## we set the global Graphviz attributes
    graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
    graphAttrs$cluster <- NULL

    ## graphAttrs$graph$splines <- FALSE
    ## set the node shape
    graphAttrs$node$shape <- nodeShape.type

    ## set the fontsize for the nodes labels
    graphAttrs$node$fontsize <- '14'
    ## graphAttrs$node$height <- '1.0'
    ## graphAttrs$node$width <- '1.5'

    ## set the local attributes lists
    nodeAttrs <- list()
    edgeAttrs <- list()

    ## try to use adaptive node size
    ## nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))

    if (is.null(nodeInfo)) {
        nodeInfo <- character(numNodes(dag))
        names(nodeInfo) <- nodes(dag)
    } else {
        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
    }
    ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
    node.names <- nodes(dag)
    if (!useFullNames)
        nodeAttrs$label <- sapply(node.names,
                                  function(x) {
                                      return(paste(substr(x, 4, nchar(node.names[1])),
                                                   nodeInfo[x], sep = ''))
                                  })
    else {
        nodeAttrs$label <- paste(node.names, nodeInfo, sep = '')
        names(nodeAttrs$label) <- node.names
    }

  ## we will change the shape and the color of the nodes that generated the dag
  if (!is.null(wantedNodes)) {
      diffNodes <- setdiff(wantedNodes, genNodes)
      if(length(diffNodes) > 0) {
          nodeAttrs$color[diffNodes] <- rep('lightblue', .ln <- length(diffNodes))
          nodeAttrs$shape[diffNodes] <- rep('circle', .ln)
          nodeAttrs$height[diffNodes] <- rep('0.45', .ln)
          ##nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
          ##nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
      }
  }

    ## we will change the shape and the color of the nodes we want back
    if(!is.null(genNodes)) {
        nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
        nodeAttrs$shape[genNodes] <- rep('box', .ln)
        ## nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)
    }

    ## we will use different fillcolors for the nodes
    if(sigNodeInd) {
        if(!is.null(oldSigNodes)) {
            old.logSigNodes <- log10(sort(oldSigNodes[nodes(dag)]))
            old.range <- range(old.logSigNodes)
            logSigNodes <- log10(sort(sigNodes))
            logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
            logSigNodes[logSigNodes > old.range[2]] <- old.range[2]

            ## debug:  old.range == range(logSigNodes)
            ## if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)){
            ##  print(old.range)
            ##  print(range(logSigNodes))
            ##  stop('some stupid error here :)')
            ## }
        }
    else
        old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))

        sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
        old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)

        mm <- max(sigColor, old.sigColor)
        sigColor <- sigColor + (mm - max(sigColor))

        colorMap <- heat.colors(mm)
        nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
    }

    if(!showEdges)
        graphAttrs$edge$color <- 'white'
    else
        ## if we want to differentiate between 'part-of' and 'is-a' edges
        if(edgeTypes)
            ##    0 for a is_a relation,  1 for a part_of relation
            ## edgeAttrs$color <- ifelse(getEdgeWeights(dag) == 0, 'black', 'red')
            edgeAttrs$color <- ifelse(hpgltools::getEdgeWeights(dag) == 0, 'black', 'black')
    ## plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)

    return(agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                  nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}


#' hpgl_GroupDensity()  A hack of topGO's groupDensity()
#'
#' This just adds a couple wrappers to avoid errors in groupDensity.
#'
#' @param object  a topGO enrichment object
#' @param whichGO  an individual ontology group to compare with
#' @param ranks default=TRUE  rank order the set of ontologies
#' @param rm.one default=FALSE  remove pvalue=1 groups
hpgl_GroupDensity = function(object, whichGO, ranks=TRUE, rm.one=FALSE) {
    groupMembers <- try(topGO::genesInTerm(object, whichGO)[[1]])
    if (class(groupMembers)[1] == 'try-error') {
        return(NULL)
    }
    allS <- topGO::geneScore(object, use.names = TRUE)
    if (rm.one) {
        allS <- allS[allS < 0.99]
    }
    xlab <- "Gene' score"
    if (ranks) {
        allS <- BiocGenerics::rank(allS, ties.method = "random")
        xlab <- "Gene's rank"
    }
    group <- as.integer(names(allS) %in% groupMembers)
    xx <- data.frame(score=allS, group = factor(group, labels=paste(c("complementary", whichGO), "  (", table(group), ")", sep="")))
    plot = lattice::densityplot( ~ score | group, data=xx, layout=c(1,2), xlab=xlab)
    return(plot)
}

# EOF
