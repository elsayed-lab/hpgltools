## Time-stamp: <Sat Jan 23 00:30:03 2016 Ashton Trey Belew (abelew@gmail.com)>

check_clusterprofiler <- function(gff='test.gff') {
    genetable_test <- try(load("geneTable.rda"))
    if (class(genetable_test) == 'try-error') {
        if (!is.null(gff)) {
            message("simple_clus(): Generating the geneTable.rda")
            hpgltools:::hpgl_Gff2GeneTable(gff)
            ##clusterProfiler:::Gff2GeneTable(gff)
        } else {
            stop("simple_clus(): requires geneTable.rda, thus a gff file.")
        }
    } else {
        rm(genetable_test)
    }
    gomapping_test <- try(load("GO2EG.rda"), silent=TRUE)
    if (class(gomapping_test) == 'try-error') {
        message("simple_clus(): Generating GO mapping data.")
        gomap <- goids
        gomap <- gomap[,c(1,2)]
        colnames(gomap) <- c("entrezgene", "go_accession")
        ## It turns out that the author of clusterprofiler reversed these fields...
        ## Column 1 must be GO ID, column 2 must be gene accession.
        gomap <- gomap[,c("go_accession","entrezgene")]
        clusterProfiler::buildGOmap(gomap)
    } else {
        message("Using GO mapping data located in GO2EG.rda")
    }
}

#' Perform a simplified clusterProfiler analysis
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param goids a file containing mappings of genes to goids in the format expected by topgo
#' @param golevel a relative level in the tree for printing p-value plots, higher is more specific
#' @param pcutoff a p-value cutoff
#' @param qcutoff a q-value cutoff
#' @param padjust a method for adjusting the p-values
#' @param fold_changes a df of fold changes for the DE genes
#' @param include_cnetplots the cnetplots are often stupid and can be left behind
#' @param showcategory how many categories to show in p-value plots
#'
#' @return a big list including the following:
#'   mf_interesting: A table of the interesting molecular function groups
#'   bp_interesting: A table of the interesting biological process groups
#'   cc_interesting: A table of the interesting cellular component groups
#'   mf_pvals: A histogram of the molecular function p-values
#'   bp_pvals: Ditto, biological process
#'   cc_pvals: And cellular component...
#'   mf_enriched: A table of the enriched molecular function groups by adjusted p-value.
#'   bp_enriched: yep, you guessed it
#'   cc_enriched: cellular component, too
#'   mf_all/bp_all/cc_all: A table of all go categories observed (mf/bp/cc respectively)
#'   mfp_plot/bpp_plot/ccp_plot: ggplot2 p-value bar plots describing the over represented groups
#'   mf_cnetplot/bp_cnetplot/cc_cnetplot: clusterProfiler cnetplots
#'   mf_group_barplot/bp_group_barplot/cc_group_barplot: The group barplots from clusterProfiler
#' @export
#' @examples
#' ## up_cluster = simple_clusterprofiler(mga2_ll_thy_top, goids=goids, gff="reference/genome/gas.gff")
#' ## > Some chattery while it runs
#' ## tail(head(up_cluster$bp_interesting, n=10), n=1)
#' ## > ID ont GeneRatio BgRatio     pvalue   p.adjust    qvalue
#' ## > 10 GO:0009311  BP     5/195 10/1262 0.01089364 0.01089364 0.1272835
#' ## >   geneID Count
#' ## >   10 M5005_Spy1632/M5005_Spy1637/M5005_Spy1635/M5005_Spy1636/M5005_Spy1638     5
#' ## >   Description
#' ## >   10 oligosaccharide metabolic process
simple_clusterprofiler <- function(de_genes, goids=NULL, golevel=4, pcutoff=0.1,
                                   qcutoff=1.0, fold_changes=NULL, include_cnetplots=TRUE,
                                   showcategory=12, universe=NULL, organism="lm", gff=NULL,
                                   wrapped_width=20, method="Walllenius", padjust="BH", ...) {

    check_clusterprofiler(gff)
    if (is.null(de_genes$ID)) {
        gene_list <- as.character(rownames(de_genes))
    } else {
        gene_list <- as.character(de_genes$ID)
    }
##    message("Testing gseGO")
##    ego2 = try(clusterProfiler::gseGO(geneList=gene_list, organism=organism, ont="GO", nPerm=100, minGSSize=2, pvalueCutoff=1, verbose=TRUE))
##    print(ego2)
    message("simple_clus(): Starting MF(molecular function) analysis")
    mf_group <- clusterProfiler::groupGO(gene_list, organism=organism, ont="MF", level=golevel, readable=TRUE)
    mf_all <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="MF",
                                       pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod="none")
    all_mf_phist <- try(hpgltools::hpgl_histogram(mf_all@result$pvalue, bins=20))
    if (class(all_mf_phist)[1] != 'try-error') {
        y_limit <- (sort(unique(table(all_mf_phist$data)), decreasing=TRUE)[2]) * 2
        all_mf_phist <- all_mf_phist + scale_y_continuous(limits=c(0, y_limit))
    }
    enriched_mf <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="MF",
                                            pvalueCutoff=pcutoff, qvalueCutoff=qcutoff, pAdjustMethod=padjust)

    message("simple_clus(): Starting BP(biological process) analysis")
    bp_group <- clusterProfiler::groupGO(gene_list, organism=organism, ont="BP", level=golevel, readable=TRUE)
    bp_all <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="BP", pvalueCutoff=1.0,
                                       qvalueCutoff=1.0, pAdjustMethod="none")
    all_bp_phist <- try(hpgltools::hpgl_histogram(bp_all@result$pvalue, bins=20))
    if (class(all_bp_phist)[1] != 'try-error') {
        y_limit <- (sort(unique(table(all_bp_phist$data)), decreasing=TRUE)[2]) * 2
        all_bp_phist <- all_bp_phist + scale_y_continuous(limits=c(0, y_limit))
    }

    enriched_bp <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="BP",
                                            pvalueCutoff=pcutoff, qvalueCutoff=qcutoff, pAdjustMethod=padjust)

    message("simple_clus(): Starting CC(cellular component) analysis")
    cc_group <- clusterProfiler::groupGO(gene_list, organism=organism, ont="CC", level=golevel, readable=TRUE)
    cc_all <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="CC", pvalueCutoff=1.0,
                                       qvalueCutoff=1.0, pAdjustMethod="none")
    enriched_cc <- hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="CC", pvalueCutoff=pcutoff,
                                            qvalueCutoff=qcutoff, pAdjustMethod=padjust)
    all_cc_phist <- try(hpgltools::hpgl_histogram(cc_all@result$pvalue, bins=20))
    ## Try and catch if there are no significant hits.
    if (class(all_cc_phist)[1] != 'try-error') {
        y_limit <- (sort(unique(table(all_cc_phist$data)), decreasing=TRUE)[2]) * 2
        all_cc_phist <- all_cc_phist + scale_y_continuous(limits=c(0, y_limit))
    }

    mf_group_barplot <- try(barplot(mf_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(mf_group_barplot)[1] != 'try-error') {
        mf_group_barplot$data$Description <- as.character(lapply(strwrap(mf_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    bp_group_barplot <- try(barplot(bp_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(bp_group_barplot)[1] != 'try-error') {
        bp_group_barplot$data$Description <- as.character(lapply(strwrap(bp_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    cc_group_barplot <- try(barplot(cc_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(cc_group_barplot)[1] != 'try-error') {
        cc_group_barplot$data$Description <- as.character(lapply(strwrap(cc_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    all_mf_barplot <- try(barplot(mf_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    enriched_mf_barplot <- try(barplot(enriched_mf, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_mf_barplot)[1] == 'try-error') {
        message("simple_clus(): No enriched MF groups were observed.")
    } else {
        enriched_mf_barplot$data$Description <- as.character(lapply(strwrap(enriched_mf_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_mf_barplot)[1] != 'try-error') {
        all_mf_barplot$data$Description <- as.character(lapply(strwrap(all_mf_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    all_bp_barplot <- try(barplot(bp_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    enriched_bp_barplot <- try(barplot(enriched_bp, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_bp_barplot)[1] == 'try-error') {
        message("simple_clus(): No enriched BP groups observed.")
    } else {
        enriched_bp_barplot$data$Description <- as.character(lapply(strwrap(enriched_bp_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_bp_barplot)[1] != 'try-error') {
        all_bp_barplot$data$Description <- as.character(lapply(strwrap(all_bp_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    all_cc_barplot <- try(barplot(cc_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    enriched_cc_barplot <- try(barplot(enriched_cc, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_cc_barplot)[1] == 'try-error') {
        message("simple_clus(): No enriched CC groups observed.")
    } else {
        enriched_cc_barplot$data$Description <- as.character(lapply(strwrap(enriched_cc_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_cc_barplot)[1] != 'try-error') {
        all_cc_barplot$data$Description <- as.character(lapply(strwrap(all_cc_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    cnetplot_mf <- cnetplot_bp <- cnetplot_cc <- NULL
    cnetplot_mfall <- cnetplot_bpall <- cnetplot_ccall <- NULL
    if (include_cnetplots == TRUE) {
        message("simple_clus(): Attempting to include the cnetplots.")
        message("simple_clus(): If they fail, set include_cnetplots to FALSE")
        cnetplot_mf <- try(clusterProfiler::cnetplot(enriched_mf, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_mf)[1] != 'try-error') {
            cnetplot_mf <- recordPlot()
        } else {
            message("simple_clus(): cnetplot failed for MF, no worries.")
        }
        cnetplot_bp <- try(clusterProfiler::cnetplot(enriched_bp, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_bp)[1] != 'try-error') {
            cnetplot_bp <- recordPlot()
        } else {
            message("simple_clus(): cnetplot failed for BP, no worries.")
        }
        cnetplot_cc <- try(clusterProfiler::cnetplot(enriched_cc, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_cc)[1] != 'try-error') {
            cnetplot_cc <- recordPlot()
        } else {
            message("simple_clus(): cnetplot failed for CC, no worries.")
        }
        cnetplot_mfall <- try(clusterProfiler::cnetplot(mf_all, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_mfall)[1] != 'try-error') {
            cnetplot_mfall <- recordPlot()
        }
        cnetplot_bpall <- try(clusterProfiler::cnetplot(bp_all, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_bpall)[1] != 'try-error') {
            cnetplot_bpall <- recordPlot()
        }
        cnetplot_ccall <- try(clusterProfiler::cnetplot(cc_all, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_ccall)[1] != 'try-error') {
            cnetplot_ccall <- recordPlot()
        }
    }

    if (!is.null(mf_all)) {
        mf_interesting <- mf_all@result
        rownames(mf_interesting) = NULL
        mf_interesting$ont <- "MF"
        mf_interesting <- mf_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue",
                                            "p.adjust","qvalue","geneID","Count","Description")]
        ## mf_interesting = subset(mf_interesting, pvalue <= 0.1)
        mf_interesting <- mf_interesting[ which(mf_interesting$pvalue <= 0.1), ]
    } else {
        mf_interesting <- NULL
    }
    if (!is.null(bp_all)) {
        bp_interesting <- bp_all@result
        rownames(bp_interesting) <- NULL
        bp_interesting$ont <- "BP"
        bp_interesting <- bp_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue",
                                            "p.adjust","qvalue","geneID","Count","Description")]
        ## bp_interesting = subset(bp_interesting, pvalue <= 0.1)
        bp_interesting <- bp_interesting[ which(bp_interesting$pvalue <= 0.1), ]
    } else {
        bp_interesting <- NULL
    }
    if (!is.null(cc_all)) {
        cc_interesting <- cc_all@result
        rownames(cc_interesting) <- NULL
        cc_interesting$ont <- "CC"
        cc_interesting <- cc_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue",
                                            "p.adjust","qvalue","geneID","Count","Description")]
        ## cc_interesting = subset(cc_interesting, pvalue <= 0.1)
        cc_interesting <- cc_interesting[ which(cc_interesting$pvalue <= 0.1), ]
    } else {
        cc_interesting <- NULL
    }

    return_information <- list(
        de_genes=de_genes,
        mf_interesting=mf_interesting, bp_interesting=bp_interesting, cc_interesting=cc_interesting,
        mf_pvals=all_mf_phist, bp_pvals=all_bp_phist, cc_pvals=all_cc_phist,
        mf_enriched=enriched_mf, bp_enriched=enriched_bp, cc_enriched=enriched_cc,
        mf_all=mf_all, bp_all=bp_all, cc_all=cc_all,
        mf_all_barplot=all_mf_barplot, bp_all_barplot=all_bp_barplot, cc_all_barplot=all_cc_barplot,
        mfp_plot=enriched_mf_barplot, bpp_plot=enriched_bp_barplot, ccp_plot=enriched_cc_barplot,
        mf_cnetplot=cnetplot_mf, bp_cnetplot=cnetplot_bp, cc_cnetplot=cnetplot_cc,
        mfall_cnetplot=cnetplot_mfall, bpall_cnetplot=cnetplot_bpall, ccall_cnetplot=cnetplot_ccall,
        mf_group=mf_group, bp_group=bp_group, cc_group=cc_group,
        mf_group_barplot=mf_group_barplot, bp_group_barplot=bp_group_barplot, cc_group_barplot=cc_group_barplot)
    return(return_information)
}

## Take clusterprofile group data and print it on a tree as topGO does
#' Make fun trees a la topgo from goseq data.
#'
#' @param godata data from cluster Profiler
#' @param goids a mapping of IDs to GO in the Ramigo expected format
#' @param sigforall Print significance on all nodes?
#'
#' @return plots! Trees! oh my!
#' @seealso \code{\link{Ramigo}}
#' @export
cluster_trees <- function(de_genes, cpdata, goid_map="reference/go/id2go.map", goids_df=NULL,
                          score_limit=0.2, overwrite=FALSE, selector="topDiffGenes", pval_column="adj.P.Val") {
    de_genes <- cpdata$de_genes
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO <- topGO::readMappings(file=goid_map)
    annotated_genes <- names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID <- make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes <- factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) <- annotated_genes

    message(paste0("Checking the de_table for a p-value column:", pval_column))
    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata <- new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
    } else {
        pvals <- as.vector(as.numeric(de_genes[[pval_column]]))
        names(pvals) <- rownames(de_genes)
        mf_GOdata <- new("topGOdata", description="MF", ontology="MF", allGenes=pvals,
                         geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", description="BP", ontology="BP", allGenes=pvals,
                         geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", description="CC", ontology="CC", allGenes=pvals,
                        geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    }

    mf_all <- cpdata$mf_all
    ## mf_enriched = cpdata$mf_enriched
    bp_all <- cpdata$bp_all
    ## bp_enriched = cpdata$bp_enriched
    cc_all <- cpdata$cc_all
    ## cc_enriched = cpdata$cc_enriched
    mf_all_ids <- mf_all@result$ID
    bp_all_ids <- bp_all@result$ID
    cc_all_ids <- cc_all@result$ID
    mf_all_scores <- mf_all@result$p.adjust
    bp_all_scores <- bp_all@result$p.adjust
    cc_all_scores <- cc_all@result$p.adjust
    names(mf_all_scores) <- mf_all_ids
    names(bp_all_scores) <- bp_all_ids
    names(cc_all_scores) <- cc_all_ids
    mf_included <- length(which(mf_all_scores <= score_limit))
    ## mf_tree_data = try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_all_scores, useInfo="all", sigForAll=TRUE, firstSigNodes=mf_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    mf_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_all_scores, useInfo="all",
                                                               sigForAll=TRUE, firstSigNodes=floor(mf_included * 1.5),
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(mf_tree_data)[1] == 'try-error') {
        mf_tree <- NULL
    } else {
        mf_tree <- recordPlot()
    }
    bp_included <- length(which(bp_all_scores <= score_limit))
    bp_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_all_scores, useInfo="all",
                                                              sigForAll=TRUE, firstSigNodes=bp_included,
                                                              useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(bp_tree_data)[1] == 'try-error') {
        bp_tree <- NULL
    } else {
        bp_tree <- recordPlot()
    }
    cc_included <- length(which(cc_all_scores <= score_limit))
    cc_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_all_scores, useInfo="all",
                                                               sigForAll=TRUE, firstSigNodes=cc_included,
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(cc_tree_data)[1] == 'try-error') {
        cc_tree <- NULL
    } else {
        cc_tree <- recordPlot()
    }
    trees <- list(MF=mf_tree, BP=bp_tree, CC=cc_tree, MFdata=mf_tree_data, BPdata=bp_tree_data, CCdata=cc_tree_data)
    return(trees)
}


#' A minor hack in the clusterProfiler function 'enrichGO'
#'
#' @param gene some differentially expressed genes
#' @param organism by default 'human'
#' @param ont by default 'MF'
#'
#' @return some clusterProfiler data
#' @seealso \code{\link{clusterProfiler}}
#' @export
hpgl_enrichGO <- function(gene, organism="human", ont="MF",
                          pvalueCutoff=0.05, pAdjustMethod="BH", universe,
                          qvalueCutoff=0.2, minGSSize=2, readable=FALSE) {
    information <- hpgl_enrich.internal(gene, organism=organism, pvalueCutoff=pvalueCutoff,
##        pAdjustMethod=pAdjustMethod, ont=ont, universe=universe,
                                        pAdjustMethod=pAdjustMethod, ont=ont,
                                        qvalueCutoff=qvalueCutoff, minGSSize=minGSSize)
##    print(summary(information))
    return(information)
}

#' A minor hack in the clusterProfiler function 'enrich.internal'
#'
#' @param gene some differentially expressed genes
#' @param organism by default 'human'
#' @param ont by default 'MF'
#'
#' @return some clusterProfiler data
#' @seealso \code{\link{clusterProfiler}}
#' @export
hpgl_enrich.internal <- function(gene, organism, pvalueCutoff=1, pAdjustMethod="BH",
                                 ont, minGSSize=2, qvalueCutoff=0.2, readable=FALSE, universe=NULL) {
    gene <- as.character(gene)
    class(gene) <- ont
    qExtID2TermID <- DOSE::EXTID2TERMID(gene, organism)
    qTermID <- unlist(qExtID2TermID)
    if (is.null(qTermID)) {
        return(NA)
    }
    ## Term ID -- query external ID association list.
    qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                   times=lapply(qExtID2TermID, length)),
                                   termID=qTermID)
    qExtID2TermID.df <- unique(qExtID2TermID.df)
    termID <- NULL ## to satisfy code tools
    qTermID2ExtID <- dlply(qExtID2TermID.df, .(termID),
                           .fun=function(i) as.character(i$extID))
    class(organism) <- ont
    extID <- DOSE::ALLEXTID(organism)
    if(!missing(universe)) {
        extID <- intersect(extID, universe)
    }
    qTermID2ExtID <- sapply(qTermID2ExtID, intersect, extID)
### The L.major ontologies are smaller, and so if the default (5)
### minGSSize is left in place, this comes up as null and therefore
### ends with the entire thing returning null.  I changed it to 2 for
### the moment.
    idx <- sapply(qTermID2ExtID, length) > minGSSize
    if (sum(idx) == 0) {
        return (NULL)
    }
    qTermID2ExtID <- qTermID2ExtID[idx]
    ## Term ID annotate query external ID
    qTermID <- unique(names(qTermID2ExtID))
    ## prepare parameter for hypergeometric test
    k <- sapply(qTermID2ExtID, length)
    k <- k[qTermID]
    class(qTermID) <- ont
    termID2ExtID <- DOSE::TERMID2EXTID(qTermID, organism)
    termID2ExtID <- sapply(termID2ExtID, intersect, extID)
    if (length(qTermID)== 1) {
        M <- nrow(termID2ExtID)
    } else {
        M <- sapply(termID2ExtID, length)
        M <- M[qTermID]
    }
    N <- rep(length(extID), length(M))
    ## n <- rep(length(gene), length(M)) ## those genes that have no annotation should drop.
    n <- rep(length(qExtID2TermID), length(M))
    args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                          numW=M,        ## White balls
                          numB=N-M,      ## Black balls
                          numDrawn=n)    ## balls drawn
    ## calcute pvalues based on hypergeometric model
    pvalues <- apply(args.df, 1, function(n)
                     phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
                     )
    ## gene ratio and background ratio
    GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )
    Over <- data.frame(ID=as.character(qTermID),
                       GeneRatio=GeneRatio,
                       BgRatio=BgRatio,
                       pvalue=pvalues)
    original_over = Over
    p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
    cat(sprintf("The minimum observed pvalue is: %f\n", min(pvalues)))
    qobj = try(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), silent=TRUE)
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       p.adjust = p.adj,
                       qvalue=qvalues,
                       geneID=geneID,
                       Count=k)
    class(qTermID) <- ont
    Description <- DOSE::TERM2NAME(qTermID, organism)

    if (length(qTermID) != length(Description)) {
        idx <- qTermID %in% names(tt)
        Over <- Over[idx,]
    }
    Over$Description <- Description
    nc <- ncol(Over)
    Over <- Over[, c(1,nc, 2:(nc-1))]
    Over <- Over[order(pvalues),]
    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    if (! any(is.na(Over$qvalue))) {
        Over <- Over[ Over$qvalue <= qvalueCutoff, ]
    }
    Over$ID <- as.character(Over$ID)
    Over$Description <- as.character(Over$Description)
    category <- as.character(Over$ID)
    ### On my computer this fails.
    ##    rownames(Over) <- category
    x <- new("enrichResult",
             result = Over,
             pvalueCutoff=pvalueCutoff,
             pAdjustMethod=pAdjustMethod,
             organism=as.character(organism),
             ontology=as.character(ont),
             gene=as.character(gene),
             geneInCategory=qTermID2ExtID[category]
             )
    if(readable) {
        x <- setReadable(x)
    }
    return (x)
}

#' A copy and paste of clusterProfiler's readGff
#' @export
##readGff <- function(gffFile, nrows = -1) {
##    cat("Reading ", gffFile, ": ", sep="")
##    gff <- read.table(gffFile, sep="\t", as.is=TRUE, quote="\"", fill=TRUE,
##                      header=FALSE, comment.char="#", nrows=nrows,
##                      colClasses=c("character", "character", "character", "integer",
##                          "integer", "character", "character", "character", "character"))
##    colnames(gff) = c("seqname", "source", "feature", "start", "end",
##                "score", "strand", "frame", "attributes")
##    cat("found", nrow(gff), "rows with classes:",
##        paste(sapply(gff, class), collapse=", "), "\n")
##    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
##    return(gff)
##}
##

## I had to hack getGffAttribution for gff files with bad encoding (yeast)
## Functions in this are not exported by clusterProfiler/topGO
hpgl_Gff2GeneTable <- function(gffFile, compress=TRUE, split="=") {
    ##gffFile="reference/gff/clbrener_8.1_complete_genes.gff"
    if (is.data.frame(gffFile)) {
        GeneID <- data.frame(GeneID = gffFile$ID)
        geneInfo <- gffFile
        geneInfo$start <- 1
        geneInfo$GeneID <- gffFile$ID
        geneInfo$GeneName <- gffFile$ID
        geneInfo$Locus <- gffFile$ID
        geneInfo$end <- geneInfo$width
        geneInfo$strand <- "+"
    } else {
        gff <- clusterProfiler:::readGff(gffFile)
        GeneID <- data.frame(GeneID=hpgl_getGffAttribution(gff$attributes, field="ID", split=split))
        geneInfo <- gff[gff$feature == "gene",]
        geneInfo <- geneInfo[, c("seqname", "start", "end", "strand", "attributes")]
        geneInfo$GeneID <- hpgl_getGffAttribution(geneInfo$attributes, field="ID", split=split)
        geneInfo$GeneName <- hpgl_getGffAttribution(geneInfo$attributes, field="gene", split=split)
        first_locus <- hpgl_getGffAttribution(geneInfo$attributes, field="locus_tag", split=split)
        first_sum <- sum(is.na(first_locus))
        second_locus <- NULL
        second_sum <- first_sum
        if (first_sum > (length(geneInfo$attributes) / 2)) { ## Using this to approximate whether locus_tag has a useful meaning.
            ## If more than 1/2 of the attributes have no locus tag, try using gene_id instead -- which is what yeast uses (btw).
            message("Trying to use gene_id insteady of locus_tag, because locus_tag is poorly defined.")
            second_locus <- hpgl_getGffAttribution(geneInfo$attributes, field="gene_id", split=split)
            second_sum <- sum(is.na(second_locus))
        }
        if (first_sum > second_sum) {
            geneInfo$Locus <- second_locus
        } else {
            geneInfo$Locus <- first_locus
        }
        geneInfo$GeneName[is.na(geneInfo$GeneName)] = "-"
        geneInfo <- geneInfo[, -5] ## drop "attributes" column.
    }
    ##GI2GeneID <- data.frame(GI=getGffAttribution(gff$attributes, field="GI"),
    ##GeneID=getGffAttribution(gff$attributes, field="GeneID"),
    ##Product=getGffAttribution(gff$attributes, field="product"))
    ##GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$GI),]
    ##GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$Gene),]
    ##geneTable <- merge(GI2GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- merge(GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- unique(geneTable)
    if (compress) {
        save(geneTable, file="geneTable.rda", compress="xz")
    } else {
        save(geneTable, file="geneTable.rda")
    }
    message("Gene Table file save in the working directory.")
}


hpgl_getGffAttribution <- function(x, field, attrsep=";", split='=') {
    s <- strsplit(x, split=attrsep, fixed=TRUE)
    sapply(s, function(atts) {
        atts <- gsub("^ ", "", x=atts, perl=TRUE)
        a <- strsplit(atts, split=split, fixed = TRUE)
        m <- match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv <- a[[m]][2]
        } else {
            b <- sapply(a, function(atts) {
                strsplit(atts[2], split=",", fixed = TRUE)
            })
            rv <- as.character(NA)
            sapply(b, function(atts) {
                secA <- strsplit(atts, split = ":", fixed = TRUE)
                m <- match(field, sapply(secA, "[", 1))
                if (!is.na(m)) {
                  rv <<- secA[[m]][2]
                }
            })
        }
        return(rv)
    })
}

mygroupgo <- function(gene, organism="human", ont="CC", level=2, readable=FALSE) {
    GOLevel <- clusterProfiler:::getGOLevel(ont, level)
    class(GOLevel) <- ont
    GO2ExtID <- TERMID2EXTID.GO(GOLevel, organism)
    geneID.list <- lapply(GO2ExtID, function(x) gene[gene %in% x])
    geneID <- sapply(geneID.list, function(i) paste(i, collapse = "/"))
    Count <- unlist(lapply(geneID.list, length))
    GeneRatio <- paste(Count, length(unique(unlist(gene))), sep = "/")
    Descriptions <- TERM2NAME(GOLevel)
    result = data.frame(ID=as.character(GOLevel), Description=Descriptions, Count=Count,
                        GeneRatio=GeneRatio, geneID=geneID)
    x <- new("groupGOResult", result=result, ontology=ont, level=level,
             organism=organism, gene=gene, geneInCategory=geneID.list)
    if (readable == TRUE) {
        x <- setReadable(x)
    }
    return(x)
}

# EOF
