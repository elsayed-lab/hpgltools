simple_reactome <- function(de_list, annotations, shortname="org.Sc.sgd") {
    reactomeEXTID2PATHID <- NULL
    require.auto("reactome.db")
    require.auto("annotations")
    xx <- as.list(reactomeEXTID2PATHID)
    ## ls("package:org.Sc.sgd.db")
    ## ok, I am starting to get it, I need to cross from sgd IDs to
    ## Entrez gene identifiers, from there I can pass to reactome IDs.
    ## Once I have that, I can pass along to reactomePA.
    ## The Sc.sgd.db package has sgdENTREZID
    shortname <- "org.Sc.sgd"
    x <- eval(parse(text=paste0(shortname, "ENTREZID")))
    mapped <- AnnotationDbi::mappedkeys(x)
    all_entrez <- as.list(x[mapped])
    entrezids <- all_entrez[(names(de_list))]
    ## Drop null elements
    entrezids <- entrezids[!unlist(lapply(entrezids, is.null))]
    ids <- as.character(entrezids)
    x <- ReactomePA::enrichPathway(gene=ids, pvalueCutoff=0.2, organism="yeast")
    summary(x)
    barplot(x, showCategory=8)
    DOSE::dotplot(x, showCategory=15)
    DOSE::enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
    DOSE::cnetplot(x, categorySize="pvalue", foldChange=de_list$log2FC)
    res <- clusterProfiler::compareCluster(ids, fun="enrichPathway")
    ## I think for this to work, I need the full set of DE values
    ## Then do the SGDID->Entrez conversion and reset the names() to that.
    ## y <- gsePathway(all_list, nPerm=100, minGSSize=2, pvalueCutoff=0.2, pAdjustMethod="BH", verbose=TRUE)
    ## res <- summary(y)
    ## head(res)
    ## enrichMap(y)
    ## gseaplot(y, geneSetID = "1280215")
    ## viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
}

## EOF
