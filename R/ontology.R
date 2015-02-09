#' Get a go term from ID
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = goterm("GO:0032559")
goterm = function(go="GO:0032559") {
    go = as.character(go)
    term = function(id) {
        value = try(as.character(Term(GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value = "not found"
        }             
        return(value)
    }
    go = mapply(term, go)
    return(go)
    ## count = 1
    ## for (go in goid) {
    ##     value = try(as.character(Term(GOTERM[go])), silent=TRUE)
    ##     if (class(value) == "try-error") {
    ##         value = "not found"
    ##     }             
    ##    goid[count] = value
    ##     count = count + 1
    ## }
    ## return(goid)
}

#' Get a go synonym from an ID
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosyn("GO:0032559")
gosyn = function(go) {
    go = as.character(go)
    syn = function(id) {
        value = try(as.character(AnnotationDbi::Synonym(GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value = "not found"
        }             
        return(value)
    }
    go = mapply(syn, go)
    return(go)
}

#' Get a go secondary ID from an id
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosec("GO:0032559")
gosec = function(go) {
    go = as.character(go)
    sec = function(id) {
        value = try(as.character(AnnotationDbi::Secondary(GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value = "not found"
        }             
        return(value)
    }
    go = mapply(sec, go)
    return(go)
}

#' Get a go long-form definition from an id
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosec("GO:0032559")
godef = function(go) {
    go = as.character(go)
    def = function(id) {
        value = try(as.character(AnnotationDbi::Definition(GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value = "not found"
        }             
        return(value)
    }
    go = mapply(def, go)
    return(go)
}

#' Get a go ontology name from an ID
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosec("GO:0032559")
goont = function(go) {
    go = as.character(go)
    ont = function(id) {
        value = try(as.character(AnnotationDbi::Ontology(GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value = "not found"
        }             
        return(value)
    }
    go = mapply(ont, go)
    return(go)
}


#' Get a go level approximation from an ID
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosec("GO:0032559")
golev = function(go) {
    go = as.character(go)
    level = 0
    go = "GO:0005874"
    while(class(try(as.character(AnnotationDbi::Ontology(GOTERM[[go]])), silent=FALSE)) != 'try-error') {
        if (ontology == "MF") {
            ancestors = GOMFANCESTOR[[go]]
        } else if (ontology == "BP") {
            ancestors = GOBPANCESTOR[[go]]
        } else if (ontology == "CC") {
            ancestors = GOCCANCESTOR[[go]]
        } else {
            ## There was an error
            message(paste("There was an error getting the ontology: ", as.character(id), sep=""))
            ancestors = "error"
        }
        print("Incrementing level")
        go = ancestors[1]
        level = level + 1
    }  ## End while
    return(level)
}
golevel = function(go) {
    mapply(golev, go)
}

#' Test a GO id to see if it is useful by my arbitrary definition of 'useful'
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#' 
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## text = gosec("GO:0032559")
gotest = function(go) {
    go = as.character(go)
    tst = function(id) {
        value = GOTERM[[go]]
        if (is.null(value)) {
            return(0)
        } else {
            return(1)
        }
    }
    go = mapply(tst, go)
    return(go)
}

#' Make a table of gene ontology information
#'
#' @param df a dataframe of ontology information.  This is intended to
#' be the output from goseq including information like
#' numbers/category, GOids, etc.
#' @param file a csv file to which to write the table
#' 
#' @return the ontology table with annotation information included
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#' 
#' @export
#' @examples
#' ## annotated_go = goseq_table(go_ids)
goseq_table = function(df, file=NULL) {
    ## Testing args:
    ## df=godata
    ## end testing
    keep_columns = c("goid", "over_pval","under_pval","numDEinCat","numInCat", "qvalue")
    colnames(df) = keep_columns
    df$good = gotest(df$goid)
##    df$good = mapply(gotest, df$goid)
    df = subset(df, good == 1)
    df = df[, keep_columns]
    df$term = goterm(df$goid)
##    df$term = mapply(goterm, df$goid)
##    df$term = as.character(df$term)
    df = subset(df, !is.null(term))
    df$syn = gosyn(df$goid)
##    df$syn = mapply(gosyn, df$goid)
##    df$syn = as.character(df$syn)
    df$sec = goseq(df$goid)
##    df$sec = mapply(gosec, df$goid)
##    df$sec = as.character(df$sec)
    df$ont = goont(df$goid)
##    df$ont = mapply(goont, df$goid)
##    df$ont = as.character(df$ont)
    df = subset(df, !is.null(ont))
######    df$level = mapply(golevel, df$goid)
    df$def = godef(df$goid)
##    df$def = mapply(godef, df$goid)
##    df$def = as.character(df$def)
    if (!is.null(file)) {
        write.csv(df, file=file)
    }
    return(df)
}

#' Perform a simplified goseq analysis
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param lengths the length of each gene with an ID in de_genes
#' @param goids a list of ontology accessions to gene accessions
#' 
#' @return a big list including: the pwd:pwf function, alldata:the godata dataframe, pvalue_histogram:p-value histograms, godata_interesting:the ontology information of the enhanced groups, term_table:the goterms with some information about them, mf_subset:a plot of the MF enhanced groups, mfp_plot:the pvalues of the MF group, bp_subset:a plot of the BP enhanced groups, bpp_plot, cc_subset, and ccp_plot
#' @seealso \code{\link{goseq}} and \code{\link{clusterProfiler}}
#' @export
simple_goseq = function(de_genes, lengths=NULL, goids=NULL, adjust=0.1, pvalue=0.1, qvalue=0.1, method="Wallenius") {
    message("simple_goseq() makes some pretty hard assumptions about the data it is fed:")
    message("It requires 2 tables, one of GOids which must have columns (gene)ID and GO(category)")
    message("The other table is of gene lengths with columns (gene)ID and (gene)width.")
    message("Other columns are fine, but ignored.")
    if (is.null(de_genes$ID)) {
        de_genes$ID = make.names(rownames(de_genes), unique=TRUE)
    }
    de_genes$DE = 1
    adjust=NULL
    de_table = de_genes[,c("ID","DE")]
    length_table = lengths[,c("ID","width")]
##    de_table = merge(de_table, length_table, by="ID")
    de_table = merge(de_table, length_table, by="ID", all.y=TRUE)    
    de_table[is.na(de_table)] = 0
    rownames(de_table) = make.names(de_table$ID, unique=TRUE)
    de_vector = as.vector(de_table$DE)
    names(de_vector) = rownames(de_table)
    width_vector = as.vector(de_table$width)
    names(width_vector) = de_table$ID
    pwf = goseq::nullp(DEgenes=de_vector, bias.data=width_vector, plot.fit=TRUE)
    pwf_plot = recordPlot()
##    godata = goseq(pwf, gene2cat=goids, method='Wallenius')
    colnames(goids) = c("ID", "GO")
    godata = goseq::goseq(pwf, gene2cat=goids, use_genes_without_cat=TRUE, method=method)
    goseq_p = try(hpgltools::hpgl_histogram(godata$over_represented_pvalue, bins=20))
    goseq_p_second = sort(unique(table(goseq_p$data)), decreasing=TRUE)[2]
    ## Set the y scale to 2x the second highest number
    ## (assuming always that the highest is a p-value of 1)
    goseq_y_limit = goseq_p_second * 2
    goseq_p = goseq_p + scale_y_continuous(limits=c(0, goseq_y_limit))
    message("Calculating q-values")
    qdata = godata$over_represented_pvalue
    qdata[qdata > 1] = 1 ## For scientific numbers which are 1.0000E+00 it might evaluate to 1.0000000000000001
    qdata = qvalue::qvalue(qdata)
    godata = cbind(godata, qdata$qvalues)
    colnames(godata) = c("category","over_represented_pvalue","under_represented_pvalue","numDEInCat","numInCat","qvalue")
    message("Filling godata table with term information, this takes a while.")
    godata$ont = goont(godata$category)
    godata$term = goterm(godata$category)
    if (!is.null(adjust)) {
        godata_interesting = subset(godata, p.adjust(godata$over_represented_pvalue, method=method) < adjust)
        adjust_method=method
        if (dim(godata_interesting)[1] == 0) {
            message(paste("There are no genes with an adjusted pvalue < ", adjust, " using method: ", method, ".", sep=""))
            message(sprintf("Providing genes with an un-adjusted pvalue < %s", pvalue))
            godata_interesting = subset(godata, godata$over_represented_pvalue < pvalue)
            adjust_method="none"
        }
    } else {
        godata_interesting = subset(godata, godata$over_represented_pvalue < pvalue)
        adjust_method="none"
    }
    ##    goterms = goseq_table(godata_interesting)
    message("Making pvalue plots for the ontologies.")
    pvalue_plots = goseq_pval_plots(godata)
    mf_subset = subset(godata, ont == "MF")
    bp_subset = subset(godata, ont == "BP")
    cc_subset = subset(godata, ont == "CC")
    mf_interesting = subset(godata_interesting, ont == "MF")
    rownames(mf_interesting) = mf_interesting$category
    mf_interesting = mf_interesting[,c("ont","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    bp_interesting = subset(godata_interesting, ont == "BP")
    rownames(bp_interesting) = bp_interesting$category
    bp_interesting = bp_interesting[,c("ont","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]    
    cc_interesting = subset(godata_interesting, ont == "CC")
    cc_interesting = cc_interesting[,c("ont","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]    
    return_list = list(pwf=pwf, pwf_plot=pwf_plot,
        alldata=godata, pvalue_histogram=goseq_p,
        godata_interesting=godata_interesting,
        mf_interesting=mf_interesting, bp_interesting=bp_interesting, cc_interesting=cc_interesting,
        adjust_method=adjust_method,
        ##term_table=goterms,
        mf_subset=mf_subset, mfp_plot=pvalue_plots$mfp_plot,
        bp_subset=bp_subset, bpp_plot=pvalue_plots$bpp_plot,
        cc_subset=cc_subset, ccp_plot=pvalue_plots$ccp_plot,
        qdata=qdata)
    return(return_list)
}

#' Make a pvalue plot from topgo data
#'
#' @param topgo_data some data from topgo!
#' 
#' @return a plot!
#' @seealso \code{\link{goseq}}
#' @export
topgo_pval_plot = function(topgo, wrapped_width=20, cutoff=0.1, n=12, type="fisher") {
    mf_newdf = topgo$tables$mf[,c("GO.ID", "Term", "Annotated","Significant",type)]
    mf_newdf$term = as.character(lapply(strwrap(mf_newdf$Term, wrapped_width, simplify=F), paste, collapse="\n"))
    mf_newdf$pvalue = as.numeric(mf_newdf[[type]])
    mf_newdf = subset(mf_newdf, get(type) < cutoff)
    mf_newdf = mf_newdf[order(mf_newdf$pvalue, mf_newdf[[type]]),]
    mf_newdf = head(mf_newdf, n=n)
    mf_newdf$score = mf_newdf$Significant / mf_newdf$Annotated    
    mf_pval_plot = pval_plot(mf_newdf, ontology="MF")

    bp_newdf = topgo$tables$bp[,c("GO.ID", "Term", "Annotated","Significant",type)]
    bp_newdf$term = as.character(lapply(strwrap(bp_newdf$Term, wrapped_width, simplify=F), paste, collapse="\n"))
    bp_newdf$pvalue = as.numeric(bp_newdf[[type]])
    bp_newdf = subset(bp_newdf, get(type) < cutoff)
    bp_newdf = bp_newdf[order(bp_newdf$pvalue, bp_newdf[[type]]),]
    bp_newdf = head(bp_newdf, n=n)
    bp_newdf$score = bp_newdf$Significant / bp_newdf$Annotated
    bp_pval_plot = pval_plot(bp_newdf, ontology="MF")

    cc_newdf = topgo$tables$cc[,c("GO.ID", "Term", "Annotated","Significant",type)]
    cc_newdf$term = as.character(lapply(strwrap(cc_newdf$Term, wrapped_width, simplify=F), paste, collapse="\n"))
    cc_newdf$pvalue = as.numeric(cc_newdf[[type]])
    cc_newdf = subset(cc_newdf, get(type) < cutoff)
    cc_newdf = cc_newdf[order(cc_newdf$pvalue, cc_newdf[[type]]),]
    cc_newdf = head(cc_newdf, n=n)
    cc_newdf$score = cc_newdf$Significant / cc_newdf$Annotated    
    cc_pval_plot = pval_plot(cc_newdf, ontology="CC")
    
    pval_plots = list(MF=mf_pval_plot, BP=bp_pval_plot, CC=cc_pval_plot)
    return(pval_plots)
}

#' Make a pvalue plot from goseq data
#'
#' @param topgo_data some data from goseq!
#' 
#' @return a plot!
#' @seealso \code{\link{goseq}}
#' @export
goseq_pval_plots = function(goterms, wrapped_width=20, cutoff=0.1, n=10) {
    plotting_mf = subset(goterms, complete.cases(goterms))
    plotting_mf$score = plotting_mf$numDEInCat / plotting_mf$numInCat    
    plotting_mf = subset(plotting_mf, ont == "MF")
    plotting_mf = subset(plotting_mf, term != "NULL")
    plotting_mf = subset(plotting_mf, over_represented_pvalue <= 0.1)
    plotting_mf = subset(plotting_mf, numInCat > 10)    
    plotting_mf = plotting_mf[order(plotting_mf$over_represented_pvalue),]
    plotting_mf = head(plotting_mf, n=n)
    plotting_mf = plotting_mf[,c("term","over_represented_pvalue","score")]
    colnames(plotting_mf) = c("term","pvalue","score")
    mf_pval_plot = pval_plot(plotting_mf, ontology="MF")

    plotting_bp = subset(goterms, complete.cases(goterms))
    plotting_bp$score = plotting_bp$numDEInCat / plotting_bp$numInCat
    plotting_bp = subset(plotting_bp, ont == "BP")
    plotting_bp = subset(plotting_bp, term != "NULL")
    plotting_bp = subset(plotting_bp, over_represented_pvalue <= 0.1)
    plotting_bp = subset(plotting_bp, numInCat > 10)    
    plotting_bp = plotting_bp[order(plotting_bp$over_represented_pvalue),]
    plotting_bp = head(plotting_bp, n=n)
    plotting_bp = plotting_bp[,c("term","over_represented_pvalue","score")]
    colnames(plotting_bp) = c("term","pvalue","score")
    bp_pval_plot = pval_plot(plotting_bp, ontology="BP")

    plotting_cc = subset(goterms, complete.cases(goterms))
    plotting_cc$score = plotting_cc$numDEInCat / plotting_cc$numInCat
    plotting_cc = subset(plotting_cc, ont == "CC")
    plotting_cc = subset(plotting_cc, term != "NULL")
    plotting_cc = subset(plotting_cc, over_represented_pvalue <= 0.1)
    plotting_cc = subset(plotting_cc, numInCat > 10)    
    plotting_cc = plotting_cc[order(plotting_cc$over_represented_pvalue),]
    plotting_cc = head(plotting_cc, n=n)
    plotting_cc = plotting_cc[,c("term","over_represented_pvalue","score")]
    colnames(plotting_cc) = c("term","pvalue","score")
    cc_pval_plot = pval_plot(plotting_cc, ontology="CC")
    
    pval_plots = list(mfp_plot=mf_pval_plot, bpp_plot=bp_pval_plot, ccp_plot=cc_pval_plot,
                      mf_subset=plotting_mf, bp_subset=plotting_bp, cc_subset=plotting_cc)
    return(pval_plots)    
}

#' Make a pvalue plot from a df of IDs, scores, and p-values
#'
#' @param df some data from topgo/goseq/clusterprofiler
#' 
#' @return a plot!
#' @seealso \code{\link{goseq}}
#' @export
pval_plot = function(df, ontology="MF") {
    y_name = paste("Enriched ", ontology, " categories.", sep="")
    pvalue_plot = ggplot2::ggplot(df, aes(term, score)) +
        geom_bar(stat="identity") +
        coord_flip() +
        scale_x_discrete(name=y_name) +
        aes(fill=pvalue) +
        scale_fill_continuous(low="red", high="blue") +
        theme(text=element_text(size=10))
    return(pvalue_plot)
}


hpgl_topdiffgenes = function(scores, df=de_genes, direction="up") {
    ## Testing parameters
    ##scores = pvals
    ##df = epi_cl14clbr_high
    ## Here is the original topDiffGenes
    ## topDiffGenes <- function(allScore) {
    ##   return(allScore < 0.01)
    ##}
    ## my version of this will expect a limma result table from which I will extract the entries with low p-values
    ## and logFCs which are high or low
    quartiles = summary(df)
}

#' Perform ontology searches of the output from limma
#'
#' @param limma_out a list of topTables comprising limma outputs
#' @param n a number of genes at the top/bottom to search
#' @param z a number of standard deviations to search
#'
#' @export
limma_ontology = function(limma_out, gene_lengths=NULL, goids=NULL, n=NULL, z=NULL, overwrite=FALSE, goid_map="reference/go/id2go.map", goids_df=NULL, do_goseq=TRUE, do_cluster=TRUE, do_topgo=TRUE, do_trees=FALSE, workbook="excel/ontology.xls", csv=TRUE, excel=FALSE) {
    message("This function expects a list of limma contrast tables and some annotation information.")
    message("The annotation information would be gene lengths and ontology ids")
    if (is.null(n) & is.null(z)) {
        z = 1
    }

    testdir = dirname(workbook)
    if (isTRUE(excel) | isTRUE(csv)) {
        if (!file.exists(testdir)) {
            dir.create(testdir)
            message(paste("Creating directory: ", testdir, "for writing excel/csv data.", sep=""))
        }
    }
    
    output = list()
    for (c in 1:length(limma_out)) {
        datum = limma_out[[c]]
        comparison = names(limma_out[c])
        if (is.null(n)) {
            out_summary = summary(datum$logFC)
            out_mad = mad(datum$logFC, na.rm=TRUE)
            up_median_dist = out_summary["Median"] + (out_mad * z)
            down_median_dist = out_summary["Median"] - (out_mad * z)
            up_genes = subset(datum, logFC >= up_median_dist)
            down_genes = subset(datum, logFC <= down_median_dist)
        } else if (is.null(z)) {
            upranked = datum[order(datum$logFC, decreasing=TRUE),]
            up_genes = head(upranked, n=n)
            down_genes = tail(upranked, n=n)
        }
        goseq_up_ontology = goseq_up_trees = goseq_down_ontology = goseq_down_trees = NULL
        cluster_up_ontology = cluster_up_trees = cluster_down_ontology = cluster_down_trees = NULL
        topgo_up_ontology = topgo_up_trees = topgo_down_ontology = topgo_down_trees = NULL        
        if (isTRUE(do_goseq)) {
            goseq_up_ontology = simple_goseq(up_genes, lengths=gene_lengths, goids=goids)
            goseq_down_ontology = simple_goseq(down_genes, lengths=gene_lengths, goids=goids)
            if (isTRUE(do_trees)) {
                goseq_up_trees = try(goseq_trees(up_genes, goseq_up_ontology, goid_map=goid_map, goids_df=goids, overwrite=overwrite))
                goseq_down_trees = try(goseq_trees(down_genes, goseq_down_ontology, goid_map=goid_map, goids_df=goids))
            }
            if (isTRUE(excel)) {
                sheetname = paste(comparison, "_goseq_mf_up", sep="")
                try(write_xls(data=goseq_up_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_goseq_bp_up", sep="")
                try(write_xls(data=goseq_up_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_goseq_cc_up", sep="")
                try(write_xls(data=goseq_up_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_goseq_mf_down", sep="")
                try(write_xls(data=goseq_down_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_goseq_bp_down", sep="")
                try(write_xls(data=goseq_down_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_goseq_cc_down", sep="")
                try(write_xls(data=goseq_down_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
            }
            if (isTRUE(csv)) {
                csv_base = gsub(".xls$", "", workbook)
                write.csv(goseq_up_ontology$mf_interesting, file=paste(comparison, csv_base, "_goseq_mf_up.csv", sep=""))
                write.csv(goseq_up_ontology$bp_interesting, file=paste(comparison, csv_base, "_goseq_bp_up.csv", sep=""))
                write.csv(goseq_up_ontology$cc_interesting, file=paste(comparison, csv_base, "_goseq_cc_up.csv", sep=""))
                write.csv(goseq_down_ontology$mf_interesting, file=paste(comparison, csv_base, "_goseq_mf_down.csv", sep=""))
                write.csv(goseq_down_ontology$bp_interesting, file=paste(comparison, csv_base, "_goseq_bp_down.csv", sep=""))
                write.csv(goseq_down_ontology$cc_interesting, file=paste(comparison, csv_base, "_goseq_cc_down.csv", sep=""))
            }
        }
        if (isTRUE(do_cluster)) {
            cluster_up_ontology = simple_clusterprofiler(up_genes, goids=goids, gff=goids)
            cluster_down_ontology = simple_clusterprofiler(down_genes, goids=goids, gff=goids)
            if (isTRUE(do_trees)) {
                cluster_up_trees = try(cluster_trees(up_genes, cluster_up_ontology, goid_map=goid_map, goids_df=goids))
                cluster_down_trees = try(cluster_trees(down_genes, cluster_down_ontology, goid_map=goid_map, goids_df=goids))
            }
            if (isTRUE(excel)) {
                sheetname = paste(comparison, "_cluster_mf_up", sep="")
                try(write_xls(data=cluster_up_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_cluster_bp_up", sep="")
                try(write_xls(data=cluster_up_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_cluster_cc_up", sep="")
                try(write_xls(data=cluster_up_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_cluster_mf_down", sep="")
                try(write_xls(data=cluster_down_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_cluster_bp_down", sep="")
                try(write_xls(data=cluster_down_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_cluster_cc_down", sep="")
                try(write_xls(data=cluster_down_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
            }
            if (isTRUE(csv)) {
                csv_base = gsub(".xls$", "", workbook)
                write.csv(cluster_up_ontology$mf_interesting, file=paste(comparison, csv_base, "_cluster_mf_up.csv", sep=""))
                write.csv(cluster_up_ontology$bp_interesting, file=paste(comparison, csv_base, "_cluster_bp_up.csv", sep=""))
                write.csv(cluster_up_ontology$cc_interesting, file=paste(comparison, csv_base, "_cluster_cc_up.csv", sep=""))
                write.csv(cluster_down_ontology$mf_interesting, file=paste(comparison, csv_base, "_cluster_mf_down.csv", sep=""))
                write.csv(cluster_down_ontology$bp_interesting, file=paste(comparison, csv_base, "_cluster_bp_down.csv", sep=""))
                write.csv(cluster_down_ontology$cc_interesting, file=paste(comparison, csv_base, "_cluster_cc_down.csv", sep=""))
            }
        }
        if (isTRUE(do_topgo)) {
            topgo_up_ontology = simple_topgo(up_genes, goid_map=goid_map, goids_df=goids)
            topgo_down_ontology = simple_topgo(down_genes, goid_map=goid_map, goids_df=goids)
            if (isTRUE(do_trees)) {
                topgo_up_trees = try(topgo_trees(topgo_up_ontology))
                topgo_down_trees = try(topgo_trees(topgo_down_ontology))
            }
            if (isTRUE(excel)) {
                sheetname = paste(comparison, "_topgo_mf_up", sep="")
                try(write_xls(data=topgo_up_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_bp_up", sep="")
                try(write_xls(data=topgo_up_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_cc_up", sep="")
                try(write_xls(data=topgo_up_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_mf_down", sep="")
                try(write_xls(data=topgo_down_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_bp_down", sep="")
                try(write_xls(data=topgo_down_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_cc_down", sep="")
                try(write_xls(data=topgo_down_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
            }
            if (isTRUE(csv)) {
                csv_base = gsub(".xls$", "", workbook)
                write.csv(topgo_up_ontology$mf_interesting, file=paste(comparison, csv_base, "_topgo_mf_up.csv", sep=""))
                write.csv(topgo_up_ontology$bp_interesting, file=paste(comparison, csv_base, "_topgo_bp_up.csv", sep=""))
                write.csv(topgo_up_ontology$cc_interesting, file=paste(comparison, csv_base, "_topgo_cc_up.csv", sep=""))
                write.csv(topgo_down_ontology$mf_interesting, file=paste(comparison, csv_base, "_topgo_mf_down.csv", sep=""))
                write.csv(topgo_down_ontology$bp_interesting, file=paste(comparison, csv_base, "_topgo_bp_down.csv", sep=""))
                write.csv(topgo_down_ontology$cc_interesting, file=paste(comparison, csv_base, "_topgo_cc_down.csv", sep=""))
            }
        }
        c_data = list(up_goseq=goseq_up_ontology, down_goseq=goseq_down_ontology,
            up_cluster=cluster_up_ontology, down_cluster=cluster_down_ontology,
            up_topgo=topgo_up_ontology, down_topgo=topgo_down_ontology,
            up_goseqtrees=goseq_up_trees, down_goseqtrees=goseq_down_trees,
            up_clustertrees=cluster_up_trees, down_clustertrees=cluster_down_trees,
            up_topgotrees=topgo_up_trees, down_topgotrees=topgo_down_trees)
            output[[c]] = c_data
    }
    names(output) = names(limma_out)
    return(output)
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
topDiffGenes <- function(allScore) { return(allScore < 0.01) }

#' Perform a simplified topgo analysis
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param goid_map a file containing mappings of genes to goids in the format expected by topgo
#' 
#' @return a big list including the various outputs from topgo
#' @export
simple_topgo = function(de_genes, goid_map="reference/go/id2go.map", goids_df=NULL, pvals=NULL, limitby="fisher", limit=0.1, signodes=100, sigforall=TRUE, numchar=300, selector="topDiffGenes", overwrite=FALSE) {
### Some neat ideas from the topGO documentation:
### geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")
### A variant of these operations make it possible to give topGO scores so that
### a larger array of tests may be performed
### x <- topDiffGenes(geneList)
### sum(x) ## the number of selected genes
### If we do something like above to give scores to all the 'DEgenes', then we set up the GOdata object like this:
### mf_GOdata = new("topGOdata", description="something", ontology="BP", allGenes = entire_geneList, geneSel=topDiffGenes, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=2)
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO = topGO::readMappings(file=goid_map)
    annotated_genes = names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID = make.names(rownames(de_genes), unique=TRUE)
    }
    ##    interesting_genes = factor(as.integer(annotated_genes %in% de_genes$ID))
    interesting_genes = factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) = annotated_genes
    if (is.null(pvals)) {
        mf_GOdata = new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
    } else {
        mf_GOdata = new("topGOdata", description="MF", ontology="MF", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", description="BP", ontology="BP", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", description="CC", ontology="CC", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    }
    test_stat = new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
    mf_fisher_result = topGO::getSigGroups(mf_GOdata, test_stat)
    bp_fisher_result = topGO::getSigGroups(bp_GOdata, test_stat)
    cc_fisher_result = topGO::getSigGroups(cc_GOdata, test_stat)
    test_stat = new("classicScore", testStatistic=GOKSTest, name="KS tests")
    mf_ks_result = topGO::getSigGroups(mf_GOdata, test_stat)
    bp_ks_result = topGO::getSigGroups(bp_GOdata, test_stat)
    cc_ks_result = topGO::getSigGroups(cc_GOdata, test_stat)
    test_stat = new("elimScore", testStatistic=GOKSTest, name="Fisher test", cutOff=0.01)
    mf_el_result = topGO::getSigGroups(mf_GOdata, test_stat)
    bp_el_result = topGO::getSigGroups(bp_GOdata, test_stat)
    cc_el_result = topGO::getSigGroups(cc_GOdata, test_stat)
    test_stat = new("weightCount", testStatistic=GOFisherTest, name="Fisher test", sigRatio="ratio")
    mf_weight_result = topGO::getSigGroups(mf_GOdata, test_stat)
    bp_weight_result = topGO::getSigGroups(bp_GOdata, test_stat)
    cc_weight_result = topGO::getSigGroups(cc_GOdata, test_stat)


    mf_fisher_pdist = try(hpgltools::hpgl_histogram(mf_fisher_result@score, bins=20))
    mf_ks_pdist = try(hpgltools::hpgl_histogram(mf_ks_result@score, bins=20))
    mf_el_pdist = try(hpgltools::hpgl_histogram(mf_el_result@score, bins=20))
    mf_weight_pdist = try(hpgltools::hpgl_histogram(mf_weight_result@score, bins=20))
    bp_fisher_pdist = try(hpgltools::hpgl_histogram(bp_fisher_result@score, bins=20))
    bp_ks_pdist = try(hpgltools::hpgl_histogram(bp_ks_result@score, bins=20))
    bp_el_pdist = try(hpgltools::hpgl_histogram(bp_el_result@score, bins=20))
    bp_weight_pdist = try(hpgltools::hpgl_histogram(bp_weight_result@score, bins=20))
    cc_fisher_pdist = try(hpgltools::hpgl_histogram(cc_fisher_result@score, bins=20))
    cc_ks_pdist = try(hpgltools::hpgl_histogram(cc_ks_result@score, bins=20))
    cc_el_pdist = try(hpgltools::hpgl_histogram(cc_el_result@score, bins=20))
    cc_weight_pdist = try(hpgltools::hpgl_histogram(cc_weight_result@score, bins=20))
    p_dists = list(mf_fisher=mf_fisher_pdist, bp_fisher=bp_fisher_pdist, cc_fisher=cc_fisher_pdist,
        mf_ks=mf_ks_pdist, bp_ks=bp_ks_pdist, cc_ks=cc_ks_pdist,
        mf_el=mf_el_pdist, bp_el=bp_el_pdist, cc_el=cc_el_pdist,
        mf_weight=mf_weight_pdist, bp_weight=bp_weight_pdist, cc_weight=cc_weight_pdist)
        
    results = list(mf_godata=mf_GOdata, bp_godata=bp_GOdata, cc_godata=cc_GOdata,
        mf_fisher=mf_fisher_result, bp_fisher=bp_fisher_result, cc_fisher=cc_fisher_result,
        mf_ks=mf_ks_result, bp_ks=bp_ks_result, cc_ks=cc_ks_result,
        mf_el=mf_el_result, bp_el=bp_el_result, cc_el=cc_el_result,
        mf_weight=mf_weight_result, bp_weight=bp_weight_result, cc_weight=cc_weight_result)

    tables = try(topgo_tables(results, limitby=limitby, limit=limit))

    mf_first_density = bp_first_density = cc_first_density = NULL
    if (class(tables$mf) != 'try-error') {
        mf_first_group = tables$mf[1, "GO.ID"]
        mf_first_density = try(hpgl_GroupDensity(mf_GOdata, mf_first_group, ranks=TRUE))
    }
    if (class(tables$bp) != 'try-error') {
        bp_first_group = tables$bp[1, "GO.ID"]
        bp_first_density = try(hpgl_GroupDensity(bp_GOdata, bp_first_group, ranks=TRUE))
    }
    if(class(tables$cc) != 'try-error') {
        cc_first_group = tables$cc[1, "GO.ID"]
        cc_first_density = try(hpgl_GroupDensity(cc_GOdata, cc_first_group, ranks=TRUE))
    }
    first_densities = list(mf=mf_first_density, bp=bp_first_density, cc=cc_first_density)
    
    information = list(
        mf_godata=mf_GOdata, bp_godata=bp_GOdata, cc_godata=cc_GOdata,
        results=results, tables=tables, first_densities=first_densities,
        pdists=p_dists)
    return(information)
}

#' Make tables out of topGO data
#' @export
topgo_tables = function(result, limit=0.01, limitby="fisher", numchar=300, orderby="classic", ranksof="classic") {
    ## The following if statement could be replaced by get(limitby)
    ## But I am leaving it as a way to ensure that no shenanigans ensue
    mf_allRes = bp_allRes = cc_allRes = mf_interesting = bp_interesting = cc_interesting = NULL
    if (limitby == "fisher") {
        mf_siglist = names(which(result$mf_fisher@score <= limit))        
        bp_siglist = names(which(result$bp_fisher@score <= limit))
        cc_siglist = names(which(result$bp_fisher@score <= limit))
    } else if (limitby == "KS") {
        mf_siglist = names(which(result$mf_ks@score <= limit))
        bp_siglist = names(which(result$bp_ks@score <= limit))
        cc_siglist = names(which(result$bp_ks@score <= limit))
    } else if (limitby == "EL") {
        mf_siglist = names(which(result$mf_el@score <= limit))
        bp_siglist = names(which(result$bp_el@score <= limit))
        cc_siglist = names(which(result$bp_el@score <= limit))
    } else if (limitby == "weight") {
        mf_siglist = names(which(result$mf_weight@score <= limit))
        bp_siglist = names(which(result$bp_weight@score <= limit))
        cc_siglist = names(which(result$bp_weight@score <= limit))
    } else {
        stop("I can only limit by: fisher, KS, EL, or weight.")
    }
    mf_topnodes = length(mf_siglist)
    if (mf_topnodes > 0) {
        mf_allRes = try(topGO::GenTable(result$mf_godata, classic=result$mf_fisher, KS=result$mf_ks,
            EL=result$mf_el, weight=result$mf_weight, orderBy=orderby,
            ranksOf=ranksof, topNodes=mf_topnodes, numChar=numchar))
        if (class(mf_allRes) != 'try-error') {
            mf_qvalues = as.data.frame(qvalue::qvalue(topGO::score(result$mf_fisher))$qvalues)
            mf_allRes = merge(mf_allRes, mf_qvalues, by.x="GO.ID", by.y="row.names")
            mf_allRes$classic = as.numeric(mf_allRes$classic)
            mf_allRes = mf_allRes[with(mf_allRes, order(classic)), ]
            colnames(mf_allRes) = c("GO.ID","Term","Annotated","Significant","Expected","fisher","KS","EL","weight","qvalue")
            mf_interesting = subset(mf_allRes, get(limitby) <= limit)
            rownames(mf_interesting) = NULL
            mf_interesting$ont = "MF"
            mf_interesting = mf_interesting[,c("GO.ID","ont","Annotated","Significant","Expected","fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    bp_topnodes = length(bp_siglist)
    if (bp_topnodes > 0) {
        bp_allRes = try(topGO::GenTable(result$bp_godata, classic=result$bp_fisher, KS=result$bp_ks,
            EL=result$bp_el, weight=result$bp_weight, orderBy=orderby,
            ranksOf=ranksof, topNodes=bp_topnodes, numChar=numchar))
        if (class(bp_allRes) != 'try-error') {
            bp_qvalues = as.data.frame(qvalue::qvalue(topGO::score(result$bp_fisher))$qvalues)
            bp_allRes = merge(bp_allRes, bp_qvalues, by.x="GO.ID", by.y="row.names", all.x=TRUE)
            bp_allRes$classic = as.numeric(bp_allRes$classic)
            bp_allRes = bp_allRes[with(bp_allRes, order(classic)), ]
            colnames(bp_allRes) = c("GO.ID","Term","Annotated","Significant","Expected","fisher","KS","EL","weight","qvalue")
            bp_interesting = subset(bp_allRes, get(limitby) <= limit)
            rownames(bp_interesting) = NULL
            bp_interesting$ont = "BP"
            bp_interesting = bp_interesting[,c("GO.ID","ont","Annotated","Significant","Expected","fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    cc_topnodes = length(cc_siglist)
    if (cc_topnodes > 0) {
        cc_allRes = try(topGO::GenTable(result$cc_godata, classic=result$cc_fisher, KS=result$cc_ks,
            EL=result$cc_el, weight=result$cc_weight, orderBy=orderby,
            ranksOf=ranksof, topNodes=cc_topnodes, numChar=numchar))
        if (class(cc_allRes) != 'try-error') {
            cc_qvalues = as.data.frame(qvalue::qvalue(topGO::score(result$cc_fisher))$qvalues)
            cc_allRes = merge(cc_allRes, cc_qvalues, by.x="GO.ID", by.y="row.names")
            cc_allRes$classic = as.numeric(cc_allRes$classic)
            cc_allRes = cc_allRes[with(cc_allRes, order(classic)), ]
            colnames(cc_allRes) = c("GO.ID","Term","Annotated","Significant","Expected","fisher","KS","EL","weight","qvalue")
            cc_interesting = subset(cc_allRes, get(limitby) <= limit)
            rownames(cc_interesting) = NULL
            cc_interesting$ont = "CC"
            cc_interesting = cc_interesting[,c("GO.ID","ont","Annotated","Significant","Expected","fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    tables = list(mf=mf_allRes, bp=bp_allRes, cc=cc_allRes, mf_interesting=mf_interesting, bp_interesting=bp_interesting, cc_interesting=cc_interesting)
    return(tables)
}

#' Print trees from topGO
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param goid_map a file containing mappings of genes to goids in the format expected by topgo
#' 
#' @return a big list including the various outputs from topgo
#' @export
topgo_trees = function(tg, score_limit=0.01, sigforall=TRUE, do_mf_fisher_tree=TRUE, do_bp_fisher_tree=TRUE, do_cc_fisher_tree=TRUE, do_mf_ks_tree=FALSE, do_bp_ks_tree=FALSE, do_cc_ks_tree=FALSE, do_mf_el_tree=FALSE, do_bp_el_tree=FALSE, do_cc_el_tree=FALSE, do_mf_weight_tree=FALSE, do_bp_weight_tree=FALSE, do_cc_weight_tree=FALSE) {
    mf_fisher_nodes = mf_fisher_tree = NULL
    if (do_mf_fisher_tree) {
        included = length(which(topGO::score(tg$results$mf_fisher) <= score_limit))
        mf_fisher_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_fisher), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_fisher_nodes)[1] != 'try-error') {
            mf_fisher_tree = try(recordPlot())
        }
    }
    bp_fisher_nodes = bp_fisher_tree = NULL
    if (do_bp_fisher_tree) {
        included = length(which(topGO::score(tg$results$bp_fisher) <= score_limit))
        bp_fisher_nodes = try(suppressWarnings(showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_fisher), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_fisher_nodes)[1] != 'try-error') {
            bp_fisher_tree = try(recordPlot())
        }
    }
    cc_fisher_nodes = cc_fisher_tree = NULL
    if (do_cc_fisher_tree) {
        included = length(which(topGO::score(tg$results$cc_fisher) <= score_limit))        
        cc_fisher_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_fisher), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_fisher_nodes)[1] != 'try-error') {
            cc_fisher_tree = try(recordPlot())
        }
    }
    mf_ks_nodes = mf_ks_tree = NULL
    if (do_mf_ks_tree) {
        included = length(which(topGO::score(tg$results$mf_ks) <= score_limit))        
        mf_ks_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_ks), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_ks_nodes)[1] != 'try-error') {
            mf_ks_tree = try(recordPlot())
        }
    }
    bp_ks_nodes = bp_ks_tree = NULL
    if (do_bp_ks_tree) {
        included = length(which(topGO::score(tg$results$bp_ks) <= score_limit))
        bp_ks_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_ks), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_ks_nodes)[1] != 'try-error') {
            bp_ks_tree = try(recordPlot())
        }
    }
    cc_ks_nodes = cc_ks_tree = NULL
    if (do_cc_ks_tree) {
        included = length(which(topGO::score(tg$results$cc_ks) <= score_limit))        
        cc_ks_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_ks), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_ks_nodes)[1] != 'try-error') {
            cc_ks_tree = try(recordPlot())
        }
    }
    mf_el_nodes = mf_el_tree = NULL
    if (do_mf_el_tree) {
        included = length(which(topGO::score(tg$results$mf_el) <= score_limit))        
        mf_el_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_el), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_el_nodes)[1] != 'try-error') {
            mf_el_tree = try(recordPlot())
        }
    }
    bp_el_nodes = bp_el_tree = NULL
    if (do_bp_el_tree) {
        included = length(which(topGO::score(tg$results$bp_el) <= score_limit))                
        bp_el_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_el), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_el_nodes)[1] != 'try-error') {
            bp_el_tree = try(recordPlot())
        }
    }
    cc_el_nodes = cc_el_tree = NULL
    if (do_cc_el_tree) {
        included = length(which(topGO::score(tg$results$cc_el) <= score_limit))                
        cc_el_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_el), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_el_nodes)[1] != 'try-error') {
            cc_el_tree = try(recordPlot())
        }
    }
    mf_weight_nodes = mf_weight_tree = NULL
    if (do_mf_weight_tree) {
        included = length(which(topGO::score(tg$results$mf_weight) <= score_limit))                
        mf_weight_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$mf_godata, topGO::score(tg$results$mf_weight), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(mf_weight_nodes)[1] != 'try-error') {
            mf_weight_tree = try(recordPlot())
        }
    }
    bp_weight_nodes = bp_weight_tree = NULL
    if (do_bp_weight_tree) {
        included = length(which(topGO::score(tg$results$bp_weight) <= score_limit))                
        bp_weight_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$bp_godata, topGO::score(tg$results$bp_weight), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(bp_weight_nodes)[1] != 'try-error') {
            bp_weight_tree = try(recordPlot())
        }
    }
    cc_weight_nodes = cc_weight_tree = NULL
    if (do_cc_weight_tree) {
        included = length(which(topGO::score(tg$results$cc_weight) <= score_limit))                
        cc_weight_nodes = try(suppressWarnings(topGO::showSigOfNodes(tg$cc_godata, topGO::score(tg$results$cc_weight), useInfo="all", sigForAll=sigforall, firstSigNodes=included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
        if (class(cc_weight_nodes)[1] != 'try-error') {
            cc_weight_tree = try(recordPlot())
        }
    }
    
    trees = list(
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
#' @return a big list including the various outputs from topgo
#' @export
simple_clusterprofiler = function(de_genes, goids=NULL, golevel=4, pcutoff=0.1,
    qcutoff=1.0, fold_changes=NULL, include_cnetplots=TRUE,
    showcategory=12, universe=NULL, organism="lm", gff=NULL,
    wrapped_width=20, method="Walllenius", padjust="BH") {
    genetable_test = try(load("geneTable.rda"))
    if (class(genetable_test) == 'try-error') {
        if (!is.null(gff)) {
            message("Generating the geneTable.rda")
            ## clusterProfiler::Gff2GeneTable(gff)
            hpgltools::Gff2GeneTable(gff)            
        } else {
            stop("cluster Profiler requires a geneTable.rda, which requires a gff file to read.")
        }
    } else {
        rm(genetable_test)
    }

    if (is.null(de_genes$ID)) {
        gene_list = as.character(rownames(de_genes))
    } else {
        gene_list = as.character(de_genes$ID)
    }
    gomapping_test = try(load("GO2EG.rda"))
    if (class(gomapping_test) == 'try-error') {
        message("Generating GO mapping data for cluster profiler from the goids data.")
        gomap = goids
        colnames(gomap) = c("entrezgene", "go_accession")
        clusterProfiler::buildGOmap(gomap)
    } else {
        message("Using GO mapping data located in GO2EG.rda")
    }
    message("Testing gseGO")
    ego2 = try(clusterProfiler::gseGO(geneList=gene_list, organism=organism, ont="GO", nPerm=100, minGSSize=2, pvalueCutoff=1, verbose=TRUE))
    print(ego2)
    message("Starting MF(molecular function) analysis")
    mf_group = clusterProfiler::groupGO(gene_list, organism=organism, ont="MF", level=golevel, readable=TRUE)
    mf_all = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="MF", pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod="none")
    all_mf_phist = try(hpgltools::hpgl_histogram(mf_all@result$pvalue, bins=20))
    if (class(all_mf_phist) != 'try-error') {
        y_limit = (sort(unique(table(all_mf_phist$data)), decreasing=TRUE)[2]) * 2
        all_mf_phist = all_mf_phist + scale_y_continuous(limits=c(0, y_limit))
    }
    enriched_mf = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="MF", pvalueCutoff=pcutoff, qvalueCutoff=qcutoff, pAdjustMethod=padjust)
    
    message("Starting BP(biological process) analysis")
    bp_group = clusterProfiler::groupGO(gene_list, organism=organism, ont="BP", level=golevel, readable=TRUE)
    bp_all = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="BP", pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod="none")
    all_bp_phist = try(hpgltools::hpgl_histogram(bp_all@result$pvalue, bins=20))
    if (class(all_bp_phist) != 'try-error') {
        y_limit = (sort(unique(table(all_bp_phist$data)), decreasing=TRUE)[2]) * 2
        all_bp_phist = all_bp_phist + scale_y_continuous(limits=c(0, y_limit))
    }
    
    enriched_bp = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="BP", pvalueCutoff=pcutoff, qvalueCutoff=qcutoff, pAdjustMethod=padjust)

    message("Starting CC(cellular component) analysis")
    cc_group = clusterProfiler::groupGO(gene_list, organism=organism, ont="CC", level=golevel, readable=TRUE)
    cc_all = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="CC", pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod="none")
    enriched_cc = hpgltools::hpgl_enrichGO(gene_list, organism=organism, ont="CC", pvalueCutoff=pcutoff, qvalueCutoff=qcutoff, pAdjustMethod=padjust)
    all_cc_phist = try(hpgltools::hpgl_histogram(cc_all@result$pvalue, bins=20))
    ## Try and catch if there are no significant hits.
    if (class(all_cc_phist) != 'try-error') {
        y_limit = (sort(unique(table(all_cc_phist$data)), decreasing=TRUE)[2]) * 2
        all_cc_phist = all_cc_phist + scale_y_continuous(limits=c(0, y_limit))
    }
    
    mf_group_barplot = try(barplot(mf_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(mf_group_barplot)[1] != 'try-error') {
        mf_group_barplot$data$Description = as.character(lapply(strwrap(mf_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    
    bp_group_barplot = try(barplot(bp_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(bp_group_barplot)[1] != 'try-error') {
        bp_group_barplot$data$Description = as.character(lapply(strwrap(bp_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    
    cc_group_barplot = try(barplot(cc_group, drop=TRUE, showCategory=showcategory), silent=TRUE)
    if (class(cc_group_barplot)[1] != 'try-error') {
        cc_group_barplot$data$Description = as.character(lapply(strwrap(cc_group_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    all_mf_barplot = try(barplot(mf_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)    
    enriched_mf_barplot = try(barplot(enriched_mf, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_mf_barplot)[1] == 'try-error') {
        message("No enriched MF groups were observed.")
    } else {
        enriched_mf_barplot$data$Description = as.character(lapply(strwrap(enriched_mf_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_mf_barplot)[1] != 'try-error') {
        all_mf_barplot$data$Description = as.character(lapply(strwrap(all_mf_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    all_bp_barplot = try(barplot(bp_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)        
    enriched_bp_barplot = try(barplot(enriched_bp, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_bp_barplot)[1] == 'try-error') {
        message("No enriched BP groups observed.")
    } else {
        enriched_bp_barplot$data$Description = as.character(lapply(strwrap(enriched_bp_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_bp_barplot)[1] != 'try-error') {
        all_bp_barplot$data$Description = as.character(lapply(strwrap(all_bp_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }

    all_cc_barplot = try(barplot(cc_all, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    enriched_cc_barplot = try(barplot(enriched_cc, categorySize="pvalue", showCategory=showcategory), silent=TRUE)
    if (class(enriched_cc_barplot)[1] == 'try-error') {
        message("No enriched CC groups observed.")
    } else {
        enriched_cc_barplot$data$Description = as.character(lapply(strwrap(enriched_cc_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    if (class(all_cc_barplot)[1] != 'try-error') {
        all_cc_barplot$data$Description = as.character(lapply(strwrap(all_cc_barplot$data$Description, wrapped_width, simplify=F),paste,collapse="\n"))
    }
    
    if (include_cnetplots == TRUE) {
        message("Attempting to include the cnetplots from clusterProfiler.")
        message("They fail often, if this is causing errors, set:")
        message("include_cnetplots to FALSE")
        cnetplot_mf = try(clusterProfiler::cnetplot(enriched_mf, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_mf)[1] != 'try-error') {
            cnetplot_mf = recordPlot()
        } else {
            message("cnetplot just failed for the MF ontology.  Do not be concerned with the previous error.")
        }
        cnetplot_bp = try(clusterProfiler::cnetplot(enriched_bp, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_bp)[1] != 'try-error') {
            cnetplot_bp = recordPlot()
        } else {
            message("cnetplot just failed for the BP ontology.  Do not be concerned with the previous error.")
        }            
        cnetplot_cc = try(clusterProfiler::cnetplot(enriched_cc, categorySize="pvalue", foldChange=fold_changes))
        if (class(cnetplot_cc)[1] != 'try-error') {
            cnetplot_cc = recordPlot()
        } else {
            message("cnetplot just failed for the CC ontology.  Do not be concerned with the previous error.")
        }
    }
    
    if (!is.null(mf_all)) {
        mf_interesting = mf_all@result
        rownames(mf_interesting) = NULL
        mf_interesting$ont = "MF"
        mf_interesting = mf_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Description")]    
        mf_interesting = subset(mf_interesting, pvalue <= 0.1)
    } else {
        mf_interesting = NULL
    }
    if (!is.null(bp_all)) {
        bp_interesting = bp_all@result
        rownames(bp_interesting) = NULL
        bp_interesting$ont = "BP"
        bp_interesting = bp_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Description")]    
        bp_interesting = subset(bp_interesting, pvalue <= 0.1)
    } else {
        bp_interesting = NULL
    }
    if (!is.null(cc_all)) {
        cc_interesting = cc_all@result
        rownames(cc_interesting) = NULL
        cc_interesting$ont = "CC"
        cc_interesting = cc_interesting[,c("ID","ont","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Description")]
        cc_interesting = subset(cc_interesting, pvalue <= 0.1)
    } else {
        cc_interesting = NULL
    }
    
    return_information = list(
        mf_interesting=mf_interesting, bp_interesting=bp_interesting, cc_interesting=cc_interesting,
        mf_pvals=all_mf_phist, bp_pvals=all_bp_phist, cc_pvals=all_cc_phist,        
        mf_enriched=enriched_mf, bp_enriched=enriched_bp, cc_enriched=enriched_cc,
        mf_all=mf_all, bp_all=bp_all, cc_all=cc_all,
        mf_all_barplot=all_mf_barplot, bp_all_barplot=all_bp_barplot, cc_all_barplot=all_cc_barplot,
        mf_enriched_barplot=enriched_mf_barplot, bp_enriched_barplot=enriched_bp_barplot, cc_enriched_barplot=enriched_cc_barplot,
        mf_cnetplot=cnetplot_mf, bp_cnetplot=cnetplot_bp, cc_cnetplot=cnetplot_cc,
        mf_group=mf_group, bp_group=bp_group, cc_group=cc_group,
        mf_group_barplot=mf_group_barplot, bp_group_barplot=bp_group_barplot, cc_group_barplot=cc_group_barplot)
    return(return_information)        
}

make_id2gomap = function(goid_map="reference/go/id2go.map", goids_df=NULL, overwrite=FALSE) {
    id2go_test = file.info(goid_map)
    goids_dir = dirname(goid_map)
    if (!file.exists(goids_dir)) {
        dir.create(goids_dir, recursive=TRUE)
    }

    if (isTRUE(overwrite)) {
        if (is.null(goids_df)) {
            stop("There is neither a id2go file nor a data frame of goids.")
        } else {
            message("Attempting to generate a id2go file in the format expected by topGO.")
            new_go = plyr::ddply(goids_df, .(ID), summarise, GO=paste(unique(GO), collapse=','))
            write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
            rm(id2go_test)
        }
    } else {
        if (is.na(id2go_test$size)) {
            if (is.null(goids_df)) {
                stop("There is neither a id2go file nor a data frame of goids.")
            } else {
                message("Attempting to generate a id2go file in the format expected by topGO.")
                new_go = plyr::ddply(goids_df, .(ID), summarise, GO=paste(unique(GO), collapse=','))
                write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
                rm(id2go_test)
            }
        }
    }
}

#' Make fun trees a la topgo from goseq data.
#'
#' @param de_genes some differentially expressed genes
#' @param godata data from goseq
#' @param goids a mapping of IDs to GO in the Ramigo expected format
#' @param sigforall Print significance on all nodes?
#' 
#' @return a plot!
#' @seealso \code{\link{Ramigo}}
#' @export
goseq_trees = function(de_genes, godata, goid_map="reference/go/id2go.map", score_limit=0.01, goids_df=NULL, overwrite=FALSE, selector="topDiffGenes", pval_column="adj.P.Val") {
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO = topGO::readMappings(file=goid_map)
    annotated_genes = names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID = make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes = factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) = annotated_genes    

    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata = new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)        
    } else {
        pvals = as.vector(de_genes[[pval_column]])
        names(pvals) = rownames(de_genes)
        mf_GOdata = new("topGOdata", description="MF", ontology="MF", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", description="BP", ontology="BP", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", description="CC", ontology="CC", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    }
        
    enriched_ids = godata$alldata$category
    enriched_scores = godata$alldata$over_represented_pvalue
    names(enriched_scores) = enriched_ids
    
    mf_avail_nodes = as.list(mf_GOdata@graph@nodes)
    names(mf_avail_nodes) = mf_GOdata@graph@nodes
    mf_nodes = enriched_scores[names(enriched_scores) %in% names(mf_avail_nodes)]
    mf_included = length(which(mf_nodes <= score_limit))
    mf_tree_data = try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=mf_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(mf_tree_data) == 'try-error') {
        message("There was an error generating the MF tree.")
        mf_tree = NULL
    } else {
        mf_tree = recordPlot()
    }

    bp_avail_nodes = as.list(bp_GOdata@graph@nodes)
    names(bp_avail_nodes) = bp_GOdata@graph@nodes
    bp_nodes = enriched_scores[names(enriched_scores) %in% names(bp_avail_nodes)]
    bp_included = length(which(bp_nodes <= score_limit))
    bp_tree_data = try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=bp_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(bp_tree_data) == 'try-error') {
        message("There was an error generating the BP tree.")
        bp_tree = NULL
    } else {
        bp_tree = recordPlot()
    }    

    cc_avail_nodes = as.list(cc_GOdata@graph@nodes)
    names(cc_avail_nodes) = cc_GOdata@graph@nodes
    cc_nodes = enriched_scores[names(enriched_scores) %in% names(cc_avail_nodes)]
    cc_included = length(which(cc_nodes <= score_limit))
    cc_tree_data = try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_nodes, useInfo="all", sigForAll=TRUE, firstSigNodes=cc_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(cc_tree_data) == 'try-error') {
        message("There was an error generating the CC tree.")
        cc_tree = NULL
    } else {
        cc_tree = recordPlot()
    }
    trees = list(MF=mf_tree, BP=bp_tree, CC=cc_tree, MFdata=mf_tree_data, BPdata=bp_tree_data, CCdata=cc_tree_data)
    return(trees)
}

## Take clusterprofile group data and print it on a tree as topGO does
#' Make fun trees a la topgo from goseq data.
#'
#' @param de_genes some differentially expressed genes
#' @param godata data from cluster Profiler
#' @param goids a mapping of IDs to GO in the Ramigo expected format
#' @param sigforall Print significance on all nodes?
#' 
#' @return a plot!
#' @seealso \code{\link{Ramigo}}
#' @export
cluster_trees = function(de_genes, cpdata, goid_map="reference/go/id2go.map", goids_df=NULL, score_limit=0.1, overwrite=FALSE, selector="topDiffGenes", pval_column="adj.P.Value") {
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO = topGO::readMappings(file=goid_map)
    annotated_genes = names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID = make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes = factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) = annotated_genes

    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata = new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=annFUN.gene2GO, gene2GO=geneID2GO)        
    } else {
        pvals = as.vector(de_genes[[pval_column]])
        names(pvals) = rownames(de_genes)
        mf_GOdata = new("topGOdata", description="MF", ontology="MF", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata = new("topGOdata", description="BP", ontology="BP", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata = new("topGOdata", description="CC", ontology="CC", allGenes=pvals, geneSel=get(selector), annot=annFUN.gene2GO, gene2GO=geneID2GO)
    }

    mf_all = cpdata$mf_all
    mf_enriched = cpdata$mf_enriched
    bp_all = cpdata$bp_all
    bp_enriched = cpdata$bp_enriched
    cc_all = cpdata$cc_all
    cc_enriched = cpdata$cc_enriched
    mf_all_ids = mf_all@result$ID
    bp_all_ids = bp_all@result$ID
    cc_all_ids = cc_all@result$ID
    mf_all_scores = mf_all@result$p.adjust
    bp_all_scores = bp_all@result$p.adjust
    cc_all_scores = cc_all@result$p.adjust
    names(mf_all_scores) = mf_all_ids
    names(bp_all_scores) = bp_all_ids
    names(cc_all_scores) = cc_all_ids
    mf_included = length(which(mf_all_scores <= score_limit))
    mf_tree_data = try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_all_scores, useInfo="all", sigForAll=TRUE, firstSigNodes=mf_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(mf_tree_data)[1] == 'try-error') {
        mf_tree = NULL
    } else {
        mf_tree = recordPlot()
    }
    bp_included = length(which(bp_all_scores <= score_limit))
    bp_tree_data = try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_all_scores, useInfo="all", sigForAll=TRUE, firstSigNodes=bp_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(bp_tree_data)[1] == 'try-error') {
        bp_tree = NULL
    } else {
        bp_tree = recordPlot()
    }
    cc_included = length(which(cc_all_scores <= score_limit))
    cc_tree_data = try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_all_scores, useInfo="all", sigForAll=TRUE, firstSigNodes=cc_included, useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(cc_tree_data)[1] == 'try-error') {
        cc_tree = NULL
    } else {
        cc_tree = recordPlot()
    }
    trees = list(MF=mf_tree, BP=bp_tree, CC=cc_tree, MFdata=mf_tree_data, BPdata=bp_tree_data, CCdata=cc_tree_data)
    return(trees)
}

## Please note that the KGML parser fails if other XML parsers are loaded into R
#' Print some data onto KEGG pathways
#'
#' @param de_genes some differentially expressed genes
#' @param godata data from cluster Profiler
#' @param goids a mapping of IDs to GO in the Ramigo expected format
#' @param sigforall Print significance on all nodes?
#' 
#' @return a plot!
#' @seealso \code{\link{Ramigo}}
#' @export
hpgl_pathview = function(path_data, indir="pathview_in", outdir="pathview", pathway="all", species="lma", string_from="LmjF", string_to="LMJF", suffix="_colored", second_from=NULL, second_to=NULL, verbose=FALSE) {
    ##eh = new.env(hash=TRUE, size=NA)
    ## There is a weird namespace conflict when using pathview, so I will reload it here
    try(detach("package:Rgraphviz", unload=TRUE))
    try(detach("package:topGO", unload=TRUE))
    try(detach("package:pathview", unload=TRUE))
    try(detach("package:KEGGgraph", unload=TRUE))
    try(detach("package:RamiGO", unload=TRUE))
    try(detach("package:graph", unload=TRUE))
    library("pathview")
    ## Testing parameters
    ##path_data = kegg_list
    ##indir="pathview_in"
    ##outdir="pathview_epi_high"
    ##pathway="all"
    ##species="tcr"
    ##string_from="TcCLB."
    ##string_to=""
    ##suffix="_epi_high"
    ## End testing parameters    
#    environment(eh)
    ## Massage the names to KEGG compatible names
    tmp_names = names(path_data)
    tmp_names = gsub(string_from, string_to, tmp_names)
    if (!is.null(second_from)) {
        tmp_names = gsub(second_from, second_to, tmp_names)
    }
##    tmp_names = gsub("\\.","_", tmp_names)
    names(path_data) = tmp_names
    rm(tmp_names)
    
    ## First check that the input pathview directory exists
    if (!file.exists(indir)) {
        dir.create(indir)
    }
    if (!file.exists(outdir)){
            dir.create(outdir)
        }
    paths = list()
    if (pathway == "all") {
        all_pathways = unique(KEGGREST::keggLink("pathway", species))
        paths = all_pathways
        paths = gsub("path:", "", paths)
        all_modules = unique(KEGGREST::keggLink("module", species))
    } else if (class(pathway) == "list") {
        paths = pathway
    } else {
        paths[1] = pathway
    }
    return_list = list()
    for (count in 1:length(paths)) {
        path = paths[count]
        limits=c(min(path_data, na.rm=TRUE), max(path_data, na.rm=TRUE))
        if (isTRUE(verbose)) {
            pv = try(pathview::pathview(gene.data=path_data, kegg.dir=indir, pathway.id=path, species=species, limit=list(gene=limits, cpd=limits), map.null=TRUE, gene.idtype="KEGG", out.suffix=suffix, split.group=TRUE, expand.node=TRUE, kegg.native=TRUE, map.symbol=TRUE, same.layer=FALSE, res=1200, new.signature=FALSE, cex=0.05, key.pos="topright"))
        } else {
            pv = suppressMessages(try(pathview::pathview(gene.data=path_data, kegg.dir=indir, pathway.id=path, species=species, limit=list(gene=limits, cpd=limits), map.null=TRUE, gene.idtype="KEGG", out.suffix=suffix, split.group=TRUE, expand.node=TRUE, kegg.native=TRUE, map.symbol=TRUE, same.layer=FALSE, res=1200, new.signature=FALSE, cex=0.05, key.pos="topright")))
        }
        if (class(pv) == "numeric") {
            colored_genes = NULL
            newfile = NULL
            up = NULL
            down = NULL
        } else {
            colored_genes = dim(pv$plot.data.gene)[1]
            ##        "lma04070._proeff.png"
            oldfile = paste(path, ".", suffix, ".png", sep="")
            newfile = paste(outdir,"/", path, suffix, ".png", sep="")
            if (isTRUE(verbose)) {
                message(paste("Moving file to: ", newfile, sep=""))
            }
            file.rename(from=oldfile, to=newfile)
            data_low = summary(path_data)[2]
            data_high = summary(path_data)[3]
            numbers_in_plot = as.numeric(pv$plot.data.gene$mol.data)
            up = sum(numbers_in_plot > data_high, na.rm=TRUE)
            down = sum(numbers_in_plot < data_low, na.rm=TRUE)
        }
        return_list[[path]]$file = newfile
        return_list[[path]]$genes = colored_genes
        return_list[[path]]$up = up
        return_list[[path]]$down = down
    }

    retdf = data.frame(rep(NA, length(names(return_list))))
    rownames(retdf) = names(return_list)
    retdf$genes = NA
    retdf$up = NA
    retdf$down = NA
    colnames(retdf) = c("file","genes","up","down")
    for (path in names(return_list)) {
        if (is.null(return_list[[path]]$genes)) {
            retdf[path,]$genes = 0
        } else {
            retdf[path,]$genes = as.numeric(return_list[[path]]$genes)
        }
        retdf[path,]$file = as.character(return_list[[path]]$file)
        retdf[path,]$up = as.numeric(return_list[[path]]$up)
        retdf[path,]$down = as.numeric(return_list[[path]]$down)        
    }
    retdf = retdf[with(retdf, order(up, down)), ]
    return(retdf)
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
hpgl_enrichGO = function(gene, organism="human", ont="MF",
    pvalueCutoff=0.05, pAdjustMethod="BH", universe,
    qvalueCutoff=0.2, minGSSize=2, readable=FALSE) {
    ## Testing parameters
    ##gene=gene_list
    ##organism="lm"
    ##ont="BP"
    ##minGSSize=2
    ## End testing parameters
    information = hpgl_enrich.internal(gene, organism=organism, pvalueCutoff=pvalueCutoff,
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
hpgl_enrich.internal = function(gene, organism, pvalueCutoff=1, pAdjustMethod="BH",
    ont, minGSSize=2, qvalueCutoff=0.2, readable=FALSE, universe=NULL) {
    gene <- as.character(gene)
    class(gene) <- ont
    qExtID2TermID = DOSE::EXTID2TERMID(gene, organism)
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
    cat(sprintf("The minimum observed adjusted pvalue is: %f\n", min(p.adj)))
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


## this function will plot the GO DAG or parts of it
## sigNodes:     a named vector of terms p-values, the names are the GO terms
## wantedNodes:  the nodes that we want to find, we will plot this nodes with
##               a different color. The vector contains the names pf the nodes
## oldSigNodes:  used to plot the (new) sigNodes in the same collor range
##               as the old ones
## export.to.dot.file: is a global variable given the name of the output .dot file

#' @export
getEdgeWeights <- function (graph) {  
  weightsList <- graph::edgeWeights(graph)
  to <- lapply(weightsList, names)
  from <- nodes(graph)

  if (any(is.na(unlist(to))) || any(is.na(from))) 
    stop("Edge names do not match node names.")

  edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
  edge.weights <- unlist(weightsList)
  names(edge.weights) <- edge.names

  return(edge.weights)
}

#' A minor hack in the topGO GOplot function
#'
#' @export
hpgl_GOplot <- function(dag, sigNodes, dag.name = 'GO terms', edgeTypes = T,
                     nodeShape.type = c('box', 'circle', 'ellipse', 'plaintext')[3],
                     genNodes = NULL, wantedNodes = NULL, showEdges = T, useFullNames = T,
                     oldSigNodes = NULL, nodeInfo=nodeInfo, maxchars=30) {
    
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
      ##nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
      ##nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
    }
  }

  ## we will change the shape and the color of the nodes we want back
  if(!is.null(genNodes)) {
    nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
    nodeAttrs$shape[genNodes] <- rep('box', .ln)
    #nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)    
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
      #if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)){
      #  print(old.range)
      #  print(range(logSigNodes))
      #  stop('some stupid error here :)')
      #}
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
        ##edgeAttrs$color <- ifelse(getEdgeWeights(dag) == 0, 'black', 'red')
        edgeAttrs$color <- ifelse(hpgltools::getEdgeWeights(dag) == 0, 'black', 'black')
  

  ##plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)

  return(agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}


GOplot.orig <- function(dag, sigNodes, dag.name = 'GO terms', edgeTypes = T,
                   nodeShape.type = c('box', 'circle', 'ellipse', 'plaintext')[3],
                   genNodes = NULL, wantedNodes = NULL, showEdges = T, useFullNames = F,
                   oldSigNodes = NULL, nodeInfo = NULL) {
    
  if(!missing(sigNodes))
    sigNodeInd = TRUE
  else
    sigNodeInd = FALSE
  
  ## we set the global Graphviz attributes
  graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL

  #graphAttrs$graph$splines <- FALSE
  
  ## set the node shape
  graphAttrs$node$shape <- nodeShape.type

  ## set the fontsize for the nodes labels
  graphAttrs$node$fontsize <- '14'
  #graphAttrs$node$height <- '1.0'
  #graphAttrs$node$width <- '1.5'

  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()

  ## try to use adaptive node size
  #nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))
  
  if(is.null(nodeInfo)) {
    nodeInfo <- character(numNodes(dag))
    names(nodeInfo) <- nodes(dag)
  }
  else
    nodeInfo <- paste('\\\n', nodeInfo, sep = '')
  
  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
  node.names <- nodes(dag)
  if(!useFullNames)
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
  if(!is.null(wantedNodes)) {
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
    #nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)    
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
      #if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)){
      #  print(old.range)
      #  print(range(logSigNodes))
      #  stop('some stupid error here :)')
      #}
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
  

  ##plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)

  return(agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}


hpgl_GroupDensity = function(object, whichGO, ranks=TRUE, rm.one=FALSE) {
    ## Testing parameters
    ##object = mf_GOdata
    ##whichGO = mf_first_group
    ## End testing parameters
    groupMembers <- topGO::genesInTerm(object, whichGO)[[1]]
    allS <- topGO::geneScore(object, use.names = TRUE)
    if(rm.one) {
        allS <- allS[allS < 0.99]
    }
    xlab <- "Gene' score"
    if(ranks) {
        allS <- BiocGenerics::rank(allS, ties.method = "random")
        xlab <- "Gene's rank" 
    }
    group <- as.integer(names(allS) %in% groupMembers)
    xx <- data.frame(score=allS, group = factor(group, labels=paste(c("complementary", whichGO), "  (", table(group), ")", sep="")))
    plot = lattice::densityplot( ~ score | group, data=xx, layout=c(1,2), xlab=xlab)
    return(plot)
}

## Functions in this are not exported by topGO
Gff2GeneTable <- function(gffFile, compress=TRUE) {
    ##gffFile="reference/gff/clbrener_8.1_complete_genes.gff"
    if (is.data.frame(gffFile)) {
        GeneID = data.frame(GeneID = gffFile$ID)
        geneInfo = gffFile
        geneInfo$start = 1
        geneInfo$GeneID = gffFile$ID
        geneInfo$GeneName = gffFile$ID
        geneInfo$Locus = gffFile$ID        
        geneInfo$end = geneInfo$width
        geneInfo$strand = "+"
    } else {
        gff <- readGff(gffFile)
        GeneID <- data.frame(GeneID=getGffAttribution(gff$attributes, field="ID"))
        geneInfo <- gff[gff$feature == "gene",]
        geneInfo <- geneInfo[, c("seqname", "start", "end", "strand", "attributes")]
        geneInfo$GeneID <- getGffAttribution(geneInfo$attributes, field="ID")
        geneInfo$GeneName <- getGffAttribution(geneInfo$attributes, field="Name")
        geneInfo$Locus <- getGffAttribution(geneInfo$attributes, field="locus_tag")
        geneInfo$GeneName[is.na(geneInfo$GeneName)] <- "-"
        geneInfo <- geneInfo[, -5] ## abondom "attributes" column.
    }
            ## GI2GeneID <- data.frame(GI=getGffAttribution(gff$attributes, field="GI"),
    ##                        GeneID=getGffAttribution(gff$attributes, field="GeneID")
    ##                                    #,
    ##                                    #Product=getGffAttribution(gff$attributes, field="product")
    ##                        )
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$GI),]
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$Gene),]

        

    ## geneTable <- merge(GI2GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- merge(GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- unique(geneTable)
    if (compress) {
        save(geneTable, file="geneTable.rda", compress="xz")
    } else {
        save(geneTable, file="geneTable.rda")
    }
    message("Gene Table file save in the working directory.")
}

parseKGML2Graph2 <-function (file, ...) {
    pathway <- parseKGML2(file)
    gR <- KEGGpathway2Graph2(pathway, ...)
    return(gR)
}

hpgl_base_pathview = function (gene.data = NULL, cpd.data = NULL, xml.file = NULL, 
    pathway.id, species = "hsa", kegg.dir = ".", cpd.idtype = "kegg", 
    gene.idtype = "entrez", gene.annotpkg = NULL, min.nnodes = 3, 
    kegg.native = TRUE, map.null = TRUE, expand.node = FALSE, 
    split.group = FALSE, map.symbol = TRUE, map.cpdname = TRUE, 
    node.sum = "sum", discrete = list(gene = FALSE, cpd = FALSE), 
    limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
    both.dirs = list(gene = T, cpd = T), trans.fun = list(gene = NULL, 
        cpd = NULL), low = list(gene = "green", cpd = "blue"), 
    mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", 
        cpd = "yellow"), na.col = "transparent", ...) {
    if (is.character(gene.data)) {
        gd.names = gene.data
        gene.data = rep(1, length(gene.data))
        names(gene.data) = gd.names
        both.dirs$gene = FALSE
        ng = length(gene.data)
        nsamp.g = 1
    }
    else if (!is.null(gene.data)) {
        if (length(dim(gene.data)) == 2) {
            gd.names = rownames(gene.data)
            ng = nrow(gene.data)
            nsamp.g = 2
        }
        else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
            gd.names = names(gene.data)
            ng = length(gene.data)
            nsamp.g = 1
        }
        else stop("wrong gene.data format!")
    }
    else if (is.null(cpd.data)) {
        stop("gene.data and cpd.data are both NULL!")
    }
    gene.idtype = toupper(gene.idtype)
    data(bods)
    data(gene.idtype.list)
    if (species != "ko") {
        species.data = pathview::kegg.species.code(species, na.rm = T, 
            code.only = FALSE)
    }
    else {
        species.data = c(kegg.code = "ko", entrez.gnodes = "0", 
            kegg.geneid = "K01488", ncbi.geneid = "")
        gene.idtype = "KEGG"
        msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message(msg)
    }
    if (length(dim(species.data)) == 2) {
        message("More than two valide species!")
        species.data = species.data[1, ]
    }
    species = species.data["kegg.code"]
    entrez.gnodes = species.data["entrez.gnodes"] == 1
    if (is.na(species.data["ncbi.geneid"])) {
        if (!is.na(species.data["kegg.geneid"])) {
            msg.fmt = "Only native KEGG gene ID is supported for this species,\nmake sure it looks like \"%s\"!"
            msg = sprintf(msg.fmt, species.data["kegg.geneid"])
            message(msg)
        }
        else {
            stop("This species is not annotated in KEGG!")
        }
    }
    if (is.null(gene.annotpkg)) 
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
    if (length(grep("ENTREZ|KEGG", gene.idtype)) < 1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg)) 
            stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.list) 
            stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype, 
            pkg.name = gene.annotpkg)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
    }
    if (gene.idtype == "ENTREZ" & !entrez.gnodes & !is.null(gene.data)) {
        message("Getting gene ID data from KEGG...")
        gene.idmap = keggConv("ncbi-geneid", species)
        message("Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        ncbi.ids = gsub("ncbi-geneid:", "", gene.idmap)
        gene.idmap = cbind(ncbi.ids, kegg.ids)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "KEGG"
    }
    if (is.character(cpd.data)) {
        cpdd.names = cpd.data
        cpd.data = rep(1, length(cpd.data))
        names(cpd.data) = cpdd.names
        both.dirs$cpd = FALSE
        ncpd = length(cpd.data)
    }
    else if (!is.null(cpd.data)) {
        if (length(dim(cpd.data)) == 2) {
            cpdd.names = rownames(cpd.data)
            ncpd = nrow(cpd.data)
        }
        else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
            cpdd.names = names(cpd.data)
            ncpd = length(cpd.data)
        }
        else stop("wrong cpd.data format!")
    }
    if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
        data(rn.list)
        cpd.types = c(names(rn.list), "name")
        cpd.types = tolower(cpd.types)
        cpd.types = cpd.types[-grep("kegg", cpd.types)]
        if (!tolower(cpd.idtype) %in% cpd.types) 
            stop("Wrong input cpd ID type!")
        cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
        cpd.data = mol.sum(cpd.data, cpd.idmap)
    }
    warn.fmt = "Parsing %s file failed, please check the file!"
    if (length(grep(species, pathway.id)) > 0) {
        pathway.name = pathway.id
        pathway.id = gsub(species, "", pathway.id)
    }
    else pathway.name = paste(species, pathway.id, sep = "")
    kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
    tfiles = paste(pathway.name, c("xml", "png"), sep = ".")
    if (!all(tfiles %in% kfiles)) {
        dstatus = download.kegg(pathway.id = pathway.id, species = species, 
            kegg.dir = kegg.dir)
        if (dstatus == "failed") {
            warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
            warn.msg = sprintf(warn.fmt, pathway.name)
            message(warn.msg)
            return(invisible(0))
        }
    }
    if (missing(xml.file)) 
        xml.file <- paste(kegg.dir, "/", pathway.name, ".xml", 
            sep = "")
    if (kegg.native) {
        node.data = try(node.info(xml.file), silent = T)
        if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file)
            message(warn.msg)
            return(invisible(0))
        }
        node.type = c("gene", "enzyme", "compound", "ortholog")
        sel.idx = node.data$type %in% node.type
        nna.idx = !is.na(node.data$x + node.data$y + node.data$width + 
            node.data$height)
        sel.idx = sel.idx & nna.idx
        if (sum(sel.idx) < min.nnodes) {
            warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
            warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name)
            message(warn.msg)
            return(invisible(0))
        }
        node.data = lapply(node.data, "[", sel.idx)
    }
    else {
        gR1 = try(parseKGML2Graph2(xml.file, genes = F, expand = expand.node, 
            split.group = split.group), silent = T)
        node.data = try(node.info(gR1), silent = T)
        if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file)
            message(warn.msg)
            return(invisible(0))
        }
    }
    if (species == "ko") 
        gene.node.type = "ortholog"
    else gene.node.type = "gene"
    if ((!is.null(gene.data) | map.null) & sum(node.data$type == 
        gene.node.type) > 1) {
        plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type, 
            node.sum = node.sum, entrez.gnodes = entrez.gnodes)
        kng = plot.data.gene$kegg.names
        kng.char = gsub("[0-9]", "", unlist(kng))
        if (any(kng.char > "")) 
            entrez.gnodes = FALSE
        if (map.symbol & species != "ko" & entrez.gnodes) {
            if (is.na(gene.annotpkg)) {
                warm.fmt = "No annotation package for the species %s, gene symbols not mapped!"
                warm.msg = sprintf(warm.fmt, species)
                message(warm.msg)
            }
            else {
                plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names), 
                  category = "SYMBOL", pkg.name = gene.annotpkg)[, 
                  2]
                mapped.gnodes = rownames(plot.data.gene)
                node.data$labels[mapped.gnodes] = plot.data.gene$labels
            }
        }
        cols.ts.gene = node.color(plot.data.gene, limit$gene, 
            bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene, 
            discrete = discrete$gene, low = low$gene, mid = mid$gene, 
            high = high$gene, na.col = na.col)
    }
    else plot.data.gene = cols.ts.gene = NULL
    if ((!is.null(cpd.data) | map.null) & sum(node.data$type == 
        "compound") > 1) {
        plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound", 
            node.sum = node.sum)
        if (map.cpdname & !kegg.native) {
            plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[, 
                2]
            mapped.cnodes = rownames(plot.data.cpd)
            node.data$labels[mapped.cnodes] = plot.data.cpd$labels
        }
        cols.ts.cpd = node.color(plot.data.cpd, limit$cpd, bins$cpd, 
            both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd, 
            discrete = discrete$cpd, low = low$cpd, mid = mid$cpd, 
            high = high$cpd, na.col = na.col)
    }
    else plot.data.cpd = cols.ts.cpd = NULL
    if (kegg.native) {
        pv.pars = keggview.native(plot.data.gene = plot.data.gene, 
            cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
            cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
            pathway.name = pathway.name, kegg.dir = kegg.dir, 
            limit = limit, bins = bins, both.dirs = both.dirs, 
            discrete = discrete, low = low, mid = mid, high = high, 
            na.col = na.col, ...)
    }
    else {
        pv.pars = keggview.graph(plot.data.gene = plot.data.gene, 
            cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
            cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
            path.graph = gR1, pathway.name = pathway.name, map.cpdname = map.cpdname, 
            split.group = split.group, limit = limit, bins = bins, 
            both.dirs = both.dirs, discrete = discrete, low = low, 
            mid = mid, high = high, na.col = na.col, ...)
    }
    plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
    if (!is.null(plot.data.gene)) {
        cnames = colnames(plot.data.gene)[-(1:7)]
        nsamp = length(cnames)/2
        if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                1):(2 * nsamp)], "col", sep = ".")
        }
        else cnames[2] = "mol.col"
        colnames(plot.data.gene)[-(1:7)] = cnames
    }
    plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
    if (!is.null(plot.data.cpd)) {
        cnames = colnames(plot.data.cpd)[-(1:7)]
        nsamp = length(cnames)/2
        if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                1):(2 * nsamp)], "col", sep = ".")
        }
        else cnames[2] = "mol.col"
        colnames(plot.data.cpd)[-(1:7)] = cnames
    }
    return(invisible(list(plot.data.gene = plot.data.gene, plot.data.cpd = plot.data.cpd)))
}
