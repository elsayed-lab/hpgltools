## Time-stamp: <Tue Jun 23 15:08:27 2015 Ashton Trey Belew (abelew@gmail.com)>
## Most of the functions in here probably shouldn't be exported...

#' Extract more easily readable information from a GOTERM datum
#'
#' The output from the GOTERM/GO.db functions is inconsistent, to put it nicely.
#' This attempts to extract from that heterogeneous datatype something easily readable.
#' Example:  Synonym() might return any of the following:
#' NA, NULL, "NA", "NULL", c("NA",NA,"GO:00001"), "GO:00002", c("Some text",NA, NULL, "GO:00003")
#' This function will boil that down to 'not found', '', 'GO:00004', or "GO:0001,  some text, GO:00004"
#'
#' @param The result of try(as.character(somefunction(GOTERM[id])), silent=TRUE)
#'   somefunction would be 'Synonym' 'Secondary' 'Ontology', etc...
#'
#' @return something more sane (hopefully)
#' @export
deparse_go_value = function(value) {
    result = ""
    if (class(value) == "try-error") {
        result = "Not found"
    } else {  ## Not an error
        if (is.null(value)) {
            result = ""
        } else if (is.na(value)) {
            result = ""
        } else if (value == "NULL") {
            result = ""
        } else if (grepl('^c\\(', as.character(value))) {
            value = eval(parse(text=as.character(value)))
            if (class(value) == "logical") {
                result = ""
            } else {
                value = as.character(value[which(complete.cases(value))])
                result = value
            }
        } else {  ## Just a string "GO:00023409"
            result = value
        }
    }
    return(result)
}

#' Get a go term from ID
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
#' @examples
#' ## goterm("GO:0032559")
#' ## > GO:0032559
#' ## > "adenyl ribonucleotide binding"
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
#' I think I will need to do similar parsing of the output for this function as per gosec()
#' In some cases this also returns stuff like c("some text", "GO:someID")
#' versus "some other text"  versus NULL versus NA
#'
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
#' @examples
#' ## text =  gosyn("GO:0000001")
#' ## text
#' ## > GO:000001
#' ## > "mitochondrial inheritance"
gosyn = function(go) {
    go = as.character(go)
    go = mapply(gosn, go)
    return(go)
}
gosn = function(go) {
    go = as.character(go)
    result = ""
    value = try(as.character(AnnotationDbi::Synonym(GOTERM[id])), silent=TRUE)
    result = deparse_go_value(value)
    return(result)
}

#' gosc does the real work for gosec()
#'
#' @param go A go id
#'
#' @return One of the following:
#'   "Not found" : for when the goID does not exist
#'   "" : when there is no secondary id
#'   or a character list of secondary IDs
gosc = function(go) {
    go = as.character(go)
    result = ""
    value = try(as.character(AnnotationDbi::Secondary(GOTERM[go])), silent=TRUE)
    result = deparse_go_value(value)
    return(result)
}

#' Get a go secondary ID from an id
#'
#' Unfortunately, GOTERM's returns for secondary IDs are not consistent, so this function
#' has to have a whole bunch of logic to handle the various outputs.
#' @param id A go ID -- this may be a character or list(assuming the elements, not names, are goids)
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
#' @examples
#' ## gosec("GO:0032432")
#' ## > GO:0032432
#' ## > "GO:0000141" "GO:0030482"
gosec = function(go) {
    go = mapply(gosc, go)
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
#' ## godef("GO:0032432")
#' ## > GO:0032432
#' ## > "An assembly of actin filaments that are on the same axis but may be oriented with the same or opposite polarities and may be packed with different levels of tightness."
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
#' ## goont(c("GO:0032432", "GO:0032433"))
#' ## > GO:0032432 GO:0032433
#' ## > "CC" "CC"
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
#' ## golev("GO:0032559")
#' ## > 3
golev = function(go, verbose=FALSE) {
    go = as.character(go)
    level = 0
    while(class(try(as.character(AnnotationDbi::Ontology(GOTERM[[go]])), silent=FALSE)) != 'try-error') {
        if(isTRUE(verbose)) {
            print(paste("Restarting while loop, level: ", level, go, sep=" "))
        }
        ontology = as.character(Ontology(GOTERM[[go]]))
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
        if(isTRUE(verbose)) {
            print("Incrementing level")
        }
        go = ancestors[1]
        level = level + 1
        if (go == "all") {
            return(level)
        }
    }  ## End while
    return(level)
}

#' Get a go level approximation from a set of IDs
#' This just wraps golev() in mapply.
#' @param id a character list of IDs
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
#' @examples
#' ## golevel(c("GO:0032559", "GO:0000001")
#' ## > 3 4
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
#' ## gotest("GO:0032559")
#' ## > 1
#' ## gotest("GO:0923429034823904")
#' ## > 0
gotst = function(go) {
    go = as.character(go)
    value = try(GOTERM[[go]])
    if (class(value) == 'try-error') {
        return(0)
    }
    if (is.null(value)) {
        return(0)
    } else {
        return(1)
    }
}

#' Test GO ids to see if they are useful
#' This just wraps gotst in mapply.
#'
#' @param id go IDs
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
gotest = function(go) {
    mapply(gotst, go)
}


gather_genes = function(goseq_data, ont='MF') {
    categories = NULL
    if (ont == 'MF') {
        categories = goseq_data$mf_subset
    } else if (ont == 'BP') {
        categories = goseq_data$bp_subset
    } else if (ont == 'CC') {
        categories = goseq_data$cc_subset
    } else {
        print('the ont argument didnt make sense, using mf.')
        categories = goseq_data$mf_subset
    }
    input = goseq_data$input
    categories = subset(categories, over_represented_pvalue < 0.05)
    cats = categories$category

    load("GO2ALLEG.rda")
    genes_per_ont = function(cat) {
        all_entries = GO2ALLEG[[cat]]
        entries_in_input = input[rownames(input) %in% all_entries,]
        names = as.character(rownames(entries_in_input))
        return(names)
    }
    allgenes_per_ont = function(cat) {
        all_entries = GO2ALLEG[[cat]]
        return(all_entries)
    }
    categories$gene_list = mapply(genes_per_ont, cats)
    categories$all_genes = mapply(allgenes_per_ont, cats)
    return(categories)
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
        theme(text=element_text(size=10)) + theme_bw()
    return(pvalue_plot)
}

### #' get_genelengths() Extract gene lengths from a gff file
### #'
### #' @param gff The file to extract
### #' @param ID Note the field to cross reference against to extract the genes
###get_genelengths = function(gff, ID="Note") {
###    annotations = BiocGenerics::as.data.frame(rtracklayer::import(gff, asRangedData=FALSE))
###    genes = annotations[annotations$type=="gene",]
###    genes$ID = unlist(genes[,ID])
###    genes = genes[,c("ID","width")]
###    return(genes)
###}

#' Perform ontology searches of the output from limma
#'
#' @param limma_out a list of topTables comprising limma outputs
#' @param n a number of genes at the top/bottom to search
#' @param z a number of standard deviations to search
#'
#' @export
limma_ontology = function(limma_out, gene_lengths=NULL, goids=NULL, n=NULL, z=NULL, overwrite=FALSE, goid_map="reference/go/id2go.map", gff_file=NULL, goids_df=NULL, do_goseq=TRUE, do_cluster=TRUE, do_topgo=TRUE, do_trees=FALSE, workbook="excel/ontology.xls", csv=TRUE, excel=FALSE) {
    message("This function expects a list of limma contrast tables and some annotation information.")
    message("The annotation information would be gene lengths and ontology ids")
    if (isTRUE(do_goseq) & is.null(gene_lengths)) {
        stop("Performing a goseq search requires a data frame of gene lengths.")
    }
##    if (isTRUE(do_cluster) & is.null(gff_file)) {
##        stop("Performing a clusterprofiler search requires a gff file.")
##    }
    if (is.null(n) & is.null(z)) {
        z = 1
    }

    if (!is.null(limma_out$all_tables)) {
        limma_out = limma_out$all_tables
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
        if (!is.null(datum$Row.names)) {
            rownames(datum) = datum$Row.names
            datum = datum[-1]
        }
        comparison = names(limma_out[c])
        message(paste("Performing ontology search of:", comparison, sep=""))
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
                testdir = dirname(paste(comparison, csv_base, sep=""))
                if (!file.exists(testdir)) {
                    dir.create(testdir)
                }
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


golevel_df = function(ont="MF", savefile="ontlevel.rda") {
    savefile = paste0(ont, "_", savefile)
##    if (file.exists(savefile)) {
##        load(savefile)
##        return(golevels)
##    } else {
        level = 0
        continue = 1
        golevels = data.frame(GO=NULL,level=NULL)
        while (continue == 1) {
            level = level + 1
            GO = try(clusterProfiler:::getGOLevel(ont, level), silent=TRUE)
            if (class(GO) != 'character') {
                golevels$level = as.numeric(golevels$level)
                save(golevels, file=savefile, compress="xz")
                return (golevels)
            } else {
                tmpdf = as.data.frame(cbind(GO, level))
                ## This (hopefully) ensures that each GO id is added only once, and gets the highest possible level.
                new_go = tmpdf[unique(tmpdf$GO, golevels$GO),]
                golevels = rbind(golevels, new_go)
            }
        }
 ##   }
}
