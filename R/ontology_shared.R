## Time-stamp: <Wed Sep  2 16:26:49 2015 Ashton Trey Belew (abelew@gmail.com)>
## Most of the functions in here probably shouldn't be exported...

#' deparse_go_value()  Extract more easily readable information from a GOTERM datum.
#'
#' The output from the GOTERM/GO.db functions is inconsistent, to put it nicely.
#' This attempts to extract from that heterogeneous datatype something easily readable.
#' Example:  Synonym() might return any of the following:
#' NA, NULL, "NA", "NULL", c("NA",NA,"GO:00001"), "GO:00002", c("Some text",NA, NULL, "GO:00003")
#' This function will boil that down to 'not found', '', 'GO:00004', or "GO:0001,  some text, GO:00004"
#'
#' @param value  the result of try(as.character(somefunction(GOTERM[id])), silent=TRUE)
#'   somefunction would be 'Synonym' 'Secondary' 'Ontology', etc...
#'
#' @return something more sane (hopefully)
#' @export
#' @examples
#' ## goterms = GOTERM[ids]
#' ## sane_goterms = deparse_go_value(goterms)
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

#' goterm()  Get a go term from ID.
#'
#' @param go default='GO:0032559'  a go ID or list thereof
#' this may be a character or list(assuming the elements, not names, are goids)
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

#' gosyn()  Get a go synonym from an ID.
#'
#' I think I will need to do similar parsing of the output for this function as per gosec()
#' In some cases this also returns stuff like c("some text", "GO:someID")
#' versus "some other text"  versus NULL versus NA
#'
#' This function just goes a mapply(gosn, go).
#'
#' @param go  a go ID, this may be a character or list(assuming the elements are goids).
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
#' gosn()  Does the actual go synonym extraction for gosyn()
#'
#' This merely calls as.character(Synonym(GOTERM[go]))
#'
#' @param go  a go ID.
#' @return some text
gosn = function(go) {
    go = as.character(go)
    result = ""
    value = try(as.character(AnnotationDbi::Synonym(GOTERM[go])), silent=TRUE)
    result = paste(deparse_go_value(value), collapse="; ")
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
#' gosc()  does the real work for gosec()
#'
#' @param go  a go id.
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

#' godef()  Get a go long-form definition from an id.
#'
#' @param go  a go ID, this may be a character or list (assuming the elements are goids).
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

#' goont()  Get a go ontology name from an ID.
#'
#' @param go  a go ID, this may be a character or list (assuming the elements are goids).
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

#' golev()  Get a go level approximation from an ID.
#'
#' @param go  a go ID, this may be a character or list (assuming the elements are goids).
#' @param verbose default=FALSE  print some information as it recurses.
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
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
#' golevel()  Get a go level approximation from a set of IDs.
#' This just wraps golev() in mapply.
#' @param go  a character list of IDs.
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

#' gotst()  Test a GO id to see if it is useful by my arbitrary definition of 'useful'.
#'
#' @param go  a go ID, this may be a character or list (assuming the elements are goids).
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
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
#' gotest()  Test GO ids to see if they are useful.
#' This just wraps gotst in mapply.
#'
#' @param go  go IDs as characters.
#'
#' @return Some text
#' @seealso \code{\link{GOTERM}}, \code{\link{GO.db}},
#'
#' @export
gotest = function(go) {
    mapply(gotst, go)
}

#' gather_genes()  Given a set of goseq data from simple_goseq(), make a list of genes represented in each ontology.
#'
#' This function uses the GO2ALLEG data structure to reverse map ontology categories to a list of genes represented.
#' It therefore assumes that the GO2ALLEG.rda data structure has been deposited in pwd().  This in turn may be generated
#' by clusterProfilers buildGOmap() function if it doesn't exist.  For some species it may also be auto-generated.
#' With little work this can be made much more generic, and it probably should.
#'
#' @param goseq_data  a list of goseq specific results as generated by simple_goseq()
#' @param ont default='MF'  an ontology to search
#' @param pval default=0.05  a maximum accepted pvalue to include in the list of categories to cross reference.
#'
#' @return a data frame of categories/genes.
#' @seealso \code{\link{simple_goseq}}, \code{\link{buildGOmap}},
#'
#' @export
#' @examples
#' ## data = simple_goseq(de_genes=limma_output, lengths=annotation_df, goids=goids_df)
#' ## genes_in_cats = gather_genes(data, ont='BP')
gather_genes = function(goseq_data, ontology='MF', pval=0.05, include_all=FALSE) {
    categories = NULL
    if (ontology == 'MF') {
        categories = goseq_data$mf_subset
    } else if (ontology == 'BP') {
        categories = goseq_data$bp_subset
    } else if (ontology == 'CC') {
        categories = goseq_data$cc_subset
    } else {
        print('the ont argument didnt make sense, using mf.')
        categories = goseq_data$mf_subset
    }
    input = goseq_data$input
    categories = subset(categories, over_represented_pvalue <= pval)
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
    gene_list = mapply(genes_per_ont, cats)
    all_genes = mapply(allgenes_per_ont, cats)
    require.auto("stringr")
    categories$gene_list = mapply(str_c, gene_list, collapse=' ')
    if (isTRUE(include_all)) {
        categories$all_genes = mapply(str_c, all_genes, collapse=' ')
    }
    return(categories)
}

#' pval_plot()  Make a pvalue plot from a df of IDs, scores, and p-values.
#'
#' This function seeks to make generating pretty pvalue plots as shown by clusterprofiler easier.
#' @param df  some data from topgo/goseq/clusterprofiler.
#' @param ontology default='MF'  an ontology to plot (MF,BP,CC).
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

#' all_ontology_searches()  Perform ontology searches of the output from limma.
#'
#' This passes a set of limma results to (optionally) goseq, clusterprofiler, topgo, and gostats,
#' collects the outputs, and provides them as a list.  This function needs a species argument,
#' as I recently made the simple_() functions able to automatically use the various supported organisms.
#' 
#' @param de_out  a list of topTables comprising limma/deseq/edger outputs.
#' @param n default=NULL  a number of genes at the top/bottom to search.
#' @param z default=NULL  a number of standard deviations to search. (if this and n are null, it assumes 1z)
#' @param gene_lengths default=NULL  a data frame of gene lengths for goseq.
#' @param goids default=NULL  a data frame of goids and genes.
#' @param overwrite default=FALSE  a boolean of whether to overwrite the id2go mapping for goseq.
#' @param goid_map default='reference/go/id2go.map'  a map file used by topGO, if it does not exist then provide goids_df to make it.
#' @param gff_file default=NULL  a gff file containing the annotations used by gff2genetable from clusterprofiler, which I hacked to make faster.
#' @param goids_df default=NULL  FIXME! a dataframe of genes and goids which I am relatively certain is no longer needed and superseded by goids.
#' @param do_goseq default=TRUE  perform simple_goseq()?
#' @param do_cluster default=TRUE  perform simple_clusterprofiler()?
#' @param do_topgo default=TRUE  perform simple_topgo()?
#' @param do_gostats default=TRUE  perform simple_gostats()?
#' @param do_trees default=FALSE  make topGO trees from the data?
#' @param workbook default='excel/ontology.xls'  generate an excel workbook of the ontology data?
#' @param csv default=TRUE  generate csv files using excel/ontology.csv as a basename
#' @param excel default=FALSE  generate the excel workbook?
#'
#' @return a list of up/down ontology results from goseq/clusterprofiler/topgo/gostats, and associated trees, all optionally.
#' @export
#' @examples
#' ## many_comparisons = limma_pairwise(expt=an_expt)
#' ## tables = many_comparisons$limma_result
#' ## this_takes_forever = limma_ontology(tables, gene_lengths=lengthdb, goids=goids_df, z=1.5, gff_file='length_db.gff')
all_ontology_searches = function(de_out, gene_lengths=NULL, goids=NULL, n=NULL, z=NULL, fc=NULL, overwrite=FALSE, goid_map="reference/go/id2go.map", gff_file=NULL, goids_df=NULL, do_goseq=TRUE, do_cluster=TRUE, do_topgo=TRUE, do_gostats=TRUE, do_trees=FALSE, workbook="excel/ontology.xls", csv=TRUE, excel=FALSE) {
    message("This function expects a list of de contrast tables and some annotation information.")
    message("The annotation information would be gene lengths and ontology ids")
    if (isTRUE(do_goseq) & is.null(gene_lengths)) {
        stop("Performing a goseq search requires a data frame of gene lengths.")
    }
    ## if (isTRUE(do_cluster) & is.null(gff_file)) {
    ##  stop("Performing a clusterprofiler search requires a gff file.")
    ## }

    goid_map = get0('goid_map')
    if (is.null(goid_map)) {
        goid_map = "reference/go/id2go.map"
    }
    
    ## Take a moment to detect what the data input looks like
    ## Perhaps a list of tables?
    if (!is.null(de_out$all_tables)) {
        de_out = de_out$all_tables
    }
    ## Perhaps a single data frame of logFC etc
    ## In which case, coerse it to a list of 1
    if (!is.null(de_out$logFC)) {
        tmp = de_out
        de_out = list("first"=tmp)
        rm(tmp)
    }

    testdir = dirname(workbook)
    if (isTRUE(excel) | isTRUE(csv)) {
        if (!file.exists(testdir)) {
            dir.create(testdir)
            message(paste("Creating directory: ", testdir, "for writing excel/csv data.", sep=""))
        }
    }

    output = list()
    for (c in 1:length(de_out)) {
        datum = de_out[[c]]
        if (!is.null(datum$Row.names)) {
            rownames(datum) = datum$Row.names
            datum = datum[-1]
        }
        comparison = names(de_out[c])
        message(paste("Performing ontology search of:", comparison, sep=""))
        updown_genes = get_sig_genes(datum, n=n, z=z, fc=fc)
        up_genes = updown_genes$up_genes
        down_genes = updown_genes$down_genes

        goseq_up_ontology = goseq_up_trees = goseq_down_ontology = goseq_down_trees = NULL
        cluster_up_ontology = cluster_up_trees = cluster_down_ontology = cluster_down_trees = NULL
        topgo_up_ontology = topgo_up_trees = topgo_down_ontology = topgo_down_trees = NULL
        gostats_up_ontology = gostats_up_trees = gostats_down_ontology = gostats_down_trees = NULL
        if (isTRUE(do_goseq)) {
            goseq_up_ontology = try(simple_goseq(up_genes, lengths=gene_lengths, goids=goids))
            goseq_down_ontology = try(simple_goseq(down_genes, lengths=gene_lengths, goids=goids))
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
            cluster_up_ontology = try(simple_clusterprofiler(up_genes, goids=goids, gff=goids))
            cluster_down_ontology = try(simple_clusterprofiler(down_genes, goids=goids, gff=goids))
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
            topgo_up_ontology = try(simple_topgo(up_genes, goid_map=goid_map, goids_df=goids))
            topgo_down_ontology = try(simple_topgo(down_genes, goid_map=goid_map, goids_df=goids))
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
        if (isTRUE(do_gostats)) {
            topgo_up_ontology = try(simple_gostats(up_genes, gff, goids))
            topgo_down_ontology = try(simple_gostats(down_genes, gff, goids))
            if (isTRUE(do_trees)) {
                message("gostats_trees has never been tested, this is commented out for the moment.")
                ## topgo_up_trees = try(gostats_trees(topgo_up_ontology))
                ## topgo_down_trees = try(gostats_trees(topgo_down_ontology))
            }
            if (isTRUE(excel)) {
                sheetname = paste(comparison, "_topgo_mf_up", sep="")
               try(write_xls(data=gostats_up_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_bp_up", sep="")
                try(write_xls(data=gostats_up_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_cc_up", sep="")
                try(write_xls(data=gostats_up_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_mf_down", sep="")
                try(write_xls(data=gostats_down_ontology$mf_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_bp_down", sep="")
                try(write_xls(data=gostats_down_ontology$bp_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
                sheetname = paste(comparison, "_topgo_cc_down", sep="")
                try(write_xls(data=gostats_down_ontology$cc_interesting, sheet=sheetname, file=workbook, overwrite=TRUE))
            }
            if (isTRUE(csv)) {
                csv_base = gsub(".xls$", "", workbook)
                write.csv(gostats_up_ontology$mf_interesting, file=paste(comparison, csv_base, "_topgo_mf_up.csv", sep=""))
                write.csv(gostats_up_ontology$bp_interesting, file=paste(comparison, csv_base, "_topgo_bp_up.csv", sep=""))
                write.csv(gostats_up_ontology$cc_interesting, file=paste(comparison, csv_base, "_topgo_cc_up.csv", sep=""))
                write.csv(gostats_down_ontology$mf_interesting, file=paste(comparison, csv_base, "_topgo_mf_down.csv", sep=""))
                write.csv(gostats_down_ontology$bp_interesting, file=paste(comparison, csv_base, "_topgo_bp_down.csv", sep=""))
                write.csv(gostats_down_ontology$cc_interesting, file=paste(comparison, csv_base, "_topgo_cc_down.csv", sep=""))
            }
        }        
        c_data = list(up_table=up_genes, down_table=down_genes,
                      up_goseq=goseq_up_ontology, down_goseq=goseq_down_ontology,
                      up_cluster=cluster_up_ontology, down_cluster=cluster_down_ontology,
                      up_topgo=topgo_up_ontology, down_topgo=topgo_down_ontology,
                      up_goseqtrees=goseq_up_trees, down_goseqtrees=goseq_down_trees,
                      up_clustertrees=cluster_up_trees, down_clustertrees=cluster_down_trees,
                      up_topgotrees=topgo_up_trees, down_topgotrees=topgo_down_trees,
                      up_gostats=gostats_up_ontology, down_gostats=gostats_down_ontology,
                      up_gostatstrees=gostats_up_trees, down_gostatstrees=gostats_down_trees)
        output[[c]] = c_data
    }
    names(output) = names(de_out)
    return(output)
}

#' golevel_df()  Extract a dataframe of golevels using getGOLevel() from clusterProfiler.
#'
#' This function is way faster than my previous iterative golevel function.
#' That is not to say it is very fast, so it saves the result to ontlevel.rda for future lookups.
#'
#' @param ont default='MF'  the ontology to recurse.
#' @param savefile default='ontlevel.rda'  a file to save the results for future lookups.
#' @return golevels  a dataframe of goids<->highest level
#' @export
golevel_df = function(ont="MF", savefile="ontlevel.rda") {
    savefile = paste0(ont, "_", savefile)
    golevels = NULL
    if (file.exists(savefile)) {
        load(savefile)
    } else {
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
    }
    return(golevels)
}

## EOF

