## Time-stamp: <Mon May 16 16:01:19 2016 Ashton Trey Belew (abelew@gmail.com)>
## Most of the functions in here probably shouldn't be exported...

#' Extract more easily readable information from a GOTERM datum.
#'
#' The output from the GOTERM/GO.db functions is inconsistent, to put it nicely.
#' This attempts to extract from that heterogeneous datatype something easily readable.
#' Example:  Synonym() might return any of the following:
#' NA, NULL, "NA", "NULL", c("NA",NA,"GO:00001"), "GO:00002", c("Some text",NA, NULL, "GO:00003")
#' This function will boil that down to 'not found', '', 'GO:00004', or "GO:0001,  some text, GO:00004"
#'
#' @param value Result of try(as.character(somefunction(GOTERM[id])), silent=TRUE).
#'   somefunction would be 'Synonym' 'Secondary' 'Ontology', etc...
#' @return something more sane (hopefully).
#' @examples
#' \dontrun{
#'  goterms = GOTERM[ids]
#'  sane_goterms = deparse_go_value(goterms)
#' }
#' @export
deparse_go_value <- function(value) {
    result <- ""
    if (class(value) == "try-error") {
        result <- "Not found"
    } else {  ## Not an error
        if (is.null(value)) {
            result <- ""
        } else if (is.na(value)) {
            result <- ""
        } else if (value == "NULL") {
            result <- ""
        } else if (grepl('^c\\(', as.character(value))) {
            value <- eval(parse(text=as.character(value)))
            if (class(value) == "logical") {
                result <- ""
            } else {
                value <- as.character(value[which(complete.cases(value))])
                result <- value
            }
        } else {  ## Just a string "GO:00023409"
            result <- value
        }
    }
    return(result)
}

#' Get a go term from ID.
#'
#' @param go GO id or a list thereof, this may be a character or list(assuming the elements, not
#'     names, are goids).
#' @return Some text containing the terms associated with GO id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  goterm("GO:0032559")
#' ## > GO:0032559
#' ## > "adenyl ribonucleotide binding"
#' }
#' @export
goterm <- function(go="GO:0032559") {
    go <- as.character(go)
    term <- function(id) {
        value <- try(as.character(
            AnnotationDbi::Term(GO.db::GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value <- "not found"
        }
        return(value)
    }
    go <- mapply(term, go)
    return(go)
}

#' Get a go synonym from an ID.
#'
#' I think I will need to do similar parsing of the output for this function as per gosec()
#' In some cases this also returns stuff like c("some text", "GO:someID")
#' versus "some other text"  versus NULL versus NA.
#' This function just goes a mapply(gosn, go).
#'
#' @param go GO id, this may be a character or list(assuming the elements are goids).
#' @return Some text providing the synonyms for the given id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  text =  gosyn("GO:0000001")
#'  text
#' ## > GO:000001
#' ## > "mitochondrial inheritance"
#' }
#' @export
gosyn <- function(go="GO:0000001") {
    go <- as.character(go)
    gosn <- function(go) {
        go <- as.character(go)
        result <- ""
        value <- try(as.character(
            AnnotationDbi::Synonym(GO.db::GOTERM[go])), silent=TRUE)
        result <- paste(deparse_go_value(value), collapse="; ")
        return(result)
    }
    go <- mapply(gosn, go)
    return(go)
}

#' Get a GO secondary ID from an id.
#'
#' Unfortunately, GOTERM's returns for secondary IDs are not consistent, so this function
#' has to have a whole bunch of logic to handle the various outputs.
#'
#' @param go GO ID, this may be a character or list(assuming the elements, not names, are goids).
#' @return Some text comprising the secondary GO id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  gosec("GO:0032432")
#' ## > GO:0032432
#' ## > "GO:0000141" "GO:0030482"
#' }
#' @export
gosec <- function(go="GO:0032432") {
    gosc <- function(go) {
        go <- as.character(go)
        result <- ""
        value <- try(as.character(
            AnnotationDbi::Secondary(GO.db::GOTERM[go])), silent=TRUE)
        result <- deparse_go_value(value)
        return(result)
    }
    go <- mapply(gosc, go)
    return(go)
}

#' Get a go long-form definition from an id.
#'
#' Sometimes it is nice to be able to read the full definition of some GO terms.
#'
#' @param go GO ID, this may be a character or list (assuming the elements are goids).
#' @return Some text providing the long definition of each provided GO id.
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  godef("GO:0032432")
#' ## > GO:0032432
#' ## > "An assembly of actin filaments that are on the same axis but may be oriented with the
#' ## > same or opposite polarities and may be packed with different levels of tightness."
#' }
#' @export
godef <- function(go="GO:0032432") {
    go <- as.character(go)
    def <- function(id) {
        ## This call to AnnotationDbi might be wrong
        value <- try(as.character(
            AnnotationDbi::Definition(GO.db::GOTERM[id])), silent=TRUE)
        if (class(value) == "try-error") {
            value <- "not found"
        }
        return(value)
    }
    go <- mapply(def, go)
    return(go)
}

#'  Get a go ontology name from an ID.
#'
#' @param go GO id, this may be a character or list (assuming the elements are goids).
#' @return The set of ontology IDs associated with the GO ids, thus 'MF' or 'BP' or 'CC'.
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  goont(c("GO:0032432", "GO:0032433"))
#' ## > GO:0032432 GO:0032433
#' ## > "CC" "CC"
#' }
#' @export
goont <- function(go=c("GO:0032432", "GO:0032433")) {
    go <- as.character(go)
    ont <- function(id) {
        value <- try(as.character(AnnotationDbi::Ontology(GO.db::GOTERM[id])),
                     silent=TRUE)
        if (class(value) == "try-error") {
            value <- "not found"
        }
        return(value)
    }
    go <- mapply(ont, go)
    return(go)
}

#' Get a go level approximation from an ID.
#'
#' Sometimes it is useful to know how far up/down the ontology tree a given id resides.  This
#' attmepts to answer that question.
#'
#' @param go GO id, this may be a character or list (assuming the elements are goids).
#' @return Set of numbers corresponding to approximate tree positions of the GO ids.
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  golev("GO:0032559")
#'  ## > 3
#' }
#' @export
golev <- function(go) {
    go <- as.character(go)
    level <- 0
    requireNamespace("GO.db")
    while(class(try(as.character(AnnotationDbi::Ontology(GO.db::GOTERM[[go]])),
                    silent=FALSE)) != 'try-error') {
        ontology <- as.character(AnnotationDbi::Ontology(GO.db::GOTERM[[go]]))
        if (ontology == "MF") {
            ## I am not certain if GO.db:: will work for this
            ancestors <- GO.db::GOMFANCESTOR[[go]]
        } else if (ontology == "BP") {
            ancestors <- GO.db::GOBPANCESTOR[[go]]
        } else if (ontology == "CC") {
            ancestors <- GO.db::GOCCANCESTOR[[go]]
        } else {
            ## There was an error
            message(paste("There was an error getting the ontology: ", as.character(go), sep=""))
            ancestors <- "error"
        }
        go <- ancestors[1]
        level <- level + 1
        if (go == "all") {
            return(level)
        }
    }  ## End while
    return(level)
}
#' Get a go level approximation from a set of IDs.
#'
#' This just wraps golev() in mapply.
#
#' @param go Character list of IDs.
#' @return Set pf approximate levels within the onlogy.
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  golevel(c("GO:0032559", "GO:0000001")
#' ## > 3 4
#' }
#' @export
golevel <- function(go=c("GO:0032559", "GO:0000001")) {
    mapply(golev, go)
}

#' Test GO ids to see if they are useful.
#'
#' This just wraps gotst in mapply.
#'
#' @param go  go IDs as characters.
#' @return Some text
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  gotest("GO:0032559")
#' ## > 1
#'  gotest("GO:0923429034823904")
#' ## > 0
#' }
#' @export
gotest <- function(go) {
    gotst <- function(go) {
        go <- as.character(go)
        value <- try(GO.db::GOTERM[[go]])
        if (class(value) == 'try-error') {
            return(0)
        }
        if (is.null(value)) {
            return(0)
        } else {
            return(1)
        }
    }
    mapply(gotst, go)
}

#' Given a set of goseq data from simple_goseq(),
#' make a list of genes represented in each ontology.
#'
#' This function uses the GO2ALLEG data structure to reverse map ontology
#' categories to a list of genes represented. It therefore assumes that the
#' GO2ALLEG.rda data structure has been deposited in pwd().  This in turn
#' may be generated by clusterProfilers buildGOmap() function if it doesn't
#' exist.  For some species it may also be auto-generated.
#' With little work this can be made much more generic, and it probably should.
#'
#' @param goseq_data List of goseq specific results as generated by simple_goseq().
#' @param ontology Ontology to search (MF/BP/CC).
#' @param pval Maximum accepted pvalue to include in the list of categories to cross reference.
#' @param include_all Include all genes in the ontology search?
#' @return Data frame of categories/genes.
#' @seealso \link{simple_goseq} \code{\link[clusterProfiler]{buildGOmap}},
#' @examples
#' \dontrun{
#'  data = simple_goseq(de_genes=limma_output, lengths=annotation_df, goids=goids_df)
#'  genes_in_cats = gather_genes(data, ont='BP')
#' }
#' @export
gather_genes <- function(goseq_data, ontology='MF', pval=0.05, include_all=FALSE) {
    categories <- NULL
    if (ontology == 'MF') {
        categories <- goseq_data[["mf_subset"]]
    } else if (ontology == 'BP') {
        categories <- goseq_data[["bp_subset"]]
    } else if (ontology == 'CC') {
        categories <- goseq_data[["cc_subset"]]
    } else {
        print('the ont argument didnt make sense, using mf.')
        categories <- goseq_data[["mf_subset"]]
    }
    input <- goseq_data[["input"]]
    ##categories <- subset(categories, over_represented_pvalue <= pval)
    categories <- categories[categories[["over_represented_pvalue"]] <= pval, ]
    cats <- categories[["category"]]
    load("GO2ALLEG.rda")
    GO2ALLEG <- get0("GO2ALLEG")
    genes_per_ont <- function(cat) {
        all_entries <- GO2ALLEG[[cat]]
        entries_in_input <- input[rownames(input) %in% all_entries,]
        names <- as.character(rownames(entries_in_input))
        return(names)
    }
    allgenes_per_ont <- function(cat) {
        all_entries <- GO2ALLEG[[cat]]
        return(all_entries)
    }
    gene_list <- mapply(genes_per_ont, cats)
    all_genes <- mapply(allgenes_per_ont, cats)
    ## require.auto("stringr")
    categories$gene_list <- mapply(stringr::str_c, gene_list, collapse=' ')
    if (isTRUE(include_all)) {
        categories$all_genes <- mapply(stringr::str_c, all_genes, collapse=' ')
    }
    return(categories)
}

#' Make a pvalue plot from a df of IDs, scores, and p-values.
#'
#' This function seeks to make generating pretty pvalue plots as shown by clusterprofiler easier.
#'
#' @param df Some data from topgo/goseq/clusterprofiler.
#' @param ontology  Ontology to plot (MF,BP,CC).
#' @return Ggplot2 plot of pvalues vs. ontology.
#' @seealso \link[goseq]{goseq} \pkg{ggplot2}
#' @export
plot_ontpval <- function(df, ontology="MF") {
    y_name <- paste("Enriched ", ontology, " categories.", sep="")
    pvalue_plot <- ggplot2::ggplot(df, ggplot2::aes_string(x="term", y="score", fill="pvalue")) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::coord_flip() +
        ggplot2::scale_x_discrete(name=y_name) +
##        ggplot2::aes_string(fill="pvalue") +
        ggplot2::scale_fill_continuous(low="red", high="blue") +
        ggplot2::theme(text=ggplot2::element_text(size=10)) +
        ggplot2::theme_bw()
    return(pvalue_plot)
}

#' Perform ontology searches given the results of a differential expression analysis.
#'
#' This takes a set of differential expression results, extracts a subset of up/down expressed
#' genes; passes them to goseq, clusterProfiler, topGO, GOstats, and gProfiler; collects the
#' outputs; and returns them in a (hopefully) consistent fashion.  It attempts to handle the
#' differing required annotation/GOid inputs required for each tool and/or provide supported species
#' in ways which the various tools expect.
#'
#' @param de_out List of topTables comprising limma/deseq/edger outputs.
#' @param gene_lengths Data frame of gene lengths for goseq.
#' @param goids Data frame of goids and genes.
#' @param n Number of genes at the top/bottom of the fold-changes to define 'significant.'
#' @param z Number of standard deviations from the mean fold-change used to define 'significant.'
#' @param fc Log fold-change used to define 'significant'.
#' @param p Maximum pvalue to define 'significant.'
#' @param overwrite Overwrite existing excel results file?
#' @param species Supported organism used by the tools.
#' @param goid_map Mapping file used by topGO, if it does not exist then goids_df creates it.
#' @param gff_file gff file containing the annotations used by gff2genetable from clusterprofiler.
#' @param gff_type Column to use from the gff file for the universe of genes.
#' @param goids_df FIXME! Dataframe of genes and goids which I am relatively certain is no longer needed and superseded by goids.
#' @param do_goseq Perform simple_goseq()?
#' @param do_cluster Perform simple_clusterprofiler()?
#' @param do_topgo Perform simple_topgo()?
#' @param do_gostats Perform simple_gostats()?
#' @param do_gprofiler Perform simple_gprofiler()?
#' @param do_trees make topGO trees from the data?
#' @return a list of up/down ontology results from goseq/clusterprofiler/topgo/gostats, and
#'     associated trees.
#' @examples
#' \dontrun{
#'  many_comparisons = limma_pairwise(expt=an_expt)
#'  tables = many_comparisons$limma
#'  this_takes_forever = limma_ontology(tables, gene_lengths=lengthdb,
#'                                      goids=goids_df, z=1.5, gff_file='length_db.gff')
#' }
#' @export
all_ontology_searches <- function(de_out, gene_lengths=NULL, goids=NULL, n=NULL,
                                  z=NULL, fc=NULL, p=NULL, overwrite=FALSE, species="unsupported",
                                  goid_map="reference/go/id2go.map", gff_file=NULL, gff_type="gene",
                                  goids_df=NULL, do_goseq=TRUE, do_cluster=TRUE,
                                  do_topgo=TRUE, do_gostats=TRUE, do_gprofiler=TRUE, do_trees=FALSE) {
    message("This function expects a list of de contrast tables and some annotation information.")
    message("The annotation information would be gene lengths and ontology ids")
    if (isTRUE(do_goseq) & is.null(gene_lengths)) {
        stop("Performing a goseq search requires a data frame of gene lengths.")
    }
    if (isTRUE(do_cluster) & is.null(gff_file)) {
        stop("Performing a clusterprofiler search requires a gff file.")
    }

    goid_map <- get0('goid_map')
    if (is.null(goid_map)) {
        goid_map <- "reference/go/id2go.map"
    }

    ## Take a moment to detect what the data input looks like
    ## Perhaps a list of tables?
    if (is.null(de_out$all_tables)) {
        if (sum(grepl(pattern="_vs_", x=names(de_out))) == 0) {
            stop("This assumes you are passing it a limma/deseq/edger output from limma_pairwise(), edger_pairwise(), or deseq_pairwise().")
        }
    } else {  ## In this case, you fed it something$limma rather than something$limma$all_tables
        de_out <- de_out$all_tables
    }
    ## Perhaps a single data frame of logFC etc
    ## In which case, coerse it to a list of 1
    if (!is.null(de_out$logFC)) {
        tmp <- de_out
        de_out <- list("first"=tmp)
        rm(tmp)
    }

    output <- list()
    for (c in 1:length(de_out)) {
        datum <- de_out[[c]]
        if (!is.null(datum$Row.names)) {
            rownames(datum) <- datum$Row.names
            datum <- datum[-1]
        }
        comparison <- names(de_out[c])
        message(paste("Performing ontology search of:", comparison, sep=""))
        updown_genes <- get_sig_genes(datum, n=n, z=z, fc=fc, p=p)
        up_genes <- updown_genes$up_genes
        down_genes <- updown_genes$down_genes

        goseq_up_ontology <- goseq_up_trees <- goseq_down_ontology <- goseq_down_trees <- NULL
        cluster_up_ontology <- cluster_up_trees <- cluster_down_ontology <- cluster_down_trees <- NULL
        topgo_up_ontology <- topgo_up_trees <- topgo_down_ontology <- topgo_down_trees <- NULL
        gostats_up_ontology <- gostats_up_trees <- gostats_down_ontology <- gostats_down_trees <- NULL
        gprofiler_up_ontology <- gprofiler_down_ontology <- NULL

        if (isTRUE(do_goseq)) {
            goseq_up_ontology <- try(simple_goseq(up_genes, lengths=gene_lengths, goids=goids))
            goseq_down_ontology <- try(simple_goseq(down_genes, lengths=gene_lengths, goids=goids))
            if (isTRUE(do_trees)) {
                goseq_up_trees <- try(goseq_trees(up_genes, goseq_up_ontology, goid_map=goid_map, goids_df=goids, overwrite=overwrite))
                goseq_down_trees <- try(goseq_trees(down_genes, goseq_down_ontology, goid_map=goid_map, goids_df=goids))
            }
        }

        if (isTRUE(do_cluster)) {
            cluster_up_ontology <- try(simple_clusterprofiler(up_genes, goids=goids, gff=gff_file))
            cluster_down_ontology <- try(simple_clusterprofiler(down_genes, goids=goids, gff=gff_file))
            if (isTRUE(do_trees)) {
                cluster_up_trees <- try(cluster_trees(up_genes, cluster_up_ontology, goid_map=goid_map, goids_df=goids))
                cluster_down_trees <- try(cluster_trees(down_genes, cluster_down_ontology, goid_map=goid_map, goids_df=goids))
            }
        }

        if (isTRUE(do_topgo)) {
            topgo_up_ontology <- try(simple_topgo(up_genes, goid_map=goid_map, goids_df=goids))
            topgo_down_ontology <- try(simple_topgo(down_genes, goid_map=goid_map, goids_df=goids))
            if (isTRUE(do_trees)) {
                topgo_up_trees <- try(topgo_trees(topgo_up_ontology))
                topgo_down_trees <- try(topgo_trees(topgo_down_ontology))
            }
        }

        if (isTRUE(do_gostats)) {
            gostats_up_ontology <- try(simple_gostats(up_genes, gff_file, goids, gff_type=gff_type))
            gostats_down_ontology <- try(simple_gostats(down_genes, gff_file, goids, gff_type=gff_type))
            if (isTRUE(do_trees)) {
                message("gostats_trees has never been tested, this is commented out for the moment.")
                ## topgo_up_trees = try(gostats_trees(topgo_up_ontology))
                ## topgo_down_trees = try(gostats_trees(topgo_down_ontology))
            }
        }

        if (isTRUE(do_gprofiler)) {
            gprofiler_up_ontology <- try(simple_gprofiler(up_genes, species=species))
            gprofiler_down_ontology <- try(simple_gprofiler(down_genes, species=species))
        }
        c_data <- list("up_table" = up_genes,
                       "down_table" = down_genes,
                       "up_goseq" = goseq_up_ontology,
                       "down_goseq" = goseq_down_ontology,
                       "up_cluster" = cluster_up_ontology,
                       "down_cluster" = cluster_down_ontology,
                       "up_topgo" = topgo_up_ontology,
                       "down_topgo" = topgo_down_ontology,
                       "up_gostats" = gostats_up_ontology,
                       "down_gostats" = gostats_down_ontology,
                       "up_gprofiler" = gprofiler_up_ontology,
                       "down_gprofiler" = gprofiler_down_ontology,
                       "up_goseqtrees" = goseq_up_trees,
                       "down_goseqtrees" = goseq_down_trees,
                       "up_clustertrees" = cluster_up_trees,
                       "down_clustertrees" = cluster_down_trees,
                       "up_topgotrees" = topgo_up_trees,
                       "down_topgotrees" = topgo_down_trees,
                       "up_gostatstrees" = gostats_up_trees,
                       "down_gostatstrees" = gostats_down_trees)
        output[[c]] <- c_data
    }
    names(output) <- names(de_out)
    return(output)
}

#' Perform ontology searches on up/down subsets of differential expression data.
#'
#' In the same way all_pairwise() attempts to simplify using multiple DE tools, this function seeks
#' to make it easier to extract subsets of differentially expressed data and pass them to goseq,
#' clusterProfiler, topGO, GOstats, and gProfiler.
#'
#' @param changed_counts List of changed counts as ups and downs.
#' @param doplot Include plots in the results?
#' @param ...  Extra arguments!
#' @return List of ontology search results, up and down for each contrast.
#' @export
subset_ontology_search <- function(changed_counts, doplot=TRUE, ...) {
    up_list <- changed_counts[["ups"]]
    down_list <- changed_counts[["downs"]]
    arglist <- list(...)
    up_goseq <- list()
    down_goseq <- list()
    up_cluster <- list()
    down_cluster <- list()
    up_topgo <- list()
    down_topgo <- list()
    up_gostats <- list()
    down_gostats <- list()
    up_gprofiler <- list()
    down_gprofiler <- list()
    ## goseq() requires minimally gene_lengths and goids
    lengths <- arglist[["lengths"]]
    goids <- arglist[["goids"]]
    gff <- arglist[["gff"]]
    gff_type <- arglist[["gff_type"]]
    if (!is.null(gff)) {
        go2eg <- check_clusterprofiler(gff, gomap=goids)
    }
    types_list <- c("up_goseq","down_goseq","up_cluster","down_cluster",
                    "up_topgo","down_topgo","up_gostats","down_gostats")
    names_list <- names(up_list)
    names_length <- length(names_list)
    for (cluster_count in 1:names_length) {
        name <- names_list[[cluster_count]]
        uppers <- up_list[[cluster_count]]
        downers <- down_list[[cluster_count]]
        message(paste0(cluster_count, "/", names_length, ": Starting goseq"))
        up_goseq[[name]] <- try(simple_goseq(de_genes=uppers,
                                             lengths=lengths,
                                             goids=goids,
                                             doplot=doplot,
                                             ...))
        down_goseq[[name]] <- try(simple_goseq(de_genes=downers,
                                               lengths=lengths,
                                               goids=goids,
                                               doplot=doplot,
                                               ...))
        message(paste0(cluster_count, "/", names_length, ": Starting clusterprofiler"))
        up_cluster[[name]] <- try(simple_clusterprofiler(uppers,
                                                         goids=goids,
                                                         include_cnetplots=FALSE,
                                                         gff=gff,
                                                         ...))
        down_cluster[[name]] <- try(simple_clusterprofiler(downers,
                                                           goids=goids,
                                                           include_cnetplots=FALSE,
                                                           gff=gff,
                                                           ...))
        message(paste0(cluster_count, "/", names_length, ": Starting topgo"))
        up_topgo[[name]] <- try(simple_topgo(de_genes=uppers,
                                             goids_df=goids,
                                             ...))
        down_topgo[[name]] <- try(simple_topgo(de_genes=downers,
                                               goids_df=goids,
                                               ...))
        message(paste0(cluster_count, "/", names_length, ": Starting gostats"))
        up_gostats[[name]] <- try(simple_gostats(uppers,
                                                 gff,
                                                 goids,
                                                 gff_type=gff_type,
                                                 ...))
        down_gostats[[name]] <- try(simple_gostats(downers,
                                                   gff,
                                                   goids,
                                                   gff_type=gff_type,
                                                   ...))
        up_gprofiler[[name]] <- try(simple_gprofiler(uppers,
                                                     ...))
        down_gprofiler[[name]] <- try(simple_gprofiler(downers,
                                                       ...))
    }

    ret <- list(
        "up_goseq" = up_goseq,
        "down_goseq" = down_goseq,
        "up_cluster" = up_cluster,
        "down_cluster" = down_cluster,
        "up_topgo" = up_topgo,
        "down_topgo" = down_topgo,
        "up_gostats" = up_gostats,
        "down_gostats" = down_gostats,
        "up_gprofiler" = up_gprofiler,
        "down_gprofiler" = down_gprofiler)
    if (!file.exists("savefiles")) {
        dir.create("savefiles")
    }
    try(save(list=c("ret"), file="savefiles/subset_ontology_search_result.rda"))
    return(ret)
}

#' Extract a dataframe of golevels using getGOLevel() from clusterProfiler.
#'
#' This function is way faster than my previous iterative golevel function.
#' That is not to say it is very fast, so it saves the result to ontlevel.rda for future lookups.
#'
#' @param ont   the ontology to recurse.
#' @param savefile   a file to save the results for future lookups.
#' @return golevels  a dataframe of goids<->highest level
#' @export
golevel_df <- function(ont="MF", savefile="ontlevel.rda") {
    savefile <- paste0(ont, "_", savefile)
    golevels <- NULL
    if (file.exists(savefile)) {
        load(savefile)
    } else {
        level <- 0
        continue <- 1
        golevels <- data.frame(GO=NULL,level=NULL)
        while (continue == 1) {
            level <- level + 1
            GO <- try(clusterProfiler:::getGOLevel(ont, level), silent=TRUE)
            if (class(GO) != 'character') {
                golevels$level <- as.numeric(golevels$level)
                save(golevels, file=savefile, compress="xz")
                return (golevels)
            } else {
                tmpdf <- as.data.frame(cbind(GO, level))
                ## This (hopefully) ensures that each GO id is added only once, and gets the highest possible level.
                new_go <- tmpdf[unique(tmpdf$GO, golevels$GO),]
                golevels <- rbind(golevels, new_go)
            }
        }
    }
    return(golevels)
}

#' Write gene ontology tables for excel
#'
#' Combine the results from goseq, cluster profiler, topgo, and gostats and drop them into excel
#' Hopefully with a relatively consistent look.
#'
#' @param goseq  The goseq result from simple_goseq()
#' @param cluster The result from simple_clusterprofiler()
#' @param topgo  Guess
#' @param gostats  Yep, ditto
#' @param file   the file to save the results.
#' @param dated   date the excel file
#' @param n   the number of ontology categories to include in each table.
#' @param overwritefile   overwrite an existing excel file
#' @return the list of ontology information
#' @export
write_go_xls <- function(goseq, cluster, topgo, gostats, file="excel/merged_go",
                         dated=TRUE, n=30, overwritefile=TRUE) {
    n <- get0('n')
    if (is.null(n)) {
        n <- 30
    }
    file <- get0('file')
    if (is.null(file)) {
        file <- "excel/merged_go"
    }
    excel_dir <- dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    suffix <- ".xlsx"
    file <- gsub(pattern="\\.xlsx", replacement="", file, perl=TRUE)
    file <- gsub(pattern="\\.xls", replacement="", file, perl=TRUE)
    filename <- NULL
    if (isTRUE(dated)) {
        timestamp <- format(Sys.time(), "%Y%m%d%H")
        filename <- paste0(file, "-", timestamp, suffix)
    } else {
        filename <- paste0(file, suffix)
    }

    if (file.exists(filename)) {
        if (isTRUE(overwritefile)) {
            backup_file(filename)
        }
    }

    ## Massage the goseq tables to match Najib's request
    goseq_mf <- head(goseq$mf_subset, n=n)
    goseq_bp <- head(goseq$bp_subset, n=n)
    goseq_cc <- head(goseq$cc_subset, n=n)
    goseq_mf <- goseq_mf[,c(7,1,6,2,4,5,8)]
    goseq_bp <- goseq_bp[,c(7,1,6,2,4,5,8)]
    goseq_cc <- goseq_cc[,c(7,1,6,2,4,5,8)]
    colnames(goseq_mf) <- c("Ontology","Category","Term","Over p-value", "Num. DE", "Num. in cat.", "Q-value")
    colnames(goseq_bp) <- c("Ontology","Category","Term","Over p-value", "Num. DE", "Num. in cat.", "Q-value")
    colnames(goseq_cc) <- c("Ontology","Category","Term","Over p-value", "Num. DE", "Num. in cat.", "Q-value")

    ## Massage the clusterProfiler tables similarly
    cluster_mf <- head(as.data.frame(cluster$mf_all@result), n=n)
    cluster_bp <- head(as.data.frame(cluster$bp_all@result), n=n)
    cluster_cc <- head(as.data.frame(cluster$cc_all@result), n=n)
    cluster_mf$geneID <- gsub(cluster_mf$geneID, pattern="/", replacement=" ")
    cluster_bp$geneID <- gsub(cluster_bp$geneID, pattern="/", replacement=" ")
    cluster_cc$geneID <- gsub(cluster_cc$geneID, pattern="/", replacement=" ")
    cluster_mf$ontology <- "MF"
    cluster_bp$ontology <- "BP"
    cluster_cc$ontology <- "CC"
    cluster_mf <- cluster_mf[,c(10,1,2,5,3,4,6,7,9,8)]
    cluster_bp <- cluster_bp[,c(10,1,2,5,3,4,6,7,9,8)]
    cluster_cc <- cluster_cc[,c(10,1,2,5,3,4,6,7,9,8)]
    colnames(cluster_mf) <- c("Ontology","Category","Term","Over p-value","Gene ratio",
                              "BG ratio","Adj. p-value","Q-value","Count","Genes")
    colnames(cluster_bp) <- c("Ontology","Category","Term","Over p-value","Gene ratio",
                              "BG ratio","Adj. p-value","Q-value","Count","Genes")
    colnames(cluster_cc) <- c("Ontology","Category","Term","Over p-value","Gene ratio",
                              "BG ratio","Adj. p-value","Q-value","Count","Genes")

    ## Now do the topgo data
    topgo_mf <- head(topgo$tables$mf_interesting, n=n)
    topgo_bp <- head(topgo$tables$bp_interesting, n=n)
    topgo_cc <- head(topgo$tables$cc_interesting, n=n)
    topgo_mf <- topgo_mf[,c(2,1,11,6,7,8,9,10,4,3,5)]
    topgo_bp <- topgo_bp[,c(2,1,11,6,7,8,9,10,4,3,5)]
    topgo_cc <- topgo_cc[,c(2,1,11,6,7,8,9,10,4,3,5)]
    colnames(topgo_mf) <- c("Ontology","Category","Term","Fisher p-value","Q-value","KS score",
                            "EL score","Weight score","Num. DE","Num. in cat.","Exp. in cat.")
    colnames(topgo_bp) <- c("Ontology","Category","Term","Fisher p-value","Q-value","KS score",
                            "EL score","Weight score","Num. DE","Num. in cat.","Exp. in cat.")
    colnames(topgo_cc) <- c("Ontology","Category","Term","Fisher p-value","Q-value","KS score",
                            "EL score","Weight score","Num. DE","Num. in cat.","Exp. in cat.")

    ## And the gostats data
    gostats_mf <- head(gostats$mf_over_all, n=n)
    gostats_bp <- head(gostats$bp_over_all, n=n)
    gostats_cc <- head(gostats$cc_over_all, n=n)
    gostats_mf$t <- gsub(gostats_mf$Term, pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_bp$t <- gsub(gostats_bp$Term, pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_cc$t <- gsub(gostats_cc$Term, pattern=".*\">(.*)</a>", replacement="\\1")
    gostats_mf$Term <- gsub(gostats_mf$Term, pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_bp$Term <- gsub(gostats_bp$Term, pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_cc$Term <- gsub(gostats_cc$Term, pattern="<a href=\"(.*)\">.*", replacement="\\1")
    gostats_mf$ont <- "MF"
    gostats_bp$ont <- "BP"
    gostats_cc$ont <- "CC"
    gostats_mf <- gostats_mf[,c(10,1,9,2,5,6,3,4,8,7)]
    gostats_bp <- gostats_bp[,c(10,1,9,2,5,6,3,4,8,7)]
    gostats_cc <- gostats_cc[,c(10,1,9,2,5,6,3,4,8,7)]
    colnames(gostats_mf) <- c("Ontology","Category","Term","Fisher p-value","Num. DE",
                              "Num. in cat.","Odds ratio","Exp. in cat.","Q-value","Link")
    colnames(gostats_bp) <- c("Ontology","Category","Term","Fisher p-value","Num. DE",
                              "Num. in cat.","Odds ratio","Exp. in cat.","Q-value","Link")
    colnames(gostats_cc) <- c("Ontology","Category","Term","Fisher p-value","Num. DE",
                              "Num. in cat.","Odds ratio","Exp. in cat.","Q-value","Link")

    lst <- list(goseq_mf=goseq_mf, goseq_bp=goseq_bp, goseq_cc=goseq_cc,
                cluster_mf=cluster_mf, cluster_bp=cluster_bp, cluster_cc=cluster_cc,
                topgo_mf=topgo_mf, topgo_bp=topgo_bp, topgo_cc=topgo_cc,
                gostats_mf=gostats_mf, gostats_bp=gostats_bp, gostats_cc=gostats_cc)
    ## require.auto("awalker89/openxlsx")
    wb <- openxlsx::createWorkbook(creator="atb")
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")

    ## This stanza will be repeated so I am just incrementing the new_row
    new_row <- 1
    sheet <- "goseq"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$goseq_bp, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$goseq_bp) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$goseq_mf, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$goseq_mf) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$goseq_cc, tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")

    new_row <- 1
    sheet <- "clusterProfiler"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$cluster_bp, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$cluster_bp) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$cluster_mf, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$cluster_mf) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$cluster_cc, tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

    new_row <- 1
    sheet <- "topgo"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$topgo_bp, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$topgo_bp) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$topgo_mf, tableStyle="TableStyleMedium9", startRow=new_row)
    new_row <- new_row + nrow(lst$topgo_mf) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$topgo_cc, tableStyle="TableStyleMedium9", startRow=new_row)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")

    new_row <- 1
    sheet <- "gostats"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$gostats_bp, tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst$gostats_bp$Link
    class(links) <- 'hyperlink'
    names(links) <- lst$gostats_bp$Category
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    new_row <- new_row + nrow(lst$gostats_bp) + 2
    openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$gostats_mf, tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst$gostats_mf$Link
    class(links) <- 'hyperlink'
    names(links) <- lst$gostats_mf$Category
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    new_row <- new_row + nrow(lst$gostats_mf) + 2
    openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=lst$gostats_cc, tableStyle="TableStyleMedium9", startRow=new_row)
    links <- lst$gostats_cc$Link
    class(links) <- 'hyperlink'
    names(links) <- lst$gostats_cc$Category
    openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")

    res <- openxlsx::saveWorkbook(wb, file, overwrite=TRUE)
    return(res)
}

#'   Write gene ontology tables for data subsets
#'
#' Given a set of ontology results, this attempts to write them to an excel
#' workbook in a consistent and relatively easy-to-read fashion.
#'
#' @param kept_ontology  A result from subset_ontology_search()
#' @param outfile   Workbook to which to write.
#' @param dated    Append the year-month-day-hour to the workbook.
#' @param n   How many ontology categories to write for each search
#' @param overwritefile   Overwrite an existing workbook?
#' @param add_plots   Add the various p-value plots to the end of each sheet?
#' @param table_style   The chosen table style for excel
#' @param ...  some extra parameters
#' @return a set of excel sheet/coordinates
#' @examples
#' \dontrun{
#'  all_contrasts <- all_pairwise(expt, model_batch=TRUE)
#'  keepers <- list(bob = ('numerator','denominator'))
#'  kept <- combine_de_tables(all_contrasts, keepers=keepers)
#'  changed <- extract_significant_genes(kept)
#'  kept_ontologies <- subset_ontology_search(changed, lengths=gene_lengths,
#'                                            goids=goids, gff=gff, gff_type='gene')
#'  go_writer <- write_subset_ontologies(kept_ontologies)
#' }
#' @export
write_subset_ontologies <- function(kept_ontology, outfile="excel/subset_go", dated=TRUE,
                                    n=NULL, overwritefile=TRUE,
                                    add_plots=TRUE, table_style="TableStyleMedium9", ...) {
    arglist <- list(...)
    table_style <- get0('table_style')
    if (is.null(table_style)) {
        table_style <- "TableStyleMedium9"
    }
    n <- get0('n')
    outfile <- get0('outfile')
    if (is.null(outfile)) {
        outfile <- "excel/subset_go"
    }
    excel_dir <- dirname(outfile)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    suffix <- ".xlsx"
    outfile <- gsub(pattern="\\.xlsx", replacement="", outfile, perl=TRUE)
    outfile <- gsub(pattern="\\.xls", replacement="", outfile, perl=TRUE)

    types_list <- c("up_goseq","down_goseq","up_cluster","down_cluster",
                    "up_topgo","down_topgo","up_gostats","down_gostats")
    ## names_list doesn't exist at this point, I losted it
    ## It is buried not very deep in kept_ontology I think
    names_list <- names(kept_ontology[["up_goseq"]])
    count <- 0
    for (name in names_list) {
        count <- count + 1
        up_filename <- paste0(outfile, "_up-", name)
        down_filename <- paste0(outfile, "_down-", name)
        if (isTRUE(dated)) {
            timestamp <- format(Sys.time(), "%Y%m%d%H")
            up_filename <- paste0(up_filename, "-", timestamp, suffix)
            down_filename <- paste0(down_filename, "-", timestamp, suffix)
        } else {
            up_filename <- paste0(up_filename, suffix)
            down_filename <- paste0(down_filename, suffix)
        }
        if (file.exists(up_filename)) {
            if (isTRUE(overwritefile)) {
                backup_file(up_filename)
            }
        }
        if (file.exists(down_filename)) {
            if (isTRUE(overwritefile)) {
                backup_file(down_filename)
            }
        }

        onts <- c("bp","mf","cc")
        up_stuff <- list()
        down_stuff <- list()
        for (ont in onts) {
            ONT <- toupper(ont)
            ## The goseq columns are probably wrong because I dropped one, remember that.
            varname <- paste0(ont, "_subset")
            goseq_up <- kept_ontology$up_goseq[[count]]
            goseq_up_ont <- goseq_up[[varname]]
            if (!is.null(n)) {
                goseq_up_ont <- head(goseq_up_ont, n=n)
            }
            goseq_up_ont <- goseq_up_ont[,c(7,1,6,2,4,5,8)]
            colnames(goseq_up_ont) <- c("Ontology","Category","Term","Over p-value",
                                        "Num. DE", "Num. in cat.", "Q-value")
            goseq_down <- kept_ontology$down_goseq[[count]]
            goseq_down_ont <- goseq_down[[varname]]
            if (!is.null(n)) {
                goseq_down_ont <- head(goseq_down_ont, n=n)
            }
            goseq_down_ont <- goseq_down_ont[,c(7,1,6,2,4,5,8)]
            colnames(goseq_down_ont) <- c("Ontology","Category","Term","Over p-value",
                                          "Num. DE", "Num. in cat.", "Q-value")
            element_name <- paste0("goseq_", ont)
            up_stuff[[element_name]] <- goseq_up_ont
            down_stuff[[element_name]] <- goseq_down_ont

            varname <- paste0(ont, "_all")
            cluster_up <- kept_ontology$up_cluster[[count]]
            cluster_up_ont <- as.data.frame(cluster_up[[varname]]@result)
            if (!is.null(n)) {
                cluster_up_ont <- head(cluster_up_ont, n=n)
            }
            cluster_up_ont$geneID <- gsub(cluster_up_ont$geneID, pattern="/", replacement=" ")
            cluster_up_ont$ontology <- ONT
            cluster_up_ont <- cluster_up_ont[,c(10,1,2,5,3,4,6,7,9,8)]
            colnames(cluster_up_ont) <- c("Ontology","Category","Term","Over p-value","Gene ratio",
                                         "BG ratio","Adj. p-value","Q-value","Count","Genes")
            cluster_down <- kept_ontology$down_cluster[[count]]
            cluster_down_ont <- as.data.frame(cluster_down[[varname]]@result)
            if (!is.null(n)) {
                cluster_down_ont <- head(cluster_down_ont, n=n)
            }
            cluster_down_ont$geneID <- gsub(cluster_down_ont$geneID, pattern="/", replacement=" ")
            cluster_down_ont$ontology <- ONT
            cluster_down_ont <- cluster_down_ont[,c(10,1,2,5,3,4,6,7,9,8)]
            colnames(cluster_down_ont) <- c("Ontology","Category","Term","Over p-value","Gene ratio",
                                            "BG ratio","Adj. p-value","Q-value","Count","Genes")
            element_name <- paste0("cluster_", ont)
            up_stuff[[element_name]] <- cluster_up_ont
            down_stuff[[element_name]] <- cluster_down_ont

            varname <- paste0(ont, "_interesting")
            topgo_up <- kept_ontology$up_topgo[[count]]
            topgo_up_ont <- topgo_up$tables[[varname]]
            if (!is.null(n)) {
                topgo_up_ont <- head(topgo_up_ont, n=n)
            }
            topgo_up_ont <- topgo_up_ont[,c(2,1,11,6,7,8,9,10,4,3,5)]
            colnames(topgo_up_ont) <- c("Ontology","Category","Term","Fisher p-value","Q-value","KS score",
                                        "EL score","Weight score","Num. DE","Num. in cat.","Exp. in cat.")
            topgo_down <- kept_ontology$down_topgo[[count]]
            topgo_down_ont <- topgo_down$tables[[varname]]
            if (!is.null(n)) {
                topgo_down_ont <- head(topgo_down_ont, n=n)
            }
            topgo_down_ont <- topgo_down_ont[,c(2,1,11,6,7,8,9,10,4,3,5)]
            colnames(topgo_down_ont) <- c("Ontology","Category","Term","Fisher p-value","Q-value","KS score",
                                          "EL score","Weight score","Num. DE","Num. in cat.","Exp. in cat.")
            element_name <- paste0("topgo_", ont)
            up_stuff[[element_name]] <- topgo_up_ont
            down_stuff[[element_name]] <- topgo_down_ont

            varname <- paste0(ont, "_over_all")
            gostats_up <- kept_ontology$up_gostats[[count]]
            gostats_up_ont <- gostats_up[[varname]]
            if (!is.null(n)) {
                gostats_up_ont <- head(gostats_up_ont, n=n)
            }
            gostats_up_ont$t <- gsub(gostats_up_ont$Term, pattern=".*\">(.*)</a>", replacement="\\1")
            gostats_up_ont$Term <- gsub(gostats_up_ont$Term, pattern="<a href=\"(.*)\">.*", replacement="\\1")
            gostats_up_ont$ont <- ONT
            gostats_up_ont <- gostats_up_ont[,c(10,1,9,2,5,6,3,4,8,7)]
            colnames(gostats_up_ont) <- c("Ontology","Category","Term","Fisher p-value","Num. DE",
                                          "Num. in cat.","Odds ratio","Exp. in cat.","Q-value","Link")
            gostats_down <- kept_ontology$down_gostats[[count]]
            gostats_down_ont <- gostats_down[[varname]]
            if (!is.null(n)) {
                gostats_down_ont <- head(gostats_down_ont, n=n)
            }
            gostats_down_ont$t <- gsub(gostats_down_ont$Term, pattern=".*\">(.*)</a>", replacement="\\1")
            gostats_down_ont$Term <- gsub(gostats_down_ont$Term, pattern="<a href=\"(.*)\">.*", replacement="\\1")
            gostats_down_ont$ont <- ONT
            gostats_down_ont <- gostats_down_ont[,c(10,1,9,2,5,6,3,4,8,7)]
            colnames(gostats_down_ont) <- c("Ontology","Category","Term","Fisher p-value","Num. DE",
                                            "Num. in cat.","Odds ratio","Exp. in cat.","Q-value","Link")
            element_name <- paste0("gostats_", ont)
            up_stuff[[element_name]] <- gostats_up_ont
            down_stuff[[element_name]] <- gostats_down_ont
        } ## End MF/BP/CC loop

        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
        ## This stanza will be repeated so I am just incrementing the new_row
        new_row <- 1
        sheet <- "goseq"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$goseq_bp, tableStyle=table_style, startRow=new_row)
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_goseq"]][[name]]$bpp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$goseq_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$goseq_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$goseq_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_goseq"]][[name]]$mfp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$goseq_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$goseq_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$goseq_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_goseq"]][[name]]$ccp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$goseq_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
        new_row <- 1
        sheet <- "clusterProfiler"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$cluster_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_cluster"]][[name]]$bp_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$cluster_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$cluster_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$cluster_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_cluster"]][[name]]$mf_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$cluster_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$cluster_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$cluster_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_cluster"]][[name]]$cc_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$cluster_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        new_row <- 1
        sheet <- "topgo"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$topgo_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_topgo"]][[name]]$pvalue_plots$BP
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$topgo_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$topgo_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$topgo_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_topgo"]][[name]]$pvalue_plots$MF
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$topgo_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(up_stuff$topgo_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$topgo_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_topgo"]][[name]]$pvalue_plots$CC
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$topgo_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")
        new_row <- 1
        sheet <- "gostats"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$gostats_bp, tableStyle=table_style, startRow=new_row)
        links <- up_stuff$gostats_bp$Link
        class(links) <- 'hyperlink'
        names(links) <- up_stuff$gostats_bp$Category
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_gostats"]][[name]]$pvalue_plots$bp_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$gostats_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        message("The previous line was a warning about overwriting existing data because of a link.")
        new_row <- new_row + nrow(up_stuff$gostats_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$gostats_mf, tableStyle=table_style, startRow=new_row)
        links <- up_stuff$gostats_mf$Link
        class(links) <- 'hyperlink'
        names(links) <- up_stuff$gostats_mf$Category
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_gostats"]][[name]]$pvalue_plots$mf_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$gostats_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        new_row <- new_row + nrow(up_stuff$gostats_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=up_stuff$gostats_cc, tableStyle=table_style, startRow=new_row)
        links <- up_stuff$gostats_cc$Link
        class(links) <- 'hyperlink'
        names(links) <- up_stuff$gostats_cc$Category
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["up_gostats"]][[name]]$pvalue_plots$cc_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(up_stuff$gostats_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        res <- openxlsx::saveWorkbook(wb, up_filename, overwrite=TRUE)

        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
        ## This stanza will be repeated so I am just incrementing the new_row
        new_row <- 1
        sheet <- "goseq"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$goseq_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_goseq"]][[name]]$bpp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$goseq_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$goseq_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$goseq_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_goseq"]][[name]]$mfp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$goseq_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$goseq_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$goseq_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_goseq"]][[name]]$ccp_plot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$goseq_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
        new_row <- 1
        sheet <- "clusterProfiler"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$cluster_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_cluster"]][[name]]$bp_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$cluster_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$cluster_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$cluster_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_cluster"]][[name]]$mf_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$cluster_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$cluster_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$cluster_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_cluster"]][[name]]$cc_all_barplot
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$cluster_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        new_row <- 1
        sheet <- "topgo"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$topgo_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_topgo"]][[name]]$pvalue_plots$BP
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$topgo_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$topgo_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$topgo_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_topgo"]][[name]]$pvalue_plots$MF
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$topgo_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        new_row <- new_row + nrow(down_stuff$topgo_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$topgo_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_topgo"]][[name]]$pvalue_plots$CC
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$topgo_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:11, widths="auto")
        new_row <- 1
        sheet <- "gostats"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$gostats_bp, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_gostats"]][[name]]$pvalue_plots$bp_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$gostats_bp) + 2, startRow=new_row, fileType="png", units="in")
        }
        links <- down_stuff$gostats_bp$Link
        class(links) <- 'hyperlink'
        names(links) <- down_stuff$gostats_bp$Category
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        new_row <- new_row + nrow(down_stuff$gostats_bp) + 2
        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$gostats_mf, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_gostats"]][[name]]$pvalue_plots$mf_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$gostats_mf) + 2, startRow=new_row, fileType="png", units="in")
        }
        links <- down_stuff$gostats_mf$Link
        class(links) <- 'hyperlink'
        names(links) <- down_stuff$gostats_mf$Category
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        new_row <- new_row + nrow(down_stuff$gostats_mf) + 2
        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        openxlsx::writeDataTable(wb, sheet, x=down_stuff$gostats_cc, tableStyle=table_style, startRow=new_row)
        if (isTRUE(add_plots)) {
            a_plot <- kept_ontology[["down_gostats"]][[name]]$pvalue_plots$cc_plot_over
            print(a_plot)
            openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(down_stuff$gostats_cc) + 2, startRow=new_row, fileType="png", units="in")
        }
        links <- down_stuff$gostats_cc$Link
        class(links) <- 'hyperlink'
        names(links) <- down_stuff$gostats_cc$Category
        openxlsx::writeData(wb, sheet, x=links, startRow=new_row + 1, startCol=10)
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:9, widths="auto")
        res <- openxlsx::saveWorkbook(wb, down_filename, overwrite=TRUE)
    }  ## End of name_list
}

#' Compare the results from different ontology tools
#'
#' Combine the results from goseq, cluster profiler, topgo, and gostats; poke at them
#' with a stick and see what happens.
#' The general idea is to pull the p-value data from each tool and contrast that to the
#' set of all possibile ontologies.  This allows one to do a correlation coefficient
#' between them.  In addition, take the 1-pvalue for each ontology for each tool.
#' Thus for strong p-values the score will be near 1 and so we can sum the scores
#' for all the tools.  Since topgo has 4 tools, the total possible is 7 if everything
#' has a p-value equal to 0.
#'
#' @param goseq   The goseq result from simple_goseq()
#' @param cluster   The result from simple_clusterprofiler()
#' @param topgo   Guess
#' @param gostats   Yep, ditto
#' @return a summary of the similarities of ontology searches
#' @export
compare_go_searches <- function(goseq=NULL, cluster=NULL, topgo=NULL, gostats=NULL) {
    goseq_mf_data <- goseq_bp_data <- goseq_cc_data <- NULL
    if (!is.null(goseq)) {
        goseq_mf_data <- goseq$mf_subset[,c("category","over_represented_pvalue")]
        goseq_bp_data <- goseq$bp_subset[,c("category","over_represented_pvalue")]
        goseq_cc_data <- goseq$cc_subset[,c("category","over_represented_pvalue")]
        colnames(goseq_mf_data) <- c("goseq_id","goseq_pvalue")
        colnames(goseq_bp_data) <- c("goseq_id","goseq_pvalue")
        colnames(goseq_cc_data) <- c("goseq_id","goseq_pvalue")
    }
    cluster_mf_data <- cluster_bp_data <- cluster_cc_data <- NULL
    if (!is.null(cluster)) {
        cluster_mf_data <- as.data.frame(summary(cluster$mf_all))[,c("ID","pvalue")]
        cluster_bp_data <- as.data.frame(summary(cluster$bp_all))[,c("ID","pvalue")]
        cluster_cc_data <- as.data.frame(summary(cluster$cc_all))[,c("ID","pvalue")]
        colnames(cluster_mf_data) <- c("cluster_id","cluster_pvalue")
        colnames(cluster_bp_data) <- c("cluster_id","cluster_pvalue")
        colnames(cluster_cc_data) <- c("cluster_id","cluster_pvalue")
    }
    topgo_mf_data <- topgo_bp_data <- topgo_cc_data <- NULL
    if (!is.null(topgo)) {
        topgo_mf_data <- topgo$tables$mf[,c("GO.ID","fisher","KS","EL","weight")]
        topgo_bp_data <- topgo$tables$bp[,c("GO.ID","fisher","KS","EL","weight")]
        topgo_cc_data <- topgo$tables$cc[,c("GO.ID","fisher","KS","EL","weight")]
        colnames(topgo_mf_data) <- c("topgo_id","topgo_fisher_pvalue","topgo_el_pvalue","topgo_ks_pvalue","topgo_weight_pvalue")
        colnames(topgo_bp_data) <- c("topgo_id","topgo_fisher_pvalue","topgo_el_pvalue","topgo_ks_pvalue","topgo_weight_pvalue")
        colnames(topgo_cc_data) <- c("topgo_id","topgo_fisher_pvalue","topgo_el_pvalue","topgo_ks_pvalue","topgo_weight_pvalue")
    }
    gostats_mf_data <- gostats_bp_data <- gostats_cc_data <- NULL
    if (!is.null(gostats)) {
        gostats_mf_data <- gostats$mf_over_all[,c("GOMFID","Pvalue")]
        gostats_bp_data <- gostats$bp_over_all[,c("GOBPID","Pvalue")]
        gostats_cc_data <- gostats$cc_over_all[,c("GOCCID","Pvalue")]
        colnames(gostats_mf_data) <- c("gostats_id","gostats_pvalue")
        colnames(gostats_bp_data) <- c("gostats_id","gostats_pvalue")
        colnames(gostats_cc_data) <- c("gostats_id","gostats_pvalue")
    }
    ## Now combine them...
    all_mf <- merge(goseq_mf_data, cluster_mf_data, by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
    all_mf <- merge(all_mf, topgo_mf_data, by.x="goseq_id", by.y="topgo_id", all.x=TRUE, all.y=TRUE)
    all_mf <- merge(all_mf, gostats_mf_data, by.x="goseq_id", by.y="gostats_id", all.x=TRUE, all.y=TRUE)
    rownames(all_mf) <- all_mf$goseq_id
    all_mf <- all_mf[-1]
    all_mf[is.na(all_mf)] <- 1
    all_mf$goseq_pvalue <- 1 - as.numeric(all_mf$goseq_pvalue)
    all_mf$cluster_pvalue <- 1 - as.numeric(all_mf$cluster_pvalue)
    all_mf$topgo_fisher_pvalue <- 1 - as.numeric(all_mf$topgo_fisher_pvalue)
    all_mf$topgo_el_pvalue <- 1 - as.numeric(all_mf$topgo_el_pvalue)
    all_mf$topgo_ks_pvalue <- 1 - as.numeric(all_mf$topgo_ks_pvalue)
    all_mf$topgo_weight_pvalue <- 1 - as.numeric(all_mf$topgo_weight_pvalue)
    all_mf$gostats_pvalue <- 1 - as.numeric(all_mf$gostats_pvalue)
    all_mf[is.na(all_mf)] <- 0
    message(cor(all_mf))
    mf_summary <- rowSums(all_mf)
    message(summary(mf_summary))

    all_bp <- merge(goseq_bp_data, cluster_bp_data, by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
    all_bp <- merge(all_bp, topgo_bp_data, by.x="goseq_id", by.y="topgo_id", all.x=TRUE, all.y=TRUE)
    all_bp <- merge(all_bp, gostats_bp_data, by.x="goseq_id", by.y="gostats_id", all.x=TRUE, all.y=TRUE)
    rownames(all_bp) <- all_bp$goseq_id
    all_bp <- all_bp[-1]
    all_bp[is.na(all_bp)] <- 1
    all_bp$goseq_pvalue <- 1 - as.numeric(all_bp$goseq_pvalue)
    all_bp$cluster_pvalue <- 1 - as.numeric(all_bp$cluster_pvalue)
    all_bp$topgo_fisher_pvalue <- 1 - as.numeric(all_bp$topgo_fisher_pvalue)
    all_bp$topgo_el_pvalue <- 1 - as.numeric(all_bp$topgo_el_pvalue)
    all_bp$topgo_ks_pvalue <- 1 - as.numeric(all_bp$topgo_ks_pvalue)
    all_bp$topgo_weight_pvalue <- 1 - as.numeric(all_bp$topgo_weight_pvalue)
    all_bp$gostats_pvalue <- 1 - as.numeric(all_bp$gostats_pvalue)
    all_bp[is.na(all_bp)] <- 0
    message(cor(all_bp))
    bp_summary <- rowSums(all_bp)
    message(summary(bp_summary))

    all_cc <- merge(goseq_cc_data, cluster_cc_data, by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
    all_cc <- merge(all_cc, topgo_cc_data, by.x="goseq_id", by.y="topgo_id", all.x=TRUE, all.y=TRUE)
    all_cc <- merge(all_cc, gostats_cc_data, by.x="goseq_id", by.y="gostats_id", all.x=TRUE, all.y=TRUE)
    rownames(all_cc) <- all_cc$goseq_id
    all_cc <- all_cc[-1]
    all_cc[is.na(all_cc)] <- 1
    all_cc$goseq_pvalue <- 1 - as.numeric(all_cc$goseq_pvalue)
    all_cc$cluster_pvalue <- 1 - as.numeric(all_cc$cluster_pvalue)
    all_cc$topgo_fisher_pvalue <- 1 - as.numeric(all_cc$topgo_fisher_pvalue)
    all_cc$topgo_el_pvalue <- 1 - as.numeric(all_cc$topgo_el_pvalue)
    all_cc$topgo_ks_pvalue <- 1 - as.numeric(all_cc$topgo_ks_pvalue)
    all_cc$topgo_weight_pvalue <- 1 - as.numeric(all_cc$topgo_weight_pvalue)
    all_cc$gostats_pvalue <- 1 - as.numeric(all_cc$gostats_pvalue)
    all_cc[is.na(all_cc)] <- 0
    message(cor(all_cc))
    cc_summary <- rowSums(all_cc)
    message(summary(cc_summary))
}

## EOF
