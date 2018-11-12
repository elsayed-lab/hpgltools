#' Perform a simple_ontology() on some random data.
#'
#' At the very least, the result should be less significant than the actual data!
#'
#' @param input  Some input data
#' @param method  goseq, clusterp, topgo, gostats, gprofiler.
#' @param n  how many 'genes' to analyse?
#' @param ...  Arguments passed to the method.
#' @return  An ontology result
#' @export
random_ontology <- function(input, method="goseq", n=200, ...) {
  ## Lets assume the result of *_pairwise() or combine_de_tables()
  input_table <- NULL
  if (class(input) == "expt" | class(input) == "ExpressionSet") {
    input_table <- data.frame(row.names=rownames(exprs(input)))
    input_table[["ID"]] <- rownames(input_table)
    input_table[["DE"]] <- 1
  } else if (!is.null(input[["data"]])) {
    ## Then it is from combine_de_tables

  } else if (!is.null()) {
    ## Then it came from *_pairwise

  } else {
    stop("Not sure what to do with this input.")
  }

  input_idx <- sample(x=nrow(input_table), size=n)
  input_table <- input_table[input_idx, ]
  random_result <- NULL
  switchret <- switch(
    method,
    "goseq" = {
      random_result <- simple_goseq(input_table, ...)
    },
    "clusterp" = {
      random_result <- simple_clusterprofiler(input_table, ...)
    },
    "topgo" = {
      random_result <- simple_topgo(input_table, ...)
    },
    "gostats" = {
      random_result <- simple_gostats(input_table, ...)
    },
    "gprofiler" = {
      random_result <- simple_gprofiler(input_table, ...)
    },
    {
      message("Not sure what to do with this method.")
      random_result <- NULL
    })

  return(random_result)
}

#' Take gene/exon lengths from a suitable data source (gff/TxDb/OrganismDbi)
#'
#' Primarily goseq, but also other tools on occasion require a set of gene IDs
#' and lengths. This function is resposible for pulling that data from either a
#' gff, or TxDb/OrganismDbi.
#'
#' @param db Object containing data, if it is a string then a filename is
#'   assumed to a gff file.
#' @param gene_list Set of genes to query.
#' @param type Function name used for extracting data from TxDb objects.
#' @param id Column from the resulting data structure to extract gene IDs.
#' @param possible_types Character list of types I have previously used.
#' @param ...  More arguments are passed to arglist.
#' @return Dataframe containing 2 columns: ID, length
#' @seealso \pkg{GenomicFeatures}
#' @export
extract_lengths <- function(db=NULL, gene_list=NULL,
                            type = "GenomicFeatures::transcripts", id="TXID",
                            possible_types=c("GenomicFeatures::genes",
                                             "GenomicFeatures::cds",
                                             "GenomicFeatures::transcripts"), ...) {
  arglist <- list(...)
  ## The 3 ids correspond to the columns produced by genes/cds/transcripts
  ## respectively which contain the IDs. If one is overwritten, the other
  ## should be, too.

  ## Here is the fundamental problem with extracting lengths from these
  ## databases: If we have 1 gene ID, what do we associate the length to?  A
  ## single transcript, the entire gene's length before transcription?, the
  ## length before splicing? All of these are possible given the databases of
  ## gene lengths we have available and it is not always trivial/possible to
  ## tell which is correct because we cannot be certain which ID type has been
  ## provided.  Therefore, this function will query every type and try to find
  ## the one with the best overlap against the set of gene IDs provided.
  tmpdb <- db
  metadf <- NULL
  gene_list <- gene_list[complete.cases(gene_list)]
  chosen_column <- NULL
  hits_list <- list()
  column_list <- list()
  chosen_type <- NULL
  most_hits <- 0
  ## Translating to ENTREZIDs sometimes introduces NAs which messes up the
  ## following operations. Possible_types should be a listing of the various
  ## methods we have of acquiring gene lengths. The code in the for loop
  ## should therefore invoke each of these in turn and figure out which
  ## provides the best overlap and use that.
  for (c in 1:length(possible_types)) {
    testing <- NULL
    ty <- possible_types[c]
    ## make a granges/iranges using the function in possible_types.
    test_string <- glue("testing <- {ty}(tmpdb)")
    eval(parse(text=test_string))
    ## as.data.frame is not only base, but also biocgenerics!!!
    ## Make a dataframe out of the information above and find the most
    ## likely appropriate ID column
    test_meta <- BiocGenerics::as.data.frame(testing)
    possible_columns <- colnames(test_meta)
    ## It turns out this is pretty much always the last column.
    chosen_column <- possible_columns[length(possible_columns)]
    ## Find the overlap, and if this is largest than our current best overlap...
    overlap <- gene_list %in% test_meta[[chosen_column]]
    hits_list[[ty]] <- testing
    column_list[[ty]] <- chosen_column
    message("Testing ", ty, " with column ", chosen_column, " an overlap of ",
            sum(overlap), " was observed out of ", length(gene_list), " genes.")
    ## Note it as the best type so far.
    if (sum(overlap) > most_hits) {
      chosen_type <- ty
      most_hits <- sum(overlap)
    }
  }
  message("Actually using type ", chosen_type,
          " consider one of the above if that is not good enough.")
  ## Now we have a list of all the lengths by function used to acquire them.
  testing <- hits_list[[chosen_type]]
  ## We have a second list of the columns containing the appropriate IDs.
  chosen_column <- column_list[[chosen_type]]
  ## So, bring them together here.
  meta <- BiocGenerics::as.data.frame(testing)
  if (!is.null(meta[["width"]])) {
    metadf <- as.data.frame(meta)[, c(chosen_column, "width")]
  } else if (!is.null(test_meta[["length"]])) {
    metadf <- as.data.frame(meta)[, c(chosen_column, "length")]
  } else {
    stop("This requires either length or width columns.")
  }
  colnames(metadf) <- c("ID", "length")
  rownames(metadf) <- metadf[["ID"]]
  return(metadf)
}

#' Extract a set of geneID to GOID mappings from a suitable data source.
#'
#' Like extract_lengths above, this is primarily intended to read gene ID and GO
#'  ID mappings from a OrgDb/OrganismDbi object.
#'
#' @param db Data source containing mapping information.
#' @param metadf Data frame containing extant information.
#' @param keytype Keytype used for querying
#' @return Dataframe of 2 columns: geneID and goID.
#' @seealso \pkg{AnnotationDbi}
#' @export
extract_go <- function(db, metadf=NULL, keytype="ENTREZID") {
  keytype <- toupper(keytype)
  possible_keytypes <- AnnotationDbi::keytypes(db)
  godf <- data.frame()
  success <- FALSE
  ids <- AnnotationDbi::keys(x=db, keytype=keytype)
  if ("GOID" %in% possible_keytypes) {
    godf <- sm(AnnotationDbi::select(x=db, keys=ids, keytype=keytype,
                                     columns=c("GOID")))
    godf[["ID"]] <- godf[[1]]
    godf <- godf[, c("ID", "GOID")]
    colnames(godf) <- c("ID", "GO")
  } else if ("GO" %in% possible_keytypes) {
    godf <- sm(AnnotationDbi::select(x=db, keys=ids, keytype=keytype,
                                     columns=c("GO")))
    godf[["ID"]] <- godf[[1]]
    godf <- godf[, c("ID", "GO")]
    colnames(godf) <- c("ID", "GO")
  } else if ("GOALL" %in% possible_keytypes) {
    godf <- sm(AnnotationDbi::select(x=db, keys=ids, keytype=keytype,
                                     columns=c("GOALL")))
    godf[["ID"]] <- godf[[1]]
    godf <- godf[, c("ID", "GO")]
    colnames(godf) <- c("ID", "GO")
  }
  return(godf)
}

#' Extract more easily readable information from a GOTERM datum.
#'
#' The output from the GOTERM/GO.db functions is inconsistent, to put it
#' nicely. This attempts to extract from that heterogeneous datatype something
#' easily readable. Example:  Synonym() might return any of the following:
#' NA, NULL, "NA", "NULL", c("NA",NA,"GO:00001"), "GO:00002", c("Some text",NA,
#' NULL, "GO:00003") This function will boil that down to 'not found', '',
#' 'GO:00004', or "GO:0001,  some text, GO:00004"
#'
#' @param value Result of try(as.character(somefunction(GOTERM[id])), silent=TRUE).
#'   somefunction would be 'Synonym' 'Secondary' 'Ontology', etc...
#' @return something more sane (hopefully).
#' @seealso \pkg{GO.db}
#' @examples
#' \dontrun{
#'  ## goterms = GOTERM[ids]
#'  ## sane_goterms = deparse_go_value(goterms)
#' }
#' @export
deparse_go_value <- function(value) {
  result <- ""
  if (class(value) == "try-error") {
    result <- "Not found"
  } else {
    ## Not an error
    if (is.null(value)) {
      result <- ""
    } else if (is.na(value)) {
      result <- ""
    } else if (value == "NULL") {
      result <- ""
    } else if (grepl("^c\\(", as.character(value))) {
      value <- eval(parse(text=as.character(value)))
      if (class(value) == "logical") {
        result <- ""
      } else {
        value <- as.character(value[which(complete.cases(value))])
        result <- value
      }
    } else {
      ## Just a string "GO:00023409"
      result <- value
    }
  }
  return(result)
}

#' Get a go term from ID.
#'
#' @param go GO id or a list thereof, this may be a character or list(assuming
#'   the elements, not names, are goids).
#' @return Some text containing the terms associated with GO id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  goterm("GO:0032559")
#'  ## > GO:0032559
#'  ## > "adenyl ribonucleotide binding"
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
#' I think I will need to do similar parsing of the output for this function as
#' per gosec() In some cases this also returns stuff like c("some text",
#' "GO:someID") versus "some other text"  versus NULL versus NA.
#' This function just goes a mapply(gosn, go).
#'
#' @param go GO id, this may be a character or list(assuming the elements are goids).
#' @return Some text providing the synonyms for the given id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  text =  gosyn("GO:0000001")
#'  text
#'  ## > GO:000001
#'  ## > "mitochondrial inheritance"
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
#' Unfortunately, GOTERM's returns for secondary IDs are not consistent, so this
#' function has to have a whole bunch of logic to handle the various outputs.
#'
#' @param go GO ID, this may be a character or list(assuming the elements, not
#'   names, are goids).
#' @return Some text comprising the secondary GO id(s).
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  gosec("GO:0032432")
#'  ## > GO:0032432
#'  ## > "GO:0000141" "GO:0030482"
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
#' @param go GO ID, this may be a character or list (assuming the elements are
#'   goids).
#' @return Some text providing the long definition of each provided GO id.
#' @seealso \pkg{GOTermsAnnDbBimap}
#' @examples
#' \dontrun{
#'  godef("GO:0032432")
#'  ## > GO:0032432
#'  ## > "An assembly of actin filaments that are on the same axis but may be
#'  ## > same or opposite polarities and may be packed with different levels of tightness."
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
#'  ## > GO:0032432 GO:0032433
#'  ## > "CC" "CC"
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
#' Sometimes it is useful to know how far up/down the ontology tree a given id
#' resides.  This attmepts to answer that question.
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
  while (class(try(as.character(AnnotationDbi::Ontology(GO.db::GOTERM[go])),
                   silent=TRUE)) != "try-error") {
                     term <- GO.db::GOTERM[go]
                     test <- as.character(AnnotationDbi::Ontology(term))
                     if (class(test) == "try-error") {
                       next
                     }

                     if (test == "MF") {
                       ## I am not certain if GO.db:: will work for this
                       ancestors <- GO.db::GOMFANCESTOR[go]
                     } else if (test == "BP") {
                       ancestors <- GO.db::GOBPANCESTOR[go]
                     } else if (test == "CC") {
                       ancestors <- GO.db::GOCCANCESTOR[go]
                     } else {
                       ## There was an error
                       message("There was an error getting the ontology: ", as.character(go))
                       ancestors <- "error"
                     }
                     go <- as.list(ancestors)[[1]][1]
                     level <- level + 1
                     if (go == "all") {
                       return(level)
                     }
                   }  ## End for
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
#'  ## > 3 4
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
#'  ## > 1
#'  gotest("GO:0923429034823904")
#'  ## > 0
#' }
#' @export
gotest <- function(go) {
  gotst <- function(go) {
    go <- as.character(go)
    value <- try(GO.db::GOTERM[go], silent=TRUE)
    if (class(value) == "try-error") {
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

#' Use the orgdb instances from clusterProfiler to gather annotation data for GO.
#'
#' Since clusterprofiler no longer builds gomaps, I need to start understanding
#' how to properly get information from orgDBs.
#'
#' @param goseq_data  Some data from goseq and friends.
#' @param orgdb_go  The orgDb instance with GO data.
#' @param orgdb_ensembl  The orgDb instance with ensembl data.
#' @return a go mapping
#' @seealso \pkg{clusterProfiler}
gather_genes_orgdb <- function(goseq_data, orgdb_go, orgdb_ensembl) {
  ## all_ontologies <- mappedkeys(orgdb)
  ## all_mappings <- as.list(orgdb[all_ontologies])
  my_table <- goseq_data[["godata_interesting"]]
  my_ontologies <- my_table[["category"]]
  my_genes <- goseq_data[["input"]][["ensembl_gene"]]

  my_table[["entrez_ids"]] <- ""
  my_table[["ensembl_ids"]] <- ""
  my_table[["all_ensembl_in_ontology"]] <- ""
  for (count in 1:nrow(my_table)) {
    ont <- my_table[count, "category"]
    test_map <- try(as.list(orgdb_go[ont]), silent=TRUE)
    if (class(test_map) == "list") {
      ## If it returned as a list, then we should have 1:n entrez gene IDs
      ## Sadly, we do pretty much everything as ensembl, or perhaps geneDb
      entrez_geneids <- unlist(test_map)
      my_table[count, "entrez_ids"] <- toString(entrez_geneids)
      ensembl_geneids <- as.character(orgdb_ensembl[entrez_geneids])
      my_ensembl_geneids <- my_genes[my_genes %in% ensembl_geneids]
      my_table[count, "ensembl_ids"] <- toString(my_ensembl_geneids)
      my_table[count, "all_ensembl_in_ontology"] <- toString(ensembl_geneids)
    }
  }
  return(my_table)
}

#' Perform ontology searches given the results of a differential expression analysis.
#'
#' This takes a set of differential expression results, extracts a subset of
#' up/down expressed genes; passes them to goseq, clusterProfiler, topGO,
#' GOstats, and gProfiler; collects the outputs; and returns them in a
#' (hopefully) consistent fashion.  It attempts to handle the differing required
#' annotation/GOid inputs required for each tool and/or provide supported
#' species in ways which the various tools expect.
#'
#' @param de_out List of topTables comprising limma/deseq/edger outputs.
#' @param gene_lengths Data frame of gene lengths for goseq.
#' @param goids Data frame of goids and genes.
#' @param n Number of genes at the top/bottom of the fold-changes to define 'significant.'
#' @param z Number of standard deviations from the mean fold-change used to
#'   define 'significant.'
#' @param lfc Log fold-change used to define 'significant'.
#' @param p Maximum pvalue to define 'significant.'
#' @param overwrite Overwrite existing excel results file?
#' @param species Supported organism used by the tools.
#' @param goid_map Mapping file used by topGO, if it does not exist then
#'   goids_df creates it.
#' @param gff_file gff file containing the annotations used by gff2genetable
#'   from clusterprofiler.
#' @param gff_type Column to use from the gff file for the universe of genes.
#' @param do_goseq Perform simple_goseq()?
#' @param do_cluster Perform simple_clusterprofiler()?
#' @param do_topgo Perform simple_topgo()?
#' @param do_gostats Perform simple_gostats()?
#' @param do_gprofiler Perform simple_gprofiler()?
#' @param do_trees make topGO trees from the data?
#' @param orgdb Provide an organismDbi/Orgdb to hold the various annotation
#'   data, in response to the shift of clusterprofiler and friends towards using them.
#' @param ...  Arguments to pass through in arglist.
#' @return a list of up/down ontology results from
#'   goseq/clusterprofiler/topgo/gostats, and associated trees.
#' @seealso \pkg{goseq} \pkg{clusterProfiler} \pkg{topGO} \pkg{goStats}
#'   \pkg{gProfiler} \pkg{GO.db}
#' @examples
#' \dontrun{
#'  many_comparisons = limma_pairwise(expt=an_expt)
#'  tables = many_comparisons$limma
#'  this_takes_forever = limma_ontology(tables, gene_lengths=lengthdb,
#'                                      goids=goids_df, z=1.5, gff_file='length_db.gff')
#' }
#' @export
all_ontology_searches <- function(de_out, gene_lengths=NULL, goids=NULL, n=NULL,
                                  z=NULL, lfc=NULL, p=NULL, overwrite=FALSE,
                                  species="unsupported", orgdb="org.Dm.eg.db",
                                  goid_map="reference/go/id2go.map",
                                  gff_file=NULL, gff_type="gene", do_goseq=TRUE,
                                  do_cluster=TRUE, do_topgo=TRUE,
                                  do_gostats=TRUE, do_gprofiler=TRUE,
                                  do_trees=FALSE, ...) {
  arglist <- list(...)
  message("This function expects a list of contrast tables and annotation information.")
  message("The annotation information would be gene lengths and ontology ids")
  if (isTRUE(do_goseq) & is.null(gene_lengths)) {
    stop("Performing a goseq search requires a data frame of gene lengths.")
  }
  if (isTRUE(do_cluster) & is.null(gff_file)) {
    stop("Performing a clusterprofiler search requires a gff file.")
  }
  arglist <- list(...)

  goid_map <- get0("goid_map")
  if (is.null(goid_map)) {
    goid_map <- "reference/go/id2go.map"
  }

  ## Take a moment to detect what the data input looks like
  ## Perhaps a list of tables?
  if (is.null(de_out[["all_tables"]])) {
    if (sum(grepl(pattern="_vs_", x=names(de_out))) == 0) {
      stop("This assumes you are passing it a limma/deseq/edger output from
limma_pairwise(), edger_pairwise(), or deseq_pairwise().")
    }
  } else {
    ## In this case, you fed it something$limma rather than something$limma$all_tables
    de_out <- de_out[["all_tables"]]
  }
  ## Perhaps a single data frame of logFC etc
  ## In which case, coerse it to a list of 1
  if (!is.null(de_out[["logFC"]])) {
    tmp <- de_out
    de_out <- list("first"=tmp)
    rm(tmp)
  }

  output <- list()
  for (c in 1:length(de_out)) {
    datum <- de_out[[c]]
    if (!is.null(datum[["Row.names"]])) {
      rownames(datum) <- datum[["Row.names"]]
      datum <- datum[-1]
    }
    comparison <- names(de_out[c])
    message("Performing ontology search of:", comparison)
    updown_genes <- get_sig_genes(datum, n=n, z=z, lfc=lfc, p=p)
    up_genes <- updown_genes[["up_genes"]]
    down_genes <- updown_genes[["down_genes"]]

    goseq_up_ontology <- goseq_up_trees <- NULL
    goseq_down_ontology <- goseq_down_trees <- NULL
    cluster_up_ontology <- cluster_up_trees <- NULL
    cluster_down_ontology <- cluster_down_trees <- NULL
    topgo_up_ontology <- topgo_up_trees <- NULL
    topgo_down_ontology <- topgo_down_trees <- NULL
    gostats_up_ontology <- gostats_up_trees <- NULL
    gostats_down_ontology <- gostats_down_trees <- NULL
    gprofiler_up_ontology <- gprofiler_down_ontology <- NULL

    if (isTRUE(do_goseq)) {
      goseq_up_ontology <- try(simple_goseq(up_genes, goids, gene_lengths))
      goseq_down_ontology <- try(simple_goseq(down_genes, goids, gene_lengths))
      if (isTRUE(do_trees)) {
        goseq_up_trees <- try(
          goseq_trees(goseq_up_ontology, goid_map=goid_map))
        goseq_down_trees <- try(
          goseq_trees(goseq_down_ontology, goid_map=goid_map))
      }
    }

    if (isTRUE(do_cluster)) {
      cluster_up_ontology <- try(
        simple_clusterprofiler(up_genes, datum, orgdb=orgdb, ...))
      cluster_down_ontology <- try(
        simple_clusterprofiler(down_genes, datum, orgdb=orgdb, ...))
      if (isTRUE(do_trees)) {
        cluster_up_trees <- try(cluster_trees(up_genes, cluster_up_ontology,
                                              goid_map=goid_map, go_db=goids))
        cluster_down_trees <- try(cluster_trees(down_genes, cluster_down_ontology,
                                                goid_map=goid_map, go_db=goids))
      }
    }

    if (isTRUE(do_topgo)) {
      topgo_up_ontology <- try(
        simple_topgo(up_genes, goid_map=goid_map, go_db=goids))
      topgo_down_ontology <- try(
        simple_topgo(down_genes, goid_map=goid_map, go_db=goids))
      if (isTRUE(do_trees)) {
        topgo_up_trees <- try(topgo_trees(topgo_up_ontology))
        topgo_down_trees <- try(topgo_trees(topgo_down_ontology))
      }
    }

    if (isTRUE(do_gostats)) {
      gostats_up_ontology <- try(
        simple_gostats(up_genes, gff_file, goids, gff_type=gff_type))
      gostats_down_ontology <- try(
        simple_gostats(down_genes, gff_file, goids, gff_type=gff_type))
      if (isTRUE(do_trees)) {
        message("gostats_trees has never been tested, this is commented out for the moment.")
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
#' In the same way all_pairwise() attempts to simplify using multiple DE tools,
#' this function seeks to make it easier to extract subsets of differentially
#' expressed data and pass them to goseq, clusterProfiler, topGO, GOstats, and
#' gProfiler.
#'
#' @param changed_counts List of changed counts as ups and downs.
#' @param doplot Include plots in the results?
#' @param do_goseq Perform goseq search?
#' @param do_cluster Perform clusterprofiler search?
#' @param do_topgo Perform topgo search?
#' @param do_gostats Perform gostats search?
#' @param do_gprofiler Do a gprofiler search?
#' @param according_to If results from multiple DE tools were passed, which one
#'   defines 'significant'?
#' @param ...  Extra arguments!
#' @return List of ontology search results, up and down for each contrast.
#' @seealso \pkg{goseq} \pkg{clusterProfiler} \pkg{topGO} \pkg{goStats} \pkg{gProfiler}
#' @export
subset_ontology_search <- function(changed_counts, doplot=TRUE, do_goseq=TRUE,
                                   do_cluster=TRUE, do_topgo=TRUE, do_gostats=TRUE,
                                   do_gprofiler=TRUE, according_to="limma", ...) {
  arglist <- list(...)
  up_list <- NULL
  down_list <- NULL
  if (!is.null(changed_counts[["ups"]])) {
    message("The data is not sorted by a differential expression tool.")
    up_list <- changed_counts[["ups"]]
    down_list <- changed_counts[["downs"]]
  } else if (is.null(changed_counts[[according_to]])) {
    warning(glue("Could not find the table according to: {according_to}."))
    warning("Using the first available table, which might end badly.")
    up_list <- changed_counts[[1]][["ups"]]
    down_list <- changed_counts[[2]][["downs"]]
  } else {
    up_list <- changed_counts[[according_to]][["ups"]]
    down_list <- changed_counts[[according_to]][["downs"]]
  }
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
  goids <- NULL
  if (is.null(arglist[["goids"]])) {
    goids <- arglist[["go_db"]]
  } else {
    goids <- arglist[["goids"]]
  }
  gff <- arglist[["gff"]]
  gff_type <- arglist[["gff_type"]]
  types_list <- c("up_goseq", "down_goseq", "up_cluster", "down_cluster",
                  "up_topgo", "down_topgo", "up_gostats", "down_gostats",
                  "up_gprofiler", "down_gprofiler")
  names_list <- names(up_list)
  names_length <- length(names_list)
  for (cluster_count in 1:names_length) {
    name <- names_list[[cluster_count]]
    uppers <- up_list[[cluster_count]]
    downers <- down_list[[cluster_count]]
    if (isTRUE(do_goseq)) {
      message(cluster_count, "/", names_length, ": Starting goseq")
      up_goseq[[name]] <- try(simple_goseq(de_genes=uppers, ...))
      down_goseq[[name]] <- try(simple_goseq(de_genes=downers, ...))
    }
    if (isTRUE(do_cluster)) {
      message(cluster_count, "/", names_length, ": Starting clusterprofiler")
      up_cluster[[name]] <- try(simple_clusterprofiler(uppers, ...))
      down_cluster[[name]] <- try(simple_clusterprofiler(downers, ...))
    }
    if (isTRUE(do_topgo)) {
      message(cluster_count, "/", names_length, ": Starting topgo")
      up_topgo[[name]] <- try(simple_topgo(de_genes=uppers, ...))
      down_topgo[[name]] <- try(simple_topgo(de_genes=downers, ...))
    }
    if (isTRUE(do_gostats)) {
      message(cluster_count, "/", names_length, ": Starting gostats")
      up_gostats[[name]] <- try(simple_gostats(uppers, ...))
      down_gostats[[name]] <- try(simple_gostats(downers, ...))
    }
    if (isTRUE(do_gprofiler)) {
      message(cluster_count, "/", names_length, ": starting gprofiler.")
      up_gprofiler[[name]] <- try(simple_gprofiler(uppers, ...))
      down_gprofiler[[name]] <- try(simple_gprofiler(downers, ...))
    }
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
#' That is not to say it is very fast, so it saves the result to ontlevel.rda
#' for future lookups.
#'
#' @param ont   the ontology to recurse.
#' @param savefile   a file to save the results for future lookups.
#' @return golevels  a dataframe of goids<->highest level
#' @seealso \pkg{clusterProfiler}
#' @export
golevel_df <- function(ont="MF", savefile="ontlevel.rda") {
  savefile <- glue("{ont}_{savefile}")
  golevels <- NULL
  if (file.exists(savefile)) {
    load(savefile)
  } else {
    level <- 0
    continue <- 1
    golevels <- data.frame(GO=NULL, level=NULL)
    while (continue == 1) {
      level <- level + 1
      GO <- try(getGOLevel(ont, level), silent=TRUE)
      if (class(GO) != "character") {
        golevels[["level"]] <- as.numeric(golevels[["level"]])
        save(golevels, file=savefile, compress="xz")
        return(golevels)
      } else {
        tmpdf <- as.data.frame(cbind(GO, level))
        ## This (hopefully) ensures that each GO id is added only once,
        ## and gets the highest possible level.
        new_go <- tmpdf[unique(tmpdf[["GO"]], golevels[["GO"]]), ]
        golevels <- rbind(golevels, new_go)
      }
    }
  }
  return(golevels)
}

#' Compare the results from different ontology tools
#'
#' Combine the results from goseq, cluster profiler, topgo, and gostats; poke at
#' them with a stick and see what happens. The general idea is to pull the
#' p-value data from each tool and contrast that to the set of all possibile
#' ontologies.  This allows one to do a correlation coefficient between them.
#' In addition, take the 1-pvalue for each ontology for each tool. Thus for
#' strong p-values the score will be near 1 and so we can sum the scores for all
#' the tools.  Since topgo has 4 tools, the total possible is 7 if everything
#' has a p-value equal to 0.
#'
#' @param goseq   The goseq result from simple_goseq()
#' @param cluster   The result from simple_clusterprofiler()
#' @param topgo   Guess
#' @param gostats   Yep, ditto
#' @return a summary of the similarities of ontology searches
#' @seealso \pkg{goseq} \pkg{clusterProfiler} \pkg{topGO} \pkg{goStats}
#' @export
compare_go_searches <- function(goseq=NULL, cluster=NULL, topgo=NULL, gostats=NULL) {
  goseq_mf_data <- goseq_bp_data <- goseq_cc_data <- NULL
  if (!is.null(goseq)) {
    goseq_mf_data <- goseq[["mf_subset"]][, c("category", "over_represented_pvalue")]
    goseq_bp_data <- goseq[["bp_subset"]][, c("category", "over_represented_pvalue")]
    goseq_cc_data <- goseq[["cc_subset"]][, c("category", "over_represented_pvalue")]
    colnames(goseq_mf_data) <- c("goseq_id", "goseq_pvalue")
    colnames(goseq_bp_data) <- c("goseq_id", "goseq_pvalue")
    colnames(goseq_cc_data) <- c("goseq_id", "goseq_pvalue")
  }
  cluster_mf_data <- cluster_bp_data <- cluster_cc_data <- NULL
  if (!is.null(cluster)) {
    cluster_mf_data <- as.data.frame(summary(cluster[["mf_all"]]))[, c("ID", "pvalue")]
    cluster_bp_data <- as.data.frame(summary(cluster[["bp_all"]]))[, c("ID", "pvalue")]
    cluster_cc_data <- as.data.frame(summary(cluster[["cc_all"]]))[, c("ID", "pvalue")]
    colnames(cluster_mf_data) <- c("cluster_id", "cluster_pvalue")
    colnames(cluster_bp_data) <- c("cluster_id", "cluster_pvalue")
    colnames(cluster_cc_data) <- c("cluster_id", "cluster_pvalue")
  }
  topgo_mf_data <- topgo_bp_data <- topgo_cc_data <- NULL
  if (!is.null(topgo)) {
    topgo_mf_data <- topgo[["tables"]][["mf"]][, c("GO.ID", "fisher", "KS", "EL", "weight")]
    topgo_bp_data <- topgo[["tables"]][["bp"]][, c("GO.ID", "fisher", "KS", "EL", "weight")]
    topgo_cc_data <- topgo[["tables"]][["cc"]][, c("GO.ID", "fisher", "KS", "EL", "weight")]
    colnames(topgo_mf_data) <- c("topgo_id", "topgo_fisher_pvalue",
                                 "topgo_el_pvalue", "topgo_ks_pvalue", "topgo_weight_pvalue")
    colnames(topgo_bp_data) <- c("topgo_id", "topgo_fisher_pvalue",
                                 "topgo_el_pvalue", "topgo_ks_pvalue", "topgo_weight_pvalue")
    colnames(topgo_cc_data) <- c("topgo_id", "topgo_fisher_pvalue",
                                 "topgo_el_pvalue", "topgo_ks_pvalue", "topgo_weight_pvalue")
  }
  gostats_mf_data <- gostats_bp_data <- gostats_cc_data <- NULL
  if (!is.null(gostats)) {
    gostats_mf_data <- gostats[["mf_over_all"]][, c("GOMFID", "Pvalue")]
    gostats_bp_data <- gostats[["bp_over_all"]][, c("GOBPID", "Pvalue")]
    gostats_cc_data <- gostats[["cc_over_all"]][, c("GOCCID", "Pvalue")]
    colnames(gostats_mf_data) <- c("gostats_id", "gostats_pvalue")
    colnames(gostats_bp_data) <- c("gostats_id", "gostats_pvalue")
    colnames(gostats_cc_data) <- c("gostats_id", "gostats_pvalue")
  }
  ## Now combine them...
  all_mf <- merge(goseq_mf_data, cluster_mf_data,
                  by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
  all_mf <- merge(all_mf, topgo_mf_data,
                  by.x="goseq_id", by.y="topgo_id", all.x=TRUE, all.y=TRUE)
  all_mf <- merge(all_mf, gostats_mf_data,
                  by.x="goseq_id", by.y="gostats_id", all.x=TRUE, all.y=TRUE)
  rownames(all_mf) <- all_mf[["goseq_id"]]
  all_mf <- all_mf[-1]
  all_mf[is.na(all_mf)] <- 1
  all_mf[["goseq_pvalue"]] <- 1 - as.numeric(all_mf[["goseq_pvalue"]])
  all_mf[["cluster_pvalue"]] <- 1 - as.numeric(all_mf[["cluster_pvalue"]])
  all_mf[["topgo_fisher_pvalue"]] <- 1 - as.numeric(all_mf[["topgo_fisher_pvalue"]])
  all_mf[["topgo_el_pvalue"]] <- 1 - as.numeric(all_mf[["topgo_el_pvalue"]])
  all_mf[["topgo_ks_pvalue"]] <- 1 - as.numeric(all_mf[["topgo_ks_pvalue"]])
  all_mf[["topgo_weight_pvalue"]] <- 1 - as.numeric(all_mf[["topgo_weight_pvalue"]])
  all_mf[["gostats_pvalue"]] <- 1 - as.numeric(all_mf[["gostats_pvalue"]])
  all_mf[is.na(all_mf)] <- 0
  message(cor(all_mf))
  mf_summary <- rowSums(all_mf)
  message(summary(mf_summary))

  all_bp <- merge(goseq_bp_data, cluster_bp_data,
                  by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
  all_bp <- merge(all_bp, topgo_bp_data, by.x="goseq_id",
                  by.y="topgo_id", all.x=TRUE, all.y=TRUE)
  all_bp <- merge(all_bp, gostats_bp_data, by.x="goseq_id",
                  by.y="gostats_id", all.x=TRUE, all.y=TRUE)
  rownames(all_bp) <- all_bp[["goseq_id"]]
  all_bp <- all_bp[-1]
  all_bp[is.na(all_bp)] <- 1
  all_bp[["goseq_pvalue"]] <- 1 - as.numeric(all_bp[["goseq_pvalue"]])
  all_bp[["cluster_pvalue"]] <- 1 - as.numeric(all_bp[["cluster_pvalue"]])
  all_bp[["topgo_fisher_pvalue"]] <- 1 - as.numeric(all_bp[["topgo_fisher_pvalue"]])
  all_bp[["topgo_el_pvalue"]] <- 1 - as.numeric(all_bp[["topgo_el_pvalue"]])
  all_bp[["topgo_ks_pvalue"]] <- 1 - as.numeric(all_bp[["topgo_ks_pvalue"]])
  all_bp[["topgo_weight_pvalue"]] <- 1 - as.numeric(all_bp[["topgo_weight_pvalue"]])
  all_bp[["gostats_pvalue"]] <- 1 - as.numeric(all_bp[["gostats_pvalue"]])
  all_bp[is.na(all_bp)] <- 0
  message(cor(all_bp))
  bp_summary <- rowSums(all_bp)
  message(summary(bp_summary))

  all_cc <- merge(goseq_cc_data, cluster_cc_data,
                  by.x="goseq_id", by.y="cluster_id", all.x=TRUE, all.y=TRUE)
  all_cc <- merge(all_cc, topgo_cc_data, by.x="goseq_id",
                  by.y="topgo_id", all.x=TRUE, all.y=TRUE)
  all_cc <- merge(all_cc, gostats_cc_data, by.x="goseq_id",
                  by.y="gostats_id", all.x=TRUE, all.y=TRUE)
  rownames(all_cc) <- all_cc[["goseq_id"]]
  all_cc <- all_cc[-1]
  all_cc[is.na(all_cc)] <- 1
  all_cc[["goseq_pvalue"]] <- 1 - as.numeric(all_cc[["goseq_pvalue"]])
  all_cc[["cluster_pvalue"]] <- 1 - as.numeric(all_cc[["cluster_pvalue"]])
  all_cc[["topgo_fisher_pvalue"]] <- 1 - as.numeric(all_cc[["topgo_fisher_pvalue"]])
  all_cc[["topgo_el_pvalue"]] <- 1 - as.numeric(all_cc[["topgo_el_pvalue"]])
  all_cc[["topgo_ks_pvalue"]] <- 1 - as.numeric(all_cc[["topgo_ks_pvalue"]])
  all_cc[["topgo_weight_pvalue"]] <- 1 - as.numeric(all_cc[["topgo_weight_pvalue"]])
  all_cc[["gostats_pvalue"]] <- 1 - as.numeric(all_cc[["gostats_pvalue"]])
  all_cc[is.na(all_cc)] <- 0
  message(cor(all_cc))
  cc_summary <- rowSums(all_cc)
  message(summary(cc_summary))
}

## EOF
