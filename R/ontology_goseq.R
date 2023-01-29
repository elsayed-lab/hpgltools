## ontology_goseq.r: Some methods to simplify the usage of goseq.  goseq is the
## first GSEA tool I learned about.  It is not particularly robust, this seeks
## to amend that.

#' Filter a goseq significance search
#'
#' Given a goseq result, use some simple filters to pull out the
#' categories of likely interest.
#'
#' @param godata goseq result
#' @param expand_categories Extract GO terms from GO.db and add them
#'  to the table
#' @param pvalue Significance filter.
#' @param minimum_interesting The category should have more than this
#'  number of elements.
#' @param adjust Adjusted p-value filter.
#' @param padjust_method Method for adjusting the p-values.
extract_interesting_goseq <- function(godata, expand_categories = TRUE, pvalue = 0.05,
                                      minimum_interesting = 1, adjust = 0.05, padjust_method = "BH") {
  ## Add a little logic if we are using a non-GO input database.
  ## This is relevant if you use something like msigdb/reactome/kegg
  ## I should probably make the rest of this function generic so that this is not required.
  if (is.null(godata[["ontology"]])) {
    godata[["ontology"]] <- "MF"
  }

  godata_interesting <- godata
  if (isTRUE(expand_categories)) {
    mesg("simple_goseq(): Filling godata with terms, this is slow.")
    godata_interesting <- goseq_table(godata)
  } else {
    ## Set the 'term' category for plotting.
    godata_interesting[["term"]] <- godata_interesting[["category"]]
    godata[["term"]] <- godata[["category"]]
  }

  if (is.null(adjust)) {
    interesting_idx <- godata[["over_represented_pvalue"]] <= pvalue
    godata_interesting <- godata[godata_interesting, ]
    padjust_method <- "none"
  } else {
    ## There is a requested pvalue adjustment
    interesting_idx <- godata[["over_represented_pvalue"]] <= adjust
    godata_interesting <- godata[interesting_idx, ]
    if (dim(godata_interesting)[1] < minimum_interesting) {
      message("simple_goseq(): There are no genes with an adj.p < ", adjust, " using: ",
              padjust_method, ".")
      message("simple_goseq(): Providing genes with raw pvalue < ", pvalue, ".")
      interesting_idx <- godata[["over_represented_pvalue"]] <= pvalue
      godata_interesting <- godata_interesting[interesting_idx, ]
      padjust_method <- "none"
    }
  }

  na_idx <- is.na(godata[["ontology"]])
  godata <- godata[!na_idx, ]
  mf_idx <- godata[["ontology"]] == "MF"
  mf_subset <- godata[mf_idx, ]
  rownames(mf_subset) <- mf_subset[["category"]]
  bp_idx <- godata[["ontology"]] == "BP"
  bp_subset <- godata[bp_idx, ]
  rownames(bp_subset) <- bp_subset[["category"]]
  cc_idx <- godata[["ontology"]] == "CC"
  cc_subset <- godata[cc_idx, ]
  rownames(cc_subset) <- cc_subset[["category"]]

  na_idx <- is.na(godata_interesting[["ontology"]])
  godata_interesting <- godata_interesting[!na_idx, ]
  mf_idx <- godata_interesting[["ontology"]] == "MF"
  mf_interesting <- godata_interesting[mf_idx, ]
  rownames(mf_interesting) <- mf_interesting[["category"]]
  if (is.null(mf_interesting[["term"]])) {
    mf_interesting <- mf_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue")]
  } else {
    mf_interesting <- mf_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]
  }
  bp_idx <- godata_interesting[["ontology"]] == "BP"
  bp_interesting <- godata_interesting[bp_idx, ]
  rownames(bp_interesting) <- bp_interesting[["category"]]
  if (is.null(bp_interesting[["term"]])) {
    bp_interesting <- bp_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue")]
  } else {
    bp_interesting <- bp_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]
  }
  cc_idx <- godata_interesting[["ontology"]] == "CC"
  cc_interesting <- godata_interesting[cc_idx, ]
  rownames(cc_interesting) <- cc_interesting[["category"]]
  if (is.null(cc_interesting[["term"]])) {
    cc_interesting <- cc_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue")]
  } else {
    cc_interesting <- cc_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]
  }
  retlist <- list(
      "godata" = godata,
      "interesting" = godata_interesting,
      "mf_subset" = mf_subset,
      "bp_subset" = bp_subset,
      "cc_subset" = cc_subset,
      "MF" = mf_interesting,
      "BP" = bp_interesting,
      "CC" = cc_interesting)
  return(retlist)
}

#' Pass MSigDB categorical data to goseq and run it.
#'
#' goseq is probably the easiest method to push varying data types into.  Thus
#' it was the first thing I thought of when looking to push MSigDB data into a
#' GSEA method.
#'
#' @param sig_genes Character list of genes deemed significant.  I think in the
#'  current implementation this must be just a list of IDs as opposed to the
#'  full dataframe of interesting genes because we likely need to convert IDs.
#' @param signatures Used by load_gmt_signatures(), the signature file or set.
#' @param data_pkg Used by load_gmt_signatures().
#' @param signature_category Ibid, but the name of the signatures group.
#' @param current_id Used by convert_msig_ids(), when converting IDs, the name
#'  of the existing type.
#' @param required_id What type to convert to in convert_msig_ids().
#' @param length_db Dataframe of lengths.  It is worth noting that goseq
#'  explicitly states that one might wish to use other potentially confounding
#'  factors here, but they only examine lenghts in their paper. Starting with
#'  this parameter, everything is just passed directly to simple_goseq()
#' @param doplot Print the prior plot?
#' @param adjust passed to simple_goseq()
#' @param pvalue passed to simple_goseq()
#' @param length_keytype passed to simple_goseq()
#' @param go_keytype passed to simple_goseq()
#' @param goseq_method passed to simple_goseq()
#' @param padjust_method passed to simple_goseq()
#' @param excel passed to simple_goseq()
#' @param orgdb Ideally used to help goseq collect lengths.
#' @return Some goseq data!
#' @seealso [gsva] [goseq]
#' @export
goseq_msigdb <- function(sig_genes, signatures = "c2BroadSets", data_pkg = "GSVAdata",
                         signature_category = "c2", current_id = "ENSEMBL", required_id = "ENTREZID",
                         length_db = NULL, doplot = TRUE, adjust = 0.1, pvalue = 0.1,
                         length_keytype = "transcripts", go_keytype = "entrezid",
                         goseq_method = "Wallenius", padjust_method = "BH",
                         excel = NULL, orgdb = "org.Hs.eg.db") {
  sig_data <- load_gmt_signatures(signatures = signatures, data_pkg = data_pkg,
                                  signature_category = signature_category)
  mesg("Starting to coerce the msig data to the ontology format.")
  go_db <- data.table::data.table()
  for (i in seq_along(sig_data)) {
    gsc <- sig_data[[i]]
    gsc_id <- gsc@setName
    gsc_genes <- gsc@geneIds
    tmp_db <- data.table::data.table("ID" = gsc_genes, "GO" = rep(gsc_id, length(gsc_genes)))
    go_db <- rbind(go_db, tmp_db)
  }
  mesg("Finished coercing the msig data into a df with ", nrow(go_db), " rows.")

  new_ids <- NULL
  new_sig <- data.frame()
  if ("character" %in% class(sig_genes)) {
    new_ids <- convert_ids(sig_genes, from = current_id, to = required_id, orgdb = orgdb)
    new_sig <- new_ids
    colnames(new_sig) <- c(current_id, "ID")
    new_sig <- new_sig[, c("ID", current_id)]
    rownames(new_sig) <- make.names(new_sig[["ID"]], unique = TRUE)
  } else if ("data.frame" %in% class(sig_genes)) {
    new_ids <- convert_ids(rownames(sig_genes), from = current_id,
                           to = required_id, orgdb = orgdb)
    new_sig <- merge(new_ids, sig_genes, by.x = current_id,
                     by.y = "row.names")
    new_sig[["ID"]] <- new_sig[[required_id]]
  } else {
    stop("I do not understand this input data format for sig_genes.")
  }

  new_lids <- convert_ids(rownames(length_db), from = current_id, to = required_id, orgdb = orgdb)
  new_length <- merge(new_lids, length_db, by.x = current_id, by.y = "row.names")
  new_length[["ID"]] <- NULL
  new_length[["ID"]] <- new_length[[required_id]]
  if (is.null(new_length[["length"]])) {
    new_length <- new_length[, c("ID", "width")]
    colnames(new_length) <- c("ID", "length")
  } else {
    new_length <- new_length[, c("ID", "length")]
  }

  go_result <- simple_goseq(new_sig, go_db = go_db, length_db = new_length,
                            doplot = TRUE, adjust = 0.1, pvalue = 0.1,
                            length_keytype = "transcripts", go_keytype = "entrezid",
                            goseq_method = "Wallenius", padjust_method = "BH",
                            plot_title = "Enriched MSIG categories",
                            expand_categories = FALSE,
                            excel = excel, add_trees = FALSE, gather_genes = FALSE, width = 20)
  return(go_result)
}

#' Enhance the goseq table of gene ontology information.
#'
#' While goseq has some nice functionality, the table of outputs it provides is
#' somewhat lacking. This attempts to increase that with some extra helpful data
#' like ontology categories, definitions, etc.
#'
#' @param df Dataframe of ontology information.  This is intended to be the
#'  output from goseq including information like numbers/category, GOids, etc.
#'  It requires a column 'category'
#'  which contains: GO:000001 and such.
#' @param file Csv file to which to write the table.
#' @return Ontology table with annotation information included.
#' @seealso [goseq] [GO.db]
#' @examples
#' \dontrun{
#'  annotated_go = goseq_table(go_ids)
#'  head(annotated_go, n = 1)
#'  ## >        category numDEInCat numInCat over_represented_pvalue
#'  ## > 571  GO:0006364          9       26            4.655108e-08
#'  ## >      under_represented_pvalue       qvalue ontology
#'  ## > 571                 1.0000000 6.731286e-05       BP
#'  ## >                                term
#'  ## > 571                 rRNA processing
#'  ## >                               synonym
#'  ## > 571        "35S primary transcript processing, GO:0006365"
#'  ## >        secondary    definition
#'  ## > 571    GO:0006365   Any process involved in the conversion of a primary ribosomal
#'  ##          RNA (rRNA) transcript into one or more mature rRNA molecules.
#' }
#' @export
goseq_table <- function(df, file = NULL) {
  df[["term"]] <- goterm(df[["category"]])
  df[["ontology"]] <- goont(df[["category"]])
  message("Testing that go categories are defined.")
  df[["good"]] <- gotest(df[["category"]])
  message("Removing undefined categories.")
  good_idx <- df[["good"]] == 1
  df <- df[good_idx, ]
  na_idx <- is.na(df[["category"]])
  df <- df[!na_idx, ]
  na_idx <- is.na(df[["term"]])
  df <- df[!na_idx, ]
  message("Gathering synonyms.")
  df[["synonym"]] <- gosyn(df[["category"]])
  message("Gathering category definitions.")
  df[["definition"]] <- godef(df[["category"]])
  df <- df[, c("category", "numDEInCat", "numInCat", "over_represented_pvalue",
               "under_represented_pvalue", "qvalue", "ontology", "term",
               "synonym", "definition")]
  if (!is.null(file)) {
    write.csv(df, file = file)
  }
  return(df)
}

#' Perform a simplified goseq analysis.
#'
#' goseq can be pretty difficult to get set up for non-supported organisms.
#' This attempts to make that process a bit simpler as well as give some
#' standard outputs which should be similar to those returned by
#' clusterprofiler/topgo/gostats/gprofiler.
#'
#' @param sig_genes Data frame of differentially expressed genes, containing IDs etc.
#' @param go_db Database of go to gene mappings (OrgDb/OrganismDb)
#' @param length_db Database of gene lengths (gff/TxDb)
#' @param doplot Include pwf plots?
#' @param adjust Minimum adjusted pvalue for 'significant.'
#' @param pvalue Minimum pvalue for 'significant.'
#' @param plot_title Set a title for the pvalue plots.
#' @param length_keytype Keytype to provide to extract lengths
#' @param go_keytype Keytype to provide to extract go IDs
#' @param goseq_method Statistical test for goseq to use.
#' @param padjust_method Which method to use to adjust the pvalues.
#' @param expand_categories Expand the GO categories to make the results more readable?
#' @param excel Print the results to an excel file?
#' @param ... Extra parameters which I do not recall
#' @return Big list including:
#'  the pwd:pwf function,
#'  alldata:the godata dataframe,
#'  pvalue_histogram:p-value histograms,
#'  godata_interesting:the ontology information of the enhanced groups,
#'  term_table:the goterms with some information about them,
#'  mf_subset:a plot of the MF enhanced groups,
#'  mfp_plot:the pvalues of the MF group,
#'  bp_subset:a plot of the BP enhanced groups,
#'  bpp_plot,
#'  cc_subset,
#'  and ccp_plot
#' @seealso [goseq] [GO.db] [GenomicFeatures] [stats::p.adjust()]
#' @examples
#' \dontrun{
#'  lotsotables <- simple_goseq(gene_list, godb, lengthdb)
#' }
#' @export
simple_goseq <- function(sig_genes, go_db = NULL, length_db = NULL, doplot = TRUE,
                         adjust = 0.1, pvalue = 0.1, plot_title = NULL,
                         length_keytype = "transcripts", go_keytype = "entrezid",
                         goseq_method = "Wallenius", padjust_method = "BH",
                         expand_categories = TRUE, excel = NULL,
                         ...) {
  arglist <- list(...)

  minimum_interesting <- 1
  if (!is.null(arglist[["minimum_interesting"]])) {
    minimum_interesting <- arglist[["minimum_interesting"]]
  }
  length_df <- data.frame()
  length_vector <- vector()
  de_vector <- vector()
  gene_ids <- NULL
  final_keytype <- NULL
  goids_df <- NULL
  ## sig_genes may be a list, character list, or data frame.
  gene_list <- NULL
  if (class(sig_genes)[1] == "data.table") {
    sig_genes <- as.data.frame(sig_genes)
  }
  if (class(sig_genes) == "character") {
    ## Then this is a character list of gene ids
    gene_list <- sig_genes
  } else if (class(sig_genes) == "list") {
    gene_list <- names(sig_genes)
  } else if (class(sig_genes) == "data.frame") {
    if (is.null(rownames(sig_genes)) & is.null(sig_genes[["ID"]])) {
      stop("This requires a set of gene IDs from the rownames or a column named 'ID'.")
    } else if (!is.null(sig_genes[["ID"]])) {
      ## Use a column named 'ID' first because a bunch of annotation databases
      ## use ENTREZ IDs which are just integers, which of course is not allowed
      ## by data frame row names.
      mesg("Using the ID column from your table rather than the row names.")
      gene_list <- sig_genes[["ID"]]
    } else if (!is.null(rownames(sig_genes))) {
      mesg("Using the row names of your table.")
      gene_list <- rownames(sig_genes)
    } else {
      gene_list <- sig_genes[["ID"]]
    }
  } else {
    stop("Not sure how to handle your set of significant gene ids.")
  }
  ## At this point I should have a character list of gene ids named 'gene_list'
  de_genelist <- as.data.frame(gene_list)
  de_genelist[["DE"]] <- 1
  colnames(de_genelist) <- c("ID", "DE")

  ## Database of lengths may be a gff file, TxDb, or OrganismDb
  metadf <- NULL
  if (class(length_db)[1] == "data.table") {
    length_db <- as.data.frame(length_db)
  }
  if (class(length_db)[[1]] == "character")  {
    ## Then this should be either a gff file or species name.
    if (grepl(pattern = "\\.gff", x = length_db, perl = TRUE) |
        grepl(pattern = "\\.gtf", x = length_db, perl = TRUE)) {
      ## gff file
      txdb <- GenomicFeatures::makeTxDbFromGFF(length_db)
      metadf <- extract_lengths(db = txdb, gene_list = gene_list)
    } else {
      ## Then species name
      message("A species name.")
    }
  } else if (class(length_db)[[1]] == "AnnotationDbi") {
    stop("This currently requires an actual OrganismDb, not AnnotationDbi.")
  } else if (class(length_db)[[1]] == "OrgDb") {
    stop("OrgDb objects contain links to other databases, but no gene lengths.")
  } else if (class(length_db)[[1]] == "OrganismDb" |
             class(length_db)[[1]] == "AnnotationDbi") {
    ##metadf <- extract_lengths(db = length_db, gene_list = gene_list)
    metadf <- sm(extract_lengths(db = length_db, gene_list = gene_list, ...))
  } else if (class(length_db)[[1]] == "TxDb") {
    metadf <- sm(extract_lengths(db = length_db, gene_list = gene_list, ...))
  } else if (class(length_db)[[1]] == "data.frame") {
    metadf <- length_db
  } else {
    stop("Extracting lengths requires the name of a supported species or orgdb instance.")
  }

  ## Sometimes the column with gene lengths is named 'width'
  ## In that case, fix it.
  if (is.null(metadf[["width"]]) & is.null(metadf[["length"]])) {
    stop("The length db needs to have a length or width column.")
  } else if (is.null(metadf[["length"]])) {
    ## Then it is named 'width' and I want to rename it to length
    colnames(metadf) <- gsub(x = colnames(metadf), pattern = "width", replacement = "length")
  }
  ## Now I should have the gene list and gene lengths

  godf <- data.frame()
  if (class(go_db)[[1]] == "character") {
    ## A text table or species name
    if (grepl(pattern = "\\.csv", x = go_db, perl = TRUE) |
        grepl(pattern = "\\.tab", x = go_db, perl = TRUE)) {
      ## table
      godf <- read.table(go_db, ...)
      colnames(godf) <- c("ID", "GO")
    } else {
      ## Assume species name
      supported <- TRUE
      species <- go_db
    }
  } else if (class(go_db)[[1]] == "OrganismDb") {
    godf <- extract_go(go_db, keytype = go_keytype)
  } else if (class(go_db)[[1]] == "OrgDb") {
    godf <- extract_go(go_db, keytype = go_keytype)
  } else if (class(go_db)[[1]] == "data.table" || class(go_db)[[1]] == "tbl_df") {
    godf <- as.data.frame(go_db)
  } else if (class(go_db)[[1]] == "data.frame") {
    godf <- go_db
    if (!is.null(godf[["ID"]])) {
      godf <- godf[, c("ID", "GO")]
    } else if (!is.null(godf[["GID"]])) {
      godf <- godf[, c("GID", "GO")]
      colnames(godf) <- c("ID", "GO")
    } else {
      stop("Unable to read the gene ID/ GO columns from the go dataframe.")
    }
  } else {
    stop("Unable to determine the input for creating a go dataframe.")
  }

  ## entrez IDs are numeric.  This is a problem when doing the pwf function
  ## because it sets the rownames to the IDs.  As a result, we need to call
  ## make.names() on them.
  colnames(godf) <- c("ID","GO")
  godf[["ID"]] <- make.names(godf[["ID"]])
  metadf[["ID"]] <- make.names(metadf[["ID"]])
  de_genelist[["ID"]] <- make.names(de_genelist[["ID"]])
  ## Ok, now I have a df of GOids, all gene lengths, and DE gene list. That is
  ## everything I am supposed to need for goseq.

  ## See how many entries from the godb are in the list of genes.
  id_xref <- de_genelist[["ID"]] %in% godf[["ID"]]
  meta_xref <- de_genelist[["ID"]] %in% metadf[["ID"]]
  message("Found ", sum(id_xref), " go_db genes and ", sum(meta_xref),
          " length_db genes out of ", nrow(de_genelist), ".")
  ## So lets merge the de genes and gene lengths to ensure that they are
  ## consistent. Then make the vectors expected by goseq.
  merged_ids_lengths <- metadf
  ## The following line was done in order to avoid
  ## "'unimplemented type 'list' in 'orderVector1'"
  merged_ids_lengths[["ID"]] <- as.character(merged_ids_lengths[["ID"]])
  merged_ids_lengths <- merge(merged_ids_lengths,
                              de_genelist, by.x = "ID", by.y = "ID", all.x = TRUE)
  merged_ids_lengths[["length"]] <- suppressWarnings(
      as.numeric(merged_ids_lengths[["length"]]))
  merged_ids_lengths[is.na(merged_ids_lengths)] <- 0
  ## Not casing the next lines as character/numeric causes weird errors like 'names' attribute
  ## must be the same length as the vector
  de_vector <- as.vector(as.numeric(merged_ids_lengths[["DE"]]))
  names(de_vector) <- make.names(as.character(
      merged_ids_lengths[["ID"]]), unique = TRUE)
  length_vector <- as.vector(as.numeric(merged_ids_lengths[["length"]]))
  names(length_vector) <- make.names(as.character(
      merged_ids_lengths[["ID"]]), unique = TRUE)

  pwf_plot <- NULL
  tmp_file <- tempfile(pattern = "goseq", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  pwf <- sm(suppressWarnings(goseq::nullp(DEgenes = de_vector, bias.data = length_vector,
                                          plot.fit = doplot)))
  if (isTRUE(doplot)) {
    pwf_plot <- recordPlot()
  }
  dev.off()
  removed <- file.remove(tmp_file)

  godata <- sm(goseq::goseq(pwf, gene2cat = godf, use_genes_without_cat = TRUE,
                            method = goseq_method))
  ## I want to limit the y-axis, but I think this is not the best way.
  goseq_p <- try(plot_histogram(godata[["over_represented_pvalue"]], bins = 50))
  goseq_p_nearzero <- table(goseq_p[["data"]])[[1]]
  goseq_y_limit <- goseq_p_nearzero * 2
  goseq_p <- goseq_p +
    ggplot2::scale_y_continuous(limits = c(0, goseq_y_limit))
  godata[["qvalue"]] <- stats::p.adjust(godata[["over_represented_pvalue"]],
                                        method = padjust_method)

  ## Subset the result for 'interesting' categories, which is defined
  ## simply as the set of categories with full annotations.
  interesting <- extract_interesting_goseq(godata, expand_categories = expand_categories,
                                           pvalue = pvalue, adjust = adjust,
                                           minimum_interesting = minimum_interesting,
                                           padjust_method = padjust_method)

  mesg("simple_goseq(): Making pvalue plots for the ontologies.")
  pvalue_plots <- plot_goseq_pval(godata, plot_title = plot_title,
                                  x_column = "over_represented_pvalue",
                                  ...)
  pval_plots <- list(
      "bpp_plot_over" = pvalue_plots[["bpp_plot_over"]],
      "mfp_plot_over" = pvalue_plots[["mfp_plot_over"]],
      "ccp_plot_over" = pvalue_plots[["ccp_plot_over"]])

  retlist <- list(
      "input" = sig_genes,
      "pwf" = pwf,
      "pwf_plot" = pwf_plot,
      "all_data" = interesting[["godata"]],
      "go_db" = godf,
      "godata" = godata,
      "pvalue_histogram" = goseq_p,
      "godata_interesting" = interesting[["interesting"]],
      "mf_interesting" = interesting[["MF"]],
      "bp_interesting" = interesting[["BP"]],
      "cc_interesting" = interesting[["CC"]],
      "goadjust_method" = goseq_method,
      "adjust_method" = padjust_method,
      "mf_subset" = interesting[["mf_subset"]],
      "bp_subset" = interesting[["bp_subset"]],
      "cc_subset" = interesting[["cc_subset"]],
      "pvalue_plots" = pval_plots)
  class(retlist) <- c("goseq_result", "list")
  retlist[["mf_enrich"]] <- goseq2enrich(retlist, ontology = "MF",
                                         cutoff = pvalue, padjust_method = padjust_method)
  retlist[["bp_enrich"]] <- goseq2enrich(retlist, ontology = "BP",
                                         cutoff = pvalue, padjust_method = padjust_method)
  retlist[["cc_enrich"]] <- goseq2enrich(retlist, ontology = "CC",
                                         cutoff = pvalue, padjust_method = padjust_method)
  if (!is.null(excel)) {
    excel_result <- write_goseq_data(retlist, excel = excel, ...)
  }
  return(retlist)
}

## EOF
