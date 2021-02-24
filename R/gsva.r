#' Use AnnotationDbi to translate geneIDs from type x to type y.
#'
#' This is intended to convert all the IDs in a geneSet from one ID type to
#' another and giving back the geneSet with the new IDs.
#' FIXME: This should use convert_ids() to simplify itself
#'
#' One caveat: this will collapse redundant IDs via unique().
#'
#' @param gsc geneSetCollection with IDs of a type one wishes to change.
#' @param orgdb  Annotation object containing the various IDs.
#' @param from_type  Name of the ID which your gsc is using.  This can probably
#'   be automagically detected...
#' @param to_type  Name of the ID you wish to use.
#' @return Fresh gene set collection replete with new names.
#' @export
convert_gsc_ids <- function(gsc, orgdb = "org.Hs.eg.db", from_type = NULL, to_type = "ENTREZID") {
  message("Converting the rownames() of the expressionset to ", to_type, ".")
  ##tt <- sm(library(orgdb, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  orgdb <- get0(orgdb)
  gsc_lst <- as.list(gsc)
  new_gsc <- list()
  show_progress <- interactive() && is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (g in 1:length(gsc)) {
    if (isTRUE(show_progress)) {
      pct_done <- g / length(gsc_lst)
      setTxtProgressBar(bar, pct_done)
    }
    gs <- gsc[[g]]
    old_ids <- GSEABase::geneIds(gs)
    if (is.null(from_type)) {
      from_type <- guess_orgdb_keytype(old_ids, orgdb)
    }
    new_ids <- sm(AnnotationDbi::select(x = orgdb,
                                        keys = old_ids,
                                        keytype = from_type,
                                        columns = c(to_type)))
    new_ids <- new_ids[[to_type]]
    GSEABase::geneIds(gs) <- unique(new_ids)
    ## gsc_lst[[g]] <- gs
    new_gsc[[g]] <- gs
  }
  if (isTRUE(show_progress)) {
    close(bar)
  }
  gsc <- GSEABase::GeneSetCollection(new_gsc)
  return(gsc)
}

#' Change gene IDs to the format expected by gsva using an orgdb.
#'
#' Though it is possible to use gsva without ENTREZ IDs, it is not trivial.
#' This function attempts to ensure that the IDs in one's expressionset are
#' therefore entrez IDs. It is possible that this function is at least partially
#' redundant with other functions in this package and should be replaced.
#'
#' @param ids Vector of IDS to modify.
#' @param from Change from this format.
#' @param to Change to this format.
#' @param orgdb Using this orgdb instance.
#' @return New vector of ENTREZ IDs.
convert_ids <- function(ids, from = "ENSEMBL", to = "ENTREZID", orgdb = "org.Hs.eg.db") {
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  new_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                      keys = ids,
                                      keytype = current_id,
                                      columns = c(required_id)))
  new_idx <- complete.cases(new_ids)
  new_ids <- new_ids[new_idx, ]
  message("Before conversion, the expressionset has ", length(ids),
          " entries.")
  message("After conversion, the expressionset has ",
          length(rownames(new_ids)),
          " entries.")
  return(new_ids)
}

#' Extract the GeneSets corresponding to the provided name(s).
#'
#' Many of the likely GSCs contain far more gene sets than one actually wants to
#' deal with.  This will subset them according to a the desired 'requests'.
#'
#' @param sig_data  The pile of GeneSets, probably from GSVAdata.
#' @param requests  Character list of sources to keep.
#' @return  Whatever GeneSets remain.
#' @export
get_gsvadb_names <- function(sig_data, requests = NULL) {
  requests <- toupper(requests)
  categories <- names(sig_data)
  prefixes <- gsub(pattern = "^(.*?)_.*", replacement = "\\1",
                   x = categories)
  tab <- table(prefixes)
  if (is.null(requests)) {
    new_order <- order(as.numeric(tab), decreasing = TRUE)
    tab <- tab[new_order]
    return(tab)
  }

  kept_requests <- requests %in% names(tab)
  requests <- requests[kept_requests]
  message("Subsetting the datasets to extract the: ", toString(requests), " data.")

  kept_idx <- rep(FALSE, length(names(sig_data)))
  for (kept in requests) {
    keepers <- grepl(x = names(sig_data), pattern = glue("^{kept}"))
    kept_idx <- kept_idx | keepers
  }

  message("After subsetting, ", length(remaining), " entries remain.")
  remaining <- sig_data[kept_idx]
  return(remaining)
}

#' Create a metadata dataframe of msigdb data, this hopefully will be usable to
#' fill the fData slot of a gsva returned expressionset.
#'
#' @param gsva_result Some data from GSVA to modify.
#' @param msig_xml msig XML file downloaded from broad.
#' @return list containing 2 data frames: all metadata from broad, and the set
#'  matching the sig_data GeneSets.
#' @export
get_msigdb_metadata <- function(gsva_result = NULL, msig_xml = "msigdb_v6.2.xml",
                                wanted_meta = c("ORGANISM", "DESCRIPTION_BRIEF", "AUTHORS", "PMID")) {
  msig_result <- xml2::read_xml(x = msig_xml)

  db_data <- rvest::xml_nodes(x = msig_result, xpath = "//MSIGDB")
  db_name <- rvest::html_attr(x = db_data, name = "NAME")
  db_ver <- rvest::html_attr(x = db_data, name = "VERSION")
  db_date <- rvest::html_attr(x = db_data, name = "BUILD_DATE")

  genesets <- rvest::xml_nodes(x = msig_result, "GENESET")
  row_names <- rvest::html_attr(x = genesets, name = "STANDARD_NAME")
  column_names <- names(rvest::html_attrs(x = genesets[[1]]))
  all_data <- data.frame(row.names = row_names)
  for (i in 2:length(column_names)) {
    c_name <- column_names[i]
    c_data <- rvest::html_attr(x = genesets, name = c_name)
    all_data[[c_name]] <- c_data
  }

  ## all_data is my dataframe of xml annotations, lets extract the wanted columns before
  ## adding it to the gsva result expressionset.
  sub_data <- all_data
  if (!is.null(wanted_meta)) {
    sub_data <- sub_data[, wanted_meta]
  }
  ## The full set of xml annotations has waay more IDs than any single category,
  ## so this should drop it from 30,000+ to <10,000.
  found_idx <- rownames(sub_data) %in% rownames(exprs(gsva_result))
  sub_data <- sub_data[found_idx, ]

  retlist <- list(
      "all_data" = all_data,
      "sub_data" = sub_data)
  if (!is.null(gsva_result)) {
    current_fdata <- fData(gsva_result)
    new_fdata <- merge(current_fdata, sub_data, by = "row.names")
    rownames(new_fdata) <- new_fdata[["Row.names"]]
    new_fdata[["Row.names"]] <- NULL
    fData(gsva_result) <- new_fdata
    retlist[["gsva_result"]] <- gsva_result
  }
  return(retlist)
}

#' Attempt to score the results from simple_gsva()
#'
#' This function uses a couple of methods to try to get an idea of whether the
#' results from gsva are actually interesting.  It does so via the following
#' methods:
#'   1.  Use limma on the expressionset returned by simple_gsva(), this might
#' provide an idea of if there are changing signatures among the sample types.
#'   2.  Perform a simplified likelihood estimate to get a sense of the
#' significant categories.
#'
#' @param gsva_result Result from simple_gsva()
#' @param cutoff Significance cutoff
#' @param excel Excel file to write the results.
#' @param model_batch Add batch to limma's model.
#' @param factor_column When extracting significance information, use this
#'  metadata factor.
#' @param factor Use this metadata factor as the reference.
#' @param label_size Used to make the category names easier to read at the expense
#'  of dropping some.
#' @param col_margin Attempt to make heatmaps fit better on the screen with this and...
#' @param row_margin this parameter
#' @export
get_sig_gsva_categories <- function(gsva_result, cutoff = 0.95, excel = "excel/gsva_subset.xlsx",
                                    model_batch = FALSE, factor_column = "condition", factor = NULL,
                                    label_size = NULL, col_margin = 6, row_margin = 12) {

  gsva_scores <- gsva_result[["expt"]]

  ## Use limma on the gsva result
  gsva_limma <- limma_pairwise(gsva_scores, model_batch = model_batch,
                               which_voom = "none")

  gsva_eset <- gsva_scores[["expressionset"]]
  ## Go from highest to lowest score, using the first sample as a guide.
  values <- as.data.frame(exprs(gsva_eset))
  annot <- fData(gsva_eset)
  meta <- pData(gsva_eset)

  ## Choose the reference factor
  clevels <- levels(as.factor(meta[[factor_column]]))
  fact <- factor
  if (is.null(factor)) {
    fact <- clevels[1]
  }

  ## Copy the gsva expressionset and use that to pull the 'significant' entries.
  subset_eset <- gsva_eset
  gl <- score_gsva_likelihoods(gsva_result, factor = fact)
  likelihoods <- gl[["likelihoods"]]
  keep_idx <- likelihoods[[fact]] >= cutoff
  subset_eset <- subset_eset[keep_idx, ]
  jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  scored_ht <- NULL
  if (is.null(label_size)) {
    scored_ht <- heatmap.3(exprs(subset_eset), trace = "none", col = jet_colors,
                           margins = c(col_margin, row_margin))
  } else {
    scored_ht <- heatmap.3(exprs(subset_eset), trace = "none", col = jet_colors,
                           margins = c(col_margin, row_margin),
                           cexCol = label_size, cexRow = label_size)
  }
  scored_ht_plot <- grDevices::recordPlot()

  order_column <- colnames(values)[1]
  gsva_table <- merge(annot, values, by = "row.names")
  rownames(gsva_table) <- gsva_table[["Row.names"]]
  gsva_table[["Row.names"]] <- NULL
  reordered_gsva_idx <- order(gsva_table[[order_column]], decreasing = TRUE)
  gsva_table <- gsva_table[reordered_gsva_idx, ]

  ## Here is a bizarre little fact: the rownames of fData(subset_mtrx)
  ## are not the same as the rownames of exprs(subset_mtrx)
  ## which AFAIK should not be possible, but clearly I was wrong.
  ## At least when I create an expressionset, the fData and exprs have the
  ## same rownames from beginning to end.
  ## I guess this does not really matter, since we can use the full annotation table.
  subset_tbl <- as.data.frame(exprs(subset_eset))
  order_column <- colnames(subset_tbl)[1]
  subset_table <- merge(annot, subset_tbl, by = "row.names")
  rownames(subset_table) <- subset_table[["Row.names"]]
  subset_table[["Row.names"]] <- NULL
  reordered_subset_idx <- order(subset_table[[order_column]], decreasing = TRUE)
  subset_table <- subset_table[reordered_subset_idx, ]

  gl_tbl <- as.data.frame(gl[["likelihoods"]])
  order_column <- colnames(gl_tbl)[1]
  likelihood_table <- merge(annot, gl_tbl, by = "row.names")
  rownames(likelihood_table) <- likelihood_table[["Row.names"]]
  likelihood_table[["Row.names"]] <- NULL
  likelihood_table_idx <- order(likelihood_table[[order_column]], decreasing = TRUE)
  likelihood_table <- likelihood_table[likelihood_table_idx, ]

  retlist <- list(
      ## Everything provided by simple_gsva()
      "input" = gsva_result,
      ## The table from simple_gsva merged with the annotations.
      "gsva_table" = gsva_table,
      ## Heatmap of gsva result.
      "raw_plot" = gl[["raw_plot"]],
      ## The result from gsva_likelihoods, which compares condition vs. others.
      "likelihood_table" = likelihood_table,
      ## Corresponding plot from gsva_likelihoods
      "score_plot" = gl[["likelihood_plot"]],
      ## The subset of gsva scores deemed 'significant' by gsva_likelihoods.
      "subset_table" = subset_table,
      ## The corresponding plot for the subset.
      "subset_plot" = scored_ht_plot,
      "scores" = gl,
      "score_pca" = gl[["pca"]][["plot"]],
      "subset_expt" = subset_eset,
      "gsva_limma" = gsva_limma)

  if (!is.null(excel)) {
    retlist[["excel"]] <- write_gsva(retlist, excel)
  }
  return(retlist)
}

#' Take a result from simple_gsva(), a list of gene IDs, and intersect them.
#'
#' Najib is curious about the relationship of genes in sets, the sets, and the
#' genes that comprise those sets.  This is pushing gsva towards a oroborous-ish
#' state.
#'
#' @param gsva_result  Result from simple_gsva().
#' @param lst  List of genes of interest.
#' @param freq_cutoff  Minimum number of observations to be counted.
#' @param sig_weights  When making venn diagrams, weight them?
#' @param gene_weights  When venning genes, weight them?
#' @return  List containing some venns, lists, and such.
#' @export
intersect_signatures <- function(gsva_result, lst, freq_cutoff = 2,
                                 sig_weights = TRUE, gene_weights = TRUE) {
  sig_venn <- Vennerable::Venn(Sets = lst)
  Vennerable::plot(sig_venn, doWeights = sig_weights)
  sig_plot <- grDevices::recordPlot()
  sig_int <- sig_venn@IntersectionSets
  annot <- fData(gsva_result[["expt"]])
  sig_genes <- list()
  gene_venn_lst <- list()
  venn_names <- list()
  ## Skip the non-existant set of 00, thus 2:length()
  ## Top level loop iterates through the observed intersections/unions from Vennerable.
  for (i in 2:length(sig_int)) {
    name <- names(sig_int)[i]
    ## Make a human readable version of the venn names.
    name_chars <- strsplit(x = name, split = "")[[1]]
    venn_name <- ""
    for (c in 1:length(name_chars)) {
      char <- name_chars[[c]]
      if (char == "1") {
        venn_name <- glue("{venn_name}_{names(lst)[c]}")
      }
    }
    venn_name <- gsub(pattern = "^_", replacement = "", x = venn_name)

    sigs <- sig_int[[i]]
    sig_annot <- annot[sigs, ]
    gene_ids <- sig_annot[["ids"]]
    internal_ret <- list()
    ## This loop iterates through the set of observed gene IDs in each intersection
    for (j in 1:length(gene_ids)) {
      id_lst <- gene_ids[j]
      ids <- strsplit(id_lst, ", ")[[1]]
      ## Finally, we count how many times each id is observed in each signature
      for (k in 1:length(ids)) {
        element <- ids[k]
        if (is.null(internal_ret[[element]])) {
          internal_ret[[element]] <- 1
        } else {
          internal_ret[[element]] <- internal_ret[[element]] + 1
        }
      }
    }
    internal_ret <- internal_ret[order(as.numeric(internal_ret), decreasing = TRUE)]
    ret <- as.numeric(internal_ret)
    names(ret) <- names(internal_ret)
    sig_genes[[venn_name]] <- ret
    gene_venn_lst[[venn_name]] <- names(ret[ret >= freq_cutoff])
  }

  gene_venn <- Vennerable::Venn(Sets = gene_venn_lst)
  Vennerable::plot(gene_venn, doWeights = gene_weights)
  gene_venn_plot <- grDevices::recordPlot()
  gene_int <- gene_venn@IntersectionSets

  retlst <- list(
      "signature_venn" = sig_venn,
      "signature_intersection" = sig_int,
      "signature_venn_plot" = sig_plot,
      "signature_genes" = sig_genes,
      "gene_venn" = gene_venn,
      "gene_intersection" = gene_int,
      "gene_venn_plot" = gene_venn_plot)
  return(retlst)
}

#' Load signatures from either a gmt file, xml file, or directly from the
#' GSVAdata data set in R.
#'
#' There are a bunch of places from which to acquire signature data.  This
#' function attempts to provide a single place to load them.  The easiest way to
#' get up to date signatures is to download them from msigdb and set the
#' signatures parameter to the downloaded filename.
#'
#' @param signatures Either the filename downloaded or the variable's name as
#'  found in the environment created by data_pkg.
#' @param data_pkg Used when signatures is not a filename to load a data
#'  package, presumably GSVAdata.
#' @param signature_category Probably not needed unless you download a signature
#'  file containing lots of different categories.
#' @return signature dataset which may be used by gsva()
#' @export
load_gmt_signatures <- function(signatures = "c2BroadSets", data_pkg = "GSVAdata",
                                signature_category = "c2") {
  sig_data <- NULL
  if (class(signatures)[1] == "character" && grepl(pattern = "\\.gmt$", x = signatures)) {
    sig_data <- GSEABase::getGmt(signatures,
                                 collectionType = GSEABase::BroadCollection(category = signature_category),
                                 geneIdType = GSEABase::EntrezIdentifier())
  } else if (class(signatures)[1] == "character" && grepl(pattern = "\\.xml$", x = signatures)) {
    gsc <- GSEABase::getBroadSets(signatures)
    types <- sapply(gsc, function(elt) GSEABase::bcCategory(GSEABase::collectionType(elt)))
    sig_data <- gsc[types == signature_category]
  } else if (class(signatures)[1] == "character") {
    lib_result <- sm(requireNamespace(data_pkg))
    att_result <- sm(try(attachNamespace(data_pkg), silent = TRUE))
    lst <- list("list" = signatures, "package" = data_pkg)
    test <- do.call("data", as.list(signatures, lst))
    sig_data <- get0(signatures)
  } else if (class(signatures)[1] != "GeneSetCollection") {
    stop("The data must be a GeneSetCollection.")
  } else {
    sig_data <- signatures
  }
  return(sig_data)
}

#' Create a gene set collection from a set of arbitrary IDs.
#'
#' This function attempts to simplify the creation of a gsva compatible
#' GeneSet.  Some important caveats when working with gsva, notably the gene IDs
#' we use are not usually compatible with the gene IDs used by gsva, thus the
#' primary logic in this function is intended to bridge these IDs.
#'
#' @param first_ids The required IDs for a single set.
#' @param second_ids Potentially null optionally used for a second, presumably
#'  contrasting set.
#' @param orgdb Orgdb annotation, used to translate IDs to the required type.
#' @param researcher_name Prefix of the name for the generated set(s).
#' @param study_name Second element in the name of the generated set(s).
#' @param category_name Third element in the name of the generated set(s).
#' @param phenotype_name Optional phenotype data for the generated set(s).
#' @param pair_names The suffix of the generated set(s).
#' @param current_id What type of ID is the data currently using?
#' @param required_id What type of ID should the use?
#' @return Small list comprised of the created gene set collection(s).
#' @export
make_gsc_from_ids <- function(first_ids, second_ids = NULL, orgdb = "org.Hs.eg.db",
                              researcher_name = "elsayed", study_name = "macrophage",
                              category_name = "infection", phenotype_name = NULL,
                              pair_names = "up", current_id = "ENSEMBL", required_id = "ENTREZID") {
  first <- NULL
  second <- NULL
  if (is.null(current_id) | is.null(required_id)) {
    first <- first_ids
    current_id <- ""
    second <- second_ids
    required_id <- ""
  }
  if (current_id == required_id) {
    first <- first_ids
    second <- second_ids
  } else {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    ## tt <- sm(try(do.call("library", as.list(orgdb)), silent = TRUE))
    lib_result <- sm(requireNamespace(orgdb))
    att_restul <- sm(try(attachNamespace(orgdb), silent = TRUE))
    first_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                          keys = first_ids,
                                          keytype = current_id,
                                          columns = c(required_id)))
    first_idx <- complete.cases(first_ids)
    first_ids <- first_ids[first_idx, ]
    first <- first_ids[[required_id]]
    if (!is.null(second_ids)) {
      second_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                             keys = second_ids,
                                             keytype = current_id,
                                             columns = c(required_id)))
      second_idx <- complete.cases(second_ids)
      second_ids <- second_ids[second_idx, ]
      second <- second_ids[[required_id]]
    } else {
      second <- NULL
    }
  }

  all_colored <- NULL
  sec_gsc <- NULL
  fst <- data.frame(row.names = unique(first))
  if (is.null(phenotype_name)) {
    phenotype_name <- "unknown"
  }
  fst[["direction"]] <- pair_names[1]
  fst[["phenotype"]] <- phenotype_name

  set_prefix <- glue("{researcher_name}_{study_name}_{category_name}")
  fst_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))
  fst_gsc <- GSEABase::GeneSet(
                           GSEABase::EntrezIdentifier(),
                           setName = fst_name,
                           geneIds = as.character(rownames(fst)))
  if (!is.null(second)) {
    sec <- data.frame(row.names = unique(second))
    if (is.null(phenotype_name)) {
      phenotype_name <- "unknown"
    }
    if (is.na(pair_names[2]) & pair_names[1] == "up") {
      message("Setting the second pair_names to 'down'.")
      ## Then we can assume it is down.
      pair_names <- c("up", "down")
    }
    sec_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
    sec[["direction"]] <- pair_names[2]
    sec[["phenotype"]] <- phenotype_name
    both <- rbind(fst, sec)
    color_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    sec_gsc <- GSEABase::GeneSet(
                             GSEABase::EntrezIdentifier(),
                             setName = sec_name,
                             geneIds = as.character(rownames(sec)))
    all_colored = GSEABase::GeneColorSet(
                                GSEABase::EntrezIdentifier(),
                                setName = color_name,
                                geneIds = rownames(both),
                                phenotype = phenotype_name,
                                geneColor = as.factor(both[["direction"]]),
                                phenotypeColor = as.factor(both[["phenotype"]]))
  }
  retlst <- list()
  retlst[[fst_name]] <- fst_gsc
  if (!is.na(pair_names[2])) {
    ## Then there is a colored/down
    retlst[[sec_name]] <- sec_gsc
    retlst[["colored"]] <- all_colored
  }
  return(retlst)
}

#' Given a pairwise result, make a gene set collection.
#'
#' If I want to play with gsva and friends, then I need GeneSetCollections!  To
#' that end, this function uses extract_significant_genes() in order to gather
#' sets of genes deemed 'significant'.  It then passes these sets to
#' make_gsc_from_ids().
#'
#' @param pairwise A pairwise result, or combined de result, or extracted genes.
#' @param according_to When getting significant genes, use this method.
#' @param orgdb Annotation dataset.
#' @param pair_names Describe the contrasts of the GSC: up vs. down, high vs. low, etc.
#' @param category_name What category does the GSC describe?
#' @param phenotype_name When making color sets, use this phenotype name.
#' @param set_name A name for the created gene set.
#' @param color Make a colorSet?
#' @param current_id Usually we use ensembl IDs, but that does not _need_ to be the case.
#' @param required_id gsva uses entrezids by default.
#' @param ... Extra arguments for extract_significant_genes().
#' @return List containing 3 GSCs, one containing both the ups/downs called
#'  'colored', one of the ups, and one of the downs.
#' @export
make_gsc_from_pairwise <- function(pairwise, according_to = "deseq", orgdb = "org.Hs.eg.db",
                                   pair_names = c("ups", "downs"), category_name = "infection",
                                   phenotype_name = "parasite", set_name = "elsayed_macrophage",
                                   color = TRUE, current_id = "ENSEMBL", required_id = "ENTREZID",
                                   ...) {
  ups <- list()
  downs <- list()
  if (class(pairwise)[1] == "data.frame") {
    ups <- pairwise
  } else if (class(pairwise)[1] == "all_pairwise") {
    message("Invoking combine_de_tables().")
    combined <- sm(combine_de_tables(pairwise, ...))
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(
        combined,
        according_to = according_to, ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(
        pairwise,
        according_to = according_to,
        ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "sig_genes") {
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, orgdb = orgdb,
                             pair_names = pair_names, category_name = category_name,
                             phenotype_name = phenotype_name, current_id = current_id,
                             required_id = required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  ## tt <- sm(library(orgdb, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  up_lst <- list()
  down_lst <- list()
  colored_lst <- list()
  for (c in 1:length(ups)) {
    name <- names(ups)[c]
    up <- ups[[name]]
    down <- downs[[name]]
    if (class(up) == "character") {
      up_ids <- up
      down_ids <- down
    } else if (class(up) == "data.frame") {
      up_ids <- rownames(up)
      down_ids <- rownames(down)
    }
    if (current_id == required_id) {
      up[[required_id]] <- rownames(up)
      down[[required_id]] <- rownames(down)
    } else {
      message("Converting the rownames() of the expressionset to ENTREZID.")
      up_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                         keys = up_ids,
                                         keytype = current_id,
                                         columns = c(required_id)))
      up_idx <- complete.cases(up_ids)
      up_ids <- up_ids[up_idx, ]
      up <- merge(up, up_ids, by.x = "row.names", by.y = current_id)
      if (!is.null(down_ids)) {
        down_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                             keys = down_ids,
                                             keytype = current_id,
                                             columns = c(required_id)))
        down_idx <- complete.cases(down_ids)
        down_ids <- down_ids[down_idx, ]
        down <- merge(down, down_ids, by.x = "row.names", by.y = current_id)
      } else {
        down <- NULL
      }
    }
    up[["direction"]] <- pair_names[1]
    up[["phenotype"]] <- phenotype_name
    down[["direction"]] <- pair_names[2]
    down[["phenotype"]] <- phenotype_name

    both <- rbind(up, down)
    shared_ids <- up[[required_id]] %in% down[[required_id]]
    if (sum(shared_ids) > 0) {
      warning("There are ", sum(shared_ids), " shared IDs in the up and down lists.")
      shared_id <- up[shared_ids, ][[required_id]]
      shared_rows <- both[[required_id]] == shared_id
      both <- both[!shared_rows, ]
    }

    dup_elements <- duplicated(up[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the up data.")
      up <- up[!dup_elements, ]
    }
    dup_elements <- duplicated(down[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the down data.")
      down <- down[!dup_elements, ]
    }
    dup_elements <- duplicated(both[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the shared.")
      both <- both[!dup_elements, ]
    }

    set_prefix <- glue("{set_name}_{category_name}")
    color_set_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    up_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))
    colored_gsc <- GSEABase::GeneColorSet(
                                 GSEABase::EntrezIdentifier(),
                                 setName = color_set_name,
                                 geneIds = as.character(both[[required_id]]),
                                 phenotype = phenotype_name,
                                 geneColor = as.factor(both[["direction"]]),
                                 phenotypeColor = as.factor(both[["phenotype"]]))
    colored_lst[[name]] <- colored_gsc
    up_gsc <- GSEABase::GeneSet(
                            GSEABase::EntrezIdentifier(),
                            setName = up_name,
                            geneIds = as.character(up[[required_id]]))
    up_lst[[name]] <- up_gsc
    down_gsc <- NULL
    down_lst[[name]] <- down_gsc
    if (!is.null(pair_names[2])) {
      down_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
      down_gsc <- GSEABase::GeneSet(
                                GSEABase::EntrezIdentifier(),
                                setName = down_name,
                                geneIds = as.character(down[[required_id]]))
      down_lst[[name]] <- down_gsc
    }
  } ## End of the for loop.

  retlst <- list(
      "colored" = colored_lst,
      "up" = up_lst,
      "down" = down_lst)
  return(retlst)
}

#' Given a pairwise result, make a gene set collection.
#'
#' If I want to play with gsva and friends, then I need GeneSetCollections!
#' Much like make_gsc_from_significant(), this function extract the genes deemed
#' 'abundant' and generates gene sets accordingly.
#'
#' @param pairwise A pairwise result, or combined de result, or extracted genes.
#' @param according_to When getting significant genes, use this method.
#' @param orgdb Annotation dataset.
#' @param researcher_name Prefix of the name for the generated set(s).
#' @param study_name Second element in the name of the generated set(s).
#' @param category_name Third element in the name of the generated set(s).
#' @param phenotype_name Optional phenotype data for the generated set(s).
#' @param pair_names The suffix of the generated set(s).
#' @param current_id What type of ID is the data currently using?
#' @param required_id What type of ID should the use?
#' @param ... Extra arguments for extract_abundant_genes().
#' @return List containing 3 GSCs, one containing both the highs/lows called
#'  'colored', one of the highs, and one of the lows.
#' @export
make_gsc_from_abundant <- function(pairwise, according_to = "deseq", orgdb = "org.Hs.eg.db",
                                   researcher_name = "elsayed", study_name = "macrophage",
                                   category_name = "infection", phenotype_name = NULL,
                                   pair_names = "high", current_id = "ENSEMBL",
                                   required_id = "ENTREZID", ...) {
  highs <- list()
  lows <- list()
  if (class(pairwise)[1] == "data.frame") {
    highs <- pairwise
  } else if (class(pairwise)[1] == "all_pairwise") {
    message("Invoking combine_de_tables().")
    combined <- sm(combine_de_tables(pairwise, ...))
    message("Invoking extract_significant_genes().")
    highs <- sm(extract_abundant_genes(
        combined, according_to = according_to, ...)[[according_to]])
    lows <- sm(extract_abundant_genes(
        combined, according_to = according_to, least = TRUE, ...)[[according_to]])
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    highs <- sm(extract_abundant_genes(
        pairwise, according_to = according_to, ...)[[according_to]])
    lows <- sm(extract_abundant_genes(
        pairwise, according_to = according_to, least = TRUE, ...)[[according_to]])
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, orgdb = orgdb,
                             pair_names = pair_names, category_name = category_name,
                             phenotype_name = phenotype_name,
                             current_id = current_id, required_id = required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  ## tt <- sm(library(orgdb, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  high_lst <- list()
  low_lst <- list()
  colored_lst <- list()
  for (c in 1:length(high_lst[["abundances"]])) {
    name <- names(high_lst[["abundances"]])[c]
    high <- high_lst[["abundances"]][[name]]
    low <- low_lst[["abundances"]][[name]]
    if (class(high) == "character") {
      high_ids <- high
      low_ids <- low
    } else if (class(high) == "data.frame") {
      high_ids <- rownames(high)
      low_ids <- rownames(low)
    }
    if (current_id == required_id) {
      high[[required_id]] <- rownames(high)
      low[[required_id]] <- rownames(low)
    } else {
      message("Converting the rownames() of the expressionset to ENTREZID.")
      high_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                           keys = high_ids,
                                           keytype = current_id,
                                           columns = c(required_id)))
      high_idx <- complete.cases(high_ids)
      high_ids <- high_ids[high_idx, ]
      high <- merge(high, high_ids, by.x = "row.names", by.y = current_id)
      if (!is.null(low_ids)) {
        low_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                            keys = low_ids,
                                            keytype = current_id,
                                            columns = c(required_id)))
        low_idx <- complete.cases(low_ids)
        low_ids <- low_ids[low_idx, ]
        low <- merge(low, low_ids, by.x = "row.names", by.y = current_id)
      } else {
        low <- NULL
      }
    }
    high[["direction"]] <- pair_names[1]
    high[["phenotype"]] <- phenotype_name
    low[["direction"]] <- pair_names[2]
    low[["phenotype"]] <- phenotype_name

    both <- rbind(high, low)
    shared_ids <- high[[required_id]] %in% low[[required_id]]
    if (sum(shared_ids) > 0) {
      warning("There are ", sum(shared_ids), " shared IDs in the high and low lists.")
      shared_id <- high[shared_ids, ][[required_id]]
      shared_rows <- both[[required_id]] == shared_id
      both <- both[!shared_rows, ]
    }

    duplicated_elements <- duplicated(high[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the high data.")
      high <- high[!duplicated_elements, ]
    }
    duplicated_elements <- duplicated(low[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the low data.")
      low <- low[!duplicated_elements, ]
    }
    duplicated_elements <- duplicated(both[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the shared.")
      both <- both[!duplicated_elements, ]
    }

    set_prefix <- glue("{set_name}_{category_name}")
    color_set_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    high_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))
    colored_gsc <- GSEABase::GeneColorSet(
                                 GSEABase::EntrezIdentifier(),
                                 setName = color_set_name,
                                 geneIds = as.character(both[[required_id]]),
                                 phenotype = phenotype_name,
                                 geneColor = as.factor(both[["direction"]]),
                                 phenotypeColor = as.factor(both[["phenotype"]]))
    colored_lst[[name]] <- colored_gsc
    high_gsc <- GSEABase::GeneSet(
                              GSEABase::EntrezIdentifier(),
                              setName = high_name,
                              geneIds = as.character(high[[required_id]]))
    high_lst[[name]] <- high_gsc
    low_gsc <- NULL
    low_lst[[name]] <- low_gsc
    if (!is.null(pair_names[2])) {
      low_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
      low_gsc <- GSEABase::GeneSet(
                               GSEABase::EntrezIdentifier(),
                               setName = low_name,
                               geneIds = as.character(low[[required_id]]))
      low_lst[[name]] <- low_gsc
    }
  } ## End of the for loop.

  retlst <- list(
      "colored" = colored_lst,
      "high" = high_lst,
      "low" = low_lst)
  return(retlst)
}

#' Score the results from simple_gsva().
#'
#' Yeah, this is a bit meta, but the scores from gsva seem a bit meaningless to
#' me, so I decided to look at the distribution of observed scores in some of my
#' data; I quickly realized that they follow a nicely normal distribution.
#' Therefore, I thought to calculate some scores of gsva() using that
#' information.
#'
#' The nicest thing in this, I think, is that it provides its scoring metric(s)
#' according to a few different possibilities, including:
#'   * the mean of samples found in an experimental factor
#'   * All provided scores against the distribution of observed scores as
#'     z-scores.
#'   * A single score against all scores.
#'   * Rows (gene sets) against the set of all gene sets.
#'
#' @param gsva_result Input result from simple_gsva()
#' @param score What type of scoring to perform, against a value, column, row?
#' @param category What category to use as baseline?
#' @param factor Which experimental factor to compare against?
#' @param sample Which sample to compare against?
#' @param factor_column When comparing against an experimental factor, which design
#'  column to use to find it?
#' @param method mean or median when when bringing together values?
#' @param label_size By default, enlarge the labels to readable at the cost of losing some.
#' @param col_margin Attempt to make heatmaps fit better on the screen with this and...
#' @param row_margin this parameter
#' @param cutoff Highlight only the categories deemed more significant than this.
#' @return The scores according to the provided category, factor, sample, or
#'  score(s).
#' @export
score_gsva_likelihoods <- function(gsva_result, score = NULL, category = NULL,
                                   factor = NULL, sample = NULL, factor_column = "condition",
                                   method = "mean", label_size = NULL,
                                   col_margin = 6, row_margin = 12, cutoff = 0.95) {
  values <- exprs(gsva_result[["expt"]])
  design <- pData(gsva_result[["expt"]])
  gsva_pca <- plot_pca(gsva_result[["expt"]])

  ## Start off with a plot of the gsva return values.
  color_range <- c("#00007F", "blue", "#007FFF", "cyan",
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  jet_colors <- grDevices::colorRampPalette(color_range)
  starting_ht <- NULL
  if (is.null(label_size)) {
    starting_ht <- heatmap.3(values, trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin))
  } else {
    starting_ht <- heatmap.3(values, trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin),
                             cexCol = label_size, cexRow = label_size)
  }
  starting_ht_plot <- grDevices::recordPlot()

  tests <- test_values <- against_values <- NULL
  choice <- NULL
  if (is.null(score) & is.null(category) & is.null(sample) & is.null(factor)) {
    message("Nothing was requested, examining the first column scores: ",
            colnames(values)[1], ".")
    sample <- 1
    tests <- as.numeric(values[, sample])
    choice <- "column"
  } else if (!is.null(factor_column)) {
    choice <- "against"
  } else if (!is.null(score)) {
    message("Examining the score: ", score, " against the data.")
    tests <- score
    choice <- "value"
  } else if (!is.null(category)) {
    message("Examining row: ", category, " against the data.")
    tests <- as.numeric(values[category, ])
    choice <- "row"
  } else {
    tests <- as.numeric(values[, sample])
    choice <- "column"
  }

  population_values <- values
  ## ok, this function will give 'better' scores for higher
  ## gsva scores, so that is good, I will therefore just say 'higher is better'.
  cheesy_likelihood <- function(test) {
    gsva_mean <- mean(population_values)
    gsva_sd <- sd(population_values)
    num_values <- length(population_values)
    pop_sd <- gsva_sd * sqrt((num_values - 1) / num_values)
    notz <- (test - gsva_mean) / pop_sd
    likelihood <- pnorm(notz)
    return(likelihood)
  }

  results <- test_results <- against_results <- t_vs_a_results <- score_plot <- NULL
  if (choice == "against") {
    message("Testing each factor against the others.")
    fact_lvls <- levels(as.factor(design[[factor_column]]))
    result_df <- data.frame()
    for (f in 1:length(fact_lvls)) {
      fact <- fact_lvls[f]
      message("Scoring ", fact, " against everything else.")
      sample_idx <- design[[factor_column]] == fact
      if (sum(sample_idx) == 0) {
        ## Just in case the factor has empty levels.
        next
      }

      test_values <- values[, sample_idx]
      against_values <- values[, !sample_idx]
      population_values <- against_values
      if (method == "mean") {
        if (sum(sample_idx) == 1) {
          tests <- test_values
        } else {
          tests <- rowMeans(test_values)
        }
        againsts <- rowMeans(against_values)
      } else {
        if (sum(sample_idx) == 1) {
          tests <- test_values
        } else {
          tests <- Biobase::rowMedians(test_values)
        }
        againsts <- Biobase::rowMedians(against_values)
      }
      a_column <- sapply(X = tests, FUN = cheesy_likelihood)
      if (f == 1) {
        result_df <- as.data.frame(a_column)
      } else {
        result_df <- cbind(result_df, a_column)
      }
    } ## End iterating over every level in the chosen factor.
    colnames(result_df) <- fact_lvls
    heat_colors <- grDevices::colorRampPalette(c("white", "black"))
    ht_result <- heatmap.3(as.matrix(result_df), trace = "none", col = heat_colors,
                           margins = c(col_margin, row_margin), Colv = FALSE,
                           dendrogram = "row", cexCol = label_size, cexRow = label_size)
    score_plot <- grDevices::recordPlot()
    test_results <- result_df
  } else if (choice == "column") {
    test_results <- sapply(X = tests, FUN = cheesy_likelihood)
    names(test_results) <- rownames(values)
    score_plot <- plot_histogram(test_results)
  } else if (choice == "row") {
    test_results <- sapply(X = tests, FUN = cheesy_likelihood)
    names(test_results) <- colnames(values)
    score_plot <- plot_histogram(test_results)
  } else {
    names(test_results) <- "value"
    score_plot <- plot_histogram(test_results)
  }

  retlist <- list(
      "pca" = gsva_pca,
      "raw_plot" = starting_ht_plot,
      "likelihoods" = test_results,
      "likelihood_plot" = score_plot)
  return(retlist)
}

#' Provide some defaults and guidance when attempting to use gsva.
#'
#' gsva seems to hold a tremendous amount of potential.  Unfortunately, it is
#' somewhat opaque and its requirements are difficult to pin down.  This
#' function will hopefully provide some of the requisite defaults and do some
#' sanity checking to make it more likely that a gsva analysis will succeed.
#'
#' @param expt Expt object to be analyzed.
#' @param signatures Provide an alternate set of signatures (GeneSetCollections)
#' @param data_pkg What package contains the requisite dataset?
#' @param signature_category Specify a subset category to extract from the signatures database.
#' @param cores How many CPUs to use?
#' @param current_id Where did the IDs of the genes come from?
#' @param required_id gsva (I assume) always requires ENTREZ IDs, but just in
#'  case this is a parameter.
#' @param min_catsize Minimum category size to consider interesting (passed to gsva()).
#' @param orgdb What is the data source for the rownames()?
#' @param method Which gsva method to use? Changed this from gsva to ssgsea
#'  because it was throwing segmentation faults.
#' @param kcdf Options for the gsva methods.
#' @param ranking another gsva option.
#' @return List containing three elements: first a modified expressionset using
#'  the result of gsva in place of the original expression data; second the
#'  result from gsva, and third a data frame of the annotation data for the
#'  gene sets in the expressionset.  This seems a bit redundant, perhaps I
#'  should revisit it?
#' @export
simple_gsva <- function(expt, signatures = "c2BroadSets", data_pkg = "GSVAdata",
                        signature_category = "c2", cores = 1, current_id = "ENSEMBL",
                        required_id = "ENTREZID", min_catsize = 5, orgdb = "org.Hs.eg.db",
                        method = "ssgsea", kcdf = NULL, ranking = FALSE, msig_xml = NULL,
                        wanted_meta = c("ORGANISM", "DESCRIPTION_BRIEF", "AUTHORS", "PMID")) {
  if (is.null(kcdf)) {
    if (expt[["state"]][["transform"]] == "raw") {
      kcdf <- "Poisson"
    } else {
      kcdf <- "Gaussian"
    }
  }

  if (!is.null(msig_xml) & !file.exists(msig_xml)) {
    stop("The msig_xml parameter was defined, but the file does not exist.")
  }

  ## Make sure some data is loaded.  I will no longer assume anything here.
  ## Here is how I will decide:
  ## 1.  If signatures is a (string)filename ending in '.gmt', then extract the genesetlists and use it.
  ## 2.  If signatures is a string, then load the data_pkg, presumably GSVAdata.
  ## 3.  If signatures is not a string, assume it is a genesetlist/geneset and use that.

  ## Assume the desired category is c2 unless specified.
  if (is.null(signature_category)) {
    signature_category <- "c2"
  }
  sig_data <- load_gmt_signatures(signatures = signatures, data_pkg = data_pkg,
                                  signature_category = signature_category)

  ## The expressionset must have the annotation field filled in for gsva to
  ## work.
  eset <- expt[["expressionset"]]
  eset_annotation <- annotation(eset)
  eset_pattern <- grepl(pattern = "Fill me in", x = annotation(eset))
  if (length(eset_annotation) == 0 | isTRUE(eset_pattern)) {
    message("gsva requires the annotation field to be filled in.")
    annotation(eset) <- orgdb
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  if (current_id != required_id | !is.integer(grep("ENSG", rownames(exprs(eset))))) {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    ##tt <- sm(library(orgdb, character.only = TRUE))
    lib_result <- sm(requireNamespace(orgdb))
    att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
    old_ids <- rownames(exprs(eset))
    new_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                        keys = old_ids,
                                        keytype = current_id,
                                        columns = c(required_id)))
    new_idx <- complete.cases(new_ids)
    if (!all(new_idx)) {
      message(sum(new_idx == FALSE),
              " ENSEMBL ID's didn't have a matching ENTEREZ ID. Dropping them now.")
    }

    new_ids <- new_ids[new_idx, ]
    message("Before conversion, the expressionset has ", length(rownames(eset)),
            " entries.")
    converted_eset <- eset[new_ids[[current_id]], ]
    rownames(converted_eset) <- make.names(new_ids[[required_id]], unique = TRUE)
    message("After conversion, the expressionset has ",
            length(rownames(converted_eset)),
            " entries.")
    rownames(converted_eset) <- gsub(pattern = "^X", replacement = "",
                                     x = rownames(converted_eset))
    eset <- converted_eset
    fData(eset)[[required_id]] <- rownames(fData(eset))
  }

  ## Responding to Theresa: cores is defined in the function definition.
  ## Sadly, some versions of gsva crash if one sets it to > 1, so for the moment
  ## it is set to 1 and gsva is not running in parallel, but I wanted to keep the
  ## possibility of speeding it up, ergo the cores option.
  gsva_result <- GSVA::gsva(eset, sig_data, verbose = TRUE, method = method,
                            min.sz = min_catsize, kcdf = kcdf, abs.ranking = ranking,
                            parallel.sz = cores)
  fdata_df <- data.frame(row.names = rownames(exprs(gsva_result)))
  fdata_df[["description"]] <- ""
  fdata_df[["ids"]] <- ""
  for (i in 1:length(sig_data)) {
    fdata_df[i, "description"] <- description(sig_data[[i]])
    fdata_df[i, "ids"] <- toString(GSEABase::geneIds(sig_data[[i]]))
  }

  fData(gsva_result) <- fdata_df
  new_expt <- expt
  new_expt[["expressionset"]] <- gsva_result
  if (!is.null(msig_xml)) {
    message("Adding annotations from ", msig_xml, ".")
    improved <- get_msigdb_metadata(msig_xml = msig_xml, wanted_meta = wanted_meta,
                                    gsva_result = gsva_result)
    new_expt[["expressionset"]] <- improved[["gsva_result"]]
  }

  retlist <- list(
      "method" = method,
      "signatures" = signatures,
      "signature_category" = signature_category,
      "required_id" = required_id,
      "min_catsize" = min_catsize,
      "expt" = new_expt,
      "gsva" = gsva_result,
      "fdata" = fdata_df)
  return(retlist)
}

#' Invoke xCell and pretty-ify the result.
#'
#' I initially thought xCell might prove the best tool/method for exploring cell
#' deconvolution.  I slowly figured out its limitations, but still think it
#' seems pretty nifty for its use case.  Thus this function is intended to make
#' invoking it easier/faster.
#'
#' @param expt Expressionset to query.
#' @param signatures Alternate set of signatures to use.
#' @param genes Subset of genes to query.
#' @param spill The xCell spill parameter.
#' @param expected_types Set of assumed types in the data.
#' @param label_size How large to make labels when printing the final heatmap.
#' @param col_margin Used by par() when printing the final heatmap.
#' @param row_margin Ibid.
#' @param cores How many CPUs to use?
#' @param ... Extra arguments when normalizing the data for use with xCell.
#' @return Small list providing the output from xCell, the set of signatures,
#'  and heatmap.
#' @export
simple_xcell <- function(expt, signatures = NULL, genes = NULL, spill = NULL,
                         expected_types = NULL, label_size = NULL, col_margin = 6,
                         row_margin = 12, sig_cutoff = 0.2, verbose = TRUE, cores = 4, ...) {
  arglist <- list(...)
  xcell_annot <- load_biomart_annotations()
  xref <- xcell_annot[["annotation"]][, c("ensembl_gene_id", "hgnc_symbol")]
  expt_state <- expt[["state"]][["conversion"]]
  xcell_eset <- NULL
  if (expt_state != "rpkm") {
    message("xCell strongly perfers rpkm values, re-normalizing now.")
    xcell_eset <- normalize_expt(expt, convert = "rpkm", ...)
  } else {
    xcell_eset <- normalize_expt(expt, norm = arglist[["norm"]], convert = arglist[["convert"]],
                                  filter = arglist[["filter"]], batch = arglist[["batch"]])
  }
  xcell_mtrx <- exprs(xcell_eset)
  xcell_na <- is.na(xcell_mtrx)
  xcell_mtrx[xcell_na] <- 0
  xcell_input <- merge(xcell_mtrx, xref, by.x = "row.names", by.y = "ensembl_gene_id")
  rownames(xcell_input) <- make.names(xcell_input[["hgnc_symbol"]], unique = TRUE)
  xcell_input[["Row.names"]] <- NULL
  xcell_input[["hgnc_symbol"]] <- NULL

  xCell.data <- NULL
  tt <- requireNamespace("xCell")
  data("xCell.data", package = "xCell")
  if (is.null(signatures)) {
    signatures <- xCell.data[["signatures"]]
  }
  if (is.null(genes)) {
    genes <- xCell.data[["genes"]]
  }
  if (is.null(spill)) {
    spill <- xCell.data[["spill"]]
  }

  xcell_result <- NULL
  if (isTRUE(verbose)) {
    xcell_result <- xCell::xCellAnalysis(expr = xcell_input, signatures = signatures,
                                         genes = genes, spill = spill,
                                         cell.types = expected_types, parallel.sz = cores)
  } else {
    xcell_result <- sm(xCell::xCellAnalysis(expr = xcell_input, signatures = signatures,
                                            genes = genes, spill = spill, parallel.sz = cores,
                                            cell.types = expected_types))
  }

  jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (is.null(label_size)) {
    ht <- heatmap.3(xcell_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin))
  } else {
    ht <- heatmap.3(xcell_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin),
                    cexCol = label_size, cexRow = label_size)
  }
  ht_plot <- grDevices::recordPlot()

  sig_idx <- Biobase::rowMax(xcell_result) >= sig_cutoff
  sig_result <- xcell_result[sig_idx, ]
  if (is.null(label_size)) {
    ht <- heatmap.3(sig_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin))
  } else {
    ht <- heatmap.3(sig_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin),
                    cexCol = label_size, cexRow = label_size)
  }
  sig_plot <- grDevices::recordPlot()

  retlist <- list(
      "xcell_input" = xcell_input,
      "xcell_result" = xcell_result,
      "signatures" = xCell.data[["signatures"]],
      "heatmap" = ht_plot,
      "sig_result" = sig_result,
      "sig_plot" = sig_plot)
  return(retlist)
}

#' Write out my various attempts at making sense of gsva.
#'
#' While I am trying to make sense of gsva, I will use this function to write
#' out the results I get so I can pass them to Najib/Maria Adelaida/Theresa to
#' see if I am making sense.
#'
#' @param retlist Result from running get_sig_gsva
#' @param excel Excel file to write
#' @param plot_dim Plot dimensions, likely needs adjustment.
#' @export
write_gsva <- function(retlist, excel, plot_dim = 6) {
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]

  methods <- list(
      "gsva" = "Hänzelmann et al, 2013",
      "ssgsea" = "Barbie et al, 2009",
      "zscore" = "Lee et al, 2008",
      "plage" = "Tomfohr et al, 2005")

  db_used <- retlist[["input"]][["signatures"]]
  if (class(db_used)[1] != "character") {
    db_used <- "user provided"
  }

  method <- retlist[["input"]][["method"]]
  if (method %in% names(methods)) {
    method <- paste0(method, ": ", methods[method])
  }

  ## Write the legend.
  legend <- data.frame(rbind(
      c("Signature database used:", db_used),
      c("Database subset used:", retlist[["input"]][["signature_category"]]),
      c("Required ID type:", retlist[["input"]][["required_id"]]),
      c("Minimum category size:", retlist[["input"]][["min_catsize"]]),
      c("GSVA method used:", method),
      c("", ""),
      c("Sheet 1: gsva_scores", "All scores as provided by gsva()."),
      c("Sheet 2: gsva_likelihoods", "All likelihood scores calculated using pnorm() of the values."),
      c("Sheet 3: factor_likelihoods", "Likelihood values for each experimental factor."),
      c("Sheet 4: subset", "GSVA scores for the categories deemed 'significant' using sheet 2/3."),
      c("Sheet 5 on:", "Limma scoring of differential signatures.")),
      stringsAsFactors = FALSE)
  colnames(legend) <- c("Term", "Definition")
  xls_result <- write_xlsx(
      wb, data = legend, sheet = "legend", rownames = FALSE,
      title = "Summary and sheets in this workbook.")
  xl_result <- openxlsx::writeData(wb = wb, sheet = "legend", x = "PCA of categories vs sample type.",
                                   startRow = 1, startCol = 8)
  try_result <- xlsx_plot_png(retlist[["score_pca"]], wb = wb, sheet = "legend",
                              start_row = 2, start_col = 8,
                              width=(plot_dim * 3/2), height = plot_dim,
                              plotname = "gsva_pca", savedir = excel_basename)

  ## Write the result from gsva()
  xls_result <- write_xlsx(data = retlist[["gsva_table"]], wb = wb, sheet = "gsva_scores")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 8
  plot_height <- 16
  try_result <- xlsx_plot_png(a_plot = retlist[["raw_plot"]], wb = wb, sheet = "gsva_scores",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  ## Write the likelihoods
  xls_result <- write_xlsx(data = retlist[["likelihood_table"]], wb = wb, sheet = "likelihood_scores")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 6
  plot_height <- 15
  try_result <- xlsx_plot_png(a_plot = retlist[["score_plot"]], wb = wb, sheet = "likelihood_scores",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  ## Write the subset
  xls_result <- write_xlsx(data = retlist[["subset_table"]], wb = wb, sheet = "subset_table")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 6
  plot_height <- 6
  try_result <- xlsx_plot_png(a_plot = retlist[["subset_plot"]], wb = wb, sheet = "subset_table",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  limma_tables <- retlist[["gsva_limma"]][["all_tables"]]
  table_names <- names(limma_tables)
  for (i in 1:length(table_names)) {
    table_name <- table_names[i]
    title <- glue::glue("Result from using limma to compare {table_name}.")
    table <- limma_tables[[table_name]]
    table_idx <- order(table[["adj.P.Val"]], decreasing = FALSE)
    table <- table[table_idx, ]
    xls_result <- write_xlsx(data = table, wb = wb, sheet = table_name, title = title)
  }

  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  if (class(save_result)[1] == "try-error") {
    message("Saving xlsx failed.")
  }
  return(save_result)
}

## EOF
