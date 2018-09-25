#' Create a gene set collection from a set of arbitrary IDs.
#'
#' @param first_ids  The required IDs for a single set.
#' @param second_ids  Potentially null optionally used for a second, presumably contrasting set.
#' @param orgdb  Orgdb annotation, used to translate IDs to the required type.
#' @param pair_names The suffix of the generated set(s).
#' @param category_name  Primarily used for colored gene sets, but describes the
#' category of data examined in this set.
#' @param phenotype_name  Also used for colored gene sets.
#' @param set_name  Primary identifier for the set.
#' @param current_id  What type of ID is the data currently using?
#' @param required_id  What type of ID should the use?
#' @export
make_gsc_from_ids <- function(first_ids, second_ids=NULL, orgdb="org.Hs.eg.db",
                              pair_names="up", category_name="infection",
                              phenotype_name=NULL, set_name="elsayed_macrophage",
                              current_id="ENSEMBL", required_id="ENTREZID") {
  first <- NULL
  second <- NULL
  if (current_id == required_id) {
    first <- first_ids
    second <- second_ids
  } else {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    first_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                       keys=first_ids,
                                       keytype=current_id,
                                       columns=c(required_id)))
    first_idx <- complete.cases(first_ids)
    first_ids <- first_ids[first_idx, ]
    first <- first_ids[[required_id]]
    if (!is.null(second_ids)) {
      second_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                             keys=second_ids,
                                             keytype=current_id,
                                             columns=c(required_id)))
      second_idx <- complete.cases(second_ids)
      second_ids <- second_ids[second_idx, ]
      second <- second_ids[[required_id]]
    } else {
      second <- NULL
    }
  }

  all_colored <- NULL
  sec_gsc <- NULL
  fst <- data.frame(row.names=unique(first))
  fst[["direction"]] <- directions[1]
  fst[["phenotype"]] <- phenotype_name

  set_prefix <- paste0(set_name, "_", category_name)
  fst_name <- paste0(set_prefix, "_", pair_names[1])
  fst_gsc <- GSEABase::GeneSet(
                         EntrezIdentifier(),
                         setName=fst_name,
                         geneIds=as.character(rownames(fst)))
  ##setName(fst_gsc) <- paste0(set_name, "_", category_name, "_", directions[1])
  if (!is.null(second)) {
    if (is.null(phenotype_name)) {
      phenotype_name <- "unknown"
    }
    sec[["direction"]] <- directions[2]
    sec[["phenotype"]] <- phenotype_name
    both <- rbind(fst, sec)
    color_name <- paste0(set_prefix, "_", phenotype_name)
    sec_name <- paste0(set_prefix, "_", pair_names[1])
    sec_gsc <- GSEABase::GeneSet(
                           EntrezIdentifier(),
                           setName=sec_name,
                           geneIds=as.character(rownames(sec)))
    all_colored = GSEABase::GeneColorSet(
                              EntrezIdentifier(),
                              setName=color_name,
                              geneIds=rownames(both),
                              phenotype=phenotype_name,
                              geneColor=as.factor(both[["direction"]]),
                              phenotypeColor=as.factor(both[["phenotype"]]))
  }
  retlst <- list(
    "first" = fst_gsc,
    "second" = sec_gsc,
    "colored" = all_colored)
  return(retlst)
}

#' Given a pairwise result, make a gene set collection.
#'
#' If I want to play with gsva and friends, then I need GeneSetCollections!
#'
#' @param pairwise  A pairwise result, or combined de result, or extracted genes.
#' @param according_to  When getting significant genes, use this method.
#' @param orgdb  Annotation dataset.
#' @param color  Make a colorSet?
#' @param current_id  Usually we use ensembl IDs, but that does not _need_ to be the case.
#' @param set_name  A name for the created gene set.
#' @param phenotype_name  When making color sets, use this phenotype name.
#' @param required_id  gsva uses entrezids by default.
#' @param ...  Extra arguments for extract_significant_genes().
#' @export
make_gsc_from_pairwise <- function(pairwise, according_to="deseq", orgdb="org.Hs.eg.db",
                                   pair_names=c("ups", "downs"), category_name="infection",
                                   phenotype_name="parasite", set_name="elsayed_macrophage",
                                   color=TRUE, current_id="ENSEMBL", required_id="ENTREZID", ...) {
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
      according_to=according_to, ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(
      pairwise,
      according_to=according_to,
      ...)[[according_to]])
    ##updown <- sm(extract_significant_genes(
    ##  pairwise,
    ##  according_to=according_to)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "sig_genes") {
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, orgdb=orgdb,
                             pair_names=pair_names, category_name=category_name,
                             phenotype_name=phenotype_name, set_name=set_name,
                             current_id=current_id, required_id=required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  tt <- sm(library(orgdb, character.only=TRUE))
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
      up_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                         keys=up_ids,
                                         keytype=current_id,
                                         columns=c(required_id)))
      up_idx <- complete.cases(up_ids)
      up_ids <- up_ids[up_idx, ]
      up <- merge(up, up_ids, by.x="row.names", by.y=current_id)
      if (!is.null(down_ids)) {
        down_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                             keys=down_ids,
                                             keytype=current_id,
                                             columns=c(required_id)))
        down_idx <- complete.cases(down_ids)
        down_ids <- down_ids[down_idx, ]
        down <- merge(down, down_ids, by.x="row.names", by.y=current_id)
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
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the up data.")
      up <- up[!dup_elements, ]
    }
    dup_elements <- duplicated(down[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the down data.")
      down <- down[!dup_elements, ]
    }
    dup_elements <- duplicated(both[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the shared.")
      both <- both[!dup_elements, ]
    }

    set_prefix <- paste0(set_name, "_", category_name)
    color_set_name <- paste0(set_prefix, "_", phenotype_name)
    up_name <- paste0(set_prefix, "_", pair_names[1])
    colored_gsc <- GSEABase::GeneColorSet(
                              EntrezIdentifier(),
                              setName=color_set_name,
                              geneIds=as.character(both[[required_id]]),
                              phenotype=phenotype_name,
                              geneColor=as.factor(both[["direction"]]),
                              phenotypeColor=as.factor(both[["phenotype"]]))
    colored_lst[[name]] <- colored_gsc
    up_gsc <- GSEABase::GeneSet(
                          EntrezIdentifier(),
                          setName=up_name,
                          geneIds=as.character(up[[required_id]]))
    up_lst[[name]] <- up_gsc
    down_gsc <- NULL
    down_lst[[name]] <- down_gsc
    if (!is.null(pair_names[2])) {
      down_name <- paste0(set_prefix, "_", pair_names[2])
      down_gsc <- GSEABase::GeneSet(
                              EntrezIdentifier(),
                              setName=down_name,
                              geneIds=as.character(down[[required_id]]))
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
#'
#' @param pairwise  A pairwise result, or combined de result, or extracted genes.
#' @param according_to  When getting significant genes, use this method.
#' @param orgdb  Annotation dataset.
#' @param color  Make a colorSet?
#' @param current_id  Usually we use ensembl IDs, but that does not _need_ to be the case.
#' @param set_name  A name for the created gene set.
#' @param phenotype_name  When making color sets, use this phenotype name.
#' @param required_id  gsva uses entrezids by default.
#' @param ...  Extra arguments for extract_significant_genes().
#' @export
make_gsc_from_abundant <- function(pairwise, according_to="deseq", orgdb="org.Hs.eg.db",
                                   pair_names=c("ups", "downs"), category_name="infection",
                                   phenotype_name="parasite", set_name="elsayed_macrophage",
                                   color=TRUE, current_id="ENSEMBL", required_id="ENTREZID", ...) {
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
      according_to=according_to, ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(
      pairwise,
      according_to=according_to,
      ...)[[according_to]])
    ##updown <- sm(extract_significant_genes(
    ##  pairwise,
    ##  according_to=according_to)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "sig_genes") {
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, orgdb=orgdb,
                             pair_names=pair_names, category_name=category_name,
                             phenotype_name=phenotype_name, set_name=set_name,
                             current_id=current_id, required_id=required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  tt <- sm(library(orgdb, character.only=TRUE))
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
      up_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                         keys=up_ids,
                                         keytype=current_id,
                                         columns=c(required_id)))
      up_idx <- complete.cases(up_ids)
      up_ids <- up_ids[up_idx, ]
      up <- merge(up, up_ids, by.x="row.names", by.y=current_id)
      if (!is.null(down_ids)) {
        down_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                             keys=down_ids,
                                             keytype=current_id,
                                             columns=c(required_id)))
        down_idx <- complete.cases(down_ids)
        down_ids <- down_ids[down_idx, ]
        down <- merge(down, down_ids, by.x="row.names", by.y=current_id)
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
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the up data.")
      up <- up[!dup_elements, ]
    }
    dup_elements <- duplicated(down[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the down data.")
      down <- down[!dup_elements, ]
    }
    dup_elements <- duplicated(both[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements), " non-unique elements in the IDs of the shared.")
      both <- both[!dup_elements, ]
    }

    set_prefix <- paste0(set_name, "_", category_name)
    color_set_name <- paste0(set_prefix, "_", phenotype_name)
    up_name <- paste0(set_prefix, "_", pair_names[1])
    colored_gsc <- GSEABase::GeneColorSet(
                              EntrezIdentifier(),
                              setName=color_set_name,
                              geneIds=as.character(both[[required_id]]),
                              phenotype=phenotype_name,
                              geneColor=as.factor(both[["direction"]]),
                              phenotypeColor=as.factor(both[["phenotype"]]))
    colored_lst[[name]] <- colored_gsc
    up_gsc <- GSEABase::GeneSet(
                          EntrezIdentifier(),
                          setName=up_name,
                          geneIds=as.character(up[[required_id]]))
    up_lst[[name]] <- up_gsc
    down_gsc <- NULL
    down_lst[[name]] <- down_gsc
    if (!is.null(pair_names[2])) {
      down_name <- paste0(set_prefix, "_", pair_names[2])
      down_gsc <- GSEABase::GeneSet(
                              EntrezIdentifier(),
                              setName=down_name,
                              geneIds=as.character(down[[required_id]]))
      down_lst[[name]] <- down_gsc
    }
  } ## End of the for loop.

  retlst <- list(
    "colored" = colored_lst,
    "up" = up_lst,
    "down" = down_lst)
  return(retlst)
}

#' Provide some defaults and guidance when attempting to use gsva.
#'
#' gsva seems to hold a tremendous amount of potential.  Unfortunately, it is
#' somewhat opaque and its requirements are difficult to pin down.  This
#' function will hopefully provide some of the requisite defaults and do some
#' sanity checking to make it more likely that a gsva analysis will succeed.
#'
#' @param expt  Expt object to be analyzed.
#' @param datasets  Name of the variable from which to acquire the gsva data, if
#'   it does not exist, then data() will be called upon it.
#' @param data_pkg  What package contains the requisite dataset?
#' @param signatures Provide an alternate set of signatures (GeneSetCollections)
#' @param current_id  Where did the IDs of the genes come from?
#' @param required_id  gsva (I assume) always requires ENTREZ IDs, but just in
#'   case this is a parameter.
#' @param orgdb  What is the data source for the rownames()?
#' @return  Something from GSVA::gsva()!
#' @export
simple_gsva <- function(expt, datasets="c2BroadSets", data_pkg="GSVAdata", signatures=NULL,
                        cores=0, current_id="ENSEMBL", required_id="ENTREZID",
                        orgdb="org.Hs.eg.db", method="gsva", kcdf=NULL, ranking=FALSE) {
  if (is.null(kcdf)) {
    if (expt[["state"]][["transform"]] == "raw") {
      kcdf <- "Poisson"
    } else {
      kcdf <- "Gaussian"
    }
  }
  ## Make sure some data is loaded.  Assume the c2BroadSets from GSVAdata.
  sig_data <- NULL
  if (!is.null(signatures)) {
    sig_data <- signatures
  } else {
    if (exists(datasets)) {
      sig_data <- datasets
      if (class(sig_data)[[1]] == "character") {
        tt <- sm(library(data_pkg, character.only=TRUE))
        lst <- list("list"=datasets, "package"=data_pkg)
        test <- do.call("data", as.list(datasets, lst))
        sig_data <- get0(datasets)
      }
      if (class(sig_data)[[1]] != "GeneSetCollection") {
        stop("The data must be a GeneSetCollection.")
      }
    } else {
      tt <- sm(library(data_pkg, character.only=TRUE))
      lst <- list("list"=datasets, "package"=data_pkg)
      test <- do.call("data", as.list(datasets, lst))
      sig_data <- get0(datasets)
    }
  }

  ## The expressionset must have the annotation field filled in for gsva to
  ## work.
  eset <- expt[["expressionset"]]
  if (length(annotation(eset)) == 0 |
      grep(pattern="Fill me in", x=annotation(eset))) {
    message("gsva requires the annotation field to be filled in.")
    annotation(eset) <- orgdb
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  if (current_id != required_id) {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    tt <- sm(library(orgdb, character.only=TRUE))
    old_ids <- rownames(exprs(eset))
    new_ids <- sm(AnnotationDbi::select(x=get0(orgdb),
                                        keys=old_ids,
                                        keytype=current_id,
                                        columns=c(required_id)))
    new_idx <- complete.cases(new_ids)
    new_ids <- new_ids[new_idx, ]
    message("Before conversion, the expressionset has ", length(rownames(eset)),
  " entries.")
    converted_eset <- eset[new_ids[[current_id]], ]
    rownames(converted_eset) <- make.names(new_ids[[required_id]], unique=TRUE)
    message("After conversion, the expressionset has ",
            length(rownames(converted_eset)),
            " entries.")
    rownames(converted_eset) <- gsub(pattern="^X", replacement="",
                                     x=rownames(converted_eset))
    eset <- converted_eset
    fData(eset)[[required_id]] <- rownames(fData(eset))
  }

  gsva_result <- GSVA::gsva(eset, sig_data, verbose=TRUE, method=method,
                            kcdf=kcdf, abs.ranking=ranking, parallel.sz=cores)
  fdata_df <- data.frame(row.names=rownames(exprs(gsva_result)))
  fdata_df[["description"]] <- ""
  fdata_df[["ids"]] <- ""
  for (i in 1:length(sig_data)) {
    fdata_df[i, "description"] <- description(sig_data[[i]])
    fdata_df[i, "ids"] <- toString(GSEABase::geneIds(sig_data[[i]]))
  }

  fData(gsva_result) <- fdata_df
  new_expt <- expt
  new_expt[["expressionset"]] <- gsva_result

  retlist <- list(
    "expt" = new_expt,
    "gsva" = gsva_result,
    "fdata" = fdata_df)
  return(retlist)
}

simple_xcell <- function(expt, label_size=NULL, col_margin=6, row_margin=12, ...) {
  arglist <- list(...)
  xcell_annot <- load_biomart_annotations()
  xref <- xcell_annot[["annotation"]][, c("ensembl_gene_id", "hgnc_symbol")]
  expt_state <- expt[["state"]][["conversion"]]
  xcell_input <- NULL
  if (expt_state != "rpkm") {
    message("xCell strongly perfers rpkm values, re-normalizing now.")
    xcell_input <- normalize_expt(expt, convert="rpkm", ...)
  } else {
    xcell_input <- normalize_expt(expt, norm=arglist[["norm"]], convert=arglist[["convert"]],
                                  filter=arglist[["filter"]], batch=arglist[["batch"]])
  }
  xcell_input <- exprs(xcell_input)
  xcell_na <- is.na(xcell_input)
  xcell_input[xcell_na] <- 0
  xcell_input <- merge(xcell_input, xref, by.x="row.names", by.y="ensembl_gene_id")
  rownames(xcell_input) <- make.names(xcell_input[["hgnc_symbol"]], unique=TRUE)
  xcell_input[["Row.names"]] <- NULL
  xcell_input[["hgnc_symbol"]] <- NULL

  xCell.data <- NULL
  library(xCell)
  data("xCell.data", package="xCell")
  xcell_result <- xCell::xCellAnalysis(xcell_input)

  jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (is.null(label_size)) {
    ht <- heatmap.3(xcell_result, trace="none", col=jet_colors, margins=c(col_margin, row_margin))
  } else {
    ht <- heatmap.3(xcell_result, trace="none", col=jet_colors, margins=c(col_margin, row_margin),
                    cexCol=label_size, cexRow=label_size)
  }

  ht_plot <- grDevices::recordPlot()

  retlist <- list(
    "xcell_result" = xcell_result,
    "signatures" = xCell.data[["signatures"]],
    "heatmap" = ht_plot)
  return(retlist)
}

#' Extract the GeneSets corresponding to the provided name(s).
#'
#' I dunno what I want to put here and I am getting tired.
#'
#' @param sig_data  The pile of GeneSets probably from GSVAdata.
#' @param requests  Character list of sources.
#' @return  Whatever GeneSets remain.
#' @export
get_gsvadb_names <- function(sig_data, requests=NULL) {
  requests <- toupper(requests)
  categories <- names(sig_data)
  prefixes <- gsub(pattern="^(.*?)_.*", replacement="\\1",
                   x=categories)
  tab <- table(prefixes)
  if (is.null(requests)) {
    new_order <- order(as.numeric(tab), decreasing=TRUE)
    tab <- tab[new_order]
    return(tab)
  }

  kept_requests <- requests %in% names(tab)
  requests <- requests[kept_requests]
  message("Subsetting the datasets to extract the: ", toString(requests), " data.")

  kept_idx <- rep(FALSE, length(names(sig_data)))
  for (kept in requests) {
    keepers <- grepl(x=names(sig_data), pattern=paste0("^", kept))
    kept_idx <- kept_idx | keepers
  }

  message("After subsetting, ", length(remaining), " entries remain.")
  remaining <- sig_data[kept_idx]
  return(remaining)
}

#' Create a metadata dataframe of msigdb data, this hopefully will be usable to
#' fill the fData slot of a gsva returned expressionset.
#'
#' @param sig_data  GeneSetCollection from the broad msigdb.
#' @param msig_xml  msig XML file downloaded from broad.
#' @param gsva_result  Some data from GSVA to modify.
#' @return list containing 2 data frames: all metadata from broad, and the set
#'   matching the sig_data GeneSets.
#' @export
get_msigdb_metadata <- function(sig_data, msig_xml="msigdb_v6.2.xml", gsva_result=NULL) {
  msig_result <- xml2::read_xml(msig_xml)

  db_data <- rvest::xml_nodes(msig_result, xpath="//MSIGDB")
  db_name <- rvest::html_attr(x=db_data, name="NAME")
  db_ver <- rvest::html_attr(x=db_data, name="VERSION")
  db_date <- rvest::html_attr(x=db_data, name="BUILD_DATE")

  genesets <- rvest::xml_nodes(msig_result, "GENESET")
  row_names <- rvest::html_attr(x=genesets, name="STANDARD_NAME")
  column_names <- names(rvest::html_attrs(x=genesets[[1]]))
  all_data <- data.frame(row.names=row_names)
  for (i in 2:length(column_names)) {
    c_name <- column_names[i]
    c_data <- rvest::html_attr(x=genesets, name=c_name)
    all_data[[c_name]] <- c_data
  }

  sig_found_idx <- rownames(all_data) %in% names(sig_data)
  ret_data <- all_data[sig_found_idx, ]
  retlist <- list(
    "all_data" = all_data,
    "sub_data" = ret_data)
  if (!is.null(gsva_result)) {
    fData(gsva_result) <- ret_data
    retlist[["gsva_result"]] <- gsva_result
  }
  return(retlist)
}

## EOF
