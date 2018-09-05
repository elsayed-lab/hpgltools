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
#' @param current_id  Where did the IDs of the genes come from?
#' @param required_id  gsva (I assume) always requires ENTREZ IDs, but just in
#'   case this is a parameter.
#' @param orgdb  What is the data source for the rownames()?
#' @return  Something from GSVA::gsva()!
#' @export
simple_gsva <- function(expt, datasets="c2BroadSets", data_pkg="GSVAdata",
                        current_id="ENSEMBL", required_id="ENTREZID",
                        orgdb="org.Hs.eg.db") {

  ## Make sure some data is loaded.  Assume the c2BroadSets from GSVAdata.
  sig_data <- NULL
  if (exists(datasets)) {
    sig_data <- datasets
    if (class(sig_data)[[1]] != "GeneSetCollection") {
      stop("The data must be a GeneSetCollection.")
    }
  } else {
    tt <- sm(library(data_pkg, character.only=TRUE))
    lst <- list("list"=datasets, "package"=data_pkg)
    test <- do.call("data", as.list(datasets, lst))
    sig_data <- get0(datasets)
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

  gsva_result <- GSVA::gsva(eset, sig_data, verbose=TRUE)
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

simple_xcell <- function(expt) {
  xcell_annot <- load_biomart_annotations()
  xref <- xcell_annot[["annotation"]][, c("ensembl_gene_id", "hgnc_symbol")]
  xcell_input <- sm(exprs(normalize_expt(expt, norm="quant", convert="cpm", transform="log2")))
  xcell_input <- merge(xcell_input, xref, by.x="row.names", by.y="ensembl_gene_id")
  rownames(xcell_input) <- make.names(xcell_input[["hgnc_symbol"]], unique=TRUE)
  xcell_input[["Row.names"]] <- NULL
  xcell_input[["hgnc_symbol"]] <- NULL

  library(xCell)
  data("xCell.data", package="xCell")
  xcell_result <- xCell::xCellAnalysis(xcell_input)

  jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  ht <- heatmap.3(xcell_result, trace="none", col=jet_colors)
  ht_plot <- grDevices::recordPlot()

  retlist <- list(
    "xcell_result" = xcell_result,
    "heatmap" = ht_plt)
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
