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

  previous_file <- "02_estimation_macrophage.Rmd"
  ver <- "20170820"
  tmp <- sm(loadme(filename=paste0(gsub(pattern="\\.Rmd", replace="",
                                        x=previous_file), "-v", ver, ".rda.xz")))
  expt <- hs_expt

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
  if (length(annotation(eset)) == 0) {
    message("gsva requires the annotation field to be filled in.")
    annotation(eset) <- orgdb
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  if (current_id != required_id) {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    tt <- sm(library(orgdb, character.only=TRUE))
    old_ids <- rownames(exprs(eset))
    new_ids <- sm(select(x=get0(orgdb),
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
  result(gsva_result)
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

## EOF
