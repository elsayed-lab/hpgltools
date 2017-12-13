#' Grab gene length/width/size from an annotation database.
#'
#' This function tries to gather an appropriate gene length column from
#' whatever annotation data source is provided.
#'
#' @param type What type of annotation data are we using?
#' @param gene_type Annotation type to use (3rd column of a gff file).
#' @param type_column Type identifier (10th column of a gff file).
#' @param key  What column has ID information?
#' @param length_names  Provide some column names which give gene length information?
#' @param ... Extra arguments likely for load_annotations()
#' @return Data frame of gene IDs and widths.
#' @seealso \pkg{rtracklayer}
#'  \code{\link{load_gff_annotations}}
#' @examples
#' \dontrun{
#'  tt = get_genesizes(gff="pa14.gff")
#'  head(tt)
#' ##           ID width
#' ## 1   YAL069W   312
#' ## 2   YAL069W   315
#' ## 3   YAL069W     3
#' ## 4 YAL068W-A   252
#' ## 5 YAL068W-A   255
#' ## 6 YAL068W-A     3
#' }
#' @export
get_genesizes <- function(type="gff", gene_type="gene", type_column="type",
                          key=NULL, length_names=NULL, ...) {
  annot <- load_annotations(type=type, ...)
  ## Pull out the rows of interest.
  if (!is.null(gene_type)) {
    desired_rows <- annot[, type_column] == gene_type
    if (sum(desired_rows) == 0) {
      message(paste0("There appear to be no genes of type ",
                     gene_type, " in the ", type_column, " column, taking them all."))
    } else {
      message(paste0("Taking only the ", sum(desired_rows), " ", gene_type, " rows."))
      annot <- annot[desired_rows, ]
    }
  }

  ## Try to get a coherent set of rownames and/or IDs.
  row_names <- NULL
  if (is.null(key)) {
    row_names <- rownames(annot)
  } else {
    row_names <- make.names(annot[[key]], unique=TRUE)
  }
  ret <- data.frame(row.names=row_names)
  if (!is.null(annot[["ID"]])) {
    ret[["ID"]] <- annot[["ID"]]
  }

  ## Now try to a column with the information of interest.
  if (is.null(length_names)) {
    length_names <- c("width", "length", "gene_size")
  }
  for (ln in length_names) {
    if (!is.null(annot[[ln]])) {
      ret[["gene_size"]] <- annot[[ln]]
    }
  }

  ## Send back the result.
  return(ret)
}

#' Use one of the load_*_annotations() functions to gather annotation data.
#'
#' We should be able to have an agnostic annotation loader which can take some
#' standard arguments and figure out where to gather data on its own.
#'
#' @param type  Explicitly state the type of annotation data to load.  If not
#'   provided, try to figure it out automagically.
#' @param ...  Arguments passed to the other load_*_annotations().
#' @return  Some annotations, hopefully.
#' @export
load_annotations <- function(type=NULL, ...) {
  annotations <- NULL
  ## FIXME: Add some logic here to figure out what search to perform.
  switchret <- switch(
    type,
    "biomart" = {
      annotations <- load_biomart_annotations(...)
    },
    "gff" = {
      annotations <- load_gff_annotations(...)
    },
    "genbank" = {
      annotations <- load_genbank_annotations(...)
    },
    "kegg" = {
      annotations <- load_kegg_annotations(...)
    },
    "microbesonline" = {
      annotations <- load_microbesonline_annotations(...)
    },
    "trinotate" = {
      annotations <- load_trinotate_annotations(...)
    },
    "uniprot" = {
      annotations <- load_uniprot_annotations(...)
    },
    {
      message(paste0("Not sure what type you chose, defaulting to biomart."))
      annotations <- load_biomart_annotations(...)
    })
  return(annotations)
}

## EOF
