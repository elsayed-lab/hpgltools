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
    "querymany" = {
      annotations <- load_querymany_annotations(...)
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
