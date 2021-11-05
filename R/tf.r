#' Invoke dorothea in an attempt to hunt down cool TFs.
#'
#' dorothea: https://github.com/saezlab/dorothea appears to provide
#' experimentally verified mappings from genes->transcription
#' factors as well as a set of functions which allow one to pass it an
#' expressionset/matrix of counts(log scale) and get back scores by
#' tf.  This function is an attempt to smooth it out and prod it for
#' usability.
#'
#' @param expt Expressionset
#' @param gene_column Column in fData containing the gene IDs.
#' @param hgnc_column fData column containing the HGNC symbols as used
#'  by dorothea/viper
#' @param transform Explicitly set the scale to log2 (TODO: improve
#'  this)
#' @param conf Vector of confidence scores to filter the data.
#' @param dorothea_options Optional configuration list.
#' @param lfc Fold-change cutoff.
#' @param p P-value cutoff.
#' @param species Either human or mouse.
#' @return list containing some information from dorothea and limma.
#' @export
simple_dorothea <- function(expt, gene_column="ensembl_gene_id",
                            hgnc_column="hgnc_symbol", transform="log2",
                            conf=c("A", "B", "C"), dorothea_options=NULL,
                            lfc = 1, p = 0.05, species="hsapiens") {
  ## dorothea requires hgnc_symbols, my expressionsets are ensembl gene IDs in general...
  ## I need to read further to see if these methods expect normalized input data.

  if (is.null(dorothea_options)) {
    dorothea_options <- list("method" = "scale",
                             "minsize" = 4,
                             "eset.filter" = FALSE,
                             "cores" = 1,
                             "verbose" = FALSE)
  }

  ## Their example data is all on the log scale, so let us ensure that here.
  mtrx <- exprs(normalize_expt(expt, transform=transform))

  ## Now recast the rownames as hgnc IDs.  Those IDs are weird and
  ## redundant, so keep that in mind.
  hgnc <- fData(expt)[, c(gene_column, hgnc_column)]
  tmp_mtrx <- merge(mtrx, hgnc, by.x="row.names", by.y=gene_column)
  rownames(tmp_mtrx) <- make.names(tmp_mtrx[[hgnc_column]],
                                   unique=TRUE)
  tmp_mtrx[["Row.names"]] <- NULL
  tmp_mtrx[["hgnc_symbol"]] <- NULL

  ## They have two datasets, one for humans and one for mice.
  dorothea_hs <- dorothea_mm <- dorothea_df <- NULL
  if (species == "hsapiens") {
    data(dorothea_hs, package="dorothea")
    dorothea_df <- dorothea_hs
  } else if (species == "mmusculus") {
    data(dorothea_mm, package="dorothea")
    dorothea_df <- dorothea_mm
  } else {
    stop("I do not know what to do with the species: ", species)
  }

  ## Their confidence scores are 'A', 'B', 'C', 'D' and maybe more.
  confidence <- NULL
  regulons = dorothea_df %>%
    filter(confidence %in% conf)

  ## viper is what does the real work as far as I can tell.
  ## "Virtual Inference of Protein-activity by Enriched Regulon analysis."
  ## Though the dorothea documentation appears (to my eyes) to send it
  ## data containing logFC/adjP values, I do not see anywhere in the
  ## viper documentation where it is able to interpret those columns
  ## in any useful way.  My reading of the documentation (I haven't
  ## started reading the source code yet) suggests to me that it is
  ## similar to the enrichment analyses performed by gostats etc with
  ## the input mapping categories set to experimentally derived
  ## mappings from a few experiments mostly from Bcells.  I have a
  ## strong feeling that we could get similar/identical results by
  ## pulling the mapping data out of it and passing it to gsva.
  dorothea_exprs <- dorothea::run_viper(tmp_mtrx, regulons,
                                        options=dorothea_options)

  ## I think what I need to do, if I wish to play with the annotations,
  ## is to reformat the entire_database so that there is one row for each
  ## transcription factor ID.
  tmp_fdata <- data.frame(row.names=rownames(dorothea_exprs))
  tmp_fdata <- merge(tmp_fdata, fData(expt), by.x = "row.names",
                     by.y = hgnc_column)
  fd_duplicated <- duplicated(tmp_fdata[["Row.names"]])
  tmp_fdata <- tmp_fdata[!fd_duplicated, ]
  rownames(tmp_fdata) <- tmp_fdata[["Row.names"]]
  tmp_fdata[["Row.names"]] <- NULL

  ## This is very much the same thing as performed by viper, they have
  ## a series of S4 dispatchers which recreate expressionsets/matrices
  ## depending on the input data.
  new_expt <- create_expt(gene_info=tmp_fdata,
                          count_dataframe=dorothea_exprs,
                          metadata=pData(expt))

  ## Since this looks just like gsva to me, I will treat it as such.
  delta_tf <- limma_pairwise(new_expt, which_voom="none")
  sig_tf_lst <- list()
  for (tbl in 1:length(delta_tf[["all_tables"]])) {
    tbl_name <- names(delta_tf[["all_tables"]])[tbl]
    tbl_df <- delta_tf[["all_tables"]][[tbl_name]]
    sig_tf_lst[[tbl_name]] <- get_sig_genes(tbl_df, lfc = lfc, p = p)
  }

  retlist <- list(
      "limma_result" = delta_tf,
      "dorothea_sig" = sig_tf_lst,
      "dorothea_result" = new_expt)
  return(retlist)
}
