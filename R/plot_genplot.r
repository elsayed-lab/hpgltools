#' Try plotting a chromosome (region)
#'
#' genoplotr is cool, I don't yet understand it though
#'
#' @param accession An accession to plot, this will download it.
#' @param start First segment to plot (doesn't quite work yet).
#' @param end Final segment to plot (doesn't quite work yet).
#' @param title Put a title on the resulting plot.
#' @return Hopefully a pretty plot of a genome
#' @seealso [genoPlotR]
#' @export
genoplot_chromosome <- function(accession = "AE009949", start = NULL, end = NULL,
                                title = "Genome plot") {
  tt <- download_gbk(accession)
  segments <- try(genoPlotR::read_dna_seg_from_file(glue("{accession}.gb")))
  if (is.null(start)) {
    start <- 1
  }
  if (is.null(end)) {
    end <- nrow(segments)
  }

  mid_pos <- genoPlotR::middle(segments)
  xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(start, end))
  tt <- sm(please_install("genoPlotR"))
  genoPlotR::plot_gene_map(dna_segs = list(segments),
                           main = title,
                           gene_type = "side_blocks",
                           dna_seg_scale = TRUE,
                           scale = FALSE)
}

## EOF
