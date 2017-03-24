#' Make an html version of an MA plot: M(log ratio of conditions) / A(mean average).
#'
#' A fun snippet from wikipedia: "In many microarray gene expression experiments, an underlying
#' assumption is that most of the genes would not see any change in their expression therefore the
#' majority of the points on the y-axis (M) would be located at 0, since Log(1) is 0. If this is not
#' the case, then a normalization method such as LOESS should be applied to the data before
#' statistical analysis. If the median line is not straight, the data should be normalized.
#'
#' @param df Data frame of counts which have been normalized counts by sample-type, which is to say the
#'  output from voom/voomMod/hpgl_voom().
#' @param tooltip_data Df of tooltip information (gene names, etc).
#' @param filename Filename to write a fancy html graph.
#' @param base_url String with a basename used for generating URLs for clicking dots on the graph.
#' @param ... more options are more options!
#' @return NULL, but along the way an html file is generated which contains a googleVis MA plot.
#'  See plot_de_ma() for details.
#' @seealso \pkg{googleVis}
#'  \code{\link{plot_ma_de}}
#' @examples
#' \dontrun{
#'  plot_gvis_ma(df, filename="html/fun_ma_plot.html",
#'               base_url="http://yeastgenome.org/accession?")
#' }
#' @export
plot_gvis_ma <- function(df, tooltip_data=NULL, filename="html/gvis_ma_plot.html", base_url="", ...) {
    gvis_chartid <- gsub("\\.html$", "", basename(filename))
    gvis_df <- data.frame("AvgExp" = df[["avg"]],
                          "LogFC" = df[["logfc"]],
                          "AdjPVal" = df[["pval"]])
    ## gvis_sig = subset(gvis_df, AdjPVal <= 0.05)
    gvis_sig <- gvis_df[ which(gvis_df[["AdjPVal"]] <= 0.05), ]
    gvis_sig <- gvis_sig[, c(1, 2)]
    gvis_sig <- merge(gvis_sig, tooltip_data, by="row.names")
    rownames(gvis_sig) <- gvis_sig[["Row.names"]]
    gvis_sig <- gvis_sig[, -1]
    colnames(gvis_sig) <- c("AvgExp", "Significant", "sig.tooltip")
    ## gvis_nonsig = subset(gvis_df, AdjPVal > 0.05)
    gvis_nonsig <- gvis_df[ which(gvis_df[["AdjPVal"]] > 0.05), ]
    gvis_nonsig <- gvis_nonsig[, c(1, 2)]
    gvis_nonsig <- merge(gvis_nonsig, tooltip_data, by="row.names")
    rownames(gvis_nonsig) <- gvis_nonsig[["Row.names"]]
    gvis_nonsig <- gvis_nonsig[, -1]
    colnames(gvis_nonsig) <- c("AvgExp", "NonSignificant", "nsig.tooltip")
    gvis_final_df <- merge(gvis_df, gvis_nonsig, by="row.names", all.x=TRUE)
    gvis_final_df <- merge(gvis_final_df, gvis_sig, by.x="Row.names", by.y="row.names", all.x=TRUE)
    rownames(gvis_final_df) <- gvis_final_df[["Row.names"]]
    gvis_final_df <- gvis_final_df[, c(2, 6, 7, 9, 10)]
    colnames(gvis_final_df) <- c("AvgExp", "NonSignificant", "nsig.tooltip",
                                 "Significant", "sig.tooltip")
    ma_json_ids <- rjson::toJSON(row.names(gvis_final_df))
    ma_jscode <- paste("
 var IDs=", ma_json_ids, "
 var sel = chart.getSelection();
 var row = sel[0].row;
 var text = IDs[row];
 window.open('", base_url, "' + text, '_blank');
")
    gvis_options <- list(pointSize=2, height=800, width=800,
                         tooltip="[{isHtml: true},{trigger:'selection'}]",
                         hAxes="[{title:'AvgExp'}]",
                         vAxes="[{title:'LogFC'}]",
                         title="MA Plot!",
                         gvis.listener.jscode=ma_jscode,
                         axisTitlesPosition="out")
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_final_df),
                                                          chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

#' Make an html version of an volcano plot.
#'
#' Volcano plots provide some visual clues regarding the success of a given contrast.  For our data,
#' it has the -log10(pvalue) on the y-axis and fold-change on the x.  Here is a neat snippet from
#' wikipedia describing them generally: "The concept of volcano plot can be generalized to other
#' applications, where the x-axis is related to a measure of the strength of a statistical signal,
#' and y-axis is related to a measure of the statistical significance of the signal."
#'
#' @param toptable_data Df of toptable() data.
#' @param fc_cutoff Fold change cutoff.
#' @param p_cutoff Maximum p value to allow.
#' @param filename Filename to write a fancy html graph.
#' @param tooltip_data Df of tooltip information.
#' @param base_url String with a basename used for generating URLs for clicking dots on the graph.
#' @param ... more options
#' @return NULL, but along the way an html file is generated which contains a googleVis volcano
#'  plot.
#' @seealso \pkg{googleVis}
#'  \code{\link{plot_volcano}}
#' @examples
#' \dontrun{
#'  plot_gvis_volcano(voomed_data, toptable_data, filename="html/fun_ma_plot.html",
#'                     base_url="http://yeastgenome.org/accession?")
#' }
#' @export
plot_gvis_volcano <- function(toptable_data, fc_cutoff=0.8, p_cutoff=0.05,
                              tooltip_data=NULL, filename="html/gvis_vol_plot.html",
                              base_url="", ...) {
    gvis_raw_df <- toptable_data[, c("logFC", "modified_p", "P.Value")]
    if (!is.null(tooltip_data)) {
        gvis_raw_df <- merge(gvis_raw_df, tooltip_data, by="row.names")
    }
    ## gvis_sig = subset(gvis_raw_df, P.Value <= p_cutoff)
    gvis_sig <- gvis_raw_df[ which(gvis_raw_df$P.Value <= p_cutoff), ]
    ## gvis_nsig = subset(gvis_raw_df, P.Value > p_cutoff)
    gvis_nsig <- gvis_raw_df[ which(gvis_raw_df$P.Value > p_cutoff), ]
    colnames(gvis_sig) <- c("Row.names", "logFCsig", "sig_modp", "sig_p", "sig.tooltip")
    colnames(gvis_nsig) <- c("Row.names", "logFCnsig", "nsig_modp", "nsig_p", "nsig.tooltip")
    gvis_sig <- gvis_sig[, c("Row.names", "sig_modp", "sig.tooltip")]
    gvis_nsig <- gvis_nsig[, c("Row.names", "nsig_modp", "nsig.tooltip")]
    gvis_df <- merge(gvis_raw_df, gvis_nsig, by.x="Row.names", by.y="Row.names", all.x=TRUE)
    gvis_df <- merge(gvis_df, gvis_sig, by.x="Row.names", by.y="Row.names", all.x=TRUE)
    rownames(gvis_df) <- gvis_df$Row.names
    gvis_df <- gvis_df[-1]
    gvis_df <- gvis_df[, c("logFC", "nsig_modp", "nsig.tooltip", "sig_modp", "sig.tooltip")]
    colnames(gvis_df) <- c("logFC", "nsig_p", "nsig.tooltip", "sig_p", "sig.tooltip")
    gvis_chartid <- gsub("\\.html$", "", basename(filename))
    vol_json_ids <- rjson::toJSON(row.names(gvis_df))
    vol_jscode <- paste("
 var IDs=", vol_json_ids, "
 var sel = chart.getSelection();
 var row = sel[0].row;
 var text = IDs[row];
 window.open('", base_url, "' + text, '_blank');
")
    gvis_options <- list(pointSize=2, height=800, width=800,
                         tooltip="[{isHtml: true},{trigger:'selection'}]",
                         hAxes="[{title:'AvgExp'}]",
                         vAxes="[{title:'LogFC'}]",
                         series="[{color:'blue'}, {color:'red'}]",
                         title="Volcano Plot!",
                         gvis.listener.jscode=vol_jscode,
                         axisTitlesPosition="out")
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_df),
                                                          chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

#' Make an html version of a scatter plot.
#'
#' Given an arbitrary scatter plot, we can make it pretty and javascript-tacular using this function.
#'
#' @param df Df of two columns to compare.
#' @param filename Filename to write a fancy html graph.
#' @param tooltip_data Df of tooltip information for gvis graphs.
#' @param base_url Url to send click events which will be suffixed with the gene name.
#' @param trendline Add a trendline?
#' @return NULL, but along the way an html file is generated which contains a googleVis scatter
#'  plot.  See plot_scatter() for details.
#' @seealso \pkg{googleVis}
#'  \code{\link[googleVis]{gvisScatterChart}}
#' @examples
#' \dontrun{
#'  gvis_scatter(a_dataframe_twocolumns, filename="html/fun_scatter_plot.html",
#'               base_url="http://yeastgenome.org/accession?")
#' }
#' @export
plot_gvis_scatter <- function(df, tooltip_data=NULL, filename="html/gvis_scatter.html",
                              base_url="", trendline=NULL) {
    gvis_df <- df
    gvis_df <- merge(gvis_df, tooltip_data, by="row.names", all.x=TRUE)
    rownames(gvis_df) <- gvis_df$Row.names
    gvis_df <- gvis_df[-1]
    json_ids <- rjson::toJSON(row.names(gvis_df))
    gvis_chartid <- gsub("\\.html$", "", basename(filename))
    scatter_jscode <- paste0(" var IDs=", json_ids, "
 var sel = chart.getSelection();
 var row = sel[0].row;
 var text = IDs[row];
 window.open('", base_url, "' + text, '_blank');
")
    if (is.null(trendline)) {
        gvis_options <- list(pointSize=2, height=800, width=800,
                             tooltip="[{isHtml: true},{trigger:'selection'}]",
                             gvis.listener.jscode=scatter_jscode,
                             axisTitlesPosition="out")
    } else {
        trendline_string <- sprintf("{0: {type: '%s', visibleInLegend: 'true', color: 'green', lineWidth: 10, opacity: 0.5}}", trendline)
        gvis_options <- list(pointSize=2, height=800, width=800,
                             tooltip="[{isHtml: true},{trigger:'selection'}]",
                             gvis.listener.jscode=scatter_jscode,
                             trendlines=trendline_string,
                             axisTitlesPosition="out")
    }
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_df),
                                                          chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

## EOF
