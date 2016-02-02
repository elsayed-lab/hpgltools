## Time-stamp: <Tue Feb  2 13:12:00 2016 Ashton Trey Belew (abelew@gmail.com)>

#' hpgl_gvis_ma_plot()  Make an html version of an MA plot.
#'
#' @param counts  df of linear-modelling, normalized counts by sample-type,
#' which is to say the output from voom/voomMod/hpgl_voom().
#' @param de_genes  df from toptable or its friends containing p-values.
#' @param adjpval_cutoff default=0.05  a cutoff defining significant from not.
#' @param filename default='html/gvis_ma_plot.html'  a filename to write a fancy html graph.
#' @param tooltip_data default=NULL  a df of tooltip information.
#' @param base_url default=''  a string with a basename used for generating URLs for clicking dots on the graph.
#'
#' @return NULL, but along the way an html file is generated which
#' contains a googleVis MA plot.  See hpgl_ma_plot() for details.
#'
#' @seealso \code{\link{hpgl_ma_plot}}
#'
#' @export
#' @examples
#' ## hpgl_gvis_ma_plot(voomed_data, toptable_data, filename="html/fun_ma_plot.html", base_url="http://yeastgenome.org/accession?")
hpgl_gvis_ma_plot <- function(counts, degenes, tooltip_data=NULL, filename="html/gvis_ma_plot.html", base_url="", ...) {
    gvis_chartid <- gsub("\\.html$", "", basename(filename))
    gvis_df <- data.frame(AvgExp=rowMeans(counts[rownames(degenes),]), LogFC=degenes$logFC, AdjPVal=degenes$adj.P.Val)
    ## gvis_sig = subset(gvis_df, AdjPVal <= 0.05)
    gvis_sig <- gvis_df[ which(gvis_df$AdjPVal <= 0.05), ]
    gvis_sig <- gvis_sig[,c(1,2)]
    gvis_sig <- merge(gvis_sig, tooltip_data, by="row.names")
    rownames(gvis_sig) <- gvis_sig$Row.names
    gvis_sig <- gvis_sig[-1]
    colnames(gvis_sig) <- c("AvgExp","Significant","sig.tooltip")
    ## gvis_nonsig = subset(gvis_df, AdjPVal > 0.05)
    gvis_nonsig <- gvis_df[ which(gvis_df$AdjPVal > 0.05), ]
    gvis_nonsig <- gvis_nonsig[,c(1,2)]
    gvis_nonsig <- merge(gvis_nonsig, tooltip_data, by="row.names")
    rownames(gvis_nonsig) <- gvis_nonsig$Row.names
    gvis_nonsig <- gvis_nonsig[-1]
    colnames(gvis_nonsig) <- c("AvgExp","NonSignificant","nsig.tooltip")
    gvis_final_df <- merge(gvis_df, gvis_nonsig, by="row.names", all.x=TRUE)
    gvis_final_df <- merge(gvis_final_df, gvis_sig, by.x="Row.names", by.y="row.names", all.x=TRUE)
    rownames(gvis_final_df) <- gvis_final_df$Row.names
    gvis_final_df <- gvis_final_df[,c(2,6,7,9,10)]
    colnames(gvis_final_df) <- c("AvgExp","NonSignificant","nsig.tooltip","Significant","sig.tooltip")
    ma_json_ids <- rjson::toJSON(row.names(gvis_final_df))
    ma_jscode <- paste("
 var IDs=", ma_json_ids, "
 var sel = chart.getSelection();
 var row = sel[0].row;
 var text = IDs[row];
 window.open('", base_url , "' + text, '_blank');
")
    gvis_options <- list(pointSize=2,height=800, width=800,
                         tooltip="[{isHtml: true},{trigger:'selection'}]",
                         hAxes="[{title:'AvgExp'}]",
                         vAxes="[{title:'LogFC'}]",
                         title="MA Plot!",
                         gvis.listener.jscode=ma_jscode,
                         axisTitlesPosition="out")
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_final_df), chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

#' hpgl_gvis_volcano_plot()  Make an html version of an volcano plot.
#'
#' @param toptable_data  df of toptable() data
#' @param fc_cutoff default=0.8  fold change cutoff.
#' @param p_cutoff default=0.05  maximum p value to allow.
#' @param filename default='html/gvis_vol_plot.html'  a filename to write a fancy html graph.
#' @param tooltip_data default=NULL  a df of tooltip information.
#' @param base_url default=''  a string with a basename used for generating URLs for clicking dots on the graph.
#'
#' @return NULL, but along the way an html file is generated which
#' contains a googleVis MA plot.  See hpgl_ma_plot() for details.
#'
#' @seealso \code{\link{hpgl_volcano_plot}}
#'
#' @export
#' @examples
#' ## hpgl_gvis_ma_plot(voomed_data, toptable_data, filename="html/fun_ma_plot.html", base_url="http://yeastgenome.org/accession?")
hpgl_gvis_volcano_plot <- function(toptable_data, fc_cutoff=0.8, p_cutoff=0.05,
                                  tooltip_data=NULL, filename="html/gvis_vol_plot.html", base_url="", ...) {
    gvis_raw_df <- toptable_data[,c("logFC", "modified_p", "P.Value")]
    if (!is.null(tooltip_data)) {
        gvis_raw_df <- merge(gvis_raw_df, tooltip_data, by="row.names")
    }
    ## gvis_sig = subset(gvis_raw_df, P.Value <= p_cutoff)
    gvis_sig <- gvis_raw_df[ which(gvis_raw_df$P.Value <= p_cutoff), ]
    ## gvis_nsig = subset(gvis_raw_df, P.Value > p_cutoff)
    gvis_nsig <- gvis_raw_df[ which(gvis_raw_df$P.Value > p_cutoff), ]
    colnames(gvis_sig) <- c("Row.names","logFCsig","sig_modp","sig_p","sig.tooltip")
    colnames(gvis_nsig) <- c("Row.names","logFCnsig","nsig_modp","nsig_p","nsig.tooltip")
    gvis_sig <- gvis_sig[,c("Row.names","sig_modp", "sig.tooltip")]
    gvis_nsig <- gvis_nsig[,c("Row.names","nsig_modp","nsig.tooltip")]
    gvis_df <- merge(gvis_raw_df, gvis_nsig, by.x="Row.names", by.y="Row.names", all.x=TRUE)
    gvis_df <- merge(gvis_df, gvis_sig, by.x="Row.names", by.y="Row.names", all.x=TRUE)
    rownames(gvis_df) <- gvis_df$Row.names
    gvis_df <- gvis_df[-1]
    gvis_df <- gvis_df[,c("logFC","nsig_modp","nsig.tooltip","sig_modp","sig.tooltip")]
    colnames(gvis_df) <- c("logFC","nsig_p","nsig.tooltip","sig_p","sig.tooltip")
    gvis_chartid <- gsub("\\.html$", "", basename(filename))
    vol_json_ids <- rjson::toJSON(row.names(gvis_df))
    vol_jscode <- paste("
 var IDs=", vol_json_ids, "
 var sel = chart.getSelection();
 var row = sel[0].row;
 var text = IDs[row];
 window.open('", base_url , "' + text, '_blank');
")
    gvis_options <- list(pointSize=2,height=800, width=800,
                         tooltip="[{isHtml: true},{trigger:'selection'}]",
                         hAxes="[{title:'AvgExp'}]",
                         vAxes="[{title:'LogFC'}]",
                         series="[{color:'blue'}, {color:'red'}]",
                         title="Volcano Plot!",
                         gvis.listener.jscode=vol_jscode,
                         axisTitlesPosition="out")
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_df), chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

#' hpgl_gvis_scatter()  Make an html version of a scatter plot.
#'
#' @param df  df of two columns to compare
#' @param filename default='html/gvis_scatter.html'  a filename to write a fancy html graph.
#' @param tooltip_data default=NULL  a df of tooltip information for gvis graphs.
#' @param base_url default=''  a url to send click events which will be suffixed with the gene name
#' @param trendline default=NULL  add a trendline?
#'
#' @return NULL, but along the way an html file is generated which
#' contains a googleVis scatter plot.  See hpgl_scatter_plot() for details.
#'
#' @seealso \code{\link{gvisScatterChart}}
#'
#' @export
#' @examples
#' ## hpgl_gvis_scatter(a_dataframe_twocolumns, filename="html/fun_scatter_plot.html", base_url="http://yeastgenome.org/accession?")
hpgl_gvis_scatter <- function(df, tooltip_data=NULL, filename="html/gvis_scatter.html", base_url="", trendline=NULL) {
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
 window.open('", base_url , "' + text, '_blank');
")
    if (is.null(trendline)) {
        gvis_options <- list(pointSize=2,height=800, width=800,
                             tooltip="[{isHtml: true},{trigger:'selection'}]",
                             gvis.listener.jscode=scatter_jscode,
                             axisTitlesPosition="out")
    } else {
        trendline_string <- sprintf("{0: {type: '%s', visibleInLegend: 'true', color: 'green', lineWidth: 10, opacity: 0.5}}", trendline)
        gvis_options <- list(pointSize=2,height=800, width=800,
                             tooltip="[{isHtml: true},{trigger:'selection'}]",
                             gvis.listener.jscode=scatter_jscode,
                             trendlines=trendline_string,
                             axisTitlesPosition="out")
    }
    hpgl_gvis_scatterchart <- googleVis::gvisScatterChart(as.data.frame(gvis_df), chartid=gvis_chartid, options=gvis_options)
    print(hpgl_gvis_scatterchart, file=filename)
}

## EOF
