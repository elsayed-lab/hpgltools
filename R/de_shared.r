# Test for infected/control/beads -- a placebo effect?
## The goal is therefore to find responses different than beads
## The null hypothesis is (H0): (infected == uninfected) || (infected == beads)
## The alt hypothesis is (HA): (infected != uninfected) && (infected != beads)
disjunct_tab <- function(contrast_fit, coef1, coef2, ...) {
    stat <- BiocGenerics::pmin(abs(contrast_fit[, coef1]), abs(contrast_fit[, coef2]))
    pval <- BiocGenerics::pmax(contrast_fit$p.val[, coef1], contrast_fit$p.val[, coef2])
}
## An F-test only does inf==uninf && inf==bead
## So the solution is to separately perform the two subtests and subset for the set of genes for which both are true.
## However, if you do that, the f-statistics are a little screwey, but there are a few ways to handle it:
## Perform the two separate tests and perform the following combination of the stat and p-value:
##    stat = min(|inf-uninf|, |inf-bead|)  (logFC)
##    ^^pval^^ = max(pval(inf-uninf), pval(inf-beads))
##    adj.pval = p.adjust(^^pval^^, method='BH')
## ReportingTools hwriter

#' Perform limma, DESeq2, EdgeR pairwise analyses.
#'
#' This takes an expt object, collects the set of all possible pairwise comparisons, sets up
#' experimental models appropriate for the differential expression analyses, and performs them.
#'
#' @param input Dataframe/vector or expt class containing count tables, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the model?  This is likely always true.
#' @param model_batch Include batch in the model?  This may be true/false/"sva" or other methods supported by get_model_adjust().
#' @param model_intercept Use an intercept model instead of cell means?
#' @param extra_contrasts Optional extra contrasts beyone the pairwise comparisons.  This can be
#'     pretty neat, lets say one has conditions A,B,C,D,E and wants to do (C/B)/A and (E/D)/A or
#'     (E/D)/(C/B) then use this with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla =
#'     (E-D)-A, de_vs_cb = (E-D)-(C-B)".
#' @param alt_model Alternate model to use rather than just condition/batch.
#' @param libsize Library size of the original data to help voom().
#' @param annot_df Annotations to add to the result tables.
#' @param ... Picks up extra arguments into arglist, currently only passed to write_limma().
#' @return A list of limma, deseq, edger results.
#' @examples
#' \dontrun{
#'  finished_comparison = eBayes(limma_output)
#'  data_list = write_limma(finished_comparison, workbook="excel/limma_output.xls")
#' }
#' @export
all_pairwise <- function(input, conditions=NULL, batches=NULL, model_cond=TRUE,
                         model_batch=TRUE, model_intercept=TRUE, extra_contrasts=NULL,
                         alt_model=NULL, libsize=NULL, annot_df=NULL, ...) {
    arglist <- list(...)
    surrogates <- 1
    if (!is.null(arglist[["surrogates"]])) {
        surrogates <- arglist[["surrogates"]]
    }
    if (is.null(model_cond)) {
        model_cond <- TRUE
    }
    if (is.null(model_batch)) {
        model_batch <- FALSE
    }
    if (is.null(model_intercept)) {
        model_intercept <- FALSE
    }
    if (class(model_batch) == 'character') {
        params <- get_model_adjust(input, estimate_type=model_batch, surrogates=surrogates)
        model_batch <- params[["model_adjust"]]
    }

    ##limma_result <- limma_pairwise(input, conditions=conditions, batches=batches, model_cond=model_cond, model_batch=model_batch, model_intercept=model_intercept, extra_contrasts=extra_contrasts, alt_model=alt_model, libsize=libsize, annot_df=annot_df)
    limma_result <- limma_pairwise(input, conditions=conditions, batches=batches,
                                   model_cond=model_cond, model_batch=model_batch,
                                   model_intercept=model_intercept, extra_contrasts=extra_contrasts,
                                   alt_model=alt_model, libsize=libsize, annot_df=annot_df, ...)
    deseq_result <- deseq2_pairwise(input, conditions=conditions, batches=batches,
                                    model_cond=model_cond, model_batch=model_batch,
                                    model_intercept=model_intercept, extra_contrasts=extra_contrasts,
                                    alt_model=alt_model, libsize=libsize, annot_df=annot_df, ...)
    edger_result <- edger_pairwise(input, conditions=conditions, batches=batches,
                                   model_cond=model_cond, model_batch=model_batch,
                                   model_intercept=model_intercept, extra_contrasts=extra_contrasts,
                                   alt_model=alt_model, libsize=libsize, annot_df=annot_df, ...)
    basic_result <- basic_pairwise(input, conditions, ...)

    result_comparison <- compare_tables(limma=limma_result, deseq=deseq_result,
                                        edger=edger_result, basic=basic_result,
                                        annot_df=annot_df, ...)
    ret <- list(
        "input" = input,
        "limma" = limma_result,
        "deseq" = deseq_result,
        "edger" = edger_result,
        "basic" = basic_result,
        "comparison" = result_comparison)
    return(ret)
}

#' See how similar are results from limma/deseq/edger.
#'
#' limma, DEseq2, and EdgeR all make somewhat different assumptions.
#' and choices about what makes a meaningful set of differentially.
#' expressed genes.  This seeks to provide a quick and dirty metric
#' describing the degree to which they (dis)agree.
#'
#' @param limma Data from limma_pairwise().
#' @param deseq Data from deseq2_pairwise().
#' @param edger Data from edger_pairwise().
#' @param basic Data from basic_pairwise().
#' @param include_basic include the basic data?
#' @param annot_df Include annotation data?
#' @param ... More options!
#' @return Heatmap showing how similar they are along with some
#'     correlations betwee the three players.
#' @seealso \code{\link{limma_pairwise}} \code{\link{edger_pairwise}} \code{\link{deseq2_pairwise}}
#' @examples
#' \dontrun{
#'  l = limma_pairwise(expt)
#'  d = deseq_pairwise(expt)
#'  e = edger_pairwise(expt)
#'  fun = compare_tables(limma=l, deseq=d, edger=e)
#' }
#' @export
compare_tables <- function(limma=NULL, deseq=NULL, edger=NULL, basic=NULL,
                           include_basic=TRUE, annot_df=NULL, ...) {
    arglist <- list(...)
    ## Fill each column/row of these with the correlation between tools for one contrast performed
    if (class(limma) == "list") {
        ## Then this was fed the raw output from limma_pairwise,
        ## lets assume the same is true for deseq/edger too and pull out the result tables.
        limma <- limma[["all_tables"]]
        deseq <- deseq[["all_tables"]]
        edger <- edger[["all_tables"]]
        basic <- basic[["all_tables"]]
    }
    len <- length(names(deseq))
    limma_vs_edger <- list()
    limma_vs_deseq <- list()
    limma_vs_basic <- list()
    edger_vs_deseq <- list()
    edger_vs_basic <- list()
    deseq_vs_basic <- list()
    limma_vs_edger_scatter <- list()
    limma_vs_deseq_scatter <- list()
    limma_vs_basic_scatter <- list()
    edger_vs_deseq_scatter <- list()
    edger_vs_basic_scatter <- list()
    deseq_vs_basic_scatter <- list()
    cc <- 0
    for (comp in names(deseq)) {
        ## assume all three have the same names() -- note that limma has more than the other two though
        cc <- cc + 1
        message(paste0("Comparing analyses ", cc, "/", len, ": ", comp))
        l <- data.frame(limma[[comp]])
        e <- data.frame(edger[[comp]])
        d <- data.frame(deseq[[comp]])
        b <- data.frame(basic[[comp]])
        le <- merge(l, e, by.x="row.names", by.y="row.names")
        le <- le[,c("logFC.x","logFC.y")]
        colnames(le) <- c("limma logFC", "edgeR logFC")
        lec <- stats::cor.test(x=le[, 1], y=le[, 2])[["estimate"]]
        les <- plot_scatter(le) + ggplot2::labs(title=paste0(comp, ": limma vs. edgeR.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        ld <- merge(l, d, by.x="row.names", by.y="row.names")
        ld <- ld[, c("logFC.x","logFC.y")]
        colnames(ld) <- c("limma logFC", "DESeq2 logFC")
        ldc <- stats::cor.test(ld[,1], ld[,2])[["estimate"]]
        lds <- plot_scatter(ld) + ggplot2::labs(title=paste0(comp, ": limma vs. DESeq2.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        lb <- merge(l, b, by.x="row.names", by.y="row.names")
        lb <- lb[, c("logFC.x","logFC.y")]
        colnames(lb) <- c("limma logFC", "basic logFC")
        lbc <- stats::cor.test(lb[,1], lb[,2])[["estimate"]]
        lbs <- plot_scatter(lb) + ggplot2::labs(title=paste0(comp, ": limma vs. basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        ed <- merge(e, d, by.x="row.names", by.y="row.names")
        ed <- ed[, c("logFC.x","logFC.y")]
        colnames(ed) <- c("edgeR logFC", "DESeq2 logFC")
        edc <- stats::cor.test(ed[,1], ed[,2])[["estimate"]]
        eds <- plot_scatter(ed) + ggplot2::labs(title=paste0(comp, ": edgeR vs. DESeq2.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        eb <- merge(e, b, by.x="row.names", by.y="row.names")
        eb <- eb[, c("logFC.x","logFC.y")]
        colnames(eb) <- c("edgeR logFC", "basic logFC")
        ebc <- stats::cor.test(eb[,1], eb[,2])[["estimate"]]
        ebs <- plot_scatter(eb) + ggplot2::labs(title=paste0(comp, ": edgeR vs. basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        db <- merge(d, b, by.x="row.names", by.y="row.names")
        db <- db[, c("logFC.x","logFC.y")]
        colnames(db) <- c("DESeq2 logFC", "basic logFC")
        dbc <- stats::cor.test(db[,1], db[,2])[["estimate"]]
        dbs <- plot_scatter(db) +
            ggplot2::labs(title=paste0(comp, ": DESeq2 vs basic.")) +
            ggplot2::geom_abline(intercept=0.0, slope=1.0, colour="blue")
        limma_vs_edger[[comp]] <- lec
        limma_vs_deseq[[comp]] <- ldc
        edger_vs_deseq[[comp]] <- edc
        limma_vs_basic[[comp]] <- lbc
        edger_vs_basic[[comp]] <- ebc
        deseq_vs_basic[[comp]] <- dbc
        limma_vs_edger_scatter[[comp]] <- les
        limma_vs_deseq_scatter[[comp]] <- lds
        edger_vs_deseq_scatter[[comp]] <- eds
        limma_vs_basic_scatter[[comp]] <- lbs
        edger_vs_basic_scatter[[comp]] <- ebs
        deseq_vs_basic_scatter[[comp]] <- dbs
    } ## End loop
    names(limma_vs_edger) <- names(deseq)
    names(limma_vs_deseq) <- names(deseq)
    names(edger_vs_deseq) <- names(deseq)
    names(limma_vs_basic) <- names(deseq)
    names(edger_vs_basic) <- names(deseq)
    names(deseq_vs_basic) <- names(deseq)
    names(limma_vs_edger_scatter) <- names(deseq)
    names(limma_vs_deseq_scatter) <- names(deseq)
    names(edger_vs_deseq_scatter) <- names(deseq)
    names(limma_vs_basic_scatter) <- names(deseq)
    names(edger_vs_basic_scatter) <- names(deseq)
    names(deseq_vs_basic_scatter) <- names(deseq)

    comparison_df <- rbind(as.numeric(limma_vs_edger), as.numeric(limma_vs_deseq))
    comparison_df <- rbind(comparison_df, as.numeric(edger_vs_deseq))
    if (isTRUE(include_basic)) {
        comparison_df <- rbind(comparison_df, as.numeric(limma_vs_basic))
        comparison_df <- rbind(comparison_df, as.numeric(edger_vs_basic))
        comparison_df <- rbind(comparison_df, as.numeric(deseq_vs_basic))
        rownames(comparison_df) <- c("le", "ld", "ed", "lb", "eb", "db")
    } else {
        rownames(comparison_df) <- c("le", "ld", "ed")
    }
    comparison_df <- as.matrix(comparison_df)
    colnames(comparison_df) <- names(deseq)
    heat_colors <- grDevices::colorRampPalette(c("white","black"))
    comparison_heatmap <- try(heatmap.3(comparison_df, scale="none", trace="none",
                                        linewidth=0.5, keysize=2, margins=c(8,8),
                                        col=heat_colors, dendrogram="none", Rowv=FALSE,
                                        Colv=FALSE, main="Compare DE tools"), silent=TRUE)
    heat <- NULL
    if (class(comparison_heatmap) != 'try-error') {
        heat <- recordPlot()
    }
    ret <- list(
        "limma_vs_edger" = limma_vs_edger,
        "limma_vs_deseq" = limma_vs_deseq,
        "limma_vs_basic" = limma_vs_basic,
        "edger_vs_deseq" = edger_vs_deseq,
        "edger_vs_basic" = edger_vs_basic,
        "deseq_vs_basic" = deseq_vs_basic,
        "limma_vs_edger_scatter" = limma_vs_edger_scatter,
        "limma_vs_deseq_scatter" = limma_vs_deseq_scatter,
        "limma_vs_basic_scatter" = limma_vs_basic_scatter,
        "edger_vs_deseq_scatter" = edger_vs_deseq_scatter,
        "edger_vs_basic_scatter" = edger_vs_basic_scatter,
        "deseq_vs_basic_scatter" = deseq_vs_basic_scatter,
        "comp" = comparison_df,
        "heat" = heat)
    return(ret)
}

#' Combine portions of deseq/limma/edger table output.
#'
#' This hopefully makes it easy to compare the outputs from
#' limma/DESeq2/EdgeR on a table-by-table basis.
#'
#' @param all_pairwise_result Output from all_pairwise().
#' @param extra_annot Add some annotation information?
#' @param excel Filename for the excel workbook, or null if not printed.
#' @param excel_title Title for the excel sheet(s).  If it has the
#'     string 'YYY', that will be replaced by the contrast name.
#' @param excel_sheet Name the excel sheet.
#' @param csv  On some computers (Edson!) printing to excel runs the machine oom for big data sets.
#' @param keepers List of reformatted table names to explicitly keep
#'     certain contrasts in specific orders and orientations.
#' @param include_basic Include my stupid basic logFC tables?
#' @param add_plots Add plots to the end of the sheets with expression values?
#' @param compare_plots  In an attempt to save memory when printing to excel, make it possible to
#'     exclude comparison plots in the summary sheet.
#' @param plot_dim Number of inches squared for the plot if added.
#' @return Table combining limma/edger/deseq outputs.
#' @seealso \code{\link{all_pairwise}}
#' @examples
#' \dontrun{
#' pretty = combine_de_tables(big_result, table='t12_vs_t0')
#' }
#' @export
combine_de_tables <- function(all_pairwise_result, extra_annot=NULL, csv=NULL,
                              excel=NULL, excel_title="Table SXXX: Combined Differential Expression of YYY",
                              excel_sheet="combined_DE", keepers="all",
                              include_basic=TRUE, add_plots=TRUE, plot_dim=6, compare_plots=TRUE) {
    ## The ontology_shared function which creates multiple sheets works a bit differently
    ## It creates all the tables, then does a createWorkbook()
    ## Does a createWorkbook() / addWorksheet()
    ## Then a writeData() / writeDataTable() / print(plot) / insertPlot() / saveWorkbook()
    ## Lets try that here.
    retlist <- NULL

    limma <- all_pairwise_result[["limma"]]
    deseq <- all_pairwise_result[["deseq"]]
    edger <- all_pairwise_result[["edger"]]
    basic <- all_pairwise_result[["basic"]]

    csv_basename <- NULL
    if (!is.null(csv)) {
        if (is.null(excel)) {
            csv_basename <- "excel/csv_export"
        } else {
            csv_basename <- excel
            csv_basename <- gsub(pattern="\\.xlsx", replacement="", x=csv_basename)
        }
    }

    wb <- NULL
    if (!is.null(excel)) {
        excel_dir <- dirname(excel)
        if (!file.exists(excel_dir)) {
            dir.create(excel_dir, recursive=TRUE)
        }
        if (file.exists(excel)) {
            message(paste0("Deleting the file ", excel, " before writing the tables."))
            file.remove(excel)
        }
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    }

    message("Writing a legend of columns.")
    legend <- data.frame(rbind(
        c("The first ~3-10 columns:", "are annotation data taken from whatever annotation source we used."),
        c("Next 6 columns", "The logFC and p-values reported by limma/edger/deseq2."),
        c("limma_logfc", "The log2 fold change reported by limma."),
        c("deseq_logfc", "The log2 fold change reported by DESeq2."),
        c("edger_logfc", "The log2 fold change reported by edgeR."),
        c("limma_adjp", "The adjusted-p value reported by limma."),
        c("deseq_adjp", "The adjusted-p value reported by DESeq2."),
        c("edger_adjp", "The adjusted-p value reported by edgeR."),
        c("The next 5 columns", "Statistics generated by limma."),
        c("limma_ave", "Average log2 expression observed by limma across _all_ samples."),
        c("limma_t", "T-statistic reported by limma given the log2FC and variances."),
        c("limma_p", "Derived from limma_t, the p-value asking 'is this logfc significant?'"),
        c("limma_b", "Use a Bayesian estimate to calculate log-odds significance instead of a student's test."),
        c("limma_q", "A q-value correction of the p-value above."),
        c("The next 5 columns", "Statistics generated by DESeq2."),
        c("deseq_basemean", "Analagous to limma's ave column, the base mean of all samples according to DESeq2."),
        c("deseq_lfcse", "The standard error observed given the log2 fold change."),
        c("deseq_stat", "T-statistic reported by DESeq2 given the log2FC and observed variances."),
        c("deseq_p", "Resulting p-value."),
        c("deseq_q", "False-positive corrected p-value."),
        c("The next 4 columns", "Statistics generated by edgeR."),
        c("edger_logcpm", "Similar to limma's ave and DESeq2's basemean, except only including the samples in the comparison."),
        c("edger_lr", "Undocumented, I am reasonably certain it is the T-statistic calculated by edgeR."),
        c("edger_p", "The observed p-value from edgeR."),
        c("edger_q", "The observed corrected p-value from edgeR."),
        c("The next 8 columns", "Statistics generated by the basic analysis written by trey."),
        c("basic_nummed", "log2 median values of the numerator for this comparison (like edgeR's basemean)."),
        c("basic_denmed", "log2 median values of the denominator for this comparison."),
        c("basic_numvar", "Variance observed in the numerator values."),
        c("basic_denvar", "Variance observed in the denominator values."),
        c("basic_logfc", "The log2 fold change observed by the basic analysis."),
        c("basic_t", "T-statistic from basic."),
        c("basic_p", "Resulting p-value."),
        c("basic_adjp", "BH correction of the p-value."),
        c("The next 5 columns", "Summaries of the limma/deseq/edger results."),
        c("fc_meta", "The mean fold-change value of limma/deseq/edger."),
        c("fc_var", "The variance between limma/deseq/edger."),
        c("fc_varbymed", "The ratio of the variance/median (closer to 0 means better agreement.)"),
        c("p_meta", "A meta-p-value of the mean p-values."),
        c("p_var", "Variance among the 3 p-values."),
        c("The following columns", "3 plots showing the expression coefficients of limma, edgeR, and DESeq2 respectively.")
        ))

    colnames(legend) <- c("column name", "column definition")
    xls_result <- write_xls(wb, data=legend, sheet="legend", rownames=FALSE,
                            title="Columns used in the following tables.")

    annot_df <- Biobase::fData(all_pairwise_result[["input"]][["expressionset"]])
    if (!is.null(extra_annot)) {
        annot_df <- merge(annot_df, extra_annot, by="row.names", all.x=TRUE)
        rownames(annot_df) <- annot_df[["Row.names"]]
        annot_df <- annot_df[-1]
    }

    combo <- list()
    limma_plots <- list()
    edger_plots <- list()
    deseq_plots <- list()
    sheet_count <- 0
    de_summaries <- data.frame()
    if (class(keepers) == "list") {
        ## Then keep specific tables in specific orientations.
        a <- 0
        keeper_len <- length(names(keepers))
        table_names <- list()
        for (name in names(keepers)) {
            a <- a + 1
            message(paste0("Working on ", a, "/", keeper_len, ": ",  name))
            sheet_count <- sheet_count + 1
            numerator <- keepers[[name]][[1]]
            denominator <- keepers[[name]][[2]]
            same_string <- paste0(numerator, "_vs_", denominator)
            inverse_string <- paste0(denominator, "_vs_", numerator)
            dat <- NULL
            plt <- NULL
            summary <- NULL
            found <- 0
            found_table <- NULL
            do_inverse <- NULL
            for (tab in names(edger[["contrast_list"]])) {
                if (tab == same_string) {
                    do_inverse <- FALSE
                    found <- found + 1
                    found_table <- same_string
                    message(paste0("Found table with ", same_string))
                } else if (tab == inverse_string) {
                    do_inverse <- TRUE
                    found <- found + 1
                    found_table <- inverse_string
                    message(paste0("Found inverse table with ", inverse_string))
                }
            }
            if (found > 0) {
                message(paste0("Running create_combined_table with inverse=", do_inverse))
                combined <- create_combined_table(limma, edger, deseq, basic, found_table, inverse=do_inverse,
                                                  annot_df=annot_df, include_basic=include_basic)
                dat <- combined[["data"]]
                summary <- combined[["summary"]]
                limma_plt <- NULL
                edger_plt <- NULL
                deseq_plt <- NULL
                if (isTRUE(do_inverse)) {
                    limma_plt <- sm(limma_coefficient_scatter(limma, x=denominator, y=numerator, gvis_filename=NULL))[["scatter"]]
                    edger_plt <- sm(edger_coefficient_scatter(edger, x=denominator, y=numerator, gvis_filename=NULL))[["scatter"]]
                    deseq_plt <- sm(deseq_coefficient_scatter(deseq, x=denominator, y=numerator, gvis_filename=NULL))[["scatter"]]
                } else {
                    limma_plt <- sm(limma_coefficient_scatter(limma, x=numerator, y=denominator, gvis_filename=NULL))[["scatter"]]
                    edger_plt <- sm(edger_coefficient_scatter(edger, x=numerator, y=denominator, gvis_filename=NULL))[["scatter"]]
                    deseq_plt <- sm(deseq_coefficient_scatter(deseq, x=numerator, y=denominator, gvis_filename=NULL))[["scatter"]]
                }
            } ## End checking that we found the numerator/denominator
            else {
                stop(paste0("Did not find either ", same_string, " nor ", inverse_string, "."))
            }
            combo[[name]] <- dat
            limma_plots[[name]] <- limma_plt
            edger_plots[[name]] <- edger_plt
            deseq_plots[[name]] <- deseq_plt
            de_summaries <- rbind(de_summaries, summary)
            table_names[[a]] <- summary[["table"]]
        }
        ## If you want all the tables in a dump
    } else if (class(keepers) == "character" & keepers == "all") {
        a <- 0
        names_length <- length(names(edger[["contrast_list"]]))
        table_names <- names(edger[["contrast_list"]])
        for (tab in names(edger[["contrast_list"]])) {
            a <- a + 1
            message(paste0("Working on table ", a, "/", names_length, ": ", tab))
            sheet_count <- sheet_count + 1
            combined <- create_combined_table(limma, edger, deseq, basic,
                                         tab, annot_df=annot_df, include_basic=include_basic)
            de_summaries <- rbind(de_summaries, combined[["summary"]])
            combo[[tab]] <- combined[["data"]]
            splitted <- strsplit(x=tab, split="_vs_")
            xname <- splitted[[1]][1]
            yname <- splitted[[1]][2]
            limma_plots[[tab]] <- sm(limma_coefficient_scatter(limma, x=xname, y=yname, gvis_filename=NULL))[["scatter"]]
            edger_plots[[tab]] <- sm(edger_coefficient_scatter(edger, x=xname, y=yname, gvis_filename=NULL))[["scatter"]]
            deseq_plots[[tab]] <- sm(deseq_coefficient_scatter(deseq, x=xname, y=yname, gvis_filename=NULL))[["scatter"]]
        }

        ## Or a single specific table
    } else if (class(keepers) == "character") {
        table <- keepers
        sheet_count <- sheet_count + 1
        if (table %in% names(edger[["contrast_list"]])) {
            message(paste0("I found ", table, " in the available contrasts."))
        } else {
            message(paste0("I did not find ", table, " in the available contrasts."))
            message(paste0("The available tables are: ", names(edger[["contrast_list"]])))
            table <- names(edger[["contrast_list"]])[[1]]
            message(paste0("Choosing the first table: ", table))
        }
        combined <- create_combined_table(limma, edger, deseq, basic,
                                          table, annot_df=annot_df, include_basic=include_basic)
        combo[[table]] <- combined[["data"]]
        splitted <- strsplit(x=tab, split="_vs_")
        de_summaries <- rbind(de_summaries, combined[["summary"]])
        table_names[[a]] <- combined[["summary"]][["table"]]
        xname <- splitted[[1]][1]
        yname <- splitted[[1]][2]
        limma_plots[[name]] <- sm(limma_coefficient_scatter(limma, x=xname, y=yname))[["scatter"]]
        edger_plots[[name]] <- sm(edger_coefficient_scatter(edger, x=xname, y=yname))[["scatter"]]
        deseq_plots[[name]] <- sm(limma_coefficient_scatter(deseq, x=xname, y=yname))[["scatter"]]
    } else {
        stop("I don't know what to do with your specification of tables to keep.")
    }

    comp <- NULL
    if (!is.null(excel)) {
        ## Starting a new counter of sheets.
        count <- 0
        for (tab in names(combo)) {
            count <- count + 1
            ddd <- combo[[count]]
            oddness = summary(ddd) ## until I did this I was getting errors I am guessing devtools::load_all() isn't clearing everything
            final_excel_title <- gsub(pattern='YYY', replacement=tab, x=excel_title)
            xls_result <- write_xls(data=ddd, wb=wb, sheet=tab, title=final_excel_title)
            if (!is.null(csv)) {
                csv_filename <- paste0(csv_basename, "_", tab, ".csv")
                write.csv(x=ddd, file=csv_filename)
            }
            if (isTRUE(add_plots)) {
                plot_column <- xls_result[["end_col"]] + 2
                message(paste0("Attempting to add a coefficient plot for ", names(combo)[[count]], " at column ", plot_column))
                openxlsx::writeData(wb, tab, x="Limma expression coefficients", startRow=1, startCol=plot_column)
                limma_plot <- limma_plots[[count]]
                print(limma_plot)
                openxlsx::insertPlot(wb, tab, width=plot_dim, height=plot_dim,
                                     startCol=plot_column, startRow=2, fileType="png", units="in")
                openxlsx::writeData(wb, tab, x="EdgeR expression coefficients", startRow=33, startCol=plot_column)
                edger_plot <- edger_plots[[count]]
                print(edger_plot)
                openxlsx::insertPlot(wb, tab, width=plot_dim, height=plot_dim,
                                     startCol=plot_column, startRow=34, fileType="png", units="in")
                openxlsx::writeData(wb, tab, x="DESeq2 expression coefficients", startRow=65, startCol=plot_column)
                deseq_plot <- deseq_plots[[count]]
                print(deseq_plot)
                openxlsx::insertPlot(wb, tab, width=plot_dim, height=plot_dim,
                                     startCol=plot_column, startRow=66, fileType="png", units="in")

            }
        }  ## End for loop
        count <- count + 1

        message("Writing summary information.")
        if (isTRUE(compare_plots)) {
            ## Add a graph on the final sheet of how similar the result types were
            comp_summary <- all_pairwise_result[["comparison"]][["comp"]]
            comp_plot <- all_pairwise_result[["comparison"]][["heat"]]
            de_summaries <- as.data.frame(de_summaries)
            rownames(de_summaries) <- table_names
            xls_result <- write_xls(wb, data=de_summaries, sheet="pairwise_summary",
                                    title="Summary of contrasts.")
            new_row <- xls_result[["end_row"]] + 2
            xls_result <- write_xls(wb, data=comp_summary, sheet="pairwise_summary",
                                    title="Pairwise correlation coefficients among differential expression tools.",
                                    start_row=new_row)
            new_row <- xls_result[["end_row"]] + 2
            message(paste0("Attempting to add the comparison plot to pairwise_summary at row: ", new_row + 1, " and column: ", 1))
            print(comp_plot)
            comp <- recordPlot()
            openxlsx::insertPlot(wb, "pairwise_summary", width=6, height=6,
                                 startRow=new_row + 1, startCol=1, fileType="png", units="in")
            logfc_comparisons <- compare_logfc_plots(combo)
            logfc_names <- names(logfc_comparisons)
            new_row <- new_row + 2
            for (c in 1:length(logfc_comparisons)) {
                new_row <- new_row + 32
                le <- logfc_comparisons[[c]][["le"]]
                ld <- logfc_comparisons[[c]][["ld"]]
                de <- logfc_comparisons[[c]][["de"]]

                tmpcol <- 1
                openxlsx::writeData(wb, "pairwise_summary", x=paste0("Comparing DE tools for the comparison of: ", logfc_names[c]),
                                    startRow=new_row - 2, startCol=tmpcol)
                openxlsx::writeData(wb, "pairwise_summary", x="Log2FC(Limma vs. EdgeR)", startRow=new_row - 1, startCol=tmpcol)
                print(le)
                openxlsx::insertPlot(wb, "pairwise_summary", width=6, height=6,
                                     startRow=new_row, startCol=tmpcol, fileType="png", units="in")

                tmpcol <- 8
                openxlsx::writeData(wb, "pairwise_summary", x="Log2FC(Limma vs. DESeq2)", startRow=new_row - 1, startCol=tmpcol)
                print(ld)
                openxlsx::insertPlot(wb, "pairwise_summary", width=6, height=6,
                                     startRow=new_row, startCol=tmpcol, fileType="png", units="in")

                tmpcol <- 15
                openxlsx::writeData(wb, "pairwise_summary", x="Log2FC(DESeq2 vs. EdgeR)", startRow=new_row - 1, startCol=tmpcol)
                print(de)
                openxlsx::insertPlot(wb, "pairwise_summary", width=6, height=6,
                                     startRow=new_row, startCol=tmpcol, fileType="png", units="in")
            }
        } ## End if compare_plots is TRUE
        message("Performing save of the workbook.")
        save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
        if (class(save_result) == "try-error") {
            message("Saving xlsx failed.  Rerunning now with arguments to save to csv files.")
            retlist <- combine_de_tables(all_pairwise_result,
                                         extra_annot=extra_annot,
                                         csv=excel,
                                         excel=NULL,
                                         excel_title=excel_title,
                                         excel_sheet=excel_sheet,
                                         keepers=keepers,
                                         include_basic=include_basic,
                                         add_plots=FALSE,
                                         compare_plots=FALSE)
        }
    } ## End if !is.null(excel)

    ret <- NULL
    if (is.null(retlist)) {
        ret <- list(
            "data" = combo,
            "limma_plots" = limma_plots,
            "edger_plots" = edger_plots,
            "deseq_plots" = deseq_plots,
            "comp_plot" = comp,
            "de_summay" = de_summaries)
    } else {
        ret <- retlist
    }
    return(ret)
}

#' Given a limma, edger, and deseq table, combine them into one.
#'
#' This combines the outputs from the various differential expression
#' tools and formalizes some column names to make them a little more
#' consistent.
#'
#' @param li Limma output table.
#' @param ed Edger output table.
#' @param de Deseq2 output table.
#' @param ba Basic output table.
#' @param table_name Name of the table to merge.
#' @param annot_df Add some annotation information?
#' @param inverse Invert the fold changes?
#' @param include_basic Include the basic table?
#' @param fc_cutoff Preferred logfoldchange cutoff.
#' @param p_cutoff Preferred pvalue cutoff.
#' @return List containing a) Dataframe containing the merged
#'     limma/edger/deseq/basic tables, and b) A summary of how many
#'     genes were observed as up/down by output table.
#' @export
create_combined_table <- function(li, ed, de, ba, table_name, annot_df=NULL, inverse=FALSE,
                                  include_basic=TRUE, fc_cutoff=1, p_cutoff=0.05) {
    li <- li[["all_tables"]][[table_name]]
    colnames(li) <- c("limma_logfc","limma_ave","limma_t","limma_p","limma_adjp","limma_b","limma_q")
    li <- li[, c("limma_logfc","limma_ave","limma_t","limma_b","limma_p","limma_adjp","limma_q")]
    de <- de[["all_tables"]][[table_name]]
    colnames(de) <- c("deseq_basemean","deseq_logfc","deseq_lfcse","deseq_stat","deseq_p","deseq_adjp","deseq_q")
    de <- de[, c("deseq_logfc","deseq_basemean","deseq_lfcse","deseq_stat","deseq_p","deseq_adjp","deseq_q")]
    ed <- ed[["all_tables"]][[table_name]]
    colnames(ed) <- c("edger_logfc","edger_logcpm","edger_lr","edger_p","edger_adjp","edger_q")
    ba <- ba[["all_tables"]][[table_name]]
    ba <- ba[, c("numerator_median","denominator_median","numerator_var","denominator_var", "logFC", "t", "p", "adjp")]
    colnames(ba) <- c("basic_nummed","basic_denmed", "basic_numvar", "basic_denvar", "basic_logfc", "basic_t", "basic_p", "basic_adjp")

    comb <- merge(li, de, by="row.names")
    comb <- merge(comb, ed, by.x="Row.names", by.y="row.names")
    if (isTRUE(include_basic)) {
        comb <- merge(comb, ba, by.x="Row.names", by.y="row.names")
    }
    rownames(comb) <- comb[["Row.names"]]
    comb <- comb[-1]
    comb[is.na(comb)] <- 0
    if (isTRUE(include_basic)) {
        comb <- comb[, c("limma_logfc","deseq_logfc","edger_logfc","limma_adjp","deseq_adjp","edger_adjp","limma_ave","limma_t","limma_p","limma_b","limma_q","deseq_basemean","deseq_lfcse","deseq_stat","deseq_p","deseq_q", "edger_logcpm","edger_lr","edger_p","edger_q","basic_nummed","basic_denmed", "basic_numvar", "basic_denvar", "basic_logfc", "basic_t", "basic_p", "basic_adjp")]
    } else {
        comb <- comb[, c("limma_logfc","deseq_logfc","edger_logfc","limma_adjp","deseq_adjp","edger_adjp","limma_ave","limma_t","limma_p","limma_b","limma_q","deseq_basemean","deseq_lfcse","deseq_stat","deseq_p","deseq_q", "edger_logcpm","edger_lr","edger_p","edger_q")]
    }
    if (isTRUE(inverse)) {
        comb[["limma_logfc"]] <- comb[["limma_logfc"]] * -1
        comb[["deseq_logfc"]] <- comb[["deseq_logfc"]] * -1
        comb[["deseq_stat"]] <- comb[["deseq_stat"]] * -1
        comb[["edger_logfc"]] <- comb[["edger_logfc"]] * -1
        if (isTRUE(include_basic)) {
            comb[["basic_logfc"]] <- comb[["basic_logfc"]] * -1
        }
    }
    ## I made an odd choice in a moment to normalize.quantils the combined fold changes
    ## This should be reevaluated
    temp_fc <- cbind(as.numeric(comb[["limma_logfc"]]),
                     as.numeric(comb[["edger_logfc"]]),
                     as.numeric(comb[["deseq_logfc"]]))
    temp_fc <- preprocessCore::normalize.quantiles(as.matrix(temp_fc))
    comb[["fc_meta"]] <- rowMeans(temp_fc, na.rm=TRUE)
    comb[["fc_var"]] <- genefilter::rowVars(temp_fc, na.rm=TRUE)
    comb[["fc_varbymed"]] <- comb$fc_var / comb$fc_meta
    temp_p <- cbind(as.numeric(comb[["limma_p"]]),
                    as.numeric(comb[["edger_p"]]),
                    as.numeric(comb[["deseq_p"]]))
    comb[["p_meta"]] <- rowMeans(temp_p, na.rm=TRUE)
    comb[["p_var"]] <- genefilter::rowVars(temp_p, na.rm=TRUE)
    comb[["fc_meta"]] <- signif(x=comb[["fc_meta"]], digits=4)
    comb[["fc_var"]] <- format(x=comb[["fc_var"]], digits=4, scientific=TRUE)
    comb[["fc_varbymed"]] <- format(x=comb[["fc_varbymed"]], digits=4, scientific=TRUE)
    comb[["p_var"]] <- format(x=comb[["p_var"]], digits=4, scientific=TRUE)
    comb$p_meta <- format(x=comb[["p_meta"]], digits=4, scientific=TRUE)
    if (!is.null(annot_df)) {
        ## colnames(annot_df) <- gsub("[[:digit:]]", "", colnames(annot_df))
        colnames(annot_df) <- gsub("[[:punct:]]", "", colnames(annot_df))
        comb <- merge(annot_df, comb, by="row.names", all.y=TRUE)
        rownames(comb) <- comb[["Row.names"]]
        comb <- comb[-1]
        colnames(comb) <- make.names(tolower(colnames(comb)), unique=TRUE)
    }

    up_fc <- fc_cutoff
    down_fc <- -1 * fc_cutoff
    message(paste0("The table is: ", table_name))
    summary_table_name <- table_name
    if (isTRUE(inverse)) {
        summary_table_name <- paste0(summary_table_name, "-inverted")
    }
    summary_lst <- list(
        "table" = summary_table_name,
        "total" = nrow(comb),
        "limma_up" = sum(comb[["limma_logfc"]] >= up_fc),
        "limma_sigup" = sum(comb[["limma_logfc"]] >= up_fc & as.numeric(comb[["limma_adjp"]]) <= p_cutoff),
        "deseq_up" = sum(comb[["deseq_logfc"]] >= up_fc),
        "deseq_sigup" = sum(comb[["deseq_logfc"]] >= up_fc & as.numeric(comb[["deseq_adjp"]]) <= p_cutoff),
        "edger_up" = sum(comb[["edger_logfc"]] >= up_fc),
        "edger_sigup" = sum(comb[["edger_logfc"]] >= up_fc & as.numeric(comb[["edger_adjp"]]) <= p_cutoff),
        "basic_up" = sum(comb[["basic_logfc"]] >= up_fc),
        "basic_sigup" = sum(comb[["basic_logfc"]] >= up_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
        "limma_down" = sum(comb[["limma_logfc"]] <= down_fc),
        "limma_sigdown" = sum(comb[["limma_logfc"]] <= down_fc & as.numeric(comb[["limma_adjp"]]) <= p_cutoff),
        "deseq_down" = sum(comb[["deseq_logfc"]] <= down_fc),
        "deseq_sigdown" = sum(comb[["deseq_logfc"]] <= down_fc & as.numeric(comb[["deseq_adjp"]]) <= p_cutoff),
        "edger_down" = sum(comb[["edger_logfc"]] <= down_fc),
        "edger_sigdown" = sum(comb[["edger_logfc"]] <= down_fc & as.numeric(comb[["edger_adjp"]]) <= p_cutoff),
        "basic_down" = sum(comb[["basic_logfc"]] <= down_fc),
        "basic_sigdown" = sum(comb[["basic_logfc"]] <= down_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
        "meta_up" = sum(comb[["fc_meta"]] >= up_fc),
        "meta_sigup" = sum(comb[["fc_meta"]] >= up_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff),
        "meta_down" = sum(comb[["fc_meta"]] <= down_fc),
        "meta_sigdown" = sum(comb[["fc_meta"]] <= down_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff)
        )

    ret <- list(
        "data" = comb,
        "summary" = summary_lst)
    return(ret)
}


extract_siggenes <- function(...) { extract_significant_genes(...) }
#' Extract the sets of genes which are significantly up/down regulated
#' from the combined tables.
#'
#' Given the output from combine_de_tables(), extract the genes in
#' which we have the greatest likely interest, either because they
#' have the largest fold changes, lowest p-values, fall outside a
#' z-score, or are at the top/bottom of the ranked list.
#'
#' @param combined Output from combine_de_tables().
#' @param according_to What tool(s) decide 'significant?'  One may use
#'     the deseq, edger, limma, basic, meta, or all.
#' @param fc Log fold change to define 'significant'.
#' @param p (Adjusted)p-value to define 'significant'.
#' @param z Z-score to define 'significant'.
#' @param n Take the top/bottom-n genes.
#' @param p_type use an adjusted p-value?
#' @param excel Write the results to this excel file, or NULL.
#' @param csv Write csv instead of xlsx when running OOM.
#' @return The set of up-genes, down-genes, and numbers therein.
#' @seealso \code{\link{combine_de_tables}}
#' @export
extract_significant_genes <- function(combined, according_to="all", fc=1.0,
                                      p=0.05, z=NULL, n=NULL, p_type="adj",
                                      excel="excel/significant_genes.xlsx", csv=NULL) {
    if (!is.null(combined[["plots"]])) {
        combined <- combined[["data"]]
    }
    trimmed_up <- list()
    trimmed_down <- list()
    up_titles <- list()
    down_titles <- list()
    title_append <- ""
    if (!is.null(fc)) {
        title_append <- paste0(title_append, " log2fc><", fc)
    }
    if (!is.null(p)) {
        title_append <- paste0(title_append, " p<", p)
    }
    if (!is.null(z)) {
        title_append <- paste0(title_append, " z><", z)
    }
    if (!is.null(n)) {
        title_append <- paste0(title_append, " top|bottom n=", n)
    }
    num_tables <- length(names(combined))
    table_count <- 0
    if (according_to == "all") {
        according_to <- c("limma","edger","deseq","basic")
    }

    wb <- openxlsx::createWorkbook(creator="hpgltools")
    ret <- list()
    summary_count <- 0
    sheet_count <- 0
    for (according in according_to) {
        summary_count <- summary_count + 1
        ret[[according]] <- list()
        change_counts_up <- list()
        change_counts_down <- list()
        for (table_name in names(combined[["data"]])) {
            table_count <- table_count + 1
            message(paste0("Writing excel data sheet ", table_count, "/", num_tables))
            table <- combined[["data"]][[table_name]]
            fc_column <- paste0(according, "_logfc")
            p_column <- paste0(according, "_adjp")
            if (p_type != "adj") {
                p_column <- paste0(according, "_p")
            }
            trimming <- get_sig_genes(table, fc=fc, p=p, z=z, n=n,
                                      column=fc_column, p_column=p_column)
            trimmed_up[[table_name]] <- trimming[["up_genes"]]
            change_counts_up[[table_name]] <- nrow(trimmed_up[[table_name]])
            trimmed_down[[table_name]] <- trimming[["down_genes"]]
            change_counts_down[[table_name]] <- nrow(trimmed_down[[table_name]])
            up_title <- paste0("Table SXXX: Genes deemed significantly up in ", table_name, " with", title_append, " according to ", according)
            up_titles[[table_name]] <- up_title
            down_title <- paste0("Table SXXX: Genes deemed significantly down in ", table_name, " with", title_append, " according to ", according)
            down_titles[[table_name]] <- down_title
        } ## End extracting significant genes for loop

        change_counts <- cbind(change_counts_up, change_counts_down)
        summary_title <- paste0("Counting the number of changed genes by contrast according to ", according, " with ", title_append)
        ## xls_result <- write_xls(data=change_counts, sheet="number_changed_genes", file=sig_table,
        ##                         title=summary_title,
        ##                         overwrite_file=TRUE, newsheet=TRUE)
        ret[[according]] <- list(
            "ups" = trimmed_up,
            "downs" = trimmed_down,
            "counts" = change_counts,
            "up_titles" = up_titles,
            "down_titles" = down_titles,
            "counts_title" = summary_title)
        if (is.null(excel)) {
            message("Not printing excel sheets for the significant genes.")
        } else {
            message(paste0("Printing significant genes to the file: ", excel))
            xlsx_ret <- print_ups_downs(ret[[according]], wb=wb, excel=excel, according=according, summary_count=summary_count, csv=csv)
            wb <- xlsx_ret[["workbook"]]
        }
    } ## End list of according_to's
    if (!is.null(excel)) {
        excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
    }
    return(ret)
}

#' Reprint the output from extract_significant_genes().
#'
#' I found myself needing to reprint these excel sheets because I
#' added some new information. This shortcuts that process for me.
#'
#' @param upsdowns Output from extract_significant_genes().
#' @param wb Workbook object to use for writing, or start a new one.
#' @param excel Filename for writing the data.
#' @param csv Write a csv instead/also?
#' @param according Use limma, deseq, or edger for defining 'significant'.
#' @param summary_count For spacing sequential tables one after another.
#' @return Return from write_xls.
#' @seealso \code{\link{combine_de_tables}}
#' @export
print_ups_downs <- function(upsdowns, wb=NULL, excel="excel/significant_genes.xlsx", csv=NULL,
                            according="limma", summary_count=1) {
    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="hpgltools")
    }
    csv_basename <- NULL
    if (!is.null(csv)) {
        if (is.null(excel)) {
            csv_basename <- "excel/csv_export"
        } else {
            csv_basename <- excel
            csv_basename <- gsub(pattern="\\.xlsx", replacement="", x=csv_basename)
        }
    }
    ups <- upsdowns[["ups"]]
    downs <- upsdowns[["downs"]]
    up_titles <- upsdowns[["up_titles"]]
    down_titles <- upsdowns[["down_titles"]]
    summary <- upsdowns[["counts"]]
    summary_title <- upsdowns[["counts_title"]]
    table_count <- 0
    summary_count <- summary_count - 1
    num_tables <- length(names(ups))
    summary_start <- ((num_tables + 2) * summary_count) + 1
    xls_summary_result <- write_xls(wb, data=summary, start_col=2, start_row=summary_start, sheet="number_changed_genes", title=summary_title)
    if (!is.null(csv)) {
        csv_filename <- paste0(csv_basename, "_num_changed.csv")
        write.csv(x=summary, file=csv_filename)
    }
    for (base_name in names(ups)) {
        table_count <- table_count + 1
        up_name <- paste0("up_", table_count, according, "_", base_name)
        down_name <- paste0("down_", table_count, according, "_", base_name)
        up_table <- ups[[table_count]]
        down_table <- downs[[table_count]]
        up_title <- up_titles[[table_count]]
        down_title <- down_titles[[table_count]]
        message(paste0(table_count, "/", num_tables, ": Writing excel data sheet ", up_name))
        xls_result <- write_xls(data=up_table, wb=wb, sheet=up_name, title=up_title)
        message(paste0(table_count, "/", num_tables, ": Writing excel data sheet ", down_name))
        xls_result <- write_xls(data=down_table, wb=wb, sheet=down_name, title=down_title)
        if (!is.null(csv)) {
            csv_filename <- paste0(csv_basename, "_", up_name, ".csv")
            write.csv(x=up_table, file=csv_filename)
            csv_filename <- paste0(csv_basename, "_", down_name, ".csv")
            write.csv(x=down_table, file=csv_filename)
        }
    } ## End for each name in ups
    message("Writing changed genes summary on last sheet.")
    return(xls_result)
}

#' Compare logFC values from limma and friends
#'
#' There are some peculiar discrepencies among these tools, what is up with that?
#'
#' @param combined_tables The combined tables from limma et al.
#' @return Some plots
#' @export
compare_logfc_plots <- function(combined_tables) {
    plots <- list()
    data <- NULL
    if (!is.null(combined_tables[["data"]])) {
        data <- combined_tables[["data"]]
    } else {
        data <- combined_tables
    }
    for (count in 1:length(data)) {
        tab <- data[[count]]
        le_data <- tab[, c("limma_logfc", "edger_logfc", "limma_adjp", "edger_adjp")]
        le <- sm(plot_linear_scatter(le_data, pretty_colors=FALSE)[["scatter"]])
        ld_data <- tab[, c("limma_logfc", "deseq_logfc", "limma_adjp", "deseq_adjp")]
        ld <- sm(plot_linear_scatter(ld_data, pretty_colors=FALSE)[["scatter"]])
        de_data <- tab[, c("deseq_logfc", "edger_logfc", "deseq_adjp", "edger_adjp")]
        de <- sm(plot_linear_scatter(de_data, pretty_colors=FALSE)[["scatter"]])
        lb_data <- tab[, c("limma_logfc", "basic_logfc", "limma_adjp", "basic_p")]
        lb <- sm(plot_linear_scatter(lb_data, pretty_colors=FALSE)[["scatter"]])
        db_data <- tab[, c("deseq_logfc", "basic_logfc", "deseq_adjp", "basic_p")]
        db <- sm(plot_linear_scatter(db_data, pretty_colors=FALSE)[["scatter"]])
        eb_data <- tab[, c("edger_logfc", "basic_logfc", "edger_adjp", "basic_p")]
        eb <- sm(plot_linear_scatter(eb_data, pretty_colors=FALSE)[["scatter"]])
        name <- names(data)[[count]]
        compared <- list(
            "le" = le, "ld" = ld, "de" = de,
            "lb" = lb, "db" = db, "eb" = eb)
        plots[[name]] <- compared
    }
    return(plots)
}

#' Small hack of limma's exampleData() to allow for arbitrary data set
#' sizes.
#'
#' exampleData has a set number of genes/samples it creates. This
#' relaxes that restriction.
#'
#' @param ngenes How many genes in the fictional data set?
#' @param columns How many samples in this data set?
#' @return Matrix of pretend counts.
#' @seealso \pkg{limma}
#' @examples
#' \dontrun{
#'  pretend = make_exampledata()
#' }
#' @export
make_exampledata <- function (ngenes=1000, columns=5) {
    q0 <- stats::rexp(ngenes, rate = 1/250)
    is_DE <- stats::runif(ngenes) < 0.3
    lfc <- stats::rnorm(ngenes, sd = 2)
    q0A <- ifelse(is_DE, q0 * 2^(lfc / 2), q0)
    q0B <- ifelse(is_DE, q0 * 2^(-lfc / 2), q0)
    ##    true_sf <- c(1, 1.3, 0.7, 0.9, 1.6)
    true_sf <- abs(stats::rnorm(columns, mean=1, sd=0.4))
    cond_types <- ceiling(sqrt(columns))
    ##    conds <- c("A", "A", "B", "B", "B")
    ##x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
    conds <- sample(LETTERS[1:cond_types], columns, replace=TRUE)
    m <- t(sapply(seq_len(ngenes),
                  function(i) sapply(1:columns,
                                     function(j) rnbinom(1,
                                                         mu = true_sf[j] * ifelse(conds[j] == "A", q0A[i], q0B[i]),
                                                         size = 1/0.2))))
    rownames(m) <- paste("gene", seq_len(ngenes), ifelse(is_DE, "T", "F"), sep = "_")
    example <- DESeq::newCountDataSet(m, conds)
    return(example)
}

#' Run makeContrasts() with all pairwise comparisons.
#'
#' In order to have uniformly consistent pairwise contrasts, I decided
#' to avoid potential human erors(sic) by having a function generate
#' all contrasts.
#'
#' @param model Model describing the conditions/batches/etc in the experiment.
#' @param conditions Factor of conditions in the experiment.
#' @param do_identities Include all the identity strings? Limma can
#'     use this information while edgeR can not.
#' @param do_pairwise Include all pairwise strings? This shouldn't
#'     need to be set to FALSE, but just in case.
#' @param extra_contrasts Optional string of extra contrasts to include.
#' @return List including the following information:
#'   all_pairwise_contrasts = the result from makeContrasts(...)
#'   identities = the string identifying each condition alone
#'   all_pairwise = the string identifying each pairwise comparison alone
#'   contrast_string = the string passed to R to call makeContrasts(...)
#'   names = the names given to the identities/contrasts
#' @seealso \link[limma]{makeContrasts}
#' @examples
#' \dontrun{
#'  pretend = make_pairwise_contrasts(model, conditions)
#' }
#' @export
make_pairwise_contrasts <- function(model, conditions, do_identities=TRUE,
                                    do_pairwise=TRUE, extra_contrasts=NULL) {
    tmpnames <- colnames(model)
    tmpnames <- gsub(pattern="data[[:punct:]]", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="-", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="+", replacement="", x=tmpnames)
    tmpnames <- gsub(pattern="conditions^(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub(pattern="conditions", replacement="", x=tmpnames)
    colnames(model) <- tmpnames
    conditions <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=conditions)
    condition_table <- table(conditions)
    identities <- list()
    contrast_string <- ""
    eval_strings <- list()
    for (c in 1:length(condition_table)) {
        identity_name <- names(condition_table[c])
        identity_string <- paste(identity_name, " = ", identity_name, ",", sep="")
        identities[identity_name] <- identity_string
    }
    ## If I also create a sample condition 'alice', and also perform a subtraction
    ## of 'alice' from 'bob', then the full makeContrasts() will be:
    ## makeContrasts(bob=bob, alice=alice, bob_vs_alice=(bob)-(alice), levels=design)
    ## The parentheses in this case are superfluous, but they remind me that I am finally
    ## doing some match, and also they remind me that we can do much more complex things like:
    ## ((bob-alice)-(jane-alice)) or whatever....
    all_pairwise <- list()
    identity_names <- names(identities)
    lenminus <- length(identities) - 1
    for (c in 1:lenminus) {
        c_name <- names(identities[c])
        nextc <- c+1
        for (d in nextc:length(identities)) {
            d_name <- names(identities[d])
            minus_string <- paste(d_name, "_vs_", c_name, sep="")
            exprs_string <- paste(minus_string, "=", d_name, "-", c_name, ",", sep="")
            all_pairwise[minus_string] <- exprs_string
        }
    }
    ## At this point, I have strings which represent the definition of every
    ## sample condition as well as strings which represent every possible
    ## B-A where B comes somewhere after A in the model matrix.
    ## The goal now is to create the variables in the R environment
    ## and add them to makeContrasts()
    if (isTRUE(do_identities)) {
        eval_strings <- append(eval_strings, identities)
    }
    if (isTRUE(do_pairwise)) {
        eval_strings <- append(eval_strings, all_pairwise)
    }
    eval_names <- names(eval_strings)
    if (!is.null(extra_contrasts)) {
        extra_eval_strings <- strsplit(extra_contrasts, "\\n")
        extra_eval_names <- extra_eval_strings
        extra_eval_names <- stringi::stri_replace_all_regex(extra_eval_strings[[1]], "^(\\s*)(\\w+)=.*$", "$2")
        eval_strings <- append(eval_strings, extra_contrasts)
    }
    ## for (f in 1:length(eval_strings)) {
    ##     eval_name = names(eval_strings[f])
    ##     message(paste("Setting ", eval_name, " with expression:<<", eval_strings[f], ">>", sep=""))
    ##     eval(parse(text=as.character(eval_strings[f])))
    ## }
    ## Now we have bob=(somestuff) in memory in R's environment
    ## Add them to makeContrasts()
    contrast_string <- paste0("all_pairwise_contrasts = limma::makeContrasts(")
    for (f in 1:length(eval_strings)) {
        ## eval_name = names(eval_strings[f])
        eval_string <- paste0(eval_strings[f])
        contrast_string <- paste(contrast_string, eval_string, sep="   ")
    }
    ## The final element of makeContrasts() is the design from voom()
    contrast_string <- paste0(contrast_string, "levels=model)")
    eval(parse(text=contrast_string))
    ## I like to change the column names of the contrasts because by default
    ## they are kind of obnoxious and too long to type

    if (!is.null(extra_contrasts)) {
        eval_names <- append(eval_names, extra_eval_names)
    }
    colnames(all_pairwise_contrasts) <- eval_names
    result <- list(
        "all_pairwise_contrasts" = all_pairwise_contrasts,
        "identities" = identities,
        "identity_names" = identity_names,
        "all_pairwise" = all_pairwise,
        "contrast_string" = contrast_string,
        "names" = eval_names)
    return(result)
}

#' Get a set of up/down differentially expressed genes.
#'
#' Take one or more criteria (fold change, rank order, (adj)p-value,
#' z-score from median FC) and use them to extract the set of genes
#' which are defined as 'differentially expressed.'  If no criteria
#' are provided, it arbitrarily chooses all genes outside of 1-z.
#'
#' @param table Table from limma/edger/deseq.
#' @param n Rank-order top/bottom number of genes to take.
#' @param z Number of z-scores >/< the median to take.
#' @param fc Fold-change cutoff.
#' @param p P-value cutoff.
#' @param fold Identifier reminding how to get the bottom portion of a
#'     fold-change (plusminus says to get the negative of the
#'     positive, otherwise 1/positive is taken).  This effectively
#'     tells me if this is a log fold change or not.
#' @param column Table's column used to distinguish top vs. bottom.
#' @param p_column Table's column containing (adjusted or not)p-values.
#' @return Subset of the up/down genes given the provided criteria.
#' @export
get_sig_genes <- function(table, n=NULL, z=NULL, fc=NULL, p=NULL,
                          column="logFC", fold="plusminus", p_column="adj.P.Val") {
    if (is.null(z) & is.null(n) & is.null(fc)) {
        message("No n, z, nor fc provided, setting z to 1.")
        z <- 1
    }
    up_genes <- table
    down_genes <- table

    if (!is.null(p)) {
        up_idx <- as.numeric(up_genes[[p_column]]) <= p
        ## Remember we have these reformatted as scientific
        up_genes <- up_genes[up_idx, ]
        down_idx <- as.numeric(down_genes[[p_column]]) <= p
        down_genes <- down_genes[down_idx, ]
        ## Going to add logic in case one does not ask for fold change
        ## In that case, a p-value assertion should still know the difference between up and down
        ## But it should also still know the difference between ratio and log changes
        if (fold == 'plusminus' | fold == 'log') {
            message(paste0("Assuming the fold changes are on the log scale and so taking >< 0"))
            ##up_idx <- up_genes[, column] > 0.0
            up_idx <- as.numeric(up_genes[[column]]) > 0.0
            up_genes <- up_genes[up_idx, ]
            down_idx <- as.numeric(down_genes[[column]]) < 0.0
            down_genes <- down_genes[down_idx, ]
        } else {
            ## plusminus refers to a positive/negative number of logfold changes from a logFC(1) = 0
            up_idx <- as.numeric(up_genes[[column]]) > 1.0
            up_genes <- up_genes[up_idx, ]
            down_idx <- as.numeric(down_genes[[column]]) < 1.0
            down_genes <- down_genes[down_idx, ]
        }
        message(paste0("After (adj)p filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After (adj)p filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(fc)) {
        up_idx <- as.numeric(up_genes[[column]]) >= fc
        up_genes <- up_genes[up_idx, ]
        if (fold == 'plusminus' | fold == 'log') {
            message(paste0("Assuming the fold changes are on the log scale and so taking -1 * fc"))
            ## plusminus refers to a positive/negative number of logfold changes from a logFC(1) = 0
            down_idx <- as.numeric(down_genes[[column]]) <= (fc * -1)
            down_genes <- down_genes[down_idx, ]
        } else {
            message(paste0("Assuming the fold changes are on a ratio scale and so taking 1/fc"))
            ## If it isn't log fold change, then values go from 0..x where 1 is unchanged
            down_idx <- as.numeric(down_genes[[column]]) <= (1 / fc)
            down_genes <- down_genes[down_idx, ]
        }
        message(paste0("After fold change filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After fold change filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(z)) {
        ## Take an arbitrary number which are >= and <= a value which is z zscores from the median.
        message(paste0("Getting the genes >= ", z, " z scores away from the median of all."))
        ## Use the entire table for the summary
        out_summary <- summary(as.numeric(table[[column]]))
        out_mad <- stats::mad(as.numeric(table[[column]]), na.rm=TRUE)
        up_median_dist <- out_summary["Median"] + (out_mad * z)
        down_median_dist <- out_summary["Median"] - (out_mad * z)
        ## But use the (potentially already trimmed) up/down tables for indexing
        up_idx <- as.numeric(up_genes[[column]]) >= up_median_dist
        up_genes <- up_genes[up_idx, ]
        down_idx <- as.numeric(down_genes[[column]]) <= down_median_dist
        down_genes <- down_genes[down_idx, ]
        message(paste0("After z filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After z filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }

    if (!is.null(n)) {
        ## Take a specific number of genes at the top/bottom of the rank ordered list.
        message(paste0("Getting the top and bottom ", n, " genes."))
        upranked <- up_genes[order(as.numeric(up_genes[[column]]), decreasing=TRUE), ]
        up_genes <- head(upranked, n=n)
        downranked <- down_genes[order(as.numeric(down_genes[[column]])), ]
        down_genes <- head(downranked, n=n)
        message(paste0("After top-n filter, the up genes table has ", dim(up_genes)[1], " genes."))
        message(paste0("After bottom-n filter, the down genes table has ", dim(down_genes)[1], " genes."))
    }
    up_genes <- up_genes[order(as.numeric(up_genes[[column]]), decreasing=TRUE), ]
    down_genes <- down_genes[order(as.numeric(down_genes[[column]]), decreasing=FALSE), ]
    ret = list(
        "up_genes" = up_genes,
        "down_genes" = down_genes)
    return(ret)
}

#' Remove multicopy genes from up/down gene expression lists.
#'
#' In our parasite data, there are a few gene types which are
#' consistently obnoxious.  Multi-gene families primarily where the
#' coding sequences are divergent, but the UTRs nearly identical.  For
#' these genes, our sequence based removal methods fail and so this
#' just excludes them by name.
#'
#' @param de_list List of sets of genes deemed significantly
#'     up/down with a column expressing approximate count numbers.
#' @param max_copies Keep only those genes with <= n putative
#'     copies.
#' @param semantic Set of strings with gene names to exclude.
#' @param semantic_column Column in the DE table used to find the
#'     semantic strings for removal.
#' @return Smaller list of up/down genes.
#' @export
semantic_copynumber_filter <- function(de_list, max_copies=2, semantic=c('mucin','sialidase','RHS','MASP','DGF'), semantic_column='1.tooltip') {
    count <- 0
    for (table in de_list[["ups"]]) {
        count <- count + 1
        tab <- de_list[["ups"]][[count]]
        table_name <- names(de_list[["ups"]])[[count]]
        message(paste0("Working on ", table_name))
        file <- paste0("singletons/gene_counts/up_", table_name, ".fasta.out.count")
        tmpdf <- try(read.table(file), silent=TRUE)
        if (class(tmpdf) == 'data.frame') {
            colnames(tmpdf) = c("ID", "members")
            tab <- merge(tab, tmpdf, by.x="row.names", by.y="ID")
            rownames(tab) <- tab$Row.names
            tab <- tab[-1]
            tab <- tab[count <= max_copies, ]
            for (string in semantic) {
                idx <- grep(pattern=string, x=tab[, semantic_column])
                tab <- tab[-idx]
            }
            de_list[["ups"]][[count]] <- tab
        }
    }
    count <- 0
    for (table in de_list[["downs"]]) {
        count <- count + 1
        tab <- de_list[["downs"]][[count]]
        table_name <- names(de_list[["downs"]])[[count]]
        message(paste0("Working on ", table_name))
        file <- paste0("singletons/gene_counts/down_", table_name, ".fasta.out.count")
        tmpdf <- try(read.table(file), silent=TRUE)
        if (class(tmpdf) == 'data.frame') {
            colnames(tmpdf) = c("ID","members")
            tab <- merge(tab, tmpdf, by.x="row.names", by.y="ID")
            rownames(tab) <- tab[["Row.names"]]
            tab <- tab[-1]
            tab <- tab[count <= max_copies, ]
            for (string in semantic) {
                ## Is this next line correct?  shouldn't it be tab[, semantic_column]?
                ## idx <- grep(pattern=string, x=tab[[, semantic_column]])
                idx <- grep(pattern=string, x=tab[, semantic_column])
                tab <- tab[-idx]
            }
            de_list[["downs"]][[count]] <- tab
        }
    }
    return(de_list)
}

#' Given a DE table with fold changes and p-values, show how 'significant' changes with changing cutoffs.
#'
#' Sometimes one might want to know how many genes are deemed significant while shifting the bars
#' which define significant.  This provides that metrics as a set of tables of numbers of
#' significant up/down genes when p-value is held constant, as well as number when fold-change is
#' held constant.
#'
#' @param table DE table to examine.
#' @param p_column Column in the DE table defining the changing p-value cutoff.
#' @param fc_column Column in the DE table defining the changing +/- log fold change.
#' @param bins Number of incremental changes in p-value/FC to examine.
#' @param constant_p When plotting changing FC, where should the p-value be held?
#' @param constant_fc When plotting changing p, where should the FC be held?
#' @return Plots and dataframes describing the changing definition of 'significant.'
#' @export
plot_num_siggenes <- function(table, p_column="limma_adjp", fc_column="limma_logfc", bins=100, constant_p=0.05, constant_fc=0) {
    fc_column="limma_logfc"
    p_column="limma_adjp"
    bins = 100
    num_genes <- nrow(table)
    min_fc <- min(table[[fc_column]])
    neutral_fc <- 0.0
    max_fc <- max(table[[fc_column]])
    min_p <- 0.0
    max_p <- 1.0
    up_increments <- max_fc / bins
    down_increments <- min_fc / bins
    p_increments <- (max_p - min_p) / bins

    constant_up_fc <- constant_fc
    constant_down_fc <- constant_fc * -1.0
    start_up <- max_fc
    start_down <- min_fc
    start_p <- 0.00
    current_up_fc <- start_up
    current_down_fc <- start_down
    current_p <- start_p
    up_nums <- data.frame()
    down_nums <- data.frame()
    p_nums <- data.frame()
    for (inc in 1:bins) {
        current_up_fc <- current_up_fc - up_increments
        current_down_fc <- current_down_fc - down_increments
        current_p <- current_p + p_increments
        num_up <- sum(table[[fc_column]] >= current_up_fc & table[[p_column]] <= constant_p)
        num_down <- sum(table[[fc_column]] <= current_down_fc & table[[p_column]] <= constant_p)
        num_pup <- sum(table[[fc_column]] >= constant_up_fc & table[[p_column]] <= current_p)
        num_pdown <- sum(table[[fc_column]] <= constant_down_fc & table[[p_column]] <= current_p)
        up_nums <- rbind(up_nums, c(current_up_fc, num_up))
        down_nums <- rbind(down_nums, c(current_down_fc, num_down))
        p_nums <- rbind(p_nums, c(current_p, num_pup, num_pdown))
    }
    colnames(p_nums) <- c("p","up","down")
    colnames(up_nums) <- c("fc","num")
    colnames(down_nums) <- c("fc", "num")

    putative_up_inflection <- inflection::findiplist(x=as.matrix(up_nums[[1]]), y=as.matrix(up_nums[[2]]), 0)
    up_point_num <- putative_up_inflection[2,1]
    up_label <- paste0("At fc=", signif(up_nums[up_point_num, ][["fc"]], 4), " and p=", constant_p, ", ", up_nums[up_point_num, ][["num"]], " genes are de.")
    up_plot <- ggplot(data=up_nums, aes_string(x="fc", y="num")) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=up_nums[[2]][[up_point_num]]) +
        ggplot2::geom_vline(xintercept=up_nums[[1]][[up_point_num]]) +
        ggplot2::geom_vline(xintercept=1.0, colour="red")

    putative_down_inflection <- inflection::findiplist(x=as.matrix(down_nums[[1]]), y=as.matrix(down_nums[[2]]), 0)
    down_point_num <- putative_down_inflection[1,2]
    down_plot <- ggplot(data=down_nums, aes_string(x="fc", y="num")) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=down_nums[[2]][[down_point_num]]) +
        ggplot2::geom_vline(xintercept=down_nums[[1]][[down_point_num]]) +
        ggplot2::geom_vline(xintercept=-1.0, colour="red")

    putative_pup_inflection <- inflection::findiplist(x=as.matrix(p_nums[[1]]), y=as.matrix(p_nums[[2]]), 1)
    pup_point_num <- putative_pup_inflection[2,1]
    putative_pdown_inflection <- inflection::findiplist(x=as.matrix(p_nums[[1]]), y=as.matrix(p_nums[[2]]), 1)
    pdown_point_num <- putative_pdown_inflection[2,1]
    p_plot <- ggplot(data=p_nums) +
        ggplot2::geom_point(aes_string(x="p", y="up"), colour="darkred") +
        ggplot2::geom_point(aes_string(x="p", y="down"), colour="darkblue") +
        ggplot2::geom_vline(xintercept=0.05, colour="red") +
        ggplot2::geom_hline(yintercept=p_nums[[2]][[pup_point_num]], colour="darkred") +
        ggplot2::geom_vline(xintercept=p_nums[[1]][[pup_point_num]], colour="black") +
        ggplot2::geom_hline(yintercept=p_nums[[3]][[pdown_point_num]], colour="darkblue")

    retlist <- list(
        "up" = up_plot,
        "down" = down_plot,
        "p" = p_plot,
        "up_data" = up_nums,
        "down_data" = down_nums,
        "p_data" = p_nums)
    return(retlist)
}

#' Try out a few experimental models and return a likely working option.
#'
#' The _pairwise family of functions all demand an experimental model.  This tries to choose a
#' consistent and useful model for all for them.  This does not try to do multi-factor, interacting,
#' nor dependent variable models, if you want those do them yourself and pass them off as alt_model.
#'
#' @param conditions Factor of conditions in the putative model.
#' @param batches Factor of batches in the putative model.
#' @param model_batch Try to include batch in the model?
#' @param model_cond Try to include condition in the model? (Yes!)
#' @param model_intercept Use an intercept model instead of cell-means?
#' @param intercept Choose an intercept for the model as opposed to 0.
#' @param reverse  Reverse condition/batch in the model?  This shouldn't/doesn't matter but I wanted
#'     to test.
#' @param alt_model Use your own model.
#' @param alt_string String describing an alternate model.
#' @return List including a model matrix and strings describing cell-means and intercept models.
choose_model <- function(conditions, batches, model_batch=TRUE,
                         model_cond=TRUE, model_intercept=TRUE,
                         alt_model=NULL, alt_string=NULL,
                         intercept=0, reverse=FALSE) {
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    ## Make a model matrix which will have one entry for
    ## each of the condition/batches
    ## It would be much smarter to generate the models in the following if() {} blocks
    ## But I have it in my head to eventually compare results using different models.
    cond_int_string <- "~ 0 + condition"
    cond_int_model <- stats::model.matrix(~ 0 + conditions)
    batch_int_string <- "~ 0 + batch"
    batch_int_model <- try(stats::model.matrix(~ 0 + batches), silent=TRUE)
    condbatch_int_string <- "~ 0 + condition + batch"
    condbatch_int_model <- try(stats::model.matrix(~ 0 + conditions + batches), silent=TRUE)
    batchcond_int_string <- "~ 0 + batch + condition"
    batchcond_int_model <- try(stats::model.matrix(~ 0 + batches + conditions), silent=TRUE)
    cond_noint_string <- "~ condition"
    cond_noint_model <- try(stats::model.matrix(~ conditions), silent=TRUE)
    batch_noint_string <- "~ batch"
    batch_noint_model <- try(stats::model.matrix(~ batches), silent=TRUE)
    condbatch_noint_string <- "~ condition + batch"
    condbatch_noint_model <- try(stats::model.matrix(~ conditions + batches), silent=TRUE)
    batchcond_noint_string <- "~ batch + condition"
    batchcond_noint_model <- try(stats::model.matrix(~ batches + conditions), silent=TRUE)
    noint_model <- NULL
    int_model <- NULL
    noint_string <- NULL
    int_string <- NULL
    including <- NULL
    if (is.null(model_batch)) {
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    } else if (isTRUE(model_cond) & isTRUE(model_batch)) {
        if (class(condbatch_int_model) == "try-error") {
            message("The condition+batch model failed.  Does your experimental design support both condition and batch?")
            message("Using only a conditional model.")
            int_model <- cond_int_model
            noint_model <- cond_noint_model
            int_string <- cond_int_string
            noint_string <- cond_noint_string
            including <- "condition"
        } else if (isTRUE(reverse)) {
            int_model <- batchcond_int_model
            noint_model <- batchcond_noint_model
            int_string <- batchcond_int_string
            noint_string <- batchcond_noint_string
            including <- "batch+condition"
        } else {
            int_model <- condbatch_int_model
            noint_model <- condbatch_noint_model
            int_string <- condbatch_int_string
            noint_string <- condbatch_noint_string
            including <- "condition+batch"
        }
    } else if (class(model_batch) == "numeric" | class(model_batch) == "matrix") {
        message("Including batch estimates from sva/ruv/pca in the model.")
        int_model <- stats::model.matrix(~ 0 + conditions + model_batch)
        noint_model <- stats::model.matrix(~ conditions + model_batch)
        int_string <- condbatch_int_string
        noint_string <- condbatch_noint_string
        including <- "condition+batchestimate"
    } else if (isTRUE(model_cond)) {
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    } else if (isTRUE(model_batch)) {
        int_model <- batch_int_model
        noint_model <- batch_noint_model
        int_string <- batch_int_string
        noint_string <- batch_noint_string
        including <- "batch"
    } else {
        ## Default to the conditional model
        int_model <- cond_int_model
        noint_model <- cond_noint_model
        int_string <- cond_int_string
        noint_string <- cond_noint_string
        including <- "condition"
    }

    tmpnames <- colnames(int_model)
    tmpnames <- gsub("data[[:punct:]]", "", tmpnames)
    tmpnames <- gsub("-", "", tmpnames)
    tmpnames <- gsub("+", "", tmpnames)
    ## The next lines ensure that conditions/batches which are all numeric will not cause weird errors for contrasts
    ## Ergo, if a condition is something like '111', now it will be 'c111'
    ## Similarly, a batch '01' will be 'b01'
    tmpnames <- gsub("^conditions(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub("^batches(\\d+)$", replacement="b\\1", x=tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    tmpnames <- gsub("batches", "", tmpnames)
    colnames(int_model) <- tmpnames

    tmpnames <- colnames(noint_model)
    tmpnames <- gsub("data[[:punct:]]", "", tmpnames)
    tmpnames <- gsub("-", "", tmpnames)
    tmpnames <- gsub("+", "", tmpnames)
    ## The next lines ensure that conditions/batches which are all numeric will not cause weird errors for contrasts
    ## Ergo, if a condition is something like '111', now it will be 'c111'
    ## Similarly, a batch '01' will be 'b01'
    tmpnames <- gsub("conditions^(\\d+)$", replacement="c\\1", x=tmpnames)
    tmpnames <- gsub("batches^(\\d+)$", replacement="b\\1", x=tmpnames)
    tmpnames <- gsub("conditions", "", tmpnames)
    tmpnames <- gsub("batches", "", tmpnames)
    colnames(noint_model) <- tmpnames

    chosen_model <- NULL
    chosen_string <- NULL
    if (isTRUE(model_intercept)) {
        message("Choosing the intercept containing model.")
        chosen_model <- int_model
        chosen_string <- int_string
    } else {
        chosen_model <- noint_model
        chosen_string <- noint_string
    }

    if (!is.null(alt_model)) {
        fun_model <- alt_model
        fun_string <- alt_string
        including <- "alt"
    }

    retlist <- list(
        "int_model" = int_model,
        "noint_model" = noint_model,
        "int_string" = int_string,
        "noint_string" = noint_string,
        "chosen_model" = chosen_model,
        "chosen_string" = chosen_string,
        "including" = including)
    return(retlist)
}

#' Choose a suitable data set for Edger/DESeq
#'
#' The _pairwise family of functions all demand data in specific formats.
#' This tries to make that consistent.
#'
#' @param input Expt input.
#' @param force Force non-standard data
#' @param ... More options for future expansion
#' @return List the data, conditions, and batches in the data.
choose_dataset <- function(input, force=FALSE, ...) {
    arglist <- list(...)
    input_class <- class(input)[1]
    ## I think I would like to make this function smarter so that it will remove the log2 from transformed data.

    if (input_class == "expt") {
        conditions <- input[["conditions"]]
        batches <- input[["batches"]]
        data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))

        ## As I understand it, EdgeR fits a binomial distribution
        ## and expects data as integer counts, not floating point nor a log2 transformation
        ## Thus, having the 'normalization' state set to something other than 'raw' is a likely
        ## violation of its stated preferred/demanded input.  There are of course ways around this
        ## but one should not take them lightly, perhaps never.
        if (!is.null(input[["state"]][["normalization"]])) {
            ## These if statements may be insufficient to check for the appropriate input for deseq.
            if (isTRUE(force)) {
                ## Setting force to TRUE allows one to round the data to fool edger into accepting it
                ## This is a pretty terrible thing to do
                warning("About to round the data, this is a")
                warning("pretty terrible thing to do.")
                warning("But if you, like me, want to see")
                warning("what happens when you put")
                warning("non-standard data into deseq,")
                warning("then here you go.")
                data <- round(data)
                data[ data < 0] <- 0
            } else if (input[["state"]][["normalization"]] != "raw" |
                       (!is.null(input[["state"]][["transform"]]) & input[["state"]][["transform"]] != "raw")) {
                ## This makes use of the fact that the order of operations in the normalization function is static.
                ## filter->normalization->convert->batch->transform.
                ## Thus, if the normalized state is not raw, we can look back either to the filtered or original data
                ## The same is true for the transformation state.
                if (input[["state"]][["filter"]] == "raw") {
                    message("EdgeR/DESeq expect raw data as input, reverting to the count filtered data.")
                    data <- input[["normalized"]][["intermediate_counts"]][["filter"]][["count_table"]]
                    if (is.null(data)) {
                        data <- input[["normalized"]][["intermediate_counts"]][["original"]]
                    }
                } else {
                    message("EdgeR/DESeq expect raw data as input, reverting to the original expressionset.")
                    data <- Biobase::exprs(input[["original_expressionset"]])
                }
            } else {
                message("The data should be suitable for EdgeR/DESeq.")
                message("If EdgeR/DESeq freaks out, check the state of the count table and ensure that it is in integer counts.")
            }
            ## End testing if normalization has been performed
        }
    } else {
        data <- as.data.frame(input)
    }
    retlist <- list(
        "conditions" = conditions,
        "batches" = batches,
        "data" = data)
    return(retlist)
}

## EOF
