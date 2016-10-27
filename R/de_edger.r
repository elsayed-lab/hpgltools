#' Set up a model matrix and set of contrasts to do pairwise comparisons using EdgeR.
#'
#' This function performs the set of possible pairwise comparisons using EdgeR.
#'
#' @param input Dataframe/vector or expt class containing data, normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the experimental model?
#' @param model_batch Include batch in the model?  In most cases this is a good thing(tm).
#' @param model_intercept Use cell means or intercept?
#' @param alt_model Alternate experimental model to use?
#' @param extra_contrasts Add some extra contrasts to add to the list of pairwise contrasts.
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df Annotation information to the data tables?
#' @param force Force edgeR to accept inputs which it should not have to deal with.
#' @param edger_method  I found a couple/few ways of doing edger in the manual, choose with this.
#' @param ... The elipsis parameter is fed to write_edger() at the end.
#' @return List including the following information:
#'   contrasts = The string representation of the contrasts performed.
#'   lrt = A list of the results from calling glmLRT(), one for each contrast.
#'   contrast_list = The list of each call to makeContrasts()
#'   I do this to avoid running into the limit on # of contrasts addressable by topTags()
#'   all_tables = a list of tables for the contrasts performed.
#' @seealso \pkg{edgeR} \code{\link[edgeR]{topTags}} \code{\link[edgeR]{glmLRT}}
#'   \code{\link{make_pairwise_contrasts}} \code{\link[edgeR]{DGEList}}
#'   \code{\link[edgeR]{calcNormFactors}} \code{\link[edgeR]{estimateTagwiseDisp}}
#'   \code{\link[edgeR]{estimateCommonDisp}} \code{\link[edgeR]{estimateGLMCommonDisp}}
#'   \code{\link[edgeR]{estimateGLMTrendedDisp}} \code{\link[edgeR]{glmFit}}
#' @examples
#' \dontrun{
#'  pretend = edger_pairwise(data, conditions, batches)
#' }
#' @export
edger_pairwise <- function(input=NULL, conditions=NULL,
                           batches=NULL, model_cond=TRUE,
                           model_batch=TRUE, model_intercept=TRUE,
                           alt_model=NULL, extra_contrasts=NULL,
                           annot_df=NULL, force=FALSE,
                           edger_method="long", ...) {
    arglist <- list(...)
    if (!is.null(arglist[["input"]])) {
        input <- arglist[["input"]]
    }
    if (!is.null(arglist[["conditions"]])) {
        conditions <- arglist[["conditions"]]
    }
    if (!is.null(arglist[["batches"]])) {
        batches <- arglist[["batches"]]
    }
    if (!is.null(arglist[["model_cond"]])) {
        model_cond <- arglist[["model_cond"]]
    }
    if (!is.null(arglist[["model_batch"]])) {
        model_batch <- arglist[["model_batch"]]
    }
    if (!is.null(arglist[["model_intercept"]])) {
        model_intercept <- arglist[["model_intercept"]]
    }
    if (!is.null(arglist[["alt_model"]])) {
        alt_model <- arglist[["alt_model"]]
    }
    if (!is.null(arglist[["extra_contrasts"]])) {
        extra_contrasts <- arglist[["extra_contrasts"]]
    }
    if (!is.null(arglist[["annot_df"]])) {
        annot_df <- arglist[["annot_df"]]
    }
    if (!is.null(arglist[["force"]])) {
        force <- arglist[["force"]]
    }
    if (!is.null(arglist[["edger_method"]])) {
        edger_method <- arglist[["edger_method"]]
    }
    test_type <- "lrt"
    if (!is.null(arglist[["test_type"]])) {
        test_type <- arglist[["test_type"]]
    }
    message("Starting edgeR pairwise comparisons.")
    input_data <- choose_binom_dataset(input, force=force)
    design <- Biobase::pData(input[["expressionset"]])
    conditions <- input_data[["conditions"]]
    batches <- input_data[["batches"]]
    data <- input_data[["data"]]

    model_choice <- choose_model(input, conditions, batches,
                                 model_batch=model_batch,
                                 model_cond=model_cond,
                                 model_intercept=model_intercept,
                                 alt_model=alt_model, ...)
    model_including <- model_choice[["including"]]
    if (class(model_choice[["model_batch"]]) == "matrix") {
        model_batch <- model_choice[["model_batch"]]
    }
    model_data <- model_choice[["chosen_model"]]

    ## I have a strong sense that the most recent version of edgeR changed its dispersion estimate code
    ## Here is a note from the user's guide, which may have been there previously and I merely did not notice:
    ## To estimate common dispersion, trended dispersions and tagwise dispersions in one run
    ## y <- estimateDisp(y, design)
    raw <- edgeR::DGEList(counts=data, group=conditions)
    message("EdgeR step 1/9: normalizing data.")
    norm <- edgeR::calcNormFactors(raw)
    final_norm <- NULL
    if (edger_method == "short") {
        message("EdgeR steps 2 through 6/9: All in one!")
        final_norm <- edgeR::estimateDisp(norm, design=model_data, robust=TRUE)
    } else {
        message("EdgeR step 2/9: Estimating the common dispersion.")
        disp_norm <- edgeR::estimateCommonDisp(norm)
        message("EdgeR step 3/9: Estimating dispersion across genes.")
        tagdisp_norm <- edgeR::estimateTagwiseDisp(disp_norm)
        message("EdgeR step 4/9: Estimating GLM Common dispersion.")
        glm_norm <- edgeR::estimateGLMCommonDisp(tagdisp_norm, model_data)
        message("EdgeR step 5/9: Estimating GLM Trended dispersion.")
        glm_trended <- edgeR::estimateGLMTrendedDisp(glm_norm, model_data)
        message("EdgeR step 6/9: Estimating GLM Tagged dispersion.")
        final_norm <- edgeR::estimateGLMTagwiseDisp(glm_trended, model_data)
    }
    cond_fit <- NULL
    if (test_type == "lrt") {
        message("EdgeR step 7/9: Running glmFit, switch to glmQLFit by changing the argument 'test_type'.")
        cond_fit <- edgeR::glmFit(final_norm, design=model_data, robust=TRUE)
    } else {
        message("EdgeR step 7/9: Running glmQLFit, switch to glmFit by changing the argument 'test_type'.")
        cond_fit <- edgeR::glmQLFit(final_norm, design=model_data, robust=TRUE)
    }

    message("EdgeR step 8/9: Making pairwise contrasts.")
    apc <- make_pairwise_contrasts(model_data, conditions,
                                   extra_contrasts=extra_contrasts,
                                   do_identities=FALSE)

    ## This section is convoluted because glmLRT only seems to take up to 7 contrasts at a time.
    ## As a result, I iterate through the set of possible contrasts one at a time and ask for each
    ## separately.
    contrast_list <- list()
    result_list <- list()
    lrt_list <- list()
    sc <- vector("list", length(apc[["names"]]))
    end <- length(apc[["names"]])
    for (con in 1:length(apc[["names"]])) {
        name <- apc[["names"]][[con]]
        message(paste0("EdgeR step 9/9: ", con, "/", end, ": Creating table: ", name, ".")) ## correct
        sc[[name]] <- gsub(pattern=",", replacement="", apc[["all_pairwise"]][[con]])
        tt <- parse(text=sc[[name]])
        ctr_string <- paste0("tt = limma::makeContrasts(", tt, ", levels=model_data)")
        eval(parse(text=ctr_string))
        contrast_list[[name]] <- tt
        lrt_list[[name]] <- NULL
        if (test_type == "lrt") {
            lrt_list[[name]] <- edgeR::glmLRT(cond_fit, contrast=contrast_list[[name]])
        } else {
            lrt_list[[name]] <- edgeR::glmQLFTest(cond_fit, contrast=contrast_list[[name]])
        }
        res <- edgeR::topTags(lrt_list[[name]], n=nrow(data), sort.by="logFC")
        res <- as.data.frame(res)
        res[["logFC"]] <- signif(x=as.numeric(res[["logFC"]]), digits=4)
        res[["logCPM"]] <- signif(x=as.numeric(res[["logCPM"]]), digits=4)
        if (!is.null(res[["LR"]])) {
            res[["LR"]] <- signif(x=as.numeric(res[["LR"]]), digits=4)
        } else if (!is.null(res[["F"]])) {
            res[["F"]] <- signif(x=as.numeric(res[["F"]]), digits=4)
        }
        res[["PValue"]] <- signif(x=as.numeric(res[["PValue"]]), digits=4)
        res[["FDR"]] <- signif(x=as.numeric(res[["FDR"]]), digits=4)
        res[["qvalue"]] <- tryCatch({
            ttmp <- as.numeric(res[["PValue"]])
            ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
            format(x=ttmp, digits=4, scientific=TRUE)
        },
        error=function(cond) {
            message(paste0("The qvalue estimation failed for ", name, "."))
            return(1)
        },
        finally={
        })
        result_list[[name]] <- res
    } ## End for loop
    final <- list(
        "model" = model_data,
        "contrasts" = apc,
        "contrasts_performed" = apc[["names"]],
        "lrt" = lrt_list,
        "contrast_list" = contrast_list,
        "all_tables" = result_list)
    return(final)
}

#' Writes out the results of a edger search using topTags()
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the results() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data Output from topTags().
#' @param adjust Pvalue adjustment chosen.
#' @param n Number of entries to report, 0 says do them all.
#' @param coef Which coefficients/contrasts to report, NULL says do them all.
#' @param workbook Excel filename into which to write the data.
#' @param excel Write an excel workbook?
#' @param csv Write out csv files of the tables?
#' @param annot_df Optional data frame including annotation information to include with the tables.
#' @return List of data frames comprising the toptable output for each coefficient, I also added a
#'     qvalue entry to these toptable() outputs.
#' @seealso \link[edger]{toptable} \link{write_xls}
#' @examples
#' \dontrun{
#'  finished_comparison = topTags(edger_output)
#'  data_list = write_edger(finished_comparison, workbook="excel/edger_output.xls")
#' }
#' @export
write_edger <- function(data, adjust="fdr", n=0, coef=NULL, workbook="excel/edger.xls",
                       excel=FALSE, csv=FALSE, annot_df=NULL) {
    testdir <- dirname(workbook)

    ## Figure out the number of genes if not provided
    if (n == 0) {
        n <- nrow(data[["coefficients"]])
    }

    ## If specific contrast(s) is/are not requested, get them all.
    if (is.null(coef)) {
        coef <- colnames(data[["contrasts"]])
    } else {
        coef <- as.character(coef)
    }
    return_data <- list()
    end <- length(coef)
    for (c in 1:end) {
        comparison <- coef[c]
        message(paste0("Edger step 6/6: ", c, "/", end, ": Creating table: ", comparison, "."))
        data_table <- edger::topTable(data, adjust=adjust, n=n, coef=comparison)
        ## Reformat the numbers so they are not so obnoxious
        ## data_table$logFC <- refnum(data_table$logFC, sci=FALSE)
        ## data_table$AveExpr <- refnum(data_table$AveExpr, sci=FALSE)
        ## data_table$t <- refnum(data_table$t, sci=FALSE)
        ## data_table$P.Value <- refnum(data_table$P.Value)
        ## data_table$adj.P.Val <- refnum(data_table$adj.P.Val)
        ## data_table$B <- refnum(data_table$B, sci=FALSE)
        data_table[["logFC"]] <- signif(x=as.numeric(data_table[["logFC"]]), digits=4)
        data_table[["AveExpr"]] <- signif(x=as.numeric(data_table[["AveExpr"]]), digits=4)
        data_table[["t"]] <- signif(x=as.numeric(data_table[["t"]]), digits=4)
        data_table[["P.Value"]] <- signif(x=as.numeric(data_table[["P.Value"]]), digits=4)
        data_table[["adj.P.Val"]] <- signif(x=as.numeric(data_table[["adj.P.Val"]]), digits=4)
        data_table[["B"]] <- signif(x=as.numeric(data_table[["B"]]), digits=4)
        data_table[["qvalue"]] <- tryCatch({
            ttmp <- as.numeric(data_table[["P.Value"]])
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)[["qvalues"]]
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            message(paste("The qvalue estimation failed for ", comparison, ".", sep=""))
            return(1)
        },
        finally={
        })
        if (!is.null(annot_df)) {
            data_table <- merge(data_table, annot_df, by.x="row.names", by.y="row.names")
            ###data_table = data_table[-1]
        }
        ## This write_xls runs out of memory annoyingly often
        if (isTRUE(excel) | isTRUE(csv)) {
            if (!file.exists(testdir)) {
                dir.create(testdir)
                message(paste0("Creating directory: ", testdir, " for writing excel/csv data."))
            }
        }
        if (isTRUE(excel)) {
            try(write_xls(data=data_table, sheet=comparison, file=workbook, overwritefile=TRUE))
        }
        ## Therefore I will write a csv of each comparison, too
        if (isTRUE(csv)) {
            csv_filename <- gsub(".xls$", "", workbook)
            csv_filename <- paste0(csv_filename, "_", comparison, ".csv")
            write.csv(data_table, file=csv_filename)
        }
        return_data[[comparison]] <- data_table
    }
    return(return_data)
}

## EOF
