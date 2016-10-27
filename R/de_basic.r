#' The simplest possible differential expression method.
#'
#' Perform a pairwise comparison among conditions which takes
#' nothing into account.  It _only_ takes the conditions, a mean value/variance among
#' them, divides by condition, and returns the result.  No fancy nomalizations, no
#' statistical models, no nothing.  It should be the very worst method possible.
#' But, it should also provide a baseline to compare the other tools against, they should
#' all do better than this, always.
#'
#' @param input Count table by sample.
#' @param design Data frame of samples and conditions.
#' @param force Force as input non-normalized data?
#' @param ... Extra options passed to arglist.
#' @return Df of pseudo-logFC, p-values, numerators, and denominators.
#' @seealso \pkg{limma} \pkg{DESeq2} \pkg{edgeR}
#' @examples
#' \dontrun{
#' stupid_de <- basic_pairwise(expt)
#' }
#' @export
basic_pairwise <- function(input=NULL, design=NULL,
                           force=FALSE, ...) {
    arglist <- list(...)
    if (!is.null(arglist[["input"]])) {
        input <- arglist[["input"]]
    }
    if (!is.null(arglist[["design"]])) {
        conditions <- arglist[["design"]]
    }
    if (!is.null(arglist[["force"]])) {
        batches <- arglist[["force"]]
    }
    message("Starting basic pairwise comparison.")
    input_data <- choose_basic_dataset(input, force=force)
    design <- Biobase::pData(input[["expressionset"]])
    conditions <- input_data[["conditions"]]
    batches <- input_data[["batches"]]
    data <- input_data[["data"]]

    conditions <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=conditions)
    batches <- gsub(pattern="^(\\d+)$", replacement="b\\1", x=batches)
    types <- levels(as.factor(conditions))
    num_conds <- length(types)
    ## These will be filled with num_conds columns and numRows(input) rows.
    median_table <- data.frame()
    variance_table <- data.frame()
    ## First use conditions to rbind a table of medians by condition.
    message("Basic step 1/3: Creating median and variance tables.")
    median_colnames <- c()
    for (c in 1:num_conds) {
        condition_name <- types[c]
        median_colnames <- append(median_colnames, condition_name)
        columns <- which(conditions == condition_name)
        if (length(columns) == 1) {
            med <- data.frame(data[, columns])
            var <- as.data.frame(matrix(NA, ncol=1, nrow=nrow(med)))
        } else {
            med_input <- data[, columns]
            med <- data.frame(Biobase::rowMedians(as.matrix(med_input)))
            colnames(med) <- c(condition_name)
            var <- as.data.frame(genefilter::rowVars(as.matrix(med_input)))
            colnames(var) <- c(condition_name)
        }
        if (c == 1) {
            median_table <- med
            variance_table <- var
        } else {
            median_table <- cbind(median_table, med)
            variance_table <- cbind(variance_table, var)
        }
    } ## end creation of median and variance tables.
    colnames(median_table) <- median_colnames
    colnames(variance_table) <- median_colnames
    rownames(median_table) <- rownames(data)
    rownames(variance_table) <- rownames(data)
    ## We have tables of the median values by condition
    ## Now perform the pairwise comparisons
    comparisons <- data.frame()
    tvalues <- data.frame()
    pvalues <- data.frame()
    lenminus <- num_conds - 1
    num_done <- 0
    column_list <- c()
    message("Basic step 2/3: Performing comparisons.")
    num_comparisons <- sum(1:lenminus)

    contrasts_performed <- c()
    for (c in 1:lenminus) {
        c_name <- types[c]
        nextc <- c + 1
        for (d in nextc:length(types)) {
            num_done <- num_done + 1
            d_name <- types[d]
            contrast <- paste0(d_name, "_vs_", c_name)
            contrasts_performed <- append(contrast, contrasts_performed)
            message(paste0("Basic step 2/3: ", num_done, "/", num_comparisons, ": Performing log2 subtraction: ", contrast, "."))
            division <- data.frame(
                median_table[, d] - median_table[, c])
            comparison_name <- paste0(d_name, "_vs_", c_name)
            column_list <- append(column_list, comparison_name)
            colnames(division) <- comparison_name
            ## Lets see if I can make a dirty p-value
            xcols <- which(conditions == c_name)
            ycols <- which(conditions == d_name)
            xdata <- as.data.frame(data[, xcols])
            ydata <- as.data.frame(data[, ycols])

            t_data <- vector("list", nrow(xdata))
            p_data <- vector("list", nrow(xdata))
            for (j in 1:nrow(xdata)) {
                test_result <- try(t.test(xdata[j, ], ydata[j, ]), silent=TRUE)
                if (class(test_result) == 'htest') {
                    t_data[[j]] <- test_result[[1]]
                    p_data[[j]] <- test_result[[3]]
                } else {
                    t_data[[j]] <- 0
                    p_data[[j]] <- 1
                }
            } ## Done calculating cheapo p-values

            if (num_done == 1) {
                comparisons <- division
                tvalues <- t_data
                pvalues <- p_data
            } else {
                comparisons <- cbind(comparisons, division)
                tvalues <- cbind(tvalues, t_data)
                pvalues <- cbind(pvalues, p_data)
            }
        } ## End for each d
    }  ## End for each c

    ## Because of the way I made tvalues/pvalues into a list
    ## If only 1 comparison was performed, the resulting data structure never gets coerced into a
    ## data frame therefore I am performing this check which, if a single comparison was done, adds
    ## a second column, performs the coercion, then strips it away.  This is a stupid way
    ## of doing what I want.
    if (num_done == 1) {
        tvalues <- cbind(tvalues, t_data)
        pvalues <- cbind(pvalues, p_data)
        tvalues <- as.data.frame(tvalues)
        pvalues <- as.data.frame(pvalues)
        tvalues <- tvalues[-1]
        pvalues <- pvalues[-1]
    }
    comparisons[is.na(comparisons)] <- 0
    tvalues[is.na(tvalues)] <- 0
    pvalues[is.na(pvalues)] <- 1
    rownames(comparisons) <- rownames(data)
    rownames(tvalues) <- rownames(data)
    rownames(pvalues) <- rownames(data)
    all_tables <- list()

    message("Basic step 3/3: Creating faux DE Tables.")
    for (e in 1:length(colnames(comparisons))) {
        colname <- colnames(comparisons)[[e]]
        fc_column <- comparisons[, e]
        t_column <- as.numeric(tvalues[, e])
        p_column <- as.numeric(pvalues[, e])
        fc_column[mapply(is.infinite, fc_column)] <- 0
        numer_denom <- strsplit(x=colname, split="_vs_")[[1]]
        numerator <- numer_denom[1]
        denominator <- numer_denom[2]
        fc_table <- data.frame(
            "numerator_median" = median_table[[numerator]],
            "denominator_median" = median_table[[denominator]],
            "numerator_var" = variance_table[[numerator]],
            "denominator_var" = variance_table[[denominator]],
            "t" = t_column,
            "p" = p_column,
            "logFC" = fc_column)
        fc_table[["adjp"]] <- stats::p.adjust(as.numeric(fc_table[["p"]]), method="BH")

        fc_table[["numerator_median"]] <- signif(x=fc_table[["numerator_median"]], digits=4)
        fc_table[["denominator_median"]] <- signif(x=fc_table[["denominator_median"]], digits=4)
        fc_table[["numerator_var"]] <- format(x=fc_table[["numerator_var"]], digits=4, scientific=TRUE)
        fc_table[["denominator_var"]] <- format(x=fc_table[["denominator_var"]], digits=4, scientific=TRUE)
        fc_table[["t"]] <- signif(x=fc_table[["t"]], digits=4)
        fc_table[["p"]] <- format(x=fc_table[["p"]], digits=4, scientific=TRUE)
        fc_table[["adjp"]] <- format(x=fc_table[["adjp"]], digits=4, scientific=TRUE)
        fc_table[["logFC"]] <- signif(x=fc_table[["logFC"]], digits=4)
        rownames(fc_table) <- rownames(data)
        all_tables[[e]] <- fc_table
    }
    message("Basic: Returning tables.")
    names(all_tables) <- colnames(comparisons)
    retlist <- list(
        "input_data" = data,
        "conditions_table" = table(conditions),
        "conditions" = conditions,
        "all_pairwise" = comparisons,
        "all_tables" = all_tables,
        "medians" = median_table,
        "contrasts_performed" = contrasts_performed,
        "variances" = variance_table)
    return(retlist)
}

choose_basic_dataset <- function(input, force=FALSE, ...) {
    arglist <- list(...)
    warn_user <- 0
    conditions <- input[["conditions"]]
    batches <- input[["batches"]]
    data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))
    tran_state <- input[["state"]][["transform"]]
    if (is.null(tran_state)) {
        tran_state <- "raw"
    }
    conv_state <- input[["state"]][["conversion"]]
    ## Note that voom takes care of this for us.
    if (is.null(conv_state)) {
        conv_state <- "raw"
    }
    norm_state <- input[["state"]][["normalization"]]
    if (is.null(norm_state)) {
        norm_state <- "raw"
    }
    filt_state <- input[["state"]][["filter"]]
    if (is.null(filt_state)) {
        filt_state <- "raw"
    }

    ready <- input
    if (isTRUE(force)) {
        message("Leaving the data alone, regardless of normalization state.")
    } else {
        if (filt_state == "raw") {
            message("Basic step 0/3: Filtering data.")
            ready <- sm(normalize_expt(ready, filter="cbcb"))
        }
        if (norm_state == "raw") {
            message("Basic step 0/3: Normalizing data.")
            ready <- sm(normalize_expt(ready, norm="quant"))
        }
        if (conv_state == "raw") {
            message("Basic step 0/3: Converting data.")
            ready <- sm(normalize_expt(ready, convert="cbcbcpm"))
        }

    }
    ## No matter what we do, it must be logged.
    if (tran_state == "raw") {
        message("Basic step 0/3: Transforming data.")
        ready <- sm(normalize_expt(ready, transform="log2"))
    }
    data <- as.data.frame(Biobase::exprs(ready[["expressionset"]]))
    rm(ready)
    retlist <- list(
        "conditions" = conditions,
        "batches" = batches,
        "data" = data)
    return(retlist)
}

#' Writes out the results of a basic search using basic_pairwise()
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the results() output for them in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param data Table from basic_pairwise().
#' @param adjust Pvalue adjustment chosen.
#' @param n Number of entries to report, 0 says do them all.
#' @param coef Which coefficients/contrasts to report, NULL says do them all.
#' @param workbook Excel filename into which to write the data.
#' @param excel Write an excel workbook?
#' @param csv Write out csv files of the tables?
#' @param annot_df Optional data frame including annotation information to include with the tables.
#' @return List of data frames comprising the toptable output for each coefficient, I also added a
#'     qvalue entry to these toptable() outputs.
#' @seealso \link[basic]{toptable} \link{write_xls}
#' @examples
#' \dontrun{
#'  finished_comparison = basic_comparison(basic_output)
#'  data_list = write_basic(finished_comparison, workbook="excel/basic_output.xls")
#' }
#' @export
write_basic <- function(data, adjust="fdr", n=0, coef=NULL, workbook="excel/basic.xls",
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
        message(paste0("Basic step 6/6: ", c, "/", end, ": Creating table: ", comparison, "."))
        data_table <- basic::topTable(data, adjust=adjust, n=n, coef=comparison)
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
        }, error=function(cond) {
            message(paste("The qvalue estimation failed for ", comparison, ".", sep=""))
            return(1)
        }, finally={
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
