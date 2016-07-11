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
basic_pairwise <- function(input, design=NULL, force=FALSE, ...) {
    message("Starting basic pairwise comparison.")
    input_class <- class(input)[1]
    arglist <- list(...)
    norm <- "quant"
    if (!is.null(arglist[["norm"]])) {
        norm <- arglist[["norm"]]
    }
    convert <- "cbcbcpm"
    if (!is.null(arglist[["convert"]])) {
        convert <- arglist[["convert"]]
    }
    transform <- "log2"
    if (!is.null(arglist[["transform"]])) {
        transform <- arglist[["transform"]]
    }
    filter <- FALSE
    if (!is.null(arglist[["filter"]])) {
        filter <- arglist[["filter"]]
    }

    if (input_class == 'expt') {
        design <- input[["design"]]
        conditions <- input[["conditions"]]
        conditions <- gsub(pattern="^(\\d+)$", replacement="c\\1", x=conditions)
        batches <- input[["batches"]]
        batches <- gsub(pattern="^(\\d+)$", replacement="b\\1", x=batches)
        data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))
        if (!is.null(input[["state"]])) {
            if (isTRUE(force)) {
                warning("This function really does require some sort of
 normalization to be applied to the data, but if you insist,
 I will not mess with it.")
            } else if (input[["state"]][["normalization"]] == "raw" &
                       input[["state"]][["conversion"]] == "raw" &
                       input[["state"]][["transform"]] == "raw") {
                message("This basic pairwise function assumes log2, converted, normalized counts, normalizing now.")
                input <- suppressMessages(normalize_expt(input, norm=norm,
                                                         convert=convert,
                                                         transform=transform,
                                                         filter=filter))
                data <- as.data.frame(Biobase::exprs(input[["expressionset"]]))
            } else if (input[["state"]][["transform"]] == "raw") {
                message("This basic pairwise function assumes log2 counts, transforming now.")
                data <- suppressMessages(transform_counts(data, transform="log2")[["count_table"]])
            } else {
                message("The data appears to have been at least transformed, leaving it alone.")
            }
        } else {  ## There is no 'state' associated with this data.
            message("Cannot tell if this data has been normalized, leaving it alone.")
        }
    } else {
        data <- as.data.frame(input)
    }
    types <- levels(as.factor(conditions))
    num_conds <- length(types)
    ## These will be filled with num_conds columns and numRows(input) rows.
    median_table <- data.frame()
    variance_table <- data.frame()
    ## First use conditions to rbind a table of medians by condition.
    message("Basic step 1/3: Creating median and variance tables.")
    for (c in 1:num_conds) {
        condition_name <- types[c]
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

    for (c in 1:lenminus) {
        c_name <- types[c]
        nextc <- c + 1
        for (d in nextc:length(types)) {
            num_done <- num_done + 1
            d_name <- types[d]
            message(paste0("Basic step 2/3: ", num_done, "/", num_comparisons, ": Performing log2 subtraction: ", d_name, "_vs_", c_name))
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
    ## If only 1 comparison was performed, the resulting data structure never gets coerced into a data frame
    ## therefore I am performing this check which, if a single comparison was done, adds a second column,
    ## performs the coercion, then strips it away.  This is probably a stupid way of doing what I want.
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
        fc_table <- data.frame(numerator_median=median_table[[numerator]],
                               denominator_median=median_table[[denominator]],
                               numerator_var=variance_table[[numerator]],
                               denominator_var=variance_table[[denominator]],
                               t=t_column,
                               p=p_column,
                               logFC=fc_column)
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
        input_data=data, conditions_table=table(conditions),
        conditions=conditions, all_pairwise=comparisons,
        all_tables=all_tables, medians=median_table,
        variances=variance_table)
    return(retlist)
}

## EOF
