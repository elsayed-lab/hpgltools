## Time-stamp: <Mon Apr 25 17:17:20 2016 Ashton Trey Belew (abelew@gmail.com)>

#' A caller for different low-count filters
#'
#' @param count_table  some counts to filter
#' @param type   Filtering method to apply (cbcb, pofa, kofa, cv right now)
#' @param p  For pofa()
#' @param A   For pofa()
#' @param k   For kofa()
#' @param cv_min   For cv()
#' @param cv_max   For cv()
#' @param thresh   Minimum threshold across samples for cbcb
#' @param min_samples   Minimum number of samples for cbcb
#' @return a data frame of lowfiltered counts
#' @seealso \pkg{genefilter}
#' @examples
#' \dontrun{
#' new <- lowfilter_counts(old)
#' }
#' @export
lowfilter_counts <- function(count_table, type='cbcb', p=0.01, A=1, k=1,
                             cv_min=0.01, cv_max=1000, thresh=2, min_samples=2) {
    if (tolower(type) == 'povera') {
        type <- 'pofa'
    } else if (tolower(type) == 'kovera') {
        type <- 'kofa'
    }
    lowfiltered_counts <- NULL
    if (type == 'cbcb') {
        lowfiltered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                                 min_samples=min_samples)
    } else if (type == 'pofa') {
        lowfiltered_counts <- genefilter_pofa_counts(count_table, p=p, A=A)
    } else if (type == 'kofa') {
        lowfiltered_counts <- genefilter_kofa_counts(count_table, k=k, A=A)
    } else if (type == 'cv') {
        lowfiltered_counts <- genefilter_cv_counts(count_table, cv_min=cv_min,
                                                   cv_max=cv_max)
    } else {
        lowfiltered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                                 min_samples=min_samples)
    }
    return(lowfiltered_counts)
}

#' Filter low-count genes from a data set.
#'
#' This was a function written by Kwame Okrah and perhaps also Laura Dillon to remove low-count genes.  It drops genes based on a threshold and number of samples.
#'
#' @param count_table  a data frame of (pseudo)counts by sample.
#' @param threshold    lower threshold of counts for each gene.
#' @param min_samples    minimum number of samples
#' @return dataframe of counts without the low-count genes
#' @seealso log2CPM which this uses to decide what to keep
#' @examples
#' \dontrun{
#' filtered_table <- cbcb_filter_counts(count_table)
#' }
#' @export
cbcb_filter_counts <- function(count_table, threshold=2, min_samples=2) {
    ## I think having a log2cpm here is kind of weird, because the next step in processing is to cpm the data.
    ##cpms = 2^log2CPM(counts, lib.size=lib.size)$y
    ## cpms = 2^hpgl_log2cpm(counts)
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        keep <- rowSums(Biobase::exprs(count_table) > threshold) >= min_samples
    } else {
        keep <- rowSums(count_table > threshold) >= min_samples
    }

    count_table <- count_table[keep,]

    message(sprintf("Removing %d low-count genes (%d remaining).",
                    num_before - nrow(count_table), nrow(count_table)))

    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using genefilter's pOverA()
#'
#' I keep thinking this function is pofa... oh well.
#'
#' @param count_table  input data frame of counts by sample
#' @param p   a minimum proportion of each gene's counts/sample to be greater than a minimum(A)
#' @param A   the minimum number of counts in the above proportion
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{pOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_pofa_counts(count_table)
#' }
#' @export
genefilter_pofa_counts <- function(count_table, p=0.01, A=100) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::pOverA(p=p, A=A)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    removed <- num_before - nrow(count_table)
    message(paste0("Removing ", removed, " low-count genes (", nrow(count_table), " remaining)."))

    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter genes from a dataset outside a range of variance
#'
#' @param count_table  input data frame of counts by sample
#' @param cv_min   a minimum coefficient of variance
#' @param cv_max   guess
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{kOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#' filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_cv_counts <- function(count_table, cv_min=0.01, cv_max=1000) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::cv(cv_min, cv_max)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    message(sprintf("Removing %d low-count genes (%d remaining).",
                    num_before - nrow(count_table), nrow(count_table)))
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using genefilter's kOverA()
#'
#' @param count_table input data frame of counts by sample
#' @param k   a minimum number of samples to have >A counts
#' @param A   the minimum number of counts for each gene's sample in kOverA()
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{kOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_kofa_counts <- function(count_table, k=1, A=1) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::kOverA(k=k, A=A)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    message(sprintf("Removing %d low-count genes (%d remaining).",
                    num_before - nrow(count_table), nrow(count_table)))
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using filterCounts() originally from cbcbSEQ
#'
#' @param count_table  input data frame of counts by sample
#' @param thresh   lower threshold of counts (default: 4)
#' @param min_samples   minimum number of samples (default: 2)
#' @param libsize the pre-quantile library sizes
#' @return dataframe of counts without the low-count genes
#' @seealso log2CPM which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = cbcb_lowfilter_counts(count_table)
#' }
#' @export
cbcb_lowfilter_counts <- function(count_table, thresh=2,
                                  min_samples=2, libsize=NULL) {
    original_dim <- dim(count_table)
    cpms <- 2 ^ hpgl_log2cpm(counts, lib.size=libsize)
    keep <- rowSums(cpms > thresh) >= min_samples
    counts_table <- as.matrix(counts[keep, ])
    following_dim <- dim(count_table)
    lost_rows <- original_dim[1] - following_dim[1]
    message(paste0("Low count filtering cost: ", lost_rows, " gene(s)."))
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

