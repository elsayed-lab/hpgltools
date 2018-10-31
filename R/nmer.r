#' Count n-mers in a given data set using Biostrings
#'
#' This just calls PDict() and vcountPDict() on a sequence database given a
#' pattern and number of mismatches.  This may be used by divide_seq()
#' normalization.
#'
#' @param genome  Sequence database, genome in this case.
#' @param pattern  Count off this string.
#' @param mismatch  How many mismatches are acceptable?
#' @return Set of counts by sequence.
count_nmer <- function(genome, pattern="ATG", mismatch=0) {
    seq_obj <- Biostrings::getSeq(genome)
    dict <- Biostrings::PDict(pattern, max.mismatch=mismatch)
    result <- Biostrings::vcountPDict(dict, seq_obj)
    return(result)
}

## EOF
