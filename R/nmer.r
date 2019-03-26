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
#' @export
count_nmer <- function(genome, pattern="ATG", mismatch=0) {
  if (class(genome)[1] == "character") {
    genome <- Rsamtools::FaFile(genome)
  }
  seq_obj <- Biostrings::getSeq(genome)
  dict <- Biostrings::PDict(pattern, max.mismatch=mismatch)
  result <- as.data.frame(Biostrings::vcountPDict(dict, seq_obj))
  rownames(result) <- pattern
  colnames(result) <- names(seq_obj)
  result <- t(result)
  return(result)
}

## EOF
