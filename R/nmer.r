
count_nmer <- function(pattern, genome, mismatch=0) {
    seq_obj <- Biostrings::getSeq(genome)
    dict <- Biostrings::PDict(pattern, max.mismatch=mismatch)
    result <- Biostrings::vcountPDict(dict, seq_obj)
    return(result)
}
