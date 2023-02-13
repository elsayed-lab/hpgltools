genomic_sequence_phylo <- function(directory, root = NULL) {
  files <- list.files(directory, pattern = "\\.fasta$")
  count_lst <- list()
  sequence_vectors <- list()
  sequence_set <- c()
  filenames <- c()
  count <- 0
  for (f in files) {
    count <- count + 1
    message("Reading ", filename)
    raw <- Rsamtools::FaFile(f)
    seq <- Biostrings::getSeq(raw)
    test <- as.vector(seq)
    single_vector <- ""
    for (v in test) {
      single_vector <- paste0(single_vector, v)
    }
    sequence_vectors[[filename]] <- single_vector
    sequence_set <- c(sequence_set, single_vector)
    filenames <- c(filenames, filename)
  }
  dnastring_test <- Biostrings::DNAStringSet(sequence_set)
  shortened <- gsub(x=filenames, pattern="\\.fasta$", replacement="")
  names(dnastring_test) <- shortened
  dnastring_dnabin <- ape::as.DNAbin(dnastring_test)
  test_dist <- kmer::kdistance(dnastring_dnabin, k=7)
  test_phy <- ape::nj(test_dist)
  if (!is.null(root)) {
    test_phy <- ape::root(test_phy, root)
  }
  test_phy <- ape::ladderize(test_phy)

  test_dnd <- as.dendrogram(test_phy)
  test_phylo <- ape::as.phylo(test_dnd) %>%
    ape::compute.brlen()

  retlist <- list(
    "dnd" = test_dnd,
    "phylo" = test_phylo,
    "phy" = test_phy)
  return(retlist)
}
