#' Use ape to generate a distance based nj tree from fasta files.
#'
#' I was thinking that a standardized version of this might be useful
#' for Theresa's recent exploration of variants in her data.
#'
#' @param directory Directory of fasta genomes.
#' @param root Species ID to place at the root of the tree.
#' @return List containing the phylogeny and some other stuff.
#' @export
genomic_sequence_phylo <- function(directory, root = NULL) {
  files <- list.files(directory, pattern = "\\.fasta$|\\.fa$|\\.fsa$")
  sequence_vectors <- list()
  sequence_set <- c()
  count <- 0
  sequence_names <- c()
  for (f in files) {
    shortened <- gsub(x = f, pattern = "\\.fasta$", replacement = "")
    sequence_names <- c(sequence_names, shortened)
    count <- count + 1
    file_path <- file.path(directory, f)
    message("Reading ", file_path)
    raw <- Rsamtools::FaFile(file_path)
    seq <- Biostrings::getSeq(raw)
    test <- as.vector(seq)
    single_vector <- ""
    for (v in test) {
      single_vector <- paste0(single_vector, v)
    }
    sequence_vectors[[shortened]] <- single_vector
    sequence_set <- c(sequence_set, single_vector)
  } ## End iterating over the fasta files.
  dnastring_test <- Biostrings::DNAStringSet(sequence_set)
  names(dnastring_test) <- sequence_names
  dnastring_dnabin <- ape::as.DNAbin(dnastring_test)
  test_dist <- kmer::kdistance(dnastring_dnabin, k = 7)
  test_phy <- ape::nj(test_dist)
  if (!is.null(root)) {
    if (root %in% sequence_names) {
      test_phy <- ape::root(test_phy, root)
    } else {
      warning("The provided root ", root, " is not in the dataset.")
      message("Here are the available samples: ")
      print(sequence_names)
    }
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
