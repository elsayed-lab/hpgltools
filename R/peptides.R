## Given a set of nucleotides/amino acids, use the Peptides
## package to create a matrix of scores describing them.

## For the moment, I am going to just use an amino acid fasta file.
## I think I would like to be able to take an expressionset, which
## would obviously only work for non-spliced organisms -- unless the
## annotations include the introns...
score_amino_acids <- function(file, sanitize_names = TRUE) {
  ## In this case, I have the cds sequences from L.major, so I will
  ## translate them.
  amino_acids <- Biostrings::readAAStringSet(file, "fasta")
  if (grepl(x = as.character(amino_acids[[1]]), pattern = "^M\\|L")) {
    ## Then let us assume this is an amino acid sequence.
    message("This appears to be an amino acid fasta file.")
  } else {
    message("This appears to be a DNA fasta file.")
    sequences <- Biostrings::readDNAStringSet(file, "fasta")
    amino_acids <- Biostrings::translate(sequences, if.fuzzy.codon="solve")
  }

  aa_names <- names(amino_acids)
  aa_sequences <- as.character(amino_acids)
  if (isTRUE(sanitize_names)) {
    aa_names <- gsub(x = aa_names, pattern = "^(\\S+)?\\s+.*$", replacement = "\\1")
  }
  names(aa_sequences) <- aa_names
  aa_sequences <- toupper(aa_sequences)
  ## These metrics, understandably, do not understand *
  aa_sequences <- gsub(x = aa_sequences, pattern = "\\*|X", replacement = "")

  ## The peptides R package provides lots of fun metrics.
  metrics <- data.frame(rownames = aa_names)

  ## Apparently there is a correlation between hydrophobicity and
  ## hydrophobic moment which is indicative of the likelihood that a
  ## segment of a protein is globular, transmembrane, or not.
  ## As a result, the following call will give back a rolling vector
  ## (11 amino acids by default) of this value.
  ## mem_data <- Peptides::membpos(aa_sequences)

  ## This might be a neat additional metrics along  with TMHMM.
  composition_lst <- Peptides::aaComp(aa_sequences)  ## Number and
  ## %molarity of tiny,
  composition_df <- t(data.frame(
      row.names = c("Tiny", "Small", "Aliphatic", "Aromatic",
                    "NonPolar", "Polar", "Charged", "Basic", "Acidic")))
  ## Agreed, this is gross, but the list data structure returned by
  ## Peptides is also a bit gross, and I wanted something which is
  ## legible.  The key observation is that each element returned by
  ## Peptides is a 2 column list of which I only want the second
  ## column, the 'Mole%'.
  for (aa_comp in composition_lst) {
    composition_df <- rbind(composition_df, t(aa_comp)[2, ])
  }
  ## small, aliphatic, aromatic, nonpolar, polar, charged, basic and
  ## acidic.
  metrics <- cbind(metrics, composition_df)

  ## Let us consider a couple of comparisons among these compositions
  ## E.g. the normalized ratio of polar/non, tiny/large, basic/acidic

  metrics[["aliphatic_index"]] <- Peptides::aIndex(aa_sequences)

  ## autoCorrelation: The Cruciani et al auto-correlation index.
  ## I do not know what this is, lets look it up!
  ## tt <- Peptides::autoCorrelation(aa_sequences, lag = 1,
  ##                                 property = AAdata$Hydrophobicity$KyteDoolittle)
  ## There is also a covariance index and crossCovariance (wtf?)
  ## tt <- Peptides::crossCovariance(aa_sequences, lag = 1, ...)

  ## I am reasonably certain this gives the blosum metric of 'rarity' for each of the
  ## various blosum matrics.  Not likely of interest for what I am doing, but neat!!
  ## tt <- Peptides::blosumIndices(aa_sequences)

  ## Calculate the bowman potential protein interaction index for each sequence.
  metrics[["bowman"]] <- Peptides::boman(aa_sequences)

  ## Calculate the charge of a protein using the Henderson-Hasselbalch equation
  ## according to Moore, D. S. from 1985.  There are some pH scales
  ## one can play with to make the result fit reality to some degree.
  metrics[["charge"]] <- Peptides::charge(aa_sequences)

  ## FASGAI vectors are 6 values calculated to reflect hydrophobicity,
  ## alpha/turn properties, bulky-ness, composition, flexibility, and
  ## electronic properties.  Neat!
  ## tt <- Peptides::fasgaiVectors(aa_sequences)

  ## Calculate the hydrophobicity moment of amino acids.
  ## This seems a bit slow.
  metrics[["hydrophobicity_moment"]] <- Peptides::hmoment(aa_sequences)
  metrics[["hydrophobicity"]] <- Peptides::hydrophobicity(aa_sequences)

  ## The instaIndex takes a long time, but is super-neat:
  ## This function calculates the instability index proposed by
  ## Guruprasad (1990). This index predicts the stability of a protein
  ## based on its amino acid composition, a protein whose instability
  ## index is smaller than 40 is predicted as stable, a value above 40
  ## predicts that the protein may be unstable.
  metrics[["instability"]] <- Peptides::instaIndex(aa_sequences)

  ## tt <- Peptides::kideraFactors(aa_sequences)
  ## This is cool, it gives the average value of the Kidera factors, which are:
  ## KF1: Helix/bend preference
  ## KF2: Side-chain size
  ## KF3: Extended structure preference (?)
  ## KF4: Hydrophobicity
  ## KF5: Double-bend preference
  ## KF6: Partial specific volume
  ## KF7: Flat extended preference
  ## KF8: Occurrence in alpha region (?meaning likelihood in alpha helices?)
  ## KF9: pK-C
  ## KF10: Surrounding hydrophobicity

  metrics[["expected_mz"]] <- Peptides::mz(aa_sequences)
  metrics[["isoelectric_point"]] <- Peptides::pI(aa_sequences)

  ## This looks like a bit of a broader pass of statistics than
  ## the FASGAI vectors, but similar in concept.
  ## tt <- Peptides::stScales(aa_sequences)

  rownames(metrics) <- metrics[["rownames"]]
  metrics[["rownames"]] <- NULL

  return(metrics)
}
