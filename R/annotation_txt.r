#' Read a csv file and make an annotation data frame.
#'
#' Yay!
#'
#' @param trinotate  CSV of trinotate annotation data.
#' @return  Dataframe of fun data.
#' @export
load_trinotate_annotations <- function(trinotate="reference/trinotate.csv") {
  big_table <- read.csv(trinotate, sep="\t", stringsAsFactors=FALSE)

  split_data <- big_table %>%
    tidyr::separate("sprot_Top_BLASTX_hit",
                    c("blastx_name", "blastx_name2", "blastx_hitlocation",
                      "blastx_identity", "blastx_evalue",
                      "blastx_recname", "blastx_taxonomy"),
                    "\\^") %>%
    tidyr::separate("sprot_Top_BLASTP_hit",
                    c("blastp_name", "blastp_name2", "blastp_hitlocation",
                      "blastp_identity", "blastp_evalue",
                      "blastp_recname", "blastp_taxonomy"),
                    "\\^") %>%
    tidyr::separate("TmHMM",
                    c("tmhmm_expaa", "tmhmm_predicted_helices", "tmhmm_topology"),
                    "\\^") %>%
    tidyr::separate("eggnog",
                    c("eggnog_id", "eggnog_description"),
                    "\\^") %>%
    tidyr::separate("blastx_hitlocation",
                    c("blastx_queryloc", "blastx_hitloc"),
                    "\\,") %>%
    tidyr::separate("blastp_hitlocation",
                    c("blastp_queryloc", "blastp_hitloc"),
                    "\\,") %>%
    tidyr::separate("RNAMMER",
                    c("rrna_subunit", "rrna_subunit_region"),
                    "\\^")

  split_data <- split_data[, c("X.gene_id", "transcript_id",
                               "blastx_name", "blastx_queryloc", "blastx_hitloc",
                               "blastx_identity", "blastx_evalue", "blastx_recname",
                               "blastx_taxonomy", "blastp_name", "blastp_queryloc",
                               "blastp_hitloc", "blastp_identity", "blastp_evalue",
                               "blastp_recname", "blastp_taxonomy",
                               "rrna_subunit", "rrna_subunit_region", "prot_id",
                               "prot_coords", "Pfam", "SignalP",
                               "tmhmm_expaa", "tmhmm_predicted_helices", "tmhmm_topology",
                               "eggnog_id", "eggnog_description", "Kegg",
                               "gene_ontology_blast", "gene_ontology_pfam", "transcript",
                               "peptide")]

  colnames(split_data) <- c("gene_id", "transcript_id",
                            "blastx_name", "blastx_queryloc", "blastx_hitloc",
                            "blastx_identity", "blastx_evalue", "blastx_recname",
                            "blastx_taxonomy", "blastp_name", "blastp_queryloc",
                            "blastp_hitloc", "blastp_identity", "blastp_evalue",
                            "blastp_recname", "blastp_taxonomy",
                            "rrna_subunit", "rrna_subunit_region", "prot_id",
                            "prot_coords", "pfam_data", "signalp_data",
                            "tmhmm_expaa", "tmhmm_helices", "tmhmm_topology",
                            "eggnog_id", "eggnog_description", "kegg_data",
                            "gene_ontology_blast", "gene_ontology_pfam", "transcript_seq",
                            "peptide")

  split_data[["gene_ontology_blast"]] <- gsub(pattern="\\^", replacement="\\,",
                                              x=split_data[["gene_ontology_blast"]])
  split_data[["gene_ontology_blast"]] <- gsub(pattern="\\`", replacement="\\; ",
                                              x=split_data[["gene_ontology_blast"]])
  split_data[["gene_ontology_blast"]] <- gsub(pattern="^\\,", replacement="",
                                              x=split_data[["gene_ontology_blast"]])
  split_data[["kegg_data"]] <- gsub(pattern="\\`", replacement="\\; ",
                                    x=split_data[["kegg_data"]])
  split_data[["tmhmm_expaa"]] <- gsub(pattern="ExpAA=", replacement="",
                                      x=split_data[["tmhmm_expaa"]])
  split_data[["tmhmm_helices"]] <- gsub(pattern="PredHel=", replacement="",
                                        x=split_data[["tmhmm_helices"]])
  split_data[["tmhmm_topology"]] <- gsub(pattern="Topology=", replacement="",
                                         x=split_data[["tmhmm_topology"]])

  split_data[["blastx_identity"]] <- gsub(pattern="%ID", replacement="",
                                          x=split_data[["blastx_identity"]])
  split_data[["blastx_evalue"]] <- gsub(pattern="E:", replacement="",
                                        x=split_data[["blastx_evalue"]])
  split_data[["blastx_queryloc"]] <- gsub(pattern="Q:", replacement="",
                                          x=split_data[["blastx_queryloc"]])
  split_data[["blastx_hitloc"]] <- gsub(pattern="H:", replacement="",
                                        x=split_data[["blastx_hitloc"]])
  split_data[["blastx_recname"]] <- gsub(pattern="RecName: ", replacement="",
                                         x=split_data[["blastx_recname"]])
  split_data[["blastx_recname"]] <- gsub(pattern="Full=", replacement="",
                                         x=split_data[["blastx_recname"]])

  split_data[["blastp_identity"]] <- gsub(pattern="%ID", replacement="",
                                          x=split_data[["blastp_identity"]])
  split_data[["blastp_evalue"]] <- gsub(pattern="E:", replacement="",
                                        x=split_data[["blastp_evalue"]])
  split_data[["blastp_queryloc"]] <- gsub(pattern="Q:", replacement="",
                                          x=split_data[["blastp_queryloc"]])
  split_data[["blastp_hitloc"]] <- gsub(pattern="H:", replacement="",
                                        x=split_data[["blastp_hitloc"]])
  split_data[["blastp_recname"]] <- gsub(pattern="RecName: ", replacement="",
                                         x=split_data[["blastp_recname"]])
  split_data[["blastp_recname"]] <- gsub(pattern="Full=", replacement="",
                                         x=split_data[["blastp_recname"]])

  ## Get rid of empty cells and cells with just '.'
  nas <- is.na(split_data)
  split_data[nas] <- ""
  dots <- split_data == "."
  split_data[dots] <- ""

  ## Now recast the numeric elements.
  split_data[["tmhmm_expaa"]] <- as.numeric(split_data[["tmhmm_expaa"]])
  na_test <- is.na(split_data[["tmhmm_expaa"]])
  split_data[na_test, "tmhmm_expaa"] <- 0
  split_data[["tmhmm_helices"]] <- as.numeric(split_data[["tmhmm_helices"]])
  na_test <- is.na(split_data[["tmhmm_helices"]])
  split_data[na_test, "tmhmm_helices"] <- 0
  split_data[["blastx_identity"]] <- as.numeric(split_data[["blastx_identity"]])
  na_test <- is.na(split_data[["blastx_identity"]])
  split_data[na_test, "blastx_identity"] <- 0
  split_data[["blastx_evalue"]] <- as.numeric(split_data[["blastx_evalue"]])
  na_test <- is.na(split_data[["blastx_evalue"]])
  split_data[na_test, "blastx_evalue"] <- 1
  split_data[["blastp_identity"]] <- as.numeric(split_data[["blastp_identity"]])
  na_test <- is.na(split_data[["blastp_identity"]])
  split_data[na_test, "blastp_identity"] <- 0
  split_data[["blastp_evalue"]] <- as.numeric(split_data[["blastp_evalue"]])
  na_test <- is.na(split_data[["blastp_evalue"]])
  split_data[na_test, "blastp_evalue"] <- 1

  rownames(split_data) <- make.names(split_data[["transcript_id"]], unique=TRUE)
  return(split_data)
}

#' Read a csv file and make a GO data frame.
#'
#' Yay!
#'
#' @param trinotate  CSV of trinotate annotation data.
#' @return  Dataframe of fun data.
#' @export
load_trinotate_go <- function(trinotate="reference/trinotate.csv") {
  big_table <- read.csv(trinotate, sep="\t", stringsAsFactors=FALSE)
  big_table[["length"]] <- stringr::str_length(as.factor(big_table[["transcript"]]))

  go_data <- big_table[, c("X.gene_id", "transcript_id", "gene_ontology_blast",
                           "gene_ontology_pfam", "length")]
  colnames(go_data) <- c("gene_id", "transcript_id", "go_blast", "go_pfam", "length")
  dots <- go_data == "."
  go_data[dots] <- ""

  expanded <- go_data %>% dplyr::mutate("GO"=strsplit(as.character(.data[["go_blast"]]), "`")) %>%
    tidyr::unnest("GO") %>%
    tidyr::separate("GO",
                    c("GO", "GO_ont", "GO_name"),
                    "\\^")

  go_data <- expanded[, c("gene_id", "transcript_id", "GO", "GO_ont", "GO_name")]

  go_table <- data.table::setDT(go_data)
  go_table <- go_table[, c("gene_id", "GO")]
  names(go_table) <- c("ID", "GO")

  length_data <- expanded[, c("gene_id", "transcript_id", "length")]
  length_table <- data.table::setDT(length_data)
  length_table <- length_table[, .(mean_gene_length = mean(length)), by=.(gene_id)]
  names(length_table) <- c("ID", "length")

  retlist <- list(
    "go_data" = go_data,
    "go_table" = go_table,
    "length_data" = length_data,
    "length_table" = length_table)
  return(retlist)
}

## EOF
