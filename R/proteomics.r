read_thermo_xlsx <- function(xlsx_file, test_row=NULL) {
  old_options <- options(java.parameters="-Xmx20G")
  message(paste0("Reading ", xlsx_file))
  result <- readxl::read_xlsx(path=xlsx_file, sheet=1, col_names=FALSE)
  group_data <- list()
  bar <- utils::txtProgressBar(style=3)
  for (r in 1:nrow(result)) {
    row <- as.data.frame(result[r, ])
    row[, is.na(row)] <- ""
    pct_done <- r / nrow(result)
    setTxtProgressBar(bar, pct_done)
    ## The following 3 stanzas handle the creation of the levels of our data structure
    ## The first defines the protein group
    if (row[, 1] == "Checked") {
      group_colnames <- as.character(row)
      group_keepers <- !grepl(pattern="^$", x=group_colnames)
      group_keepers[1] <- FALSE
      group_colnames <- group_colnames[group_keepers]
      next
    }
    ## When the 2nd column is 'Checked', then this defines a new protein in the group.
    if (row[, 2] == "Checked") {
      protein_colnames <- as.character(row)
      protein_keepers <- !grepl(pattern="^$", x=protein_colnames)
      protein_keepers[2] <- FALSE
      protein_colnames <- protein_colnames[protein_keepers]
      next
    }
    ## When the 3rd column is 'Checked', then this starts a peptide definition
    if (row[, 3] == "Checked") {
      peptide_colnames <- as.character(row)
      peptide_keepers <- !grepl(pattern="^$", x=peptide_colnames)
      peptide_keepers[3] <- FALSE
      peptide_colnames <- peptide_colnames[peptide_keepers]
      next
    }
    ## Once the column names for the data are defined, we consider how to
    ## Fill in the actual data, the protein group is probably the least interesting.
    if (row[, 1] == FALSE | row[, 1] == TRUE) {
      group_information <- row[group_keepers]
      colnames(group_information) <- group_colnames
      group_information[["ID"]] <- sub(pattern="^.* GN=(\\w+) .*$",
                                       replacement="\\1",
                                       x=group_information[["Group Description"]])
      group_accession <- group_information[["Protein Group ID"]]
      group_list <- list(
        "summary" = group_information,
        "data" = list())
      group_data[[group_accession]] <- group_list
      next
    }
    ## When the 2nd column is FALSE, then this defined a protein in the group.
    ## The protein data structure is likely the most interesting.
    if (row[, 2] == FALSE | row[, 2] == TRUE) {
      protein_information <- row[protein_keepers]
      colnames(protein_information) <- protein_colnames
      protein_information[["ID"]] <- sub(pattern="^.* GN=(\\w+) .*$",
                                         replacement="\\1",
                                         x=protein_information[["Description"]])
      protein_accession <- protein_information[["Accession"]]
      protein_list <- list(
        "summary" = protein_information,
        "data" = data.frame())
      group_data[[group_accession]][["data"]][[protein_accession]] <- protein_list
      next
    }
    ## When the 3rd group is FALSE, then this adds a peptide.
    ## The peptide data structure is the most detailed, but probably not the most interesting.
    if (row[, 3] == FALSE | row[, 3] == TRUE) {
      peptide_information <- row[peptide_keepers]
      colnames(peptide_information) <- peptide_colnames
      current <- group_data[[group_accession]][["data"]][[protein_accession]][["data"]]
      new <- rbind(current, peptide_information)
      group_data[[group_accession]][["data"]][[protein_accession]][["data"]] <- new
      next
    }
  } ## End iterating over ever row of this unholy stupid data structure.
  close(bar)
  message("Finished parsing, reorganizing the protein data.")
  protein_df <- data.frame()
  peptide_df <- data.frame()
  protein_names <- c()
  message(paste0("Starting to iterate over ", length(group_data),  " groups."))
  bar <- utils::txtProgressBar(style=3)
  for (g in 1:length(group_data)) {
    pct_done <- g / length(group_data)
    setTxtProgressBar(bar, pct_done)
    group <- as.character(names(group_data)[g])
    protein_group <- group_data[[group]][["data"]]
    protein_accessions <- names(protein_group)
    for (p in 1:length(protein_accessions)) {
      protein <- protein_accessions[p]
      protein_names <- c(protein_names, protein)
      protein_summary <- group_data[[group]][["data"]][[protein]][["summary"]]
      protein_df <- rbind(protein_df, protein_summary)
      peptide_data <- group_data[[group]][["data"]][[protein]][["data"]]
      peptide_df <- rbind(peptide_df, peptide_data)
    }
  } ## End of the for loop
  close(bar)

  current_colnames <- colnames(protein_df)
  current_colnames <- tolower(current_colnames)
  current_colnames <- gsub(pattern="%", replacement="pct", x=current_colnames)
  current_colnames <- gsub(pattern=" ", replacement="_", x=current_colnames)
  current_colnames <- gsub(pattern="_confidence", replacement="", x=current_colnames)
  current_colnames <- gsub(pattern="\\(by_search_engine\\)", replacement="", x=current_colnames)
  current_colnames <- gsub(pattern="\\[|\\]|#|:|\\.|\\/|\\,|\\-", replacement="", x=current_colnames)
  current_colnames <- gsub(pattern="^_", replacement="", x=current_colnames)
  current_colnames <- gsub(pattern="__", replacement="_", x=current_colnames)
  current_colnames <- gsub(pattern="_ht", replacement="", x=current_colnames)
  current_colnames <- gsub(pattern="_mascot_mascot", replacement="_mascot", x=current_colnames)
  current_colnames <- gsub(pattern="_sequest_sequest", replacement="_sequest", x=current_colnames)
  colnames(protein_df) <- current_colnames

  ## Now make sure the columns which should be numeric, are numeric.
  numeric_cols <- c(
    "protein_fdr_mascot", "protein_fdr_sequest", "exp_qvalue_mascot", "expt_qvalue_sequest",
    "coverage_pct", "unique_peptides", "aas", "mw_kda", "calc_pi", "score_mascot",
    "score_sequest", "peptides_mascot", "peptides_sequest")
  for (col in numeric_cols) {
    if (!is.null(protein_df[[col]])) {
      protein_df[[col]] <- as.numeric(protein_df[[col]])
    }
  }

  ## Make sure columns which are 'found_in' are factors
  for (col in colnames(protein_df)) {
    if (grepl(pattern="^found_in_", x=col)) {
      protein_df[[col]] <- as.factor(protein_df[[col]])
    }
  }

  retlist <- list(
    "names" = protein_names,
    "group_data" = group_data,
    "protein_data" = protein_df,
    "peptide_data" = peptide_df)
  return(retlist)
}
