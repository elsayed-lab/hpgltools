## annotation_uniprot.r: Some functions to simplify working with uniprot web
## services and text files.
## Recently I noticed that uniprot pretty drastically changed the ways in which one may
## access their data.  I think these changes will make it drastically easier to gather data
## from them, but I have only poked at a little thus far.

#' Download the txt uniprot data for a given accession/species.
#'
#' Uniprot is an astonishing resource, but man is it a pain to use.  Hopefully
#' this function will help.  It takes either a uniprot accession, taxonomy ID,
#' or species name and does its best to find the appropriate uniprot data.  This
#' is therefore primarily used by load_uniprot_annotations().
#'
#' @param accession Which accession to grab?
#' @param species Or perhaps species?
#' @param taxonomy Query for a specific taxonomy ID rather than species/accession?
#' @param all If there are more than 1 hit, grab them all?
#' @param first Or perhaps just grab the first hit?
#' @return A filename/accession tuple.
#' @seealso [xml2] [rvest]
#' @examples
#'  uniprot_sc_downloaded <- download_uniprot_proteome(species = "Saccharomyces cerevisiae S288c")
#'  uniprot_sc_downloaded$filename
#'  uniprot_sc_downloaded$species
#' @export
load_uniprot_annotations <- function(accession = NULL, species = "H37Rv",
                                     taxonomy = NULL, all = FALSE, first = FALSE) {
  final_species <- ""
  if (!is.null(taxonomy)) {
    request_url <- glue::glue("https://www.uniprot.org/proteomes/?query=taxonomy%3A{xml2::url_escape(taxonomy)}")
    destination <- glue("{taxonomy}.txt.gz")
    if (!file.exists(destination)) {
      ## tt <- download.file(url = request_url, destfile = destination, method = "wget", quiet = TRUE)
      tt <- download.file(url = request_url, destfile = destination, quiet = TRUE)
    }
    result <- xml2::read_html(destination)
    result_html <- rvest::html_nodes(result, "tr")
    accessions_text <- rvest::html_attr(result_html, "id")
    ## The first two elements are headers
    accessions_text <- accessions_text[3:length(accessions_text)]
    accessions <- gsub(x = accessions_text, pattern = "^(UP[0-9]+)(.*$)", replacement = "\\1")
    species_text <- rvest::html_nodes(result, "td") %>%
      rvest::html_nodes("span") %>%
      rvest::html_text()
    final_species <- species_text[species_text != ""]
    if (length(accessions) == 1) {
      accession <- accessions
    } else {
      accession <- accessions[1]
      name <- final_species[1]
    }
  } else {
    if (is.null(accession) & is.null(species)) {
      message("Defaulting to the Mycobacterium tuberculosis H37Rv strain.")
      accession <- "UP000001584"
    } else if (is.null(accession)) {
      message("Querying uniprot for the accession matching: ", species, ".")
      ##destination <- glue("{tempfile()}.txt.gz")
      ## request_url <- glue("https://www.uniprot.org/proteomes/?query={xml2::url_escape(species)}")

      ## Uniprot changed their web server in a few interesting ways which ultimately will make
      ## downloading annotations easier/better, but will also require me to rewrite a bunch
      ## of this code...
      ## The easiest way to deal with this in the short term is to use the website
      ## to 'download rest link' and get a new URL which will return TSV text which
      ## I can just dump to a tbl and play with...

      request_url <- glue("https://rest.uniprot.org/proteomes/stream?compressed=false&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query={xml2::url_escape(species)}")
      ##if (!file.exists(destination)) {
      ##  ## tt <- download.file(url = request_url, destfile = destination, method = "wget", quiet = TRUE)
      ##  tt <- download.file(url = request_url, destfile = destination)
      ##}
      ##result <- xml2::read_html(destination)
      ##result_html <- rvest::html_nodes(result, "tr")
      ##accessions_text <- rvest::html_attr(result_html, "id")
      #### The first two elements are headers
      ##accessions_text <- accessions_text[3:length(accessions_text)]
      ##accessions <- gsub(x = accessions_text, pattern = "^(UP[0-9]+)(.*$)", replacement = "\\1")
      ##species_text <- rvest::html_nodes(result, "td") %>%
      ##  rvest::html_nodes("span") %>%
      ##  rvest::html_text()
      ##final_species <- species_text[species_text != ""]
      ##removed <- file.remove(destination)
      ##if (length(accessions) == 1) {
      ##  accession <- accessions
      ##} else if (isTRUE(all)) {
      ##  for (a in 1:length(accessions)) {
      ##    name <- final_species[a]
      ##    accession <- accessions[a]
      ##    message("Downloading the proteome for ", name, ".")
      ##    tmp <- download_uniprot_proteome(accession = accession)
      ##    Sys.sleep(time = 3)
      ##  }
      ##} else if (isTRUE(first)) {
      ##  accession <- accessions[1]
      ##  name <- final_species[1]
      ##  message("Downloading the proteome for ", name, ".")
      ##  tmp <- download_uniprot_proteome(accession = accession)
      ##} else {
      ##  message("Here are the species found, please choose one and try again.")
      ##  for (a in 1:length(accessions)) {
      ##    name <- final_species[a]
      ##    accession <- accessions[a]
      ##    message(a, ") ", accession, ": ", name)
      ##  }
      ##  message(toString(final_species))
      ##  return(NULL)
      ##}
      accession_tbl <- as.data.frame(readr::read_tsv(request_url))
      message("This species provides: ", nrow(accession_tbl), " hits.")
      message("Arbitrarily downloading the first: ", accession_tbl[1, "Organism"], ".")
      accession <- as.character(accession_tbl[1, "Proteome Id"])
    }
  }
  ##request_url <- glue(
  ##    "https://www.uniprot.org/uniprot/?query=proteome:\\
  ##   {accession}&compress=yes&force=true&format=txt")


  request_url <- glue(
      "https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=accession%2Clineage%2Cvirus_hosts%2Clineage_ids%2Cgene_synonym%2Corganism_name%2Corganism_id%2Cprotein_name%2Cgene_orf%2Cgene_oln%2Cgene_names%2Cid%2Cgene_primary%2Cxref_proteomes%2Cabsorption%2Cft_act_site%2Cft_binding%2Cft_ca_bind%2Ccc_catalytic_activity%2Ccc_cofactor%2Cft_dna_bind%2Cec%2Ccc_activity_regulation%2Ccc_function%2Ckinetics%2Cft_metal%2Cft_np_bind%2Ccc_pathway%2Cph_dependence%2Credox_potential%2Crhea%2Cft_site%2Ctemp_dependence%2Cannotation_score%2Ccc_caution%2Ckeyword%2Ckeywordid%2Cprotein_existence%2Ccc_miscellaneous%2Creviewed%2Ctools%2Cuniparc_id%2Ccomment_count%2Cfeature_count%2Ccc_interaction%2Ccc_subunit%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id%2Ccc_allergen%2Ccc_biotechnology%2Ccc_disruption_phenotype%2Ccc_disease%2Cft_mutagen%2Ccc_pharmaceutical%2Ccc_toxic_dose%2Cft_intramem%2Ccc_subcellular_location%2Cft_topo_dom%2Cft_transmem%2Cft_chain%2Cft_crosslnk%2Cft_init_met%2Cft_lipid%2Ccc_ptm%2Cft_propep%2Cstructure_3d%2Cft_strand%2Cft_helix%2Cft_turn%2Clit_pubmed_id%2Cft_coiled%2Cft_compbias%2Cft_domain%2Cft_motif%2Cft_region%2Cft_repeat&format=tsv&query=proteome%3A{accession}")
  retdf <- readr::read_tsv(request_url)
  return(retdf)
}

#' Extract ontology information from a uniprot dataframe.
#'
#' @param ... Whatever args are required for load_uniprot_annotations()
#' @return Ontology dataframe
#' @seealso [load_uniprot_annotations()] [stringr] [tidyr]
#' @examples
#' \dontrun{
#'  uniprot_sc_downloaded <- download_uniprot_proteome(species = "Saccharomyces cerevisiae S288c")
#'  sc_uniprot_annot <- load_uniprot_annotations(file = uniprot_sc_downloaded$filename)
#'  sc_uniprot_go <- load_uniprot_go(sc_uniprot_annot)
#'  head(sc_uniprot_go)
#' }
#' @export
load_uniprot_go <- function(...) {
  input <- load_uniprot_annotations(...)

  kept <- input[, c("uniprot_accessions", "go", "aa_length")] %>%
    tidyr::separate_rows("uniprot_accessions")

  kept[["go"]] <- kept[["go"]] %>%
    stringr::str_extract_all(pattern = "GO:\\d+")
  kept[["go"]] <- I(kept[["go"]])
  kept[["go"]] <- as.character(kept[["go"]])
  kept <- kept %>%
    tidyr::separate_rows("go", sep = ",")
  kept[["go"]] <- gsub(pattern='"', replacement = "", x = kept[["go"]])
  kept[["go"]] <- gsub(pattern = ")", replacement = "", x = kept[["go"]])
  kept[["go"]] <- gsub(pattern = "c\\(", replacement = "", x = kept[["go"]])
  kept[["go"]] <- gsub(pattern = "\\s+", replacement = "", x = kept[["go"]])
  colnames(kept) <- c("ID", "GO", "length")
  return(kept)
}

## EOF
