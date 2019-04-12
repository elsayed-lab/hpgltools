#' Read output from mayu to get the IP/PP number corresponding to a given FDR value.
#'
#' @param file Mayu output file.
#' @param fdr Chosen fdr value to acquire.
#' @return List of two elements: the full mayu table sorted by fdr and the number
#'   corresponding to the chosen fdr value.
#' @export
extract_mayu_pps_fdr <- function(file, fdr=0.01) {
  mayu_df <- readr::read_csv(file)
  fdr_df <- mayu_df[, c("IP/PPs", "protFDR")]
  fdr_idx <- order(fdr_df[["protFDR"]])
  fdr_df <- fdr_df[fdr_idx, ]
  mayu_df <- mayu_df[fdr_idx, ]
  keepers <- fdr_df[["protFDR"]] <= fdr
  fdr_df <- fdr_df[keepers, ]
  result <- tail(fdr_df, n=1)[[1]]
  retlist <- list(
    "fdr_number" = result,
    "table" = mayu_df)
  return(retlist)
}

#' Read a mzML/mzXML file and extract from it some important metadata.
#'
#' When working with swath data, it is fundamentally important to know the
#' correct values for a bunch of the input variables.  These are not trivial
#' to acquire.  This function attempts to make this easier (but slow) by reading
#' the mzXML file and parsing out helpful data.
#'
#' @param file Filename to read.
#' @param id An id to give the result.
#' @param write_acquisitions If a filename is provided, write a tab separated
#'   table of windows.
#' @param format Either mzXML or mzML.
#' @param allow_window_overlap One may choose to foce windows to not overlap.
#' @param start_add Add a minute to the start of the windows to avoid overlaps?
#' @return List containing a table of scan and precursor data.
#' @export
extract_scan_data <- function(file, id=NULL, write_acquisitions=TRUE, format="mzXML",
                               allow_window_overlap=FALSE, start_add=0) {
  if (format == "mzML") {
    extract_mzML_scans(file, id=id, write_acquisitions=write_acqusitions,
                       allow_window_overlap=allow_window_overlap, start_add=start_add)
  } else {
    extract_mzXML_scans(file, id=id, write_acquisitions=write_acqusitions,
                        allow_window_overlap=allow_window_overlap, start_add=start_add)
  }
}

#' Parse a mzXML file and return the relevant data.
#'
#' This does the actual work for extract_scan_data().  When I wrote this
#' function, I had forgotten about the mzR library; with that in mind, this
#' seems to give a bit more information and be a bit faster than my short tests
#' with mzR (note however that my tests were to compare mzR parsing mzML files
#' vs. this function with mzXML, which is a classic apples to oranges).
#'
#' This goes a step further to pull out the windows acquired in the MS/MS scan
#' and print them in formats acceptable to TPP/OpenMS (eg. with and without
#' headers).
#'
#' @param file Input mzXML file to parse.
#' @param id Chosen ID for the given file.
#' @param write_acquisitions Write acquisition windows.
#' @param allow_window_overlap Some downstream tools cannot deal with
#'   overlapping windows. Toggle that here.
#' @param start_add Other downstream tools appear to expect some padding at the
#'   beginning of each window.  Add that here.
#' @return The list of metadata, scan data, etc from the mzXML file.
#' @export
extract_mzXML_scans <- function(file, id=NULL, write_acquisitions=TRUE,
                                allow_window_overlap=FALSE, start_add=0) {
  if (is.null(id)) {
    id <- file
  }
  message("Reading ", file)
  input <- xml2::read_html(x=file, options="NOBLANKS")
  ## peaks <- rvest::xml_nodes(input, "peaks")

  message("Extracting instrument information for ", file)
  instruments <- rvest::xml_nodes(x=input, css="msinstrument")
  instrument_data <- data.frame(row.names=(1:length(instruments)), stringsAsFactors=FALSE)
  instrument_values <- c("msmanufacturer", "msmodel", "msionisation",
                         "msmassanalyzer", "msdetector")
  for (v in instrument_values) {
    datum <- rvest::xml_nodes(instruments, v)
    instrument_data[[v]] <- datum %>% rvest::html_attr("value")
  }
  datum <- rvest::xml_nodes(instruments, "software")
  instrument_data[["software_type"]] <- datum %>% rvest::html_attr("type")
  instrument_data[["software_name"]] <- datum %>% rvest::html_attr("name")
  instrument_data[["software_version"]] <- datum %>% rvest::html_attr("version")

  message("Extracting scan information for ", file)
  scans <- rvest::xml_nodes(input, "scan")
  scan_data <- data.frame(row.names=(1:length(scans)), stringsAsFactors=FALSE)
  scan_wanted <- c("peakscount", "scantype", "centroided", "mslevel", "polarity",
                   "retentiontime", "collisionenergy", "lowmz", "highmz", "basepeakmz",
                   "basepeakintensity", "totioncurrent")
  scan_numeric <- c("peakscount", "collisionenergy", "lowmz", "highmz", "basepeakmz",
                    "basepeakintensity", "totioncurrent")
  scan_factor <- c("polarity", "mslevel", "centroided", "scantype")
  ## We will also extract the html_text() of the precursors.
  for (w in scan_wanted) {
    scan_data[[w]] <- scans %>% rvest::html_attr(w)
  }
  for (n in scan_numeric) {
    scan_data[[n]] <- as.numeric(scan_data[[n]])
  }
  for (f in scan_factor) {
    scan_data[[n]] <- as.factor(scan_data[[n]])
  }

  message("Extracting precursor information for ", file)
  precursors <- rvest::xml_nodes(scans, "precursormz")
  precursor_data <- data.frame(row.names=(1:length(precursors)), stringsAsFactors=FALSE)
  precursor_wanted <- c("precursorintensity", "activationmethod",
                        "windowwideness", "precursorscannum")
  precursor_numeric <- c("precursorintensity", "precursorscannum",
                         "window_center", "windowwideness")
  precursor_factor <- c("activationmethod")
  for (w in precursor_wanted) {
    precursor_data[[w]] <- precursors %>% rvest::html_attr(w)
  }
  precursor_data[["window_center"]] <- precursors %>% rvest::html_text()
  for (n in precursor_numeric) {
    precursor_data[[n]] <- as.numeric(precursor_data[[n]])
  }
  for (f in precursor_factor) {
    precursor_data[[f]] <- as.factor(precursor_data[[f]])
  }
  precursor_data[["window_start"]] <- precursor_data[["window_center"]] -
    (precursor_data[["windowwideness"]] / 2)
  precursor_data[["window_end"]] <- precursor_data[["window_center"]] +
    (precursor_data[["windowwideness"]] / 2)

  message("Coalescing the acquisition windows for ", file)
  acquisition_windows <- precursor_data[, c("window_start", "window_end")]
  acquisition_unique <- !duplicated(x=acquisition_windows)
  acquisition_windows <- acquisition_windows[acquisition_unique, ]
  colnames(acquisition_windows) <- c("start", "end")
  if (!isTRUE(allow_window_overlap)) {
    previous_end <- acquisition_windows[1, "end"]
    ## If we are not allowing windows to overlap, then we will add 0.01 to the previous endpoint.
    for (it in 2:nrow(acquisition_windows)) {
      new_start <- acquisition_windows[it, "start"]
      if (new_start != previous_end) {
        acquisition_windows[it, "start"] <- previous_end
      }
      previous_end <- acquisition_windows[it, "end"]
    }
  } else if (!is.null(start_add)) {
    ## Conversely, we can add a static amount of time to the start of each window.
    acquisition_windows[["start"]] <- acquisition_windows[["start"]] + start_add
  }

  ## If requested, write out the acquisition windows in a format acceptable to openswath.
  if (write_acquisitions != FALSE) {
    acq_dir <- "windows"
    acq_file <- "acquisitions.txt"
    if (isTRUE(write_acquisitions)) {
      acq_file <- paste0(gsub(pattern="\\.mzXML", replacement="", x=basename(file)), ".txt")
    } else {
      acq_dir <- dirname(write_acquisitions)
      acq_file <- write_acquisitions
    }
    if (!file.exists(acq_dir)) {
      dir.create(acq_dir, recursive=TRUE)
    }
    ## To any sane person, the following lines must look bizarre.
    ## When performing a swath analysis, one step requires windows with _no_
    ## column headers (part of making the spectral libraries), while the
    ## invocation of OpenSwathWorkFlow or whatever it is, _requires_ them.
    ## So, yeah, that is annoying, but whatever.
    pre_file <- file.path(acq_dir, acq_file)
    message("Writing acquisition file to: ", pre_file)
    no_cols <- write.table(x=acquisition_windows, file=pre_file, sep="\t", quote=FALSE,
                           row.names=FALSE, col.names=FALSE)
    osw_file <- file.path(acq_dir, glue("openswath_{acq_file}"))
    ## This is the file for openswathworkflow.
    message("Writing osw acquisitions to: ", osw_file)
    plus_cols <- write.table(x=acquisition_windows, file=osw_file,
                            sep="\t", quote=FALSE,
                            row.names=FALSE, col.names=TRUE)
  }
  retlist <- list(
    "file" = file,
    "instrument" = instrument_data,
    "scans" = scan_data,
    "precursors" = precursor_data,
    "acquisitions" = acquisition_windows)
  return(retlist)
}

#' Parse a mzML file and return the relevant data.
#'
#' This does the actual work for extract_scan_data().  This levers mzR to
#' provide the data and goes a step further to pull out the windows acquired in
#' the MS/MS scan and print them in formats acceptable to TPP/OpenMS (eg. with
#' and without headers).
#'
#' @param file Input mzML file to parse.
#' @param id Chosen ID for the given file.
#' @param write_acquisitions Write acquisition windows.
#' @param allow_window_overlap Some downstream tools cannot deal with
#'   overlapping windows. Toggle that here.
#' @param start_add Other downstream tools appear to expect some padding at the
#'   beginning of each window.  Add that here.
#' @return The list of metadata, scan data, etc from the mzXML file.
#' @export
extract_mzML_scans <- function(file, id=NULL, write_acquisitions=TRUE,
                              allow_window_overlap=FALSE, start_add=0) {
  if (is.null(id)) {
    id <- file
  }
  message("Reading ", file)
  input <- mzR::openMSfile(file)
  instrument <- mzR::instrumentInfo(input)
  info <- mzR::runInfo(input)
  scan_data <- mzR::header(input)
  closed <- mzR::close(input)
  retlist <- list(
    "instrument" = instrument,
    "info" = info,
    "scan" = scan_data)
  return(retlist)
}

#' Read a bunch of mzXML files to acquire their metadata.
#'
#' I have had difficulties getting the full set of correct parameters for a
#' DDA/DIA experiment.  After some poking, I eventually found most of these
#' required prameters in the mzXML raw files.  Ergo, this function uses them.
#' 20190310: I had forgotten about the mzR library.  I think much (all?) of this
#' is redundant with respect to it and perhaps should be removed in deference to
#' the more complete and fast implementation included in mzR.
#'
#' @param metadata Data frame describing the samples, including the mzXML
#'   filenames.
#' @param write_windows Write out SWATH window frames.
#' @param id_column What column in the sample sheet provides the ID for the samples?
#' @param allow_window_overlap What it says on the tin, some tools do not like
#'   DIA windows to overlap, if TRUE, this will make sure each annotated window
#'   starts at the end of the previous window if they overlap.
#' @param start_add Another strategy is to just add a static amount to each
#'   window.
#' @param format Currently this handles mzXML or mzML files.
#' @param parallel Perform operations using an R foreach cluster?
#' @param savefile If not null, save the resulting data structure to an rda file.
#' @param ... Extra arguments, presumably color palettes and column names and
#'   stuff like that.
#' @return List of data extracted from every sample in the MS run (DIA or DDA).
#' @export
extract_msraw_data <- function(metadata, write_windows=TRUE, id_column="sampleid",
                               allow_window_overlap=FALSE, start_add=0, format="mzXML",
                               parallel=TRUE, savefile=NULL, ...) {
  arglist <- list(...)

  ## Add a little of the code from create_expt to include some design
  ## information in the returned data structure.

  ## Since I am using this code now in two places, if I can use it without
  ## changes, I will split it into its own function so make changes in one place
  ## happen in both but until I can ensure to myself that it is functional in
  ## both contexts I will not. Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  sample_definitions <- data.frame()
  if (class(metadata) == "data.frame") {
    sample_definitions <- metadata
  } else {
      sample_definitions <- extract_metadata(metadata, ...)
      ##sample_definitions <- extract_metadata(metadata)
  }

  file_column <- "file"
  if (!is.null(arglist[["file_column"]])) {
      file_column <- arglist[["file_column"]]  ## Make it possible to have multiple count
      ## tables / sample in one sheet.
  }

  sample_column <- "sampleid"
  if (!is.null(arglist[["sample_column"]])) {
      sample_column <- arglist[["sample_column"]]
      sample_column <- tolower(sample_column)
      sample_column <- gsub(pattern="[[:punct:]]", replacement="", x=sample_column)
  }

  chosen_colors <- generate_expt_colors(sample_definitions, ...)
  ## chosen_colors <- generate_expt_colors(sample_definitions)
  meta <- sample_definitions[, c("sampleid", file_column)]
  colnames(meta) <- c("id", "file")
  existing_files <- complete.cases(meta[["file"]])
  if (sum(existing_files) != nrow(meta)) {
    warning("It appears that some files are missing in the metadata.")
  }
  meta <- meta[existing_files, ]

  ## Set up the 'cluster' and process the mzXML files.
  returns <- list()
  res <- list()
  num_files <- nrow(meta)
  if (isTRUE(parallel)) {
    tt <- sm(requireNamespace("parallel"))
    tt <- sm(requireNamespace("doParallel"))
    tt <- sm(requireNamespace("iterators"))
    tt <- sm(requireNamespace("foreach"))
    tt <- sm(try(attachNamespace("foreach"), silent=TRUE))
    ## cores <- parallel::detectCores() / 2
    cores <- 4
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
    if (isTRUE(show_progress)) {
      bar <- utils::txtProgressBar(max=num_files, style=3)
    }
    progress <- function(n) {
      setTxtProgressBar(bar, n)
    }
    pb_opts <- list()
    if (isTRUE(show_progress)) {
      pb_opts <- list("progress" = progress)
    }
    res <- foreach(i=1:num_files, .packages=c("hpgltools", "doParallel"),
                   .options.snow=pb_opts, .export=c("extract_scan_data")) %dopar% {
      file <- meta[i, "file"]
      id <- meta[i, "id"]
      file_result <- try(extract_scan_data(file, id=id, write_acquisitions=write_windows,
                                           allow_window_overlap=allow_window_overlap,
                                           format=format, start_add=start_add))
      if (class(file_result)[1] != "try-error") {
        returns[[file]] <- file_result
      }
    }
    if (isTRUE(show_progress)) {
      close(bar)
    }
    parallel::stopCluster(cl)
  } else {
    for (i in 1:num_files) {
      file <- meta[i, "file"]
      id <- meta[i, "id"]
      file_result <- try(extract_scan_data(file, id=id, write_acquisitions=write_windows,
                                           allow_window_overlap=allow_window_overlap,
                                           format=format, start_add=start_add))
      if (class(file_result)[1] != "try-error") {
        res[[file]] <- file_result
      }
    }
  }
  rownames(sample_definitions) <- make.names(sample_definitions[[id_column]], unique=TRUE)
  names(res) <- rownames(sample_definitions)

  retlist <- list(
    "colors" = chosen_colors,
    "metadata" = sample_definitions,
    "sample_data" = res)
  if (!is.null(savefile)) {
    mzxml_data <- retlist
    save_result <- try(save(list = c("mzxml_data"), file=savefile), silent=TRUE)
  }
  return(retlist)
}

#' Get some data from a peptideprophet run.
#
#' I am not sure what if any parameters this should have, but it seeks to
#' extract the useful data from a peptide prophet run.  In the situation in
#' which I wish to use it, the input command was:
#' > xinteract -dDECOY_ -OARPpd -Nfdr_library.xml comet_result.pep.xml
#' Eg. It is a peptideprophet result provided by TPP.
#' I want to read the resulting xml table and turn it into a data.table so that
#' I can plot some metrics from it.
#'
#' @param pepxml  The file resulting from the xinteract invocation.
#' @param decoy_string  What prefix do decoys have in the data.
#' @param ... Catch extra arguments passed here, currently unused.
#' @return data table of all the information I saw fit to extract
#' The columns are:
#' * protein: The name of the matching sequence (DECOYs allowed here)
#' * decoy: TRUE/FALSE, is this one of our decoys?
#' * peptide: The sequence of the matching spectrum.
#' * start_scan: The scan in which this peptide was observed
#' * end scan: Ibid
#' * index This seems to just increment
#' * precursor_neutral_mass: Calculated mass of this fragment assuming no
#'   isotope shenanigans (yeah, looking at you C13).
#' * assumed_charge: The expected charge state of this peptide.
#' * retention_time_sec: The time at which this peptide eluted during the run.
#' * peptide_prev_aa:  The amino acid before the match.
#' * peptide_next_aa:  and the following amino acid.
#' * num_tot_proteins: The number of matches not counting decoys.
#' * num_matched_ions: How many ions for this peptide matched?
#' * tot_num_ions:  How many theoretical ions are in this fragment?
#' * matched_ion_ratio: num_matched_ions / tot_num_ions, bigger is better!
#' * cal_neutral_pep_mass: This is redundant with precursor_neutral_mass, but
#'   recalculated by peptideProphet, so if there is a discrepency we should yell
#'   at someone!
#' * massdiff How far off is the observed mass vs. the calculated? (also
#'   redundant with massd later)
#' * num_tol_term: The number of peptide termini which are consistent with the
#'   cleavage (hopefully 2), but potentially 1 or even 0 if digestion was
#'   bad. (redundant with ntt later)
#' * num_missed_cleavages: How many cleavages must have failed in order for this
#'   to be a good match?
#' * num_matched_peptides: Number of alternate possible peptide matches.
#' * xcorr: cross correlation of the experimental and theoretical spectra (this
#'   is supposedly only used by sequest, but I seem to have it here...)
#' * deltacn: The normalized difference between the xcorr values for the best hit and next
#'   best hit.  Thus higher numbers suggest better matches.
#' * deltacnstar: Apparently 'important for things like phospho-searches
#'   containing homologous top-scoring peptides when analyzed by
#'   peptideprophet...' -- the comet release notes.
#' * spscore: The raw value of preliminary score from the sequest algorithm.
#' * sprank: The rank of the match in a preliminary score. 1 is good.
#' * expect: E-value of the given peptide hit.  Thus how many identifications
#'   one expect to observe by chance, lower is therefore better
#' * prophet_probability: The peptide prophet probability score, higher is
#'   better.
#' * fval: 0.6(the dot function + 0.4(the delta dot function) - (the dot bias
#'   penalty function) -- which is to say... well I dunno, but it is supposed to
#'   provide information about how similar this match is to other potential
#'   matches, so I presume higher means the match is more ambiguous.
#' * ntt: Redundant with num_tol_term above, but this time from peptide prophet.
#' * nmc: Redundant with num_missed_cleavages, except it coalesces them.
#' * massd: Redundant with massdiff
#' * isomassd: The mass difference, but taking into account stupid C13.
#' * RT: Retention time
#' * RT_score: The score of the retention time!
#' * modified_peptides: A string describing modifications in the found peptide
#' * variable_mods: A comma separated list of the variable modifications
#'   observed.
#' * static_mods: A comma separated list of the static modifications observed.
#' @export
extract_peprophet_data <- function(pepxml, decoy_string="DECOY_", ...) {
  input <- xml2::read_html(pepxml, options="NOBLANKS")

  message("Extracting spectrum queries.")
  spectrum_queries <- rvest::xml_nodes(input, "spectrum_query")
  spectra <- spectrum_queries %>% rvest::html_attr("spectrum")
  query_data <- data.frame(row.names=spectra)
  ## The interesting material at the beginning of a spectrum, these are in the
  ## <spectrum query> tag.

  message("Extracting the spectrum_query metadata.")
  toplevel_interesting <- c("start_scan", "end_scan", "precursor_neutral_mass",
                            "assumed_charge", "index", "retention_time_sec")
  for (t in toplevel_interesting) {
    query_data[[t]] <- spectrum_queries %>%
      rvest::html_attr(t)
  }

  message("Extracting the search_result metadata.")
  search_results <- rvest::xml_nodes(spectrum_queries, "search_result")
  search_hits <- search_results %>%
    rvest::html_node(xpath="search_hit")
  ## The set of fields which look interesting to me in the search_hit data.
  search_fields <- c("peptide", "peptide_prev_aa", "peptide_next_aa",
                     "protein", "num_tot_proteins",
                     "num_matched_ions", "tot_num_ions",
                     "calc_neutral_pep_mass", "massdiff", "num_tol_term",
                     "num_missed_cleavages", "num_matched_peptides")
  for (s in search_fields) {
    query_data[[s]] <- search_hits %>%
      rvest::html_attr(s)
  }
  query_data[["decoy"]] <- FALSE
  decoy_regex <- glue("^{decoy_string}")
  decoy_idx <- grepl(pattern=decoy_regex, x=query_data[["protein"]])
  query_data[decoy_idx, "decoy"] <- TRUE
  query_data[["matched_ion_ratio"]] <- as.numeric(query_data[["num_matched_ions"]]) /
    as.numeric(query_data[["tot_num_ions"]])

  ## Get modification info
  message("Extracting modification metadata.")
  query_data[["modified_peptides"]] <- search_hits %>%
    rvest::html_node(xpath="modification_info") %>%
    rvest::html_attr("modified_peptide")

  na_idx <- is.na(query_data[["modified_peptides"]])
  query_data[na_idx, "modified_peptides"] <- ""
  query_data[["variable_mods"]] <- ""
  query_data[["static_mods"]] <- ""
  modification_test <- search_hits %>%
    rvest::html_node(xpath="modification_info")

  message("Filling in modification information, this is slow.")
  show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style=3)
  }
  for (i in 1:length(modification_test)) {
    if (isTRUE(show_progress)) {
      pct_done <- i / length(modification_test)
      setTxtProgressBar(bar, pct_done)
    }
    test <- modification_test[[i]]
    if (!is.na(test)) {
      variables <- test %>%
        rvest::html_nodes(xpath="mod_aminoacid_mass") %>%
        rvest::html_attr("variable")
      statics <- test %>%
        rvest::html_nodes(xpath="mod_aminoacid_mass") %>%
        rvest::html_attr("static")
      positions <- test %>%
        rvest::html_nodes(xpath="mod_aminoacid_mass") %>%
        rvest::html_attr("position")
      masses <- test %>%
        rvest::html_nodes(xpath="mod_aminoacid_mass") %>%
        rvest::html_attr("mass")
      variable_idx <- !is.na(variables)
      if (sum(variable_idx) > 0) {
        variable_string <- toString(
          glue("position: {positions[variable_idx]} mass: {masses[variable_idx]}\\
                mod: {variables[variable_idx]}"))
        query_data[i, "variable_mods"] <- variable_string
      }
      static_idx <- !is.na(statics)
      if (sum(static_idx) > 0) {
        static_string <- toString(
          glue("position: {positions[static_idx]} mass: {masses[static_idx]}\\
                mod: {statics[static_idx]}"))
        query_data[i, "static_mods"] <- static_string
      }
    }
  }
  if (isTRUE(show_progress)) {
    close(bar)
  }

  ## Extracting the search_score tags
  message("Extracting the search_score metadata.")
  score_results <- rvest::xml_nodes(search_hits, "search_score")
  score_names <- score_results %>%
    rvest::html_attr("name")
  score_values <- score_results %>%
    rvest::html_attr("value")
  names(score_values) <- score_names
  for (v in unique(score_names)) {
    query_data[[v]] <- score_values[names(score_values) == v]
  }

  ## Get the peptideprophet_result
  message("Extracting the peptideprophet_result probabilities.")
  peptide_prophets <- rvest::xml_nodes(spectrum_queries, "peptideprophet_result")
  query_data[["prophet_probability"]] <- peptide_prophets %>%
    rvest::html_attr("probability")

  ## Get the peptideprophet parameters
  message("Extracting the search parameters.")
  parameter_results <- rvest::xml_nodes(spectrum_queries, "parameter")
  parameter_names <- parameter_results %>%
    rvest::html_attr("name")
  parameter_values <- parameter_results %>%
    rvest::html_attr("value")
  names(parameter_values) <- parameter_names
  for (p in unique(parameter_names)) {
    query_data[[p]] <- parameter_values[names(parameter_values) == p]
  }

  numeric_columns <- c("start_scan", "end_scan", "precursor_neutral_mass", "index",
                       "retention_time_sec", "calc_neutral_pep_mass", "massdiff",
                       "num_matched_peptides", "xcorr", "deltacn",
                       "deltacnstar", "spscore", "expect",
                       "prophet_probability", "fval", "massd", "RT", "RT_score")
  factor_columns <- c("assumed_charge", "num_tot_proteins", "num_matched_ions", "num_tol_term",
                      "num_missed_cleavages", "sprank", "ntt", "nmc", "isomassd")
  for (n in numeric_columns) {
    query_data[[n]] <- as.numeric(query_data[[n]])
  }
  for (f in factor_columns) {
    query_data[[f]] <- as.factor(query_data[[f]])
  }

  new_order <- c("protein", "decoy", "peptide", "start_scan", "end_scan",
                 "index", "precursor_neutral_mass", "assumed_charge",
                 "retention_time_sec", "peptide_prev_aa", "peptide_next_aa",
                 "num_tot_proteins", "num_matched_ions", "tot_num_ions",
                 "matched_ion_ratio", "calc_neutral_pep_mass", "massdiff",
                 "num_tol_term", "num_missed_cleavages", "num_matched_peptides",
                 "xcorr", "deltacn", "deltacnstar", "spscore", "sprank",
                 "expect", "prophet_probability", "fval", "ntt", "nmc", "massd",
                 "isomassd", "RT", "RT_score", "modified_peptides",
                 "variable_mods", "static_mods")
  query_data <- query_data[, new_order]
  check_masses <- testthat::expect_equal(query_data[["precursor_neutral_mass"]],
                                         query_data[["calc_neutral_pep_mass"]],
                                         tolerance=0.1)
  result <- data.table::as.data.table(query_data)
  return(result)
}

#' Read a bunch of scored swath outputs from pyprophet to acquire their metrics.
#'
#' This function is mostly cribbed from the other extract_ functions in this file.
#' With it, I hope to be able to provide some metrics of a set of openswath runs, thus
#' potentially opening the door to being able to objectively compare the same
#' run with different options and/or different runs.
#'
#' Likely columns generated by exporting OpenMS data via pyprophet include:
#' transition_group_id:  Incrementing ID of the transition in the MS(.pqp)
#' library used for matching (I am pretty sure).
#' decoy:  Is this match of a decoy peptide?
#' run_id: This is a bizarre encoding of the run, OpenMS/pyprophet re-encodes
#' the run ID from the filename to a large signed integer.
#' filename:  Which raw mzXML file provides this particular intensity value?
#' rt: Retention time in seconds for the matching peak group.
#' assay_rt: The expected retention time after normalization with the iRT. (how
#'   does the iRT change this value?)
#' delta_rt:  The difference between rt and assay_rt
#' irt: (As described in the abstract of Claudia Escher's 2012 paper: "Here we
#'   present iRT, an empirically derived dimensionless peptide-specific value that
#'   allows for highly accurate RT prediction. The iRT of a peptide is a fixed
#'   number relative to a standard set of reference iRT-peptides that can be
#'   transferred across laboratories and chromatographic systems.")
#' assay_irt: The iRT observed in the actual chromatographic run.
#' delta_irt: The difference.  I am seeing that all the delta iRTs are in the
#'   -4000 range for our actual experiment; since this is in seconds, does that
#'   mean that it is ok as long as they stay in a similar range?
#' id: unique long signed integer for the peak group.
#' sequence: The sequence of the matched peptide
#' fullunimodpeptidename: The sequence, but with unimod formatted modifications
#'   included.
#' charge:  The assumed charge of the observed peptide.
#' mz:  The m/z value of the precursor ion.
#' intensity:  The sum of all transition intensities in the peak group.
#' aggr_prec_peak_area:  Semi-colon separated list of intensities (peak areas)
#'   of the MS traces for this match.
#' aggr_prec_peak_apex:  Intensity peak apexes of the MS1 traces.
#' leftwidth: The start of the peak group in seconds.
#' rightwidth: The end of the peak group in seconds.
#' peak_group_rank: When multiple peak groups match, which one is this?
#' d_score: I think this is the score as retured by openMS (higher is better).
#' m_score: I am pretty sure this is the result of a SELECT QVALUE operation in pyprophet.
#' aggr_peak_area: The intensities of this fragment ion separated by semicolons.
#' aggr_peak_apex: The intensities of this fragment ion separated by semicolons.
#' aggr_fragment_annotation:  Annotations of the fragment ion traces by semicolon.
#' proteinname:  Name of the matching protein.
#' m_score_protein_run_specific: I am guessing the fdr for the pvalue for this run.
#' mass: Mass of the observed fragment.
#'
#' @param metadata  Data frame describing the samples, including the mzXML
#'   filenames.
#' @param pyprophet_column Which column from the metadata provides the requisite filenames?
#' @param savefile  If not null, save the data from this to the given filename.
#' @param ... Extra arguments, presumably color palettes and column names and
#'   stuff like that.
#' @return  A list of data from each sample in the pyprophet scored DIA run.
#' @export
extract_pyprophet_data <- function(metadata, pyprophet_column="diascored",
                                   savefile=NULL, ...) {
  arglist <- list(...)

  ## Add a little of the code from create_expt to include some design
  ## information in the returned data structure.

  ## Since I am using this code now in two places, if I can use it without changes, I will
  ## split it into its own function so make changes in one place happen in both
  ## but until I can ensure to myself that it is functional in both contexts I will not.
  ## Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  sample_definitions <- data.frame()
  if (class(metadata) == "data.frame") {
    sample_definitions <- metadata
  } else {
    sample_definitions <- extract_metadata(metadata, ...)
    ## sample_definitions <- extract_metadata(metadata)
  }

  if (is.null(sample_definitions[[pyprophet_column]])) {
    stop("This required a column with the tsv scored pyprophet data.")
  }

  chosen_colors <- generate_expt_colors(sample_definitions, ...)
  ## chosen_colors <- generate_expt_colors(sample_definitions)
  meta <- sample_definitions[, c("sampleid", pyprophet_column)]
  colnames(meta) <- c("id", "scored")
  existing_files <- complete.cases(meta[["scored"]])
  if (sum(existing_files) != nrow(meta)) {
    warning("It appears that some files are missing in the metadata.")
  }
  meta <- meta[existing_files, ]

  gather_masses <- function(sequence) {
    atoms <- try(BRAIN::getAtomsFromSeq(sequence), silent=TRUE)
    if (class(atoms)[1] != "try-error") {
      d <- BRAIN::useBRAIN(atoms)
      ret <- round(d[["avgMass"]])
    } else {
      ret <- 0
    }
    return(ret)
  }

  res <- list()
  num_files <- nrow(meta)
  failed_files <- c()
  for (i in 1:num_files) {
    file <- meta[i, "scored"]
    id <- meta[i, "id"]
    message("Attempting to read the tsv file for: ", id, ": ", file, ".")
    file_result <- sm(try(readr::read_tsv(file), silent=TRUE))
    if (class(file_result)[1] != "try-error") {
      colnames(file_result) <- tolower(colnames(file_result))
      file_result <- file_result %>%
        dplyr::rowwise() %>%
        dplyr::mutate(mass=gather_masses(sequence))
      res[[id]] <- file_result
    } else {
      message("Failed to read: ", id, ".")
      failed_files <- c(failed_files, file)
    }
  }
  found_ids <- names(res)
  ## Now cull the sample definitions to only include those we actually found.
  sample_idx <- sample_definitions[["sampleid"]] %in% found_ids
  sample_definitions <- sample_definitions[sample_idx, ]

  retlist <- list(
    "failed" = failed_files,
    "colors" = chosen_colors,
    "metadata" = sample_definitions,
    "sample_data" = res)
  if (!is.null(savefile)) {
    pyprophet_data <- retlist
    save_result <- try(save(list = c("pyprophet_data"), file=savefile), silent=TRUE)
  }
  return(retlist)
}

##extract_traml_data <- function(traml) {
##  message("Reading the TraML file.")
##  ## xml2 is painfully slow and annoying with these files, but I do not know why.
##  input <- xml2::read_html(x=traml, options="NOBLANKS")
##  ##children <- input %>%
##  ##  xml2::xml_children()
##  ##contents <- input %>%
##  ##  xml2::xml_contents()
##  all_nodes <- input %>%
##    rvest::html_nodes("*")
##  protein_nodes <- all_nodes %>%
##    rvest::html_nodes(xpath="//Protein")
##}

#' An attempt to address a troubling question when working with DIA data.
#'
#' My biggest concern when treating DIA data in a RNASeqish manner is the fact
#' that if a given peptide is not identified, that is not the same thing as
#' stating that it was not translated.  It is somewhat reminiscent of the often
#' mocked and repeated Donald Rumsfeld statement regarding known unknowns
#' vs. unknown unknowns.  Thus, in an RNASeq experiment, if one sees a zero, one
#' may assume that transcript was not transcribed, it may be assumed to be a
#' known zero(unknown).  In contrast, if the same thing happens in a DIA data
#' set, that represents an unknown unknown.  Perhaps it was not translated, and
#' perhaps it was not identified.
#'
#' This function therefore does the following:
#'   1.  Backfill all 0s in the matrix to NA.
#'   2.  Performs a mean across all samples which are known technical replicates
#'   of the same biological replicate.  This mean is performed using
#'   na.rm=TRUE.  Thus the entries which used to be 0 should no longer affect
#'   the result.
#'   3.  Recreate the expressionset with the modified set of samples.
#'
#' @param expt Starting expressionset to mangle.
#' @param fact Metadata factor to use when taking the mean of biological
#'   replicates.
#' @param fun Assumed to be mean, but one might want median.
#' @return new expressionset
#' @export
mean_by_bioreplicate <- function(expt, fact="bioreplicate", fun="mean") {
  ## Set all the zeros to NA so that when we do cpm and mean they will get dropped.
  exprs_set <- expt[["expressionset"]]
  mtrx <- exprs(expt)
  zero_idx <- mtrx == 0
  new <- mtrx
  new[zero_idx] <- NA
  new_libsize <- colSums(new, na.rm=TRUE)
  new <- edgeR::cpm(new, lib.size=new_libsize)
  exprs(exprs_set) <- new
  expt[["expressionset"]] <- exprs_set
  annot <- fData(expt)
  final <- median_by_factor(expt, fact=fact, fun=fun)
  current_design <- pData(expt)
  new_design <- data.frame()
  for (c in 1:length(colnames(final))) {
    colname <- colnames(final)[c]
    possible_rows <- which(current_design[[fact]] == colname)
    chosen_row_idx <- possible_rows[1]
    chosen_row <- current_design[chosen_row_idx, ]
    new_design <- rbind(new_design, chosen_row)
  }
  rownames(new_design) <- colnames(final)
  new_design[["sampleid"]] <- rownames(new_design)
  new_design[["batch"]] <- "undefined"
  new_set <- create_expt(count_dataframe=final, metadata=new_design, gene_info=annot)
  return(new_set)
}

#' Plot mzXML peak intensities with respect to m/z.
#'
#' I want to have a pretty plot of peak intensities and m/z.  The plot provided
#' by this function is interesting, but suffers from some oddities; notably that
#' it does not currently separate the MS1 and MS2 data.  Since I am stuck on
#' this forsaken plane with no hope of ever leaving, perhaps I can add that now.
#'
#' @param mzxml_data  The data structure from extract_mzxml or whatever it is.
#' @param loess  Do a loess smoothing from which to extract a function
#'   describing the data?  This is terribly slow, and in the data I have
#'   examined so far, not very helpful, so it is FALSE by default.
#' @param alpha  Make the plotted dots opaque to this degree.
#' @param ms1  Include MS1 data in the plot?
#' @param ms2  Include MS2 data in the plot?
#' @param x_scale  Plot the x-axis on a non linear scale?
#' @param y_scale  Plot the y-axis on a non linear scale?
#' @param ...  Extra arguments for the downstream functions.
#' @return  ggplot2 goodness.
#' @export
plot_intensity_mz <- function(mzxml_data, loess=FALSE, alpha=0.5, ms1=TRUE, ms2=TRUE,
                              x_scale=NULL, y_scale=NULL, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    message("Adding ", name)
    plotted_table <- sample_data[[i]][["scans"]]
    ## Caveat!  I do not have my data while on this fucking plane, so I might
    ## have forgotten the name of that column in the data.  If so, the following
    ## will fail.
    if (!isTRUE(ms1)) {
      kept_idx <- plotted_table[["level"]] != "MS1"
      plotted_table <- plotted_table[kept_idx, ]
    }
    if (!isTRUE(ms1)) {
      kept_idx <- plotted_table[["level"]] != "MS2"
      plotted_table <- plotted_table[kept_idx, ]
    }
    plotted_data <- plotted_table[, c("basepeakmz", "basepeakintensity")]
    plotted_data[["sample"]] <- name
    plotted_data <- plotted_data[, c("sample", "basepeakmz", "basepeakintensity")]
    colnames(plotted_data) <- c("sample", "mz", "intensity")
    ## Re-order the columns because I like sample first.
    plot_df <- rbind(plot_df, plotted_data)
  }

  ## Drop rows from the metadata and colors which had errors.
  metadata <- metadata[keepers, ]
  colors <- colors[keepers]

  chosen_palette <- "Dark2"
  sample_colors <- sm(
    grDevices::colorRampPalette(
                 RColorBrewer::brewer.pal(samples, chosen_palette))(samples))

  ## Randomize the rows of the df so we can see if any sample is actually overrepresented
  plot_df <- plot_df[sample(nrow(plot_df)), ]

  if (!is.null(x_scale)) {
    plot_df[["mz"]] <- check_plot_scale(plot_df[["mz"]], scale)[["data"]]
  }
  if (!is.null(y_scale)) {
    plot_df[["intensity"]] <- check_plot_scale(plot_df[["intensity"]], scale)[["data"]]
  }

  int_vs_mz <- ggplot(data=plot_df, aes_string(x="mz", y="intensity",
                                               fill="sample", colour="sample")) +
    ggplot2::geom_point(alpha=alpha, size=0.5) +
    ggplot2::scale_fill_manual(
               name="Sample", values=sample_colors,
               guide=ggplot2::guide_legend(override.aes=aes(size=3))) +
    ggplot2::scale_color_manual(
               name="Sample", values=sample_colors,
               guide=ggplot2::guide_legend(override.aes=aes(size=3))) +
    ggplot2::theme_bw(base_size=base_size)

  if (!is.null(x_scale)) {
    int_vs_mz <- int_vs_mz + ggplot2::scale_x_continuous(trans=scales::log2_trans())
  }
  if (!is.null(y_scale)) {
    int_vs_mz <- int_vs_mz + ggplot2::scale_y_continuous(trans=scales::log2_trans())
  }

  if (isTRUE(lowess)) {
    int_vs_mz <- int_vs_mz +
      ggplot2::geom_smooth(method="loess", size=1.0)
  }
  retlist <- list(
    "data" = plotted_data,
    "plot" = int_vs_mz)
  return(retlist)
}

#' Make a boxplot out of some of the various data available in the mzxml data.
#'
#' There are a few data within the mzXML raw data files which are likely
#'   candidates for simple summary via a boxplot/densityplot/whatever.  For the
#'   moment I am just doing boxplots of a few of them.  Since my metadata
#'   extractor dumps a couple of tables, one must choose a desired table and
#'   column from it to plot.
#'
#' @param mzxml_data  Provide a list of mzxml data, one element for each sample.
#' @param table  One of precursors or scans
#' @param column  One of the columns from the table; if 'scans' is chosen, then
#'   likely choices include: 'peakscount', 'basepeakmz', 'basepeakintensity'; if
#'   'precursors' is chosen, then the only likely choice for the moment is
#'   'precursorintensity'.
#' @param violin  Print the samples as violins rather than only box/whiskers?
#' @param names  Names for the x-axis of the plot.
#' @param title  Title the plot?
#' @param scale  Put the data on a specific scale?
#' @param ...  Further arguments, presumably for colors or some such.
#' @return  Boxplot describing the requested column of data in the set of mzXML files.
#' @export
plot_mzxml_boxplot <- function(mzxml_data, table="precursors", column="precursorintensity",
                               violin=FALSE, names=NULL, title=NULL, scale=NULL, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    message("Adding ", name)
    names(colors)[i] <- name
    plotted_table <- sample_data[[i]][[table]]
    plotted_data <- as.data.frame(plotted_table[[column]])
    plotted_data[["sample"]] <- name
    plotted_data[["color"]] <- colors[[name]]
    colnames(plotted_data) <- c(column, "sample", "color")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "color")]
    plot_df <- rbind(plot_df, plotted_data)
  }

  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  scale_data <- check_plot_scale(plot_df[[column]], scale)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  plot_df[[column]] <- scale_data[["data"]]
  plot_df[["color"]] <- as.factor(plot_df[["color"]])

  boxplot <- ggplot2::ggplot(data=plot_df, ggplot2::aes_string(x="sample", y=column))
  if (isTRUE(violin)) {
    boxplot <- boxplot +
      ggplot2::geom_violin(aes_string(fill="sample"),
                           width=1, scale="area", show.legend=FALSE) +
      ggplot2::geom_boxplot(na.rm=TRUE, alpha=0.3, color="black", size=0.5,
                            outlier.alpha=0.01, width=0.2)
  } else {
    boxplot <- boxplot +
      sm(ggplot2::geom_boxplot(aes_string(fill="sample"),
                               na.rm=TRUE, fill=colors, size=0.5,
                               outlier.size=1.5,
                               outlier.colour=ggplot2::alpha("black", 0.2)))
  }
  boxplot <- boxplot +
    ggplot2::scale_fill_manual(values=as.character(colors)) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
    ggplot2::xlab("Sample") + ggplot2::ylab(column)
  if (!is.null(title)) {
    boxplot <- boxplot + ggplot2::ggtitle(title)
  }
  if (!is.null(names)) {
    boxplot <- boxplot + ggplot2::scale_x_discrete(labels=names)
  }
  scale <- "log"
  if (scale == "log") {
    boxplot <- boxplot + ggplot2::scale_y_continuous(trans=scales::log2_trans())
  } else if (scale == "logdim") {
    boxplot <- boxplot + ggplot2::coord_trans(y="log2")
  } else if (isTRUE(scale)) {
    boxplot <- boxplot + ggplot2::scale_y_log10()
  }

  return(boxplot)
}

#' Make a boxplot out of some of the various data available in the pyprophet
#' data.
#'
#' This function is mostly redundant with the plot_mzxml_boxplot above.
#' Unfortunately, the two data types are subtly different enough that I felt it
#' not worth while to generalize the functions.
#'
#' @param pyprophet_data List containing the pyprophet results.
#' @param column What column of the pyprophet scored data to plot?
#' @param keep_real Do we keep the real data when plotting the data? (perhaps
#'   we only want the decoys)
#' @param keep_decoys Do we keep the decoys when plotting the data?
#' @param expt_names Names for the x-axis of the plot.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the desired column from the data.
#' @export
plot_pyprophet_distribution <- function(pyprophet_data, column="delta_rt", keep_real=TRUE,
                                        keep_decoys=TRUE, expt_names=NULL, label_chars=10,
                                        title=NULL, scale=NULL, ...) {
  arglist <- list(...)
  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)

  ## Reset the sample names if one wants a specific column from the metadata.
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      names(sample_data) <- make.names(metadata[[expt_names]], unique=TRUE)
    } else {
      names(sample_data) <- expt_names
    }
  }

  keepers <- c()
  for (i in 1:samples) {
    name <- names(sample_data)[i]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    message("Adding ", name)
    plotted_table <- sample_data[[i]]
    if (!isTRUE(keep_decoys)) {
      good_idx <- plotted_table[["decoy"]] != 1
      plotted_table <- plotted_table[good_idx, ]
    }
    if (!isTRUE(keep_real)) {
      good_idx <- plotted_table[["decoy"]] != 0
      plotted_table <- plotted_table[good_idx, ]
    }
    plotted_data <- as.data.frame(plotted_table[c("sequence", "proteinname",
                                                  "aggr_fragment_annotation", column)])
    plotted_data[["sample"]] <- name
    colnames(plotted_data) <- c("sequence", "proteinname", "fragment", column, "sample")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "sequence", "proteinname", "fragment")]
    plot_df <- rbind(plot_df, plotted_data)
  }
  ## I am not certain this is valid.
  plot_df[[column]] <- abs(plot_df[[column]])

##  testing <- data.table::as.data.table(plot_df)
##  recast_dt <- data.table::dcast.data.table(data=testing,
##                                            formula=sequence+proteinname~sample,
##                                            fun.aggregate=mean,
##                                            value.var="intensity")
##  names <- recast_dt[["proteinname"]]
##  sequences <- recast_dt[["sequence"]]
##  recast_dt[, c("proteinname", "sequence") := NULL]
##  nan_idx <- is.na(recast_dt)
##  recast_dt[nan_idx] <- 0
##  recast_norm <- log2(
##    1 + preprocessCore::normalize.quantiles.robust(as.matrix(recast_dt)))
##  remelt <- as.data.table(recast_norm)
##  remelt[["proteinname"]] <- names
##  remelt[["sequence"]] <- sequences
##  remelted <- data.table::melt(data=remelt, value.name="intensity")
##  colnames(remelted) <- c("proteinname", "sequence", "sample", "intensity")

  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  scale_data <- check_plot_scale(plot_df[[column]], scale)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  plot_df[[column]] <- scale_data[["data"]]

  if (!is.null(label_chars) & is.numeric(label_chars)) {
    plot_df[["sample"]] <- abbreviate(plot_df[["sample"]], minlength=label_chars)
  }
  boxplot <- ggplot2::ggplot(data=plot_df, ggplot2::aes_string(x="sample", y=column)) +
    sm(ggplot2::geom_boxplot(na.rm=TRUE,
                             ggplot2::aes_string(fill="sample"),
                             fill=colors,
                             size=0.5,
                             outlier.size=1.5,
                             outlier.colour=ggplot2::alpha("black", 0.2))) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
    ggplot2::xlab("Sample") + ggplot2::ylab(column)
  if (!is.null(title)) {
    boxplot <- boxplot + ggplot2::ggtitle(title)
  }
  scale <- "log"
  if (scale == "log") {
    boxplot <- boxplot + ggplot2::scale_y_continuous(trans=scales::log2_trans())
  } else if (scale == "logdim") {
    boxplot <- boxplot + ggplot2::coord_trans(y="log2")
  } else if (isTRUE(scale)) {
    boxplot <- boxplot + ggplot2::scale_y_log10()
  }

  density <- ggplot(data=plot_df, ggplot2::aes_string(x=column, colour="sample")) +
    ggplot2::geom_density(aes_string(x=column, y="..count..", fill="sample"),
                          position="identity", na.rm=TRUE) +
    ggplot2::scale_colour_manual(values=as.character(colors)) +
    ggplot2::scale_fill_manual(values=ggplot2::alpha(as.character(colors), 0.1)) +
    ggplot2::ylab("Number of genes.") +
    ggplot2::xlab("Number of hits/gene.") +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   legend.key.size=ggplot2::unit(0.3, "cm"))
  density <- directlabels::direct.label(density)

  violin <- ggplot(data=plot_df, aes_string(x="sample", y=column)) +
    ggplot2::geom_violin(aes_string(fill="sample"), width=1, scale="area") +
    ggplot2::geom_boxplot(aes_string(fill="sample"), outlier.alpha=0.01, width=0.2) +
    ggplot2::scale_fill_manual(values=as.character(colors)) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, hjust=1),
                   legend.position="none")

  dotboxplot <- boxplot +
    ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.1),
                         size=2, alpha=0.2)

  retlist <- list(
    "violin" = violin,
    "boxplot" = boxplot,
    "dotboxplot" = dotboxplot,
    "density" = density)
  return(retlist)
}

#' Read data from pyprophet and plot columns from it.
#'
#' More proteomics diagnostics!  Now that I am looking more closely, I think
#' this should be folded into plot_pyprophet_distribution().
#'
#' @param pyprophet_data Data from extract_pyprophet_data()
#' @param column Chosen column to plot.
#' @param keep_real FIXME: This should be changed to something like 'data_type'
#'   here and in plot_pyprophet_distribution.
#' @param keep_decoys Do we keep the decoys when plotting the data?
#' @param expt_names Names for the x-axis of the plot.
#' @param label_chars Maximum number of characters before abbreviating sample
#'   names.
#' @param protein chosen protein(s) to plot.
#' @param title Title the plot?
#' @param scale Put the data on a specific scale?
#' @param ... Further arguments, presumably for colors or some such.
#' @return Boxplot describing the desired column from the data.
#' @export
plot_pyprophet_protein <- function(pyprophet_data, column="intensity", keep_real=TRUE,
                                   keep_decoys=TRUE, expt_names=NULL, label_chars=10,
                                   protein=NULL, title=NULL, scale=NULL, ...) {
  arglist <- list(...)
  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)

  ## Reset the sample names if one wants a specific column from the metadata.
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      names(sample_data) <- make.names(metadata[[expt_names]], unique=TRUE)
    } else {
      names(sample_data) <- expt_names
    }
  }

  keepers <- c()
  for (i in 1:samples) {
    name <- names(sample_data)[i]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    message("Adding ", name)
    plotted_table <- sample_data[[i]]
    if (!isTRUE(keep_decoys)) {
      good_idx <- plotted_table[["decoy"]] != 1
      plotted_table <- plotted_table[good_idx, ]
    }
    if (!isTRUE(keep_real)) {
      good_idx <- plotted_table[["decoy"]] != 0
      plotted_table <- plotted_table[good_idx, ]
    }
    plotted_data <- as.data.frame(plotted_table[c("sequence", "proteinname",
                                                  "aggr_fragment_annotation", column)])
    plotted_data[["sample"]] <- name
    colnames(plotted_data) <- c("sequence", "proteinname", "fragment", column, "sample")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column, "sequence", "proteinname", "fragment")]
    plot_df <- rbind(plot_df, plotted_data)
  }

  ## Drop rows from the metadata and colors which had errors.
  if (length(keepers) > 0) {
    metadata <- metadata[keepers, ]
    colors <- colors[keepers]
  } else {
    stop("Something bad happened to the set of kept samples.")
  }

  if (is.null(protein)) {
    stop("This requires a protein ID to search.")
  } else {
    kept_prot_idx <- grepl(pattern=protein, x=plot_df[["proteinname"]])
    plot_df <- plot_df[kept_prot_idx, ]
  }
  plot_df[["sequence"]] <- as.factor(plot_df[["sequence"]])

  scale_data <- check_plot_scale(plot_df[[column]], scale, ...)
  if (is.null(scale)) {
    scale <- scale_data[["scale"]]
  }
  plot_df[[column]] <- scale_data[["data"]]

  if (!is.null(label_chars) & is.numeric(label_chars)) {
    plot_df[["sample"]] <- abbreviate(plot_df[["sample"]], minlength=label_chars)
  }
  violin <- ggplot2::ggplot(data=plot_df, ggplot2::aes_string(x="sample", y=column)) +
    ggplot2::geom_violin(aes_string(fill="sample"), width=1, scale="area") +
    ggplot2::scale_fill_manual(values=as.character(colors)) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                   axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab(column) +
    ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.1),
                         size=2, alpha=0.5)
  if (!is.null(title)) {
    violin <- violin + ggplot2::ggtitle(title)
  }
  scale <- "log"
  if (scale == "log") {
    violin <- violin + ggplot2::scale_y_continuous(
                                  labels=scales::scientific,
                                  trans=scales::log2_trans())
  } else if (scale == "logdim") {
    violin <- violin + ggplot2::coord_trans(y="log2")
  } else if (isTRUE(scale)) {
    violin <- violin + ggplot2::scale_y_log10()
  }

  return(violin)
}

#' Plot some data from the result of extract_peprophet_data()
#'
#' extract_pyprophet_data() provides a ridiculously large data table of a scored
#' openswath data after processing by pyprophet.
#'
#' @param pyprophet_data  List of pyprophet data, one element for each sample,
#'   taken from  extract_peprophet_data()
#' @param xaxis  Column to plot on the x-axis
#' @param xscale Change the scale of the x-axis?
#' @param yaxis  guess!
#' @param yscale  Change the scale of the y-axis?
#' @param alpha  How see-through to make the dots?
#' @param legend  Include a legend of samples?
#' @param size_column  Use a column for scaling the sizes of dots in the plot?
#' @param ... extra options which may be used for plotting.
#' @return a plot!
#' @export
plot_pyprophet_data <- function(pyprophet_data, xaxis="mass", xscale=NULL,
                                yaxis="leftwidth", yscale=NULL, alpha=0.4,
                                legend=TRUE, size_column="mscore", ...) {
  arglist <- list(...)

  metadata <- pyprophet_data[["metadata"]]
  colors <- pyprophet_data[["colors"]]
  sample_data <- pyprophet_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  keepers <- c()
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    if (class(sample_data[[i]])[1] == "try-error") {
      next
    }
    keepers <- c(keepers, i)
    message("Adding ", name)
    plotted_table <- sample_data[[i]]
    plotted_table[["sample"]] <- name
    plot_df <- rbind(plot_df, plotted_table)
  }

  ## Drop rows from the metadata and colors which had errors.
  metadata <- metadata[keepers, ]
  colors <- colors[keepers]

  chosen_palette <- "Dark2"
  sample_colors <- sm(
    grDevices::colorRampPalette(
                 RColorBrewer::brewer.pal(samples, chosen_palette))(samples))

  ## Randomize the rows of the df so we can see if any sample is actually overrepresented
  plot_df <- plot_df[sample(nrow(plot_df)), ]

  if (is.null(plot_df[[xaxis]])) {
    stop("The x axis data seems to be missing.")
  }
  if (is.null(plot_df[[yaxis]])) {
    stop("The y axis data seems to be missing.")
  }

  if (!is.null(xscale)) {
    plot_df[[xaxis]] <- check_plot_scale(plot_df[[xaxis]], scale)[["data"]]
  }
  if (!is.null(yscale)) {
    plot_df[[yaxis]] <- check_plot_scale(plot_df[[yaxis]], scale)[["data"]]
  }

  x_vs_y <- ggplot(data=plot_df, aes_string(x=xaxis, y=yaxis,
                                            fill="sample", colour="sample")) +
    ggplot2::geom_point(alpha=alpha, size=0.5) +
    ggplot2::scale_fill_manual(
               name="Sample", values=sample_colors,
               guide=ggplot2::guide_legend(override.aes=aes(size=3))) +
    ggplot2::scale_color_manual(
               name="Sample", values=sample_colors,
               guide=ggplot2::guide_legend(override.aes=aes(size=3))) +
    ggplot2::theme_bw(base_size=base_size)

  if (!is.null(xscale)) {
    x_vs_y <- x_vs_y + ggplot2::scale_x_continuous(trans=scales::log2_trans())
  }
  if (!is.null(yscale)) {
    x_vs_y <- x_vs_y + ggplot2::scale_y_continuous(trans=scales::log2_trans())
  }
  if (isTRUE(lowess)) {
    x_vs_y <- x_vs_y +
      ggplot2::geom_smooth(method="loess", size=1.0)
  }
  if (!isTRUE(legend)) {
    x_vs_y <- x_vs_y +
      ggplot2::theme(legend.position="none")
  }
  retlist <- list(
    "data" = plot_df,
    "plot" = x_vs_y)
  return(retlist)
}

#' Plot some data from the result of extract_peprophet_data()
#'
#' extract_peprophet_data() provides a ridiculously large data table of a comet
#' result after processing by RefreshParser and xinteract/peptideProphet.
#' This table has some 37-ish columns and I am not entirely certain which ones
#' are useful as diagnostics of the data.  I chose a few and made options to
#' pull some/most of the rest.  Lets play!
#'
#' @param table  Big honking data table from extract_peprophet_data()
#' @param xaxis  Column to plot on the x-axis
#' @param xscale Change the scale of the x-axis?
#' @param yaxis  guess!
#' @param yscale  Change the scale of the y-axis?
#' @param size_column  Use a column for scaling the sizes of dots in the plot?
#' @param ... extra options which may be used for plotting.
#' @return a plot!
#' @export
plot_peprophet_data <- function(table, xaxis="precursor_neutral_mass", xscale=NULL,
                                yaxis="num_matched_ions", yscale=NULL,
                                size_column="prophet_probability", ...) {
  arglist <- list(...)
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["chosen_palette"]])) {
    chosen_palette <- arglist[["chosen_palette"]]
  }
  color_column <- "decoy"
  if (!is.null(arglist[["color_column"]])) {
    color_column <- arglist[["color_column"]]
  }
  if (is.null(table[[color_column]])) {
    table[["color"]] <- "black"
  } else {
    table[["color"]] <- as.factor(table[[color_column]])
  }
  color_list <- NULL
  num_colors <- nlevels(as.factor(table[["color"]]))
  if (num_colors == 2) {
    color_list <- c("darkred", "darkblue")
  } else {
    color_list <- sm(
      grDevices::colorRampPalette(
                   RColorBrewer::brewer.pal(num_colors, chosen_palette))(num_colors))
  }

  if (is.null(table[[xaxis]])) {
    stop(glue("The x-axis column: {xaxis} does not appear in the data."))
  }
  if (is.null(table[[yaxis]])) {
    stop(glue("The y-axis column: {yaxis} does not appear in the data."))
  }

  table <- as.data.frame(table)
  if (is.null(table[[size_column]])) {
    table[["size"]] <- 1
  } else {
    if (class(table[[size_column]]) == "numeric") {
      ## quants <- as.numeric(quantile(unique(table[[size_column]])))
      ## size_values <- c(4, 8, 12, 16, 20)
      ## names(size_values) <- quants
      table[["size"]] <- table[[size_column]]
    } else {
      table[["size"]] <- 1
    }
  }
  ##min_val <- min(table[[size_column]])
  ##max_val <- max(table[[size_column]])
  range <- as.numeric(quantile(unique(table[["size"]])))
  table[table[["size"]] >= range[[5]], "size"] <- "06biggest"
  table[table[["size"]] >= range[[4]] &
        table[["size"]] < range[[5]], "size"] <- "05big"
  table[table[["size"]] >= range[[3]] &
        table[["size"]] < range[[4]], "size"] <- "04medium_big"
  table[table[["size"]] >= range[[2]] &
        table[["size"]] < range[[3]], "size"] <- "03medium_small"
  table[table[["size"]] >= range[[1]] &
        table[["size"]] < range[[2]], "size"] <- "02small"
  table[table[["size"]] < range[[1]], "size"] <- "01smallest"

  ## Setting the factor/vector of sizes is a bit confusing to me.
  table[["size"]] <- as.factor(table[["size"]])
  levels(table[["size"]]) <- c("01smallest", "02small", "03medium_small",
                              "04medium_big", "05big", "06biggest")
  my_sizes <- c("01smallest"=0.4, "02small"=8, "03medium_small"=1.2,
                "04medium_big"=1.6, "05big"=2.0, "06biggest"=2.4)

  scale_x_cont <- "raw"
  if (!is.null(xscale)) {
    if (is.numeric(xscale)) {
      table[[xaxis]] <- log(table[[xaxis]] + 1) / log(xscale)
    } else if (xscale == "log2") {
      scale_x_cont <- "log2"
    } else if (xscale == "log10") {
      scale_x_cont <- "log10"
    } else {
      message("I do not understand your scale.")
    }
  }
  scale_y_cont <- "raw"
  if (!is.null(yscale)) {
    if (is.numeric(yscale)) {
      table[[xaxis]] <- log(table[[yaxis]] + 1) / log(yscale)
    } else if (yscale == "log2") {
      scale_y_cont <- "log2"
    } else if (yscale == "log10") {
      scale_y_cont <- "log10"
    } else {
      message("I do not understand your scale.")
    }
  }

  table[["text"]] <- glue("{table[['protein']]}:{table[['peptide']]}")

  a_plot <- ggplot(data=table, aes_string(x=xaxis, y=yaxis, text="text",
                                          color="color", size="size")) +
    ggplot2::geom_point(alpha=0.4, aes_string(fill="color", color="color")) +
    ggplot2::scale_color_manual(name="color", values=color_list) +
    ggplot2::geom_rug() +
    ggplot2::scale_size_manual(values=c(0.2, 0.6, 1.0, 1.4, 1.8, 2.2))
  if (scale_x_cont == "log2") {
    a_plot <- ggplot2::scale_x_continuous(trans=scales::log2_trans())
  } else if (scale_x_cont == "log10") {
    a_plot <- ggplot2::scale_x_continuous(trans=scales::log10_trans())
  }
  if (scale_y_cont == "log2") {
    a_plot <- ggplot2::scale_y_continuous(trans=scales::log2_trans())
  } else if (scale_y_cont == "log10") {
    a_plot <- ggplot2::scale_y_continuous(trans=scales::log10_trans())
  }

  return(a_plot)
}

#' Parse the difficult thermo fisher xlsx file.
#'
#' The Thermo(TM) workflow has as its default a fascinatingly horrible excel
#' output.  This function parses that into a series of data frames.
#'
#' @param xlsx_file  The input xlsx file
#' @param test_row  A single row in the xlsx file to use for testing, as I have
#'   not yet seen two of these accursed files which had the same headers.
#' @return  List containing the protein names, group data, protein dataframe,
#'   and peptide dataframe.
#' @export
read_thermo_xlsx <- function(xlsx_file, test_row=NULL) {
  old_options <- options(java.parameters="-Xmx20G")
  message("Reading ", xlsx_file)
  result <- readxl::read_xlsx(path=xlsx_file, sheet=1, col_names=FALSE)
  group_data <- list()
  show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style=3)
  }
  for (r in 1:nrow(result)) {
    if (isTRUE(show_progress)) {
      pct_done <- r / nrow(result)
      setTxtProgressBar(bar, pct_done)
    }
    row <- as.data.frame(result[r, ])
    row[, is.na(row)] <- ""
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
  if (isTRUE(show_progress)) {
    close(bar)
  }
  message("Finished parsing, reorganizing the protein data.")

  protein_df <- data.frame()
  peptide_df <- data.frame()
  protein_names <- c()
  message("Starting to iterate over ", length(group_data),  " groups.")
  show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style=3)
  }
  for (g in 1:length(group_data)) {
    if (isTRUE(show_progress)) {
      pct_done <- g / length(group_data)
      setTxtProgressBar(bar, pct_done)
    }
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
  if (isTRUE(show_progress)) {
    close(bar)
  }
  current_colnames <- colnames(protein_df)
  current_colnames <- tolower(current_colnames)
  ## percent signs are stupid in columns.
  current_colnames <- gsub(
    pattern="%", replacement="pct", x=current_colnames)
  ## as are spaces.
  current_colnames <- gsub(
    pattern=" ", replacement="_", x=current_colnames)
  ## A bunch of columns have redundant adjectives.
  current_colnames <- gsub(
    pattern="_confidence", replacement="", x=current_colnames)
  ## Extra text in a column name is useless
  current_colnames <- gsub(
    pattern="\\(by_search_engine\\)", replacement="", x=current_colnames)
  ## Get rid of a bunch of doofusy punctuation.
  current_colnames <- gsub(
    pattern="\\[|\\]|#|:|\\.|\\/|\\,|\\-", replacement="", x=current_colnames)
  ## At this point we should not have any leading underscores.
  current_colnames <- gsub(
    pattern="^_", replacement="", x=current_colnames)
  ## Now should we have any double underscores.
  current_colnames <- gsub(
    pattern="__", replacement="_", x=current_colnames)
  ## Finally, because of the previous removals, there might be some duplicated
  ## terms left behind.
  current_colnames <- gsub(
    pattern="_ht", replacement="", x=current_colnames)
  current_colnames <- gsub(
    pattern="_mascot_mascot", replacement="_mascot", x=current_colnames)
  current_colnames <- gsub(
    pattern="_sequest_sequest", replacement="_sequest", x=current_colnames)
  colnames(protein_df) <- current_colnames

  ## Now make sure the columns which should be numeric, are numeric.
  numeric_cols <- c(
    "protein_fdr_mascot", "protein_fdr_sequest", "exp_qvalue_mascot",
    "expt_qvalue_sequest", "coverage_pct", "unique_peptides", "aas", "mw_kda",
    "calc_pi", "score_mascot", "score_sequest", "peptides_mascot",
    "peptides_sequest")
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

subset_pyprophet_data <- function(lst, subset=NULL, column="protein_id", operator="in") {
  data <- lst[["sample_data"]]
  ## First, I want to get the simple Rv ID from the data, this is annoyingly Tb
  ## specific and should be removed or in some way made generic for other data.
  for (c in 1:length(data)) {
    datum <- data[[c]]
    datum[["protein_id"]] <- gsub(x=datum[["proteinname"]], pattern="^1/",
                                  replacement="")
    datum[["protein_id"]] <- gsub(x=datum[["protein_id"]], pattern="_.*$",
                                  replacement="")
    data[[c]] <- datum
  }

  ## Now, separately subset the data to find the proteins of interest.
  for (c in 1:length(data)) {
    datum <- data[[c]]
    datum_idx <- rep(TRUE, nrow(datum))
    if (operator == "in") {
      datum_idx <- datum[[column]] %in% subset
    } else if (operator == "equals") {
      datum_idx <- datum[[column]] == subset
    } else if (operator == "gt") {
      datum_idx <- datum[[column]] > subset
    } else if (operator == "lt") {
      datum_idx <- datum[[column]] < subset
    } else {
      message("I do not know this operator, doing nothing.")
      return(lst)
    }
    datum <- datum[datum_idx, ]
    data[[c]] <- datum
  }
  lst[["sample_data"]] <- data

  return(lst)
}

## EOF
