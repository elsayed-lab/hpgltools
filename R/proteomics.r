#' Replace 0 with NA if not all entries for a given condition are 0.
#'
#' This will hopefully handle a troubling corner case in Volker's data:
#' He primarily wants to find proteins which are found in one condition, but
#' _not_ in another.  However, due to the unknown unknown problem in DIA
#' acquisition, answering this question is difficult.  If one uses a normal
#' expressionset or msnset or whatever, one of two things will happen:
#' either the 0/NA proteins will be entirely removed/ignored, or they will
#' lead to spurious 'significant' calls.  MSstats, to its credit, does a lot to
#' try to handle these cases; but in the case Volker is most interested, it will
#' exclude the interesting proteins entirely.
#'
#' So, here is what I am going to do: Iterate through each element of the chosen
#' experimental design factor, check if all samples for that condition are 0, if
#' so; leave them.  If not all the samples have 0 for the given condition, then
#' replace the zero entries with NA.  This should allow for stuff like
#' rowMeans(na.rm = TRUE) to provide useful information.
#'
#' Finally, this will add columns to the annotations which tell the number of
#' observations for each protein after doing this.
#'
#' @param expt Expressionset to examine.
#' @param fact Experimental design factor to use.
#' @param method Specify whether to leave the NAs as NA,
#'  or replace them with the mean of all non-NA values.
#' @return New expressionset with some, but not all, 0s replaced with NA.
#' @export
add_conditional_nas <- function(expt, fact = "condition", method = "NA") {
  exprs_set <- expt[["expressionset"]]
  mtrx <- exprs(expt)
  annotations <- fData(expt)
  if (length(fact) == 1) {
    design <- pData(expt)
    fact <- design[[fact]]
    names(fact) <- rownames(design)
  }
  types <- levels(fact)
  used_columns <- c()
  observations_df <- data.frame(row.names = rownames(mtrx))
  for (t in 1:length(types)) {
    type <- types[t]
    used_columns <- grep(pattern = type, x = fact)
    if (length(used_columns) < 1) {
      warning("The level ", type, " of the factor has no columns in the data.")
      next
    }
    sub_mtrx <- mtrx[, used_columns]
    ## Make a note of how many samples are in this factor
    total_observations <- ncol(sub_mtrx)

    ## These are the columns we will leave 0
    zero_sum_idx <- rowSums(sub_mtrx) == 0
    all_zeros <- rownames(sub_mtrx)[zero_sum_idx]
    zero_idx <- sub_mtrx == 0
    ## Now we set the zeros to NA and then reset the all-zeros to 0.
    if (method == "NA") {
      sub_mtrx[zero_idx] <- NA
    } else if (method == "mean") {
      all_mean <- mean(mtrx, na.rm = TRUE)
      sub_mtrx[zero_idx] <- all_mean
    } else if (method == "sub_mean") {
      sub_mean <- mean(sub_mtrx, na.rm = TRUE)
      sub_mtrx[zero_idx] <- sub_mean
    }
    sub_mtrx[all_zeros, ] <- 0
    message("In condition ", type, " there are ", length(all_zeros),
            " rows which are all zero.")
    ## Record the number of observations for each protein for each condition.
    observation_column <- total_observations - rowSums(zero_idx)
    observations_df <- cbind(observations_df, observation_column)
    colnames(observations_df)[t] <- glue::glue("{type}_observations")

    for (c in 1:length(used_columns)) {
      replace_col <- used_columns[c]
      mtrx[, replace_col] <- sub_mtrx[, c]
    }
  }
  ## Now put the pieces back together.
  ## I am choosing not to use cbind to ensure that the orders are not screwed up.
  ## annotations <- cbind(annotations, observations_df)
  annotations <- merge(annotations, observations_df, by = "row.names")
  rownames(annotations) <- annotations[["Row.names"]]
  annotations[["Row.names"]] <- NULL
  exprs(exprs_set) <- mtrx
  fData(exprs_set) <- annotations
  expt[["expressionset"]] <- exprs_set
  return(expt)
}

#' Read output from mayu to get the IP/PP number corresponding to a given FDR value.
#'
#' @param file Mayu output file.
#' @param fdr Chosen fdr value to acquire.
#' @return List of two elements: the full mayu table sorted by fdr and the number
#'  corresponding to the chosen fdr value.
#' @export
extract_mayu_pps_fdr <- function(file, fdr = 0.01) {
  mayu_df <- readr::read_csv(file)
  fdr_df <- mayu_df[, c("IP/PPs", "protFDR")]
  fdr_idx <- order(fdr_df[["protFDR"]])
  fdr_df <- fdr_df[fdr_idx, ]
  mayu_df <- mayu_df[fdr_idx, ]
  keepers <- fdr_df[["protFDR"]] <= fdr
  fdr_df <- fdr_df[keepers, ]
  result <- tail(fdr_df, n = 1)[[1]]
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
#'  table of windows.
#' @param format Either mzXML or mzML.
#' @param allow_window_overlap One may choose to foce windows to not overlap.
#' @param start_add Add a minute to the start of the windows to avoid overlaps?
#' @return List containing a table of scan and precursor data.
#' @export
extract_scan_data <- function(file, id = NULL, write_acquisitions = TRUE, format = "mzXML",
                              allow_window_overlap = FALSE, start_add = 0) {
  if (format == "mzML") {
    extract_mzML_scans(file, id = id, write_acquisitions = write_acquisitions,
                       allow_window_overlap = allow_window_overlap, start_add = start_add)
  } else {
    extract_mzXML_scans(file, id = id, write_acquisitions = write_acquisitions,
                        allow_window_overlap = allow_window_overlap, start_add = start_add)
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
#'  overlapping windows. Toggle that here.
#' @param start_add Other downstream tools appear to expect some padding at the
#'  beginning of each window.  Add that here.
#' @return The list of metadata, scan data, etc from the mzXML file.
#' @export
extract_mzXML_scans <- function(file, id = NULL, write_acquisitions = TRUE,
                                allow_window_overlap = FALSE, start_add = 0) {
  if (is.null(id)) {
    id <- file
  }
  message("Reading ", file)
  input <- xml2::read_html(x = file, options = "NOBLANKS")
  ## peaks <- rvest::xml_nodes(input, "peaks")

  message("Extracting instrument information for ", file)
  instruments <- rvest::xml_nodes(x = input, css = "msinstrument")
  instrument_data <- data.frame(row.names=(1:length(instruments)), stringsAsFactors = FALSE)
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
  scan_data <- data.frame(row.names=(1:length(scans)), stringsAsFactors = FALSE)
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
  precursor_data <- data.frame(row.names=(1:length(precursors)), stringsAsFactors = FALSE)
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
  acquisition_unique <- !duplicated(x = acquisition_windows)
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
      acq_file <- paste0(gsub(pattern = "\\.mzXML", replacement = "", x = basename(file)), ".txt")
    } else {
      acq_dir <- dirname(write_acquisitions)
      acq_file <- write_acquisitions
    }
    if (!file.exists(acq_dir)) {
      dir.create(acq_dir, recursive = TRUE)
    }
    ## To any sane person, the following lines must look bizarre.
    ## When performing a swath analysis, one step requires windows with _no_
    ## column headers (part of making the spectral libraries), while the
    ## invocation of OpenSwathWorkFlow or whatever it is, _requires_ them.
    ## So, yeah, that is annoying, but whatever.
    pre_file <- file.path(acq_dir, acq_file)
    message("Writing acquisition file to: ", pre_file)
    no_cols <- write.table(x = acquisition_windows, file = pre_file, sep = "\t", quote = FALSE,
                           row.names = FALSE, col.names = FALSE)
    osw_file <- file.path(acq_dir, glue("openswath_{acq_file}"))
    ## This is the file for openswathworkflow.
    message("Writing osw acquisitions to: ", osw_file)
    plus_cols <- write.table(x = acquisition_windows, file = osw_file,
                             sep = "\t", quote = FALSE,
                             row.names = FALSE, col.names = TRUE)
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
#'  overlapping windows. Toggle that here.
#' @param start_add Other downstream tools appear to expect some padding at the
#'  beginning of each window.  Add that here.
#' @return The list of metadata, scan data, etc from the mzXML file.
#' @export
extract_mzML_scans <- function(file, id = NULL, write_acquisitions = TRUE,
                               allow_window_overlap = FALSE, start_add = 0) {
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
#'  filenames.
#' @param write_windows Write out SWATH window frames.
#' @param id_column What column in the sample sheet provides the ID for the samples?
#' @param file_column Which column in the sample sheet provides the filenames?
#' @param allow_window_overlap What it says on the tin, some tools do not like
#'  DIA windows to overlap, if TRUE, this will make sure each annotated window
#'  starts at the end of the previous window if they overlap.
#' @param start_add Another strategy is to just add a static amount to each
#'  window.
#' @param format Currently this handles mzXML or mzML files.
#' @param parallel Perform operations using an R foreach cluster?
#' @param savefile If not null, save the resulting data structure to an rda file.
#' @param ... Extra arguments, presumably color palettes and column names and
#'  stuff like that.
#' @return List of data extracted from every sample in the MS run (DIA or DDA).
#' @export
extract_msraw_data <- function(metadata, write_windows = TRUE, id_column = "sampleid",
                               file_column = "raw_file", allow_window_overlap = FALSE,
                               start_add = 0, format = "mzXML",
                               parallel = TRUE, savefile = NULL, ...) {
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
  if (class(metadata)[1] == "data.frame") {
    sample_definitions <- metadata
  } else {
    sample_definitions <- extract_metadata(metadata, ...)
    ## sample_definitions <- extract_metadata(metadata)
  }

  sample_column <- "sampleid"
  if (!is.null(arglist[["sample_column"]])) {
    sample_column <- arglist[["sample_column"]]
    sample_column <- tolower(sample_column)
    sample_column <- gsub(pattern = "[[:punct:]]", replacement = "", x = sample_column)
  }

  chosen_colors <- generate_expt_colors(sample_definitions, ...)
  ## chosen_colors <- generate_expt_colors(sample_definitions)
  meta <- sample_definitions[, c(sample_column, file_column)]
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
    tt <- sm(try(attachNamespace("foreach"), silent = TRUE))
    ## cores <- parallel::detectCores() / 2
    cores <- 4
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
    if (isTRUE(show_progress)) {
      bar <- utils::txtProgressBar(max = num_files, style = 3)
    }
    progress <- function(n) {
      setTxtProgressBar(bar, n)
    }
    pb_opts <- list()
    if (isTRUE(show_progress)) {
      pb_opts <- list("progress" = progress)
    }
    res_names <- c()
    res <- foreach(i = 1:num_files, .packages = c("hpgltools", "doParallel"),
                   .options.snow = pb_opts, .export = c("extract_scan_data")) %dopar% {
                     file <- meta[i, "file"]
                     id <- meta[i, "id"]
                     file_result <- try(extract_scan_data(file, id = id, write_acquisitions = write_windows,
                                                          allow_window_overlap = allow_window_overlap,
                                                          format = format, start_add = start_add))
                     if (class(file_result)[1] == "try-error") {
                       warning("There was an error reading ", file, ".")
                     } else {
                       res_names <- c(file, res_names)
                       returns[[file]] <- file_result
                     }
                   }
    if (isTRUE(show_progress)) {
      close(bar)
    }
    parallel::stopCluster(cl)
  } else {
    res_names <- c()
    for (i in 1:num_files) {
      file <- meta[i, "file"]
      id <- meta[i, "id"]
      file_result <- try(extract_scan_data(file, id = id, write_acquisitions = write_windows,
                                           allow_window_overlap = allow_window_overlap,
                                           format = format, start_add = start_add), silent = TRUE)
      if (class(file_result)[1] == "try-error") {
        warning("There was an error reading ", file, ".")
      } else {
        res_names <- c(file, res_names)
        res[[file]] <- file_result
      }
    }
  }
  try(names(res) <- res_names, silent = TRUE)
  rownames(sample_definitions) <- make.names(sample_definitions[[id_column]], unique = TRUE)

  retlist <- list(
    "colors" = chosen_colors,
    "metadata" = sample_definitions,
    "sample_data" = res)
  if (!is.null(savefile)) {
    mzxml_data <- retlist
    save_result <- try(save(list = c("mzxml_data"), file = savefile), silent = TRUE)
  }
  return(retlist)
}

#' Get some data from a peptideprophet run.
#'
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
#'  The columns are:
#'  * protein: The name of the matching sequence (DECOYs allowed here)
#'  * decoy: TRUE/FALSE, is this one of our decoys?
#'  * peptide: The sequence of the matching spectrum.
#'  * start_scan: The scan in which this peptide was observed
#'  * end scan: Ibid
#'  * index This seems to just increment
#'  * precursor_neutral_mass: Calculated mass of this fragment assuming no
#'    isotope shenanigans (yeah, looking at you C13).
#'  * assumed_charge: The expected charge state of this peptide.
#'  * retention_time_sec: The time at which this peptide eluted during the run.
#'  * peptide_prev_aa:  The amino acid before the match.
#'  * peptide_next_aa:  and the following amino acid.
#'  * num_tot_proteins: The number of matches not counting decoys.
#'  * num_matched_ions: How many ions for this peptide matched?
#'  * tot_num_ions:  How many theoretical ions are in this fragment?
#'  * matched_ion_ratio: num_matched_ions / tot_num_ions, bigger is better!
#'  * cal_neutral_pep_mass: This is redundant with precursor_neutral_mass, but
#'    recalculated by peptideProphet, so if there is a discrepency we should yell
#'    at someone!
#'  * massdiff How far off is the observed mass vs. the calculated? (also
#'    redundant with massd later)
#'  * num_tol_term: The number of peptide termini which are consistent with the
#'    cleavage (hopefully 2), but potentially 1 or even 0 if digestion was
#'    bad. (redundant with ntt later)
#'  * num_missed_cleavages: How many cleavages must have failed in order for this
#'    to be a good match?
#'  * num_matched_peptides: Number of alternate possible peptide matches.
#'  * xcorr: cross correlation of the experimental and theoretical spectra (this
#'    is supposedly only used by sequest, but I seem to have it here...)
#'  * deltacn: The normalized difference between the xcorr values for the best hit and next
#'    best hit.  Thus higher numbers suggest better matches.
#'  * deltacnstar: Apparently 'important for things like phospho-searches
#'    containing homologous top-scoring peptides when analyzed by
#'    peptideprophet...' -- the comet release notes.
#'  * spscore: The raw value of preliminary score from the sequest algorithm.
#'  * sprank: The rank of the match in a preliminary score. 1 is good.
#'  * expect: E-value of the given peptide hit.  Thus how many identifications
#'    one expect to observe by chance, lower is therefore better
#'  * prophet_probability: The peptide prophet probability score, higher is
#'   better.
#'  * fval: 0.6(the dot function + 0.4(the delta dot function) - (the dot bias
#'    penalty function) -- which is to say... well I dunno, but it is supposed to
#'    provide information about how similar this match is to other potential
#'    matches, so I presume higher means the match is more ambiguous.
#'  * ntt: Redundant with num_tol_term above, but this time from peptide prophet.
#'  * nmc: Redundant with num_missed_cleavages, except it coalesces them.
#'  * massd: Redundant with massdiff
#'  * isomassd: The mass difference, but taking into account stupid C13.
#'  * RT: Retention time
#'  * RT_score: The score of the retention time!
#'  * modified_peptides: A string describing modifications in the found peptide
#'  * variable_mods: A comma separated list of the variable modifications
#'    observed.
#'  * static_mods: A comma separated list of the static modifications observed.
#' @export
extract_peprophet_data <- function(pepxml, decoy_string = "DECOY_", ...) {
  input <- xml2::read_html(pepxml, options = "NOBLANKS")

  message("Extracting spectrum queries.")
  spectrum_queries <- rvest::xml_nodes(input, "spectrum_query")
  spectra <- spectrum_queries %>% rvest::html_attr("spectrum")
  query_data <- data.frame(row.names = spectra)
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
    rvest::html_node(xpath = "search_hit")
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
  decoy_idx <- grepl(pattern = decoy_regex, x = query_data[["protein"]])
  query_data[decoy_idx, "decoy"] <- TRUE
  query_data[["matched_ion_ratio"]] <- as.numeric(query_data[["num_matched_ions"]]) /
    as.numeric(query_data[["tot_num_ions"]])

  ## Get modification info
  message("Extracting modification metadata.")
  query_data[["modified_peptides"]] <- search_hits %>%
    rvest::html_node(xpath = "modification_info") %>%
    rvest::html_attr("modified_peptide")

  na_idx <- is.na(query_data[["modified_peptides"]])
  query_data[na_idx, "modified_peptides"] <- ""
  query_data[["variable_mods"]] <- ""
  query_data[["static_mods"]] <- ""
  modification_test <- search_hits %>%
    rvest::html_node(xpath = "modification_info")

  message("Filling in modification information, this is slow.")
  show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (i in 1:length(modification_test)) {
    if (isTRUE(show_progress)) {
      pct_done <- i / length(modification_test)
      setTxtProgressBar(bar, pct_done)
    }
    test <- modification_test[[i]]
    if (!is.na(test)) {
      variables <- test %>%
        rvest::html_nodes(xpath = "mod_aminoacid_mass") %>%
        rvest::html_attr("variable")
      statics <- test %>%
        rvest::html_nodes(xpath = "mod_aminoacid_mass") %>%
        rvest::html_attr("static")
      positions <- test %>%
        rvest::html_nodes(xpath = "mod_aminoacid_mass") %>%
        rvest::html_attr("position")
      masses <- test %>%
        rvest::html_nodes(xpath = "mod_aminoacid_mass") %>%
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
                                         tolerance = 0.1)
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
#'  does the iRT change this value?)
#' delta_rt:  The difference between rt and assay_rt
#' irt: (As described in the abstract of Claudia Escher's 2012 paper: "Here we
#'  present iRT, an empirically derived dimensionless peptide-specific value that
#'  allows for highly accurate RT prediction. The iRT of a peptide is a fixed
#'  number relative to a standard set of reference iRT-peptides that can be
#'  transferred across laboratories and chromatographic systems.")
#' assay_irt: The iRT observed in the actual chromatographic run.
#' delta_irt: The difference.  I am seeing that all the delta iRTs are in the
#'  -4000 range for our actual experiment; since this is in seconds, does that
#'  mean that it is ok as long as they stay in a similar range?
#' id: unique long signed integer for the peak group.
#' sequence: The sequence of the matched peptide
#' fullunimodpeptidename: The sequence, but with unimod formatted modifications
#'  included.
#' charge:  The assumed charge of the observed peptide.
#' mz:  The m/z value of the precursor ion.
#' intensity:  The sum of all transition intensities in the peak group.
#' aggr_prec_peak_area:  Semi-colon separated list of intensities (peak areas)
#'  of the MS traces for this match.
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
#' @param metadata Data frame describing the samples, including the mzXML
#'  filenames.
#' @param pyprophet_column Which column from the metadata provides the requisite filenames?
#' @param savefile If not null, save the data from this to the given filename.
#' @param ... Extra arguments, presumably color palettes and column names and
#'  stuff like that.
#' @return List of data from each sample in the pyprophet scored DIA run.
#' @export
extract_pyprophet_data <- function(metadata, pyprophet_column = "diascored",
                                   savefile = NULL, ...) {
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
  if ("data.frame" %in% class(metadata)) {
    sample_definitions <- metadata
  } else {
    sample_definitions <- extract_metadata(metadata,
                                           ...)
    ## sample_definitions <- extract_metadata(metadata)
  }

  if (is.null(sample_definitions[[pyprophet_column]])) {
    stop("This required a column with the tsv scored pyprophet data.")
  }

  chosen_colors <- generate_expt_colors(sample_definitions,
                                        ...)
  ## chosen_colors <- generate_expt_colors(sample_definitions)
  meta <- sample_definitions[, c("sampleid", pyprophet_column)]
  colnames(meta) <- c("id", "scored")
  existing_files <- complete.cases(meta[["scored"]])
  if (sum(existing_files) != nrow(meta)) {
    warning("It appears that some files are missing in the metadata.")
  }
  meta <- meta[existing_files, ]

  res <- list()
  num_files <- nrow(meta)
  failed_files <- c()
  for (i in 1:num_files) {
    file <- meta[i, "scored"]
    id <- meta[i, "id"]
    message("Attempting to read the tsv file for: ", id, ": ", file, ".")
    file_result <- sm(try(readr::read_tsv(file), silent = TRUE))
    if (class(file_result)[1] != "try-error") {
      colnames(file_result) <- tolower(colnames(file_result))
      file_result <- file_result %>%
        dplyr::rowwise() %>%
        dplyr::mutate(mass = gather_masses(sequence)) %>%
        dplyr::mutate(seqlength = nchar(sequence))
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
    save_result <- try(save(list = c("pyprophet_data"), file = savefile), silent = TRUE)
  }
  class(retlist) <- c("pyprophet_tables", "list")
  return(retlist)
}

#' Use BRAIN to find the peptide mass from a sequence.
#'
#' This rounds the avgMass from BRAIN to deal with isotopes, maybe this should be changed.
#'
#' @param sequence Sequence to count.
#' @return Rounded average mass.
gather_masses <- function(sequence) {
  atoms <- try(BRAIN::getAtomsFromSeq(sequence), silent = TRUE)
  if (class(atoms)[1] != "try-error") {
    d <- BRAIN::useBRAIN(atoms)
    ret <- round(d[["avgMass"]])
  } else {
    ret <- 0
  }
  return(ret)
}

#' Impute missing values using code from DEP reworked for expressionsets.
#'
#' [impute_expt()] imputes missing values in a proteomics dataset.
#'
#' @param expt An ExpressionSet (well, expt), I think it is assumed that this should have
#'  been normalized and filtered for features which have no values across 'most' samples.
#' @param filter Use normalize_expt() to filter the data?
#' @param p When filtering with pofa, use this p parameter.
#' @param fun "bpca", "knn", "QRILC", "MLE", "MinDet",
#'  "MinProb", "man", "min", "zero", "mixed" or "nbavg",
#'  Function used for data imputation based on
#'  [MSnbase::impute-methods()]
#' @param ... Additional arguments for imputation functions.
#' @return An imputed expressionset.
#' @seealso [MSnbase]
#' @export
impute_expt <- function(expt, filter = TRUE, p = 0.5,
                        fun = c("bpca", "knn", "QRILC", "MLE",
                              "MinDet", "MinProb", "min", "zero",
                              "mixed", "nbavg"), ...) {
  ## Show error if inputs do not contain required columns
  fun <- match.arg(fun)

  ## Caveat: Imputation works only on NA values.  I reset NAs to 0,
  ## so I will need to send them back...
  found_zeros <- sum(exprs(expt) == 0)
  if (found_zeros == 0) {
    warning("No missing values in the expressionset, returning it unchanged.")
    return(expt)
  } else {
    message("Found ", found_zeros, " zeros in the data.")
  }

  if (isTRUE(filter)) {
    if (expt[["state"]][["filter"]] == "raw") {
      message("The data has not been filtered.")
    } else {
      message("The data was already filtered with: ", expt[["state"]][["filter"]], ".")
    }
    message("Filtering the data, turn off 'filter' to stop this.")
    expt <- normalize_expt(expt, filter = "pofa", p = p)
  }

  exprs_set <- expt[["expressionset"]]
  ## Annotate whether or not there are missing values and how many
  num_zeros <- apply(exprs(expt) == 0, 1, any)
  fData(exprs_set)[["imputed"]] <- num_zeros
  zero_idx <- exprs(expt) == 0
  fData(exprs_set)[["num_nas"]] <- rowSums(zero_idx)
  exprs(exprs_set)[zero_idx] <- NA

  requireNamespace("MSnbase")
  msn_data <- as(exprs_set, "MSnSet")
  starting_counts <- exprs(exprs_set)
  message("Invoking impute from MSnbase with the ", fun, " method.")

  imputed_data <- MSnbase::impute(msn_data, method = fun,
                                  ...)
  imputed_exprs <- as(imputed_data, "ExpressionSet")
  imputed_counts <- exprs(imputed_exprs)

  same <- all.equal(starting_counts, imputed_counts)
  if (isTRUE(same)) {
    message("The counts remained the same.")
  }

  expt[["expressionset"]] <- imputed_exprs
  return(expt)
}

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
#'   na.rm = TRUE.  Thus the entries which used to be 0 should no longer affect
#'   the result.
#'   3.  Recreate the expressionset with the modified set of samples.
#'
#' @param expt Starting expressionset to mangle.
#' @param fact Metadata factor to use when taking the mean of biological
#'  replicates.
#' @param fun Assumed to be mean, but one might want median.
#' @return new expressionset
#' @export
mean_by_bioreplicate <- function(expt, fact = "bioreplicate", fun = "mean") {
  ## Set all the zeros to NA so that when we do cpm and mean they will get dropped.
  exprs_set <- expt[["expressionset"]]
  mtrx <- exprs(expt)
  zero_idx <- mtrx == 0
  new <- mtrx
  new[zero_idx] <- NA
  new_libsize <- colSums(new, na.rm = TRUE)
  new <- edgeR::cpm(new, lib.size = new_libsize)
  exprs(exprs_set) <- new
  expt[["expressionset"]] <- exprs_set
  annot <- fData(expt)
  final <- median_by_factor(expt, fact = fact, fun = fun)[["medians"]]
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
  new_set <- create_expt(count_dataframe = final, metadata = new_design, gene_info = annot)
  return(new_set)
}

#' Parse the difficult thermo fisher xlsx file.
#'
#' The Thermo(TM) workflow has as its default a fascinatingly horrible excel
#' output.  This function parses that into a series of data frames.
#'
#' @param xlsx_file The input xlsx file
#' @param test_row A single row in the xlsx file to use for testing, as I have
#'  not yet seen two of these accursed files which had the same headers.
#' @return List containing the protein names, group data, protein dataframe,
#'  and peptide dataframe.
#' @export
read_thermo_xlsx <- function(xlsx_file, test_row = NULL) {
  old_options <- options(java.parameters = "-Xmx20G")
  message("Reading ", xlsx_file)
  result <- readxl::read_xlsx(path = xlsx_file, sheet = 1, col_names = FALSE)
  group_data <- list()
  show_progress <- interactive() & is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style = 3)
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
      group_keepers <- !grepl(pattern = "^$", x = group_colnames)
      group_keepers[1] <- FALSE
      group_colnames <- group_colnames[group_keepers]
      next
    }
    ## When the 2nd column is 'Checked', then this defines a new protein in the group.
    if (row[, 2] == "Checked") {
      protein_colnames <- as.character(row)
      protein_keepers <- !grepl(pattern = "^$", x = protein_colnames)
      protein_keepers[2] <- FALSE
      protein_colnames <- protein_colnames[protein_keepers]
      next
    }
    ## When the 3rd column is 'Checked', then this starts a peptide definition
    if (row[, 3] == "Checked") {
      peptide_colnames <- as.character(row)
      peptide_keepers <- !grepl(pattern = "^$", x = peptide_colnames)
      peptide_keepers[3] <- FALSE
      peptide_colnames <- peptide_colnames[peptide_keepers]
      next
    }
    ## Once the column names for the data are defined, we consider how to
    ## Fill in the actual data, the protein group is probably the least interesting.
    if (row[, 1] == FALSE | row[, 1] == TRUE) {
      group_information <- row[group_keepers]
      colnames(group_information) <- group_colnames
      group_information[["ID"]] <- sub(pattern = "^.* GN=(\\w+) .*$",
                                       replacement = "\\1",
                                       x = group_information[["Group Description"]])
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
      protein_information[["ID"]] <- sub(pattern = "^.* GN=(\\w+) .*$",
                                         replacement = "\\1",
                                         x = protein_information[["Description"]])
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
    bar <- utils::txtProgressBar(style = 3)
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
    pattern = "%", replacement = "pct", x = current_colnames)
  ## as are spaces.
  current_colnames <- gsub(
    pattern = " ", replacement = "_", x = current_colnames)
  ## A bunch of columns have redundant adjectives.
  current_colnames <- gsub(
    pattern = "_confidence", replacement = "", x = current_colnames)
  ## Extra text in a column name is useless
  current_colnames <- gsub(
    pattern = "\\(by_search_engine\\)", replacement = "", x = current_colnames)
  ## Get rid of a bunch of doofusy punctuation.
  current_colnames <- gsub(
    pattern = "\\[|\\]|#|:|\\.|\\/|\\,|\\-", replacement = "", x = current_colnames)
  ## At this point we should not have any leading underscores.
  current_colnames <- gsub(
    pattern = "^_", replacement = "", x = current_colnames)
  ## Now should we have any double underscores.
  current_colnames <- gsub(
    pattern = "__", replacement = "_", x = current_colnames)
  ## Finally, because of the previous removals, there might be some duplicated
  ## terms left behind.
  current_colnames <- gsub(
    pattern = "_ht", replacement = "", x = current_colnames)
  current_colnames <- gsub(
    pattern = "_mascot_mascot", replacement = "_mascot", x = current_colnames)
  current_colnames <- gsub(
    pattern = "_sequest_sequest", replacement = "_sequest", x = current_colnames)
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
    if (grepl(pattern = "^found_in_", x = col)) {
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

#' Gather together the various SWATH2stats filters into one place.
#'
#' There are quite a few filters available in SWATH2stats.  Reading the
#' documentation, it seems at least possible, if not appropriate, to use them
#' together when filtering DIA data before passing it to MSstats/etc.  This
#' function attempts to formalize and simplify that process.
#'
#' @param s2s_exp SWHAT2stats result from the sample_annotation()
#'  function. (s2s_exp stands for: SWATH2stats experiment)
#' @param column What column in the data contains the protein name?
#' @param pep_column What column in the data contains the peptide name (not
#'  currently used, but it should be.)
#' @param fft Ratio of false negatives to true positives, used by
#'  assess_by_fdr() and similar functions.
#' @param plot Print plots of the various rates by sample?
#' @param target_fdr When invoking mscore4assayfdr, choose an mscore which
#'  corresponds to this false discovery date.
#' @param upper_fdr Used by filter_mscore_fdr() to choose the minimum threshold
#'  of identification confidence.
#' @param mscore Mscore cutoff for the mscore filter.
#' @param percentage Cutoff for the mscore_freqobs filter.
#' @param remove_decoys Get rid of decoys in the final filter, if they were not
#'  already removed.
#' @param max_peptides A maximum number of peptides filter.
#' @param min_peptides A minimum number of peptides filter.
#' @param do_mscore Perform the mscore filter? SWATH2stats::filter_mscore()
#' @param do_freqobs Perform the mscore_freqobs filter?
#'  SWATH2stats::filter_mscore_freqobs()
#' @param do_fdr Perform the fdr filter? SWATH2stats::filter_mscore_fdr()
#' @param do_proteotypic Perform the proteotypic filter?
#'  SWATH2stats::filter_proteotypic_peptides()
#' @param do_peptide Perform the single-peptide filter?
#'  SWATH2stats::filter_all_peptides()
#' @param do_max Perform the maximum peptide filter?
#'  SWATH2stats::filter_max_peptides()
#' @param do_min Perform the minimum peptide filter?
#'  SWATH2stats::filter_min_peptides()
#' @param ... Other arguments passed down to the filters.
#' @return Smaller SWATH2stats data set.
#' @seealso [SWATH2stats]
#' @export
s2s_all_filters <- function(s2s_exp, column = "proteinname", pep_column = "fullpeptidename",
                            fft = 0.7, plot = FALSE, target_fdr = 0.02, upper_fdr = 0.05,
                            mscore = 0.01, percentage = 0.75, remove_decoys = TRUE,
                            max_peptides = 15, min_peptides = 2,
                            do_mscore = TRUE, do_freqobs = TRUE, do_fdr = TRUE,
                            do_proteotypic = TRUE, do_peptide = TRUE,
                            do_max = TRUE, do_min = TRUE, ...) {
  retlist <- list()
  retlist[["decoy_lists"]] <- SWATH2stats::assess_decoy_rate(s2s_exp)
  decoy_ratio <- as.numeric(retlist[["decoy_lists"]][["ratio"]])
  decoy_number <- grepl(pattern = "^DECOY", x = s2s_exp[["proteinname"]])
  message("There were ", sum(!decoy_number),
          " observations and ", sum(decoy_number), " decoy observations.")
  retlist[["fdr_overall"]] <- SWATH2stats::assess_fdr_overall(
                                             s2s_exp,
                                             output = "Rconsole",
                                             plot = plot)
  retlist[["byrun_fdr"]] <- SWATH2stats::assess_fdr_byrun(
                                           s2s_exp,
                                           FFT = fft,
                                           plot = plot,
                                           output = "Rconsole")
  retlist[["chosen_mscore"]] <- SWATH2stats::mscore4assayfdr(
                                               s2s_exp,
                                               FFT = fft,
                                               fdr_target = target_fdr,
                                               ...)
  retlist[["prot_score"]] <- SWATH2stats::mscore4protfdr(
                                            s2s_exp,
                                            FFT = fft,
                                            fdr_target = target_fdr,
                                            ...)
  message("Starting mscore filter.")
  retlist[["raw"]] <- s2s_exp
  filt <- s2s_exp

  if (isTRUE(do_mscore)) {
    message("Starting mscore filter.")
    filt <- try(SWATH2stats::filter_mscore(
                               s2s_exp,
                               retlist[["chosen_mscore"]],
                               ...))
    if (class(filt)[1] == "try-error") {
      warning("The mscore filter failed, reverting to the raw data.")
      filt <- s2s_exp
      retlist[["mscore_filtered"]] <- "error"
    } else {
      retlist[["mscore_filtered"]] <- filt
    }
  } else {
    message("Skipping mscore filter.")
    retlist[["mscore_filtered"]] <- NULL
  }
  filt_backup <- filt

  if (isTRUE(do_freqobs)) {
    message("Starting freqobs filter.")
    filt <- try(SWATH2stats::filter_mscore_freqobs(
                               filt,
                               mscore,
                               percentage,
                               ...))
    if (class(filt)[1] == "try-error") {
      warning("The mscore filter failed, reverting to the mscore filtered data.")
      filt <- filt_backup
      retlist[["freqobs_filtered"]] <- "error"
    } else {
      retlist[["freqobs_filtered"]] <- filt
    }
  } else {
    message("Skipping freqobs filter.")
    retlist[["freqobs_filtered"]] <- NULL
  }
  filt_backup <- filt

  if (isTRUE(do_fdr)) {
    message("Starting fdr filter.")
    ## filter_mscore_fdr should probably be modified for flexibility.
    filt <- try(SWATH2stats::filter_mscore_fdr(
                               filt,
                               FFT = fft,
                               overall_protein_fdr_target = retlist[["prot_score"]],
                               upper_overall_peptide_fdr_limit = upper_fdr,
                               ...))
    if (class(filt)[1] == "try-error") {
      warning("The fdr filter failed, reverting to the freqobs filtered data.")
      filt <- filt_backup
      retlist[["fdr_filtered"]] <- "error"
    } else {
      retlist[["fdr_filtered"]] <- filt
    }
  } else {
    message("Skipping fdr filter.")
    retlist[["fdr_filtered"]] <- NULL
  }
  filt_backup <- filt

  if (isTRUE(do_proteotypic)) {
    message("Starting proteotypic filter.")
    filt <- try(SWATH2stats::filter_proteotypic_peptides(
                               filt,
                               column = column,
                               ...))
    if (class(filt)[1] == "try-error") {
      warning("The proteotypic filter failed, reverting to the fdr filtered data.")
      filt <- filt_backup
      retlist[["proteotypic_filtered"]] <- "error"
    } else {
      retlist[["proteotypic_filtered"]] <- filt
    }
  } else {
    message("Skipping proteotypic filter.")
    retlist[["proteotypic_filtered"]] <- NULL
  }
  filt_backup <- filt

  if (isTRUE(do_peptide)) {
    message("Starting peptide filter.")
    ## Looking at this function, it just renames the peptides to remove the 1/!
    ## That is not a filter!
    filt <- try(SWATH2stats::filter_all_peptides(
                               filt,
                               column = column))
    if (class(filt)[1] == "try-error") {
      warning("The peptide filter failed, reverting to the proteotypic filtered data.")
      filt <- filt_backup
      retlist[["peptide_filtered"]] <- "error"
    } else {
      retlist[["peptide_filtered"]] <- filt
    }
  } else {
    message("Skipping peptide filter.")
    retlist[["peptide_filtered"]] <- NULL
  }
  filt_backup <- filt

  if (isTRUE(do_max)) {
    message("Starting maximum peptide filter.")
    filt <- try(SWATH2stats::filter_on_max_peptides(
                               data = filt,
                               column = column,
                               n_peptides = max_peptides,
                               ...))
    if (class(filt)[1] == "try-error") {
      warning("The maximum peptide filter failed, reverting to the proteotypic filtered data.")
      filt <- filt_backup
      retlist[["maxpeptide_filtered"]] <- "error"
    } else {
      retlist[["maxpeptide_filtered"]] <- filt
    }
  } else {
    message("Skipping max peptide filter.")
    filt_backup <- filt
  }

  if (isTRUE(do_min)) {
    message("Starting minimum peptide filter.")
    filt <- try(SWATH2stats::filter_on_min_peptides(
                               data = filt,
                               n_peptides = min_peptides,
                               column = column,
                               rm.decoy = remove_decoys))
    if (class(filt)[1] == "try-error") {
      warning("The minimum peptide filter failed, reverting to the maximum peptide filtered data.")
      filt <- filt_backup
      retlist[["minpeptide_filtered"]] <- "error"
    } else {
      retlist[["minpeptide_filtered"]] <- filt
    }
  } else {
    message("Skipping min peptide filter.")
    retlist[["minpeptide_filtered"]] <- NULL
  }
  retlist[["final"]] <- filt

  start_proteins <- length(unique(s2s_exp[[column]]))
  start_peptides <- length(unique(s2s_exp[["fullpeptidename"]]))
  end_proteins <- length(unique(retlist[["final"]][[column]]))
  end_peptides <- length(unique(retlist[["final"]][["fullpeptidename"]]))
  message("We went from ", start_proteins, "/", start_peptides, " proteins/peptides to:")
  message("             ", end_proteins, "/", end_peptides, " proteins/peptides.")
  return(retlist)
}

subset_pyprophet_data <- function(lst, subset = NULL, column = "protein_id", operator = "in") {
  data <- lst[["sample_data"]]
  ## First, I want to get the simple Rv ID from the data, this is annoyingly Tb
  ## specific and should be removed or in some way made generic for other data.
  for (c in 1:length(data)) {
    datum <- data[[c]]
    datum[["protein_id"]] <- gsub(x = datum[["proteinname"]], pattern = "^1/",
                                  replacement = "")
    datum[["protein_id"]] <- gsub(x = datum[["protein_id"]], pattern = "_.*$",
                                  replacement = "")
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
