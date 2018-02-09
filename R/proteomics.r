#' Read a mzXML file and extract from it some important metadata.
#'
#' When working with swath data, it is fundamentally important to know the
#' correct values for a bunch of the input variables.  These are not trivial
#' to acquire.  This function attempts to make this easier (but slow) by reading
#' the mzXML file and using xml2 to parse and extract some hopefully helpful data.
#'
#' @param file  Filename to read.
#' @param id  An id to give the result.
#' @param write_acquisitions  If a filename is provided, write a tab separated table of windows.
#' @return List containing a table of scan and precursor data.
#' @export
extract_scan_data <- function(file, id=NULL, write_acquisitions=TRUE) {
  if (is.null(id)) {
    id <- file
  }
  message(paste0("Reading ", file))
  input <- xml2::read_html(file, options="NOBLANKS")
  ## peaks <- rvest::xml_nodes(input, "peaks")

  message(paste0("Extracting instrument information for ", file))
  instruments <- rvest::xml_nodes(input, "msinstrument")
  instrument_data <- data.frame(row.names=(1:length(instruments)), stringsAsFactors=FALSE)
  instrument_values <- c("msmanufacturer", "msmodel", "msionisation", "msmassanalyzer",
                         "msdetector")
  for (v in instrument_values) {
    datum <- rvest::xml_nodes(instruments, v)
    instrument_data[[v]] <- datum %>% rvest::html_attr("value")
  }
  datum <- rvest::xml_nodes(instruments, "software")
  instrument_data[["software_type"]] <- datum %>% rvest::html_attr("type")
  instrument_data[["software_name"]] <- datum %>% rvest::html_attr("name")
  instrument_data[["software_version"]] <- datum %>% rvest::html_attr("version")

  message(paste0("Extracting scan information for ", file))
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

  message(paste0("Extracting precursor information for ", file))
  precursors <- rvest::xml_nodes(scans, "precursormz")
  precursor_data <- data.frame(row.names=(1:length(precursors)), stringsAsFactors=FALSE)
  precursor_wanted <- c("precursorintensity", "activationmethod",
                        "windowwideness", "precursorscannum")
  precursor_numeric <- c("precursorintensity", "precursorscannum", "window_center", "windowwideness")
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

  message(paste0("Coalescing the acquisition windows for ", file))
  acquisition_windows <- precursor_data[, c("window_start", "window_end")]
  acquisition_unique <- !duplicated(x=acquisition_windows)
  acquisition_windows <- acquisition_windows[acquisition_unique, ]
  colnames(acquisition_windows) <- c("start", "end")
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
    message(paste0("Hopefully writing acquisition file to ", pre_file))
    no_cols <- write.table(x=acquisition_windows, file=pre_file, sep="\t", quote=FALSE,
                           row.names=FALSE, col.names=FALSE)
    osw_file <- file.path(acq_dir, paste0("openswath_", acq_file))
    ## This is the file for openswathworkflow.
    message(paste0("Hopefully writing osw acquisitions to ", osw_file))
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

#' Read a bunch of mzXML files to acquire their metadata.
#'
#' I have had difficulties getting the full set of correct parameters for a
#' DDA/DIA experiment.  After some poking, I eventually found most of these
#' required prameters in the mzXML raw files.  Ergo, this function uses them.
#'
#' @param metadata  Data frame describing the samples, including the mzXML
#'   filenames.
#' @param write_windows  Write out SWATH window frames.
#' @param ... Extra arguments, presumably color palettes and column names and
#'   stuff like that.
#' @return  metadata!#'
#' @export
extract_mzxml_data <- function(metadata, write_windows=TRUE, parallel=TRUE, ...) {
  arglist <- list(...)

  ## Add a little of the code from create_expt to include some design information in the returned
  ## data structure.

  ## Since I am using this code now in two places, if I can use it without changes, I will
  ## split it into its own function so make changes in one place happen in both
  ## but until I can ensure to myself that it is functional in both contexts I will not.
  ## Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  file_column <- "file"
  if (!is.null(arglist[["file_column"]])) {
    file_column <- arglist[["file_column"]]  ## Make it possible to have multiple count
    file_column <- tolower(file_column)
    file_column <- gsub(pattern="[[:punct:]]", replacement="", x=file_column)
    ## tables / sample in one sheet.
  }

  sample_column <- "sampleid"
  if (!is.null(arglist[["sample_column"]])) {
    sample_column <- arglist[["sample_column"]]
    sample_column <- tolower(sample_column)
    sample_column <- gsub(pattern="[[:punct:]]", replacement="", x=sample_column)
  }

  sample_definitions <- extract_metadata(metadata, ...)
  ## sample_definitions <- extract_metadata(metadata)
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
    bar <- utils::txtProgressBar(max=num_files, style=3)
    progress <- function(n) {
      setTxtProgressBar(bar, n)
    }
    pb_opts <- list(progress=progress)
    res <- foreach(i=1:num_files, .packages=c("hpgltools", "doParallel"), .options.snow=pb_opts, .export=c("extract_scan_data")) %dopar% {
      file <- meta[i, "file"]
      id <- meta[i, "id"]
      returns[[file]] <- try(extract_scan_data(file, id=id))
    }
      close(bar)
      parallel::stopCluster(cl)
  } else {
    for (i in 1:num_files) {
      file <- meta[i, "file"]
      id <- meta[i, "id"]
      res[[file]] <- try(extract_scan_data(file, id=id))
    }
  }

  retlist <- list(
    "colors" = chosen_colors,
    "metadata" = sample_definitions,
    "sample_data" = res)
  return(retlist)
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

#' Plot the peak intensities with respect to m/z
#'
#' I want to have a pretty plot of peak intensities and m/z.
#'
#' @param mzxml_data  The data structure from extract_mzxml or whatever it is.
#' @param loess  Do a loess smoothing from which to extract a function
#'   describing the data?  This is terribly slow, and in the data I have
#'   examined so far, not very helpful, so it is FALSE by default.
#' @param ...  Extra arguments for the downstream functions.
#' @return  ggplot2 goodness.
#' @export
plot_intensity_mz <- function(mzxml_data, loess=FALSE, alpha=0.5, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)

  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    message(paste0("Adding ", name))
    plotted_table <- sample_data[[i]][["scans"]]
    plotted_data <- plotted_table[, c("basepeakmz", "basepeakintensity")]
    plotted_data[["sample"]] <- name
    plotted_data <- plotted_data[, c("sample", "basepeakmz", "basepeakintensity")]
    colnames(plotted_data) <- c("sample", "mz", "intensity")
    ## Re-order the columns because I like sample first.
    plot_df <- rbind(plot_df, plotted_data)
  }
  chosen_palette <- "Dark2"
  sample_colors <- sm(
    grDevices::colorRampPalette(
                 RColorBrewer::brewer.pal(samples, chosen_palette))(samples))

  ## Randomize the rows of the df so we can see if any sample is actually overrepresented
  plot_df <- plot_df[sample(nrow(plot_df)), ]

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
  if (isTRUE(lowess)) {
    int_vs_mz <- int_vs_mz +
      ggplot2::geom_smooth(method="loess", size=1.0)
  }
  return(int_vs_mz)
}

#' Make a boxplot out of some of the various data available in the mzxml data.
#'
#' There are a few data within the mzXML raw data files which are likely
#'   candidates for simple summary via a boxplot/densityplot/whatever.  For the
#'   moment I am just doing boxplots of a few of them.  Since my metadata
#'   extractor dumps a couple of tables, one must choose a desired table and
#'   column from it to plot.
#'
#' @param table  One of precursors or scans
#' @param column  One of the columns from the table; if 'scans' is chosen, then
#'   likely choices include: 'peakscount', 'basepeakmz', 'basepeakintensity'; if
#'   'precursors' is chosen, then the only likely choice for the moment is
#'   'precursorintensity'.
#' @param names  Names for the x-axis of the plot.
#' @param title  Title the plot?
#' @param scale  Put the data on a specific scale?
#' @param ...  Further arguments, presumably for colors or some such.
#' @return  Boxplot goodness!
#' @export
plot_mzxml_boxplot <- function(mzxml_data, table="precursors", column="precursorintensity",
                               names=NULL, title=NULL, scale=NULL, ...) {
  arglist <- list(...)
  metadata <- mzxml_data[["metadata"]]
  colors <- mzxml_data[["colors"]]
  sample_data <- mzxml_data[["sample_data"]]
  plot_df <- data.frame()
  samples <- length(sample_data)
  for (i in 1:samples) {
    name <- metadata[i, "sampleid"]
    message(paste0("Adding ", name))
    plotted_table <- sample_data[[i]][[table]]
    plotted_data <- as.data.frame(plotted_table[[column]])
    plotted_data[["sample"]] <- name
    colnames(plotted_data) <- c(column, "sample")
    ## Re-order the columns because I like sample first.
    plotted_data <- plotted_data[, c("sample", column)]
    plot_df <- rbind(plot_df, plotted_data)
  }

  scale_data <- check_plot_scale(plot_df[[column]], scale)
  scale <- scale_data[["scale"]]
  plot_df[[column]] <- scale_data[["data"]]

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
