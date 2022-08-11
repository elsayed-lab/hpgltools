#' Print some data onto KEGG pathways.
#'
#' KEGGREST and pathview provide neat functions for coloring molecular pathways
#' with arbitrary data. Unfortunately they are somewhat evil to use.  This
#' attempts to alleviate that.
#'
#' @param path_data Some differentially expressed genes.
#' @param indir Directory into which the unmodified kegg images will be
#'  downloaded (or already exist).
#' @param outdir Directory which will contain the colored images.
#' @param pathway Perform the coloring for a specific pathway?
#' @param species Kegg identifier for the species of interest.
#' @param from_list Regex to help in renaming KEGG categories/gene names from
#'  one format to another.
#' @param to_list Regex to help in renaming KEGG categories/gene names from one
#'  format to another.
#' @param suffix Add a suffix to the completed, colored files.
#' @param filenames Name the final files by id or name?
#' @param fc_column What is the name of the fold-change column to extract?
#' @param format Format of the resulting images, I think only png really works
#'  well.
#' @param verbose When on, this function is quite chatty.
#' @return A list of some information for every KEGG pathway
#'  downloaded/examined.  This information includes:
#'   a. The filename of the final image for each pathway.
#'   b. The number of genes which were found in each pathway image.
#'   c. The number of genes in the 'up' category
#'   d. The number of genes in the 'down' category
#' @seealso [pathview] [KEGGREST]
#' @examples
#' \dontrun{
#'  thy_el_comp2_path = hpgl_pathview(thy_el_comp2_kegg, species = "spz", indir = "pathview_in",
#'                                    outdir = "kegg_thy_el_comp2", string_from = "_Spy",
#'                                    string_to = "_Spy_", filenames = "pathname")
#' }
#' @export
simple_pathview <- function(gene_input = NULL, compound_input = NULL,
                            indir = "pathview_in", outdir = "pathview",
                            pathway = "all", species = "lma", from_list = NULL,
                            to_list = NULL, suffix = "_colored", id_column = "kegg_ids",
                            filenames = "id", fc_column = "limma_logfc",
                            format = "png", verbose = TRUE) {

  ## In the most recent iteration of using this, we were asked to create images of
  ## both compounds and genes.  Thus we will separate the input into
  ## two variables as per pathview().
  if (is.null(gene_input) & is.null(compound_input)) {
    stop("This requires at least one of gene_input or compound_input to exist.")
  }

  ## Please note that the KGML parser fails if other XML parsers are loaded into R
  ## eh = new.env(hash = TRUE, size = NA)
  ## There is a weird namespace conflict when using pathview, so I will reload it here
  tmp <- sm(library(pathview)) ## I am not sure how else to avoid the error
  ## 'unable to load 'bods''
  ## If a table from limma was passed to this, just assume that one wants logFC
  ## Similar stanzas should probably be added for deseq/edger
  ## This is added because pathview() only works with dataframes/lists with only numbers.
  ## So it is wise to pull only the column of numbers one cares about.
  if (class(gene_input) == "data.frame") {
    if (is.null(gene_input[[fc_column]])) {
      fc_column <- "logFC"
    }
    if (is.null(gene_input[[fc_column]])) {
      stop("Unable to find the fold change column.")
    }

    ## Extract the set of FC values and give it names appropriate for KEGG
    ## Oh, crap in a hat, also keep in mind that KEGG does not use the
    ## three letter prefix internally.  A fact which I annoyingly forgot.
    tmp_data <- as.vector(gene_input[[fc_column]])
    if (is.null(id_column)) {
      names(tmp_data) <- rownames(gene_input)
    } else {
      names(tmp_data) <- gene_input[[id_column]]
    }
  } else { ## If it is not a dataframe, assume it is a named vector.
    tmp_data <- gene_input
  }

  ## This little section is weird, that is because pathview does not handle non-numeric
  ## representations of the data well.
  ## So we must make certain that the FCs are numeric vectors with the names maintained.
  input_vector <- as.numeric(tmp_data)
  names(input_vector) <- names(tmp_data)
  tmp_names <- names(input_vector)

  ## Reset the names to KEGG standard names.
  if (is.null(from_list)) {
    substitutions <- get_kegg_sub(species)
    from_list <- substitutions[["patterns"]]
    to_list <- substitutions[["replaces"]]
  }
  if (!is.null(from_list)) {
    for (sub_count in seq_along(from_list)) {
      my_from <- from_list[[sub_count]]
      my_to <- to_list[[sub_count]]
      tmp_names <- gsub(pattern = my_from, replacement = my_to, x = tmp_names, perl = TRUE)
    }
    names(path_data) <- tmp_names
    rm(tmp_names)
  }

  ## First check that the input pathview directory exists
  if (!file.exists(indir)) {
    dir.create(indir)
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  paths <- list()
  if (pathway[1] == "all") {
    all_pathways <- unique(KEGGREST::keggLink("pathway", species))
    paths <- all_pathways
    paths <- gsub(pattern = "path:", replacement = "", x = paths)
  } else {
    paths <- pathway
  }
  return_list <- list()
  for (count in seq_along(paths)) {
    total_mapped <- 0
    total_pct_mapped <- 0
    unique_mapped <- 0
    unique_pct_mapped <- 0
    path <- paths[[count]]
    canonical_path <- path
    if (!grepl(pattern = glue("^{species}"), x = canonical_path)) {
      canonical_path <- glue("{species}{path}")
    }
    path_data <- try(KEGGREST::keggGet(canonical_path), silent = TRUE)
    if (class(path_data) == "try-error") {
      next
    }
    path_name <- path_data[[1]][["NAME"]]
    path_name <- gsub("(.*) - .*", "\\1", path_name)
    path_name <- tolower(path_name)
    path_name <- gsub(" ", "_", path_name)
    ## RCurl is crap and fails sometimes for no apparent reason.
    gene_examples <- try(KEGGREST::keggLink(paste("path", canonical_path, sep = ":"))[, 2])
    ## limits = c(min(path_data, na.rm = TRUE), max(path_data, na.rm = TRUE))
    limit_test <- c(abs(min(as.numeric(input_vector), na.rm = TRUE)),
                    abs(max(as.numeric(input_vector), na.rm = TRUE)))
    limit_min <- -1.0 * max(limit_test)
    limit_max <- max(limit_test)
    limits <- c(limit_min, limit_max)
    example_string <- gsub(pattern = glue("{species}:"), replacement = "",
                           x = toString(head(gene_examples)))
    if (isTRUE(verbose)) {
      if (count < 4) {
        ## Test if we have overlaps
        overlap_test <- glue::glue("{species}:{names(input_vector)}")
        num_overlap <- sum(gene_examples %in% overlap_test)
        message("Here are some path gene examples: ", example_string)
        message("Here are your genes: ", toString(head(names(input_vector))))
        message("There were ", num_overlap, " overlapping genes observed.")
      }
    }
    if (format == "png") {
      ## In this invocation, include all the possible arguments for debugging.
      pv <- suppressWarnings(
          try(pathview::pathview(gene.data = input_vector,
                                 cpd.data = compound_input,
                                 pathway.id = canonical_path,
                                 species = species,
                                 kegg.dir = indir,
                                 cpd.idtype = "kegg",
                                 gene.idtype = "KEGG",
                                 gene.annotpkg = NULL,
                                 min.nnodes = 3,
                                 kegg.native = TRUE,
                                 map.null = TRUE,
                                 expand.node = FALSE, ## Was true
                                 split.group = FALSE,
                                 map.symbol = TRUE,
                                 map.cpdname = TRUE,
                                 node.sum = "sum",
                                 discrete = list(gene = FALSE, cpd = FALSE),
                                 limit = list(gene = limits, cpd = limits),
                                 bins = list(gene = 10, cpd = 10),
                                 both.dirs = list(gene = TRUE, cpd = TRUE),
                                 trans.fun = list(gene = NULL, cpd = NULL),
                                 low = list(gene = "red", cpd = "yellow"),
                                 mid = list(gene = "gray", cpd = "gray"),
                                 high = list(gene = "green", cpd = "blue"),
                                 na.col = "transparent",
                                 out.suffix = suffix,
                                 same.layer = FALSE,
                                 res = 1200,
                                 new.signature = FALSE,
                                 cex = 0.05,
                                 key.pos = "topright")))
    } else { ## The other format is svg, and at least in the past it was rather more picky.
      pv <- suppressWarnings(
          try(pathview::pathview(gene.data = input_vector,
                                 kegg.dir = indir,
                                 pathway.id = canonical_path,
                                 species = species,
                                 limit = list(gene = limits, cpd = limits),
                                 map.null = TRUE,
                                 gene.idtype = "KEGG",
                                 out.suffix = suffix,
                                 split.group = TRUE,
                                 expand.node = TRUE,
                                 kegg.native = FALSE,
                                 map.symbol = TRUE,
                                 same.layer = FALSE,
                                 res = 1200,
                                 new.signature = FALSE,
                                 cex = 0.05,
                                 key.pos = "topright")))
    }
    if (class(pv) == "numeric") {
      warning(glue("There was a failure for: {canonical_path}."))
      colored_genes <- NULL
      newfile <- NULL
      up <- NULL
      down <- NULL
    } else {
      filetype <- ".png"
      if (format != "png") {
        filetype <- ".pdf"
      }
      colored_genes <- dim(pv[["plot.data.gene"]])[1]
      ## "lma04070._proeff.png"
      oldfile <- glue::glue("{path}.{suffix}{filetype}")
      ## An if-statement to see if the user prefers pathnames by kegg ID or pathway name
      ## Dr. McIver wants path names...
      newfile <- NULL
      if (filenames == "id") {
        newfile <- file.path(outdir, glue::glue("{path}{suffix}{filetype}"))
        ##newfile <- paste(outdir, "/", path, suffix, filetype, sep = "")
      } else {
        ## If filenames is not 'id', put in the path name...
        ##newfile <- paste(outdir, "/", path_name, suffix, filetype, sep = "")
        newfile <- file.path(outdir, glue::glue("{path_name}{suffix}{filetype}"))
      }
      rename_try <- try(file.rename(from = oldfile, to = newfile), silent = TRUE)
      if (class(rename_try)[1] == "try-error") {
        warning("There was an error renaming the png, ",
                "likely because it didn't download properly.")
        warning("It is likely easiest to just delete the pathview input directory.")
        newfile <- "undefined"
      }
      data_low <- summary(input_vector)[2]
      data_high <- summary(input_vector)[3]
      pathway_data <- as.data.frame(pv[["plot.data.gene"]])

      total_mapped <- pathway_data[["all.mapped"]] != ""
      total_pct_mapped <- signif(mean(total_mapped) * 100.0, 4)
      unique_pathway <- pathway_data
      rownames(unique_pathway) <- make.names(unique_pathway[["kegg.names"]], unique = TRUE)
      unique_path_genes <- make.names(unique(unique_pathway[["kegg.names"]]))
      unique_pathway <- unique_pathway[unique_path_genes, ]
      unique_mapped <- unique_pathway[["all.mapped"]] != ""
      unique_pct_mapped <- signif(mean(unique_mapped) * 100.0, 4)

      numbers_in_plot <- as.numeric(pv[["plot.data.gene"]][["mol.data"]])
      up <- sum(numbers_in_plot >= data_high, na.rm = TRUE)
      down <- sum(numbers_in_plot < data_low, na.rm = TRUE)

    }
    return_list[[path]][["file"]] <- newfile
    return_list[[path]][["genes"]] <- toString(colored_genes)
    return_list[[path]][["up"]] <- toString(up)
    return_list[[path]][["down"]] <- toString(down)
    return_list[[path]][["total_mapped_nodes"]] <- sum(total_mapped)
    return_list[[path]][["total_mapped_pct"]] <- total_pct_mapped
    return_list[[path]][["unique_mapped_nodes"]] <- sum(unique_mapped)
    return_list[[path]][["unique_mapped_pct"]] <- unique_pct_mapped
    if (isTRUE(verbose)) {
      message(count, "/", length(paths), ": Finished ",
              path_name, " id: ", path, " with ", total_pct_mapped, "% genes mapped(",
              sum(unique_mapped), " unique).")
    }
  } ## End for loop

  retdf <- data.frame(rep(NA, length(names(return_list))))
  rownames(retdf) <- names(return_list)
  retdf[["genes"]] <- NA
  retdf[["up"]] <- NA
  retdf[["down"]] <- NA
  retdf[["total_mapped_nodes"]] <- NA
  retdf[["total_mapped_pct"]] <- NA
  retdf[["unique_mapped_nodes"]] <- NA
  retdf[["unique_mapped_pct"]] <- NA
  colnames(retdf) <- c("file", "genes", "up", "down", "total_mapped_nodes",
                       "total_mapped_pct", "unique_mapped_nodes", "unique_mapped_pct")
  for (p in seq_along(return_list)) {
    path <- names(return_list)[p]
    if (is.null(return_list[[path]][["genes"]])) {
      retdf[path, "genes"] <- 0
    } else {
      retdf[path, "genes"] <- as.numeric(return_list[[path]][["genes"]])
    }
    if (!is.null(return_list[[path]][["file"]])) {
      retdf[path, "file"] <- as.character(return_list[[path]][["file"]])
      retdf[path, "up"] <- as.numeric(return_list[[path]][["up"]])
      retdf[path, "down"] <- as.numeric(return_list[[path]][["down"]])
      retdf[path, "total_mapped_nodes"] <- as.numeric(return_list[[path]][["total_mapped_nodes"]])
      retdf[path, "total_mapped_pct"] <- as.numeric(return_list[[path]][["total_mapped_pct"]])
      retdf[path, "unique_mapped_nodes"] <- as.numeric(return_list[[path]][["unique_mapped_nodes"]])
      retdf[path, "unique_mapped_pct"] <- as.numeric(return_list[[path]][["unique_mapped_pct"]])
      }
  }
  retdf <- retdf[with(retdf, order(up, down)), ]
  keepers <- !is.na(retdf[["total_mapped_pct"]])
  retdf <- retdf[keepers, ]
  return(retdf)
}

#' Gather all Compounds from all pathways for a given species.
#'
#' This function attempts to iterate over every pathway for a given
#' abbreviation/species and extract from them the set of compounds.
#' This was mostly copy/pasted from get_kegg_genes.
#'
#' @param pathway One or more pathways, all does what it says on the tin.
#' @param abbreviation Approximately 3 character KEGG abbreviation.
#' @param species If you do not have the abbreviation, this will try
#' to find it.
#' @param savefile Currently unused I think, but eventually should
#' make a savefile of the results.
#' @export
get_kegg_compounds <- function(pathway = "all", abbreviation = NULL,
                               species = "leishmania major", savefile = NULL) {
  if (is.null(abbreviation) & is.null(species)) {
    stop("This requires either a species or 3 letter kegg id.")
  } else if (is.null(abbreviation)) {
    ## Then the species was provided.
    abbreviation <- get_kegg_orgn(species)
    message("The abbreviation detected was: ", abbreviation)
  }

  result <- NULL
  species <- gsub(pattern = " ", replacement = "_", x = as.character(species))
  savefile <- glue("kegg_{species}.rda.xz")
  kegg_data <- NULL
  if (file.exists(savefile)) {
    message("Reading from the savefile, delete ", savefile, " to regenerate.")
    result <- new.env()
    load(savefile, envir = result)
    result <- result[["result"]]
  } else {
    paths <- list()
    if (pathway == "all") {
      all_pathways <- unique(KEGGREST::keggLink("pathway", abbreviation))
      paths <- all_pathways
      paths <- gsub("path:", "", paths)
      ## all_modules = unique(KEGGREST::keggLink("module", abbreviation))
    } else if (class(pathway) == "list") {
      paths <- pathway
    } else {
      paths[1] <- pathway
    }
    total_genes <- 0
    result <- data.frame()
    for (count in seq_along(paths)) {
      path <- paths[count]
      message("Extracting: ", path, ".")

      path_data <- KEGGREST::keggGet(path)
      kegg_class <- path_data[[1]]$CLASS
      if (is.null(kegg_class)) {
        kegg_class <- ""
      }
      kegg_description <- path_data[[1]]$DESCRIPTION
      if (is.null(kegg_description)) {
        kegg_description <- ""
      }
      kegg_name <- path_data[[1]]$NAME
      kegg_name <- gsub("(.*) - .*", "\\1", kegg_name)
      kegg_name <- tolower(kegg_name)
      kegg_name <- gsub(" ", "_", kegg_name)

      compounds <- path_data[[1]]$COMPOUND
      if (is.null(compounds)) {
        next
      }
      this_df <- as.data.frame(t(rbind(names(compounds), as.character(compounds))))
      colnames(this_df) <- c("kegg_compound_id", "compound_name")
      this_df[["kegg_pathway_id"]] <- path
      this_df[["pathway_name"]] <- kegg_name
      if (count == 1) {
        result <- this_df
      } else {
        result <- rbind(result, this_df)
      }
    } ## End iterating over pathways
  }
  return(result)
}

#' Extract the set of geneIDs matching pathways for a given species.
#'
#' This uses KEGGREST to extract the mappings for all genes for a species and
#' pathway or 'all'. Because downloading them takes a while, it will save the
#' results to kegg_species.rda.  When run interactively, it will give some
#' information regarding the number of genes observed in each pathway.
#'
#' @param pathway Either a single pathway kegg id or 'all'.
#' @param abbreviation Optional 3 letter species kegg id.
#' @param species Stringified species name used to extract the 3 letter abbreviation.
#' @param savefile Filename to which to save the relevant data.
#' @return Dataframe of the various kegg data for each pathway, 1 row/gene.
#' @seealso [KEGGREST]
#' @examples
#' \dontrun{
#'  kegg_info <- get_kegg_genes(species = "Canis familiaris")
#' }
#' @export
get_kegg_genes <- function(pathway = "all", abbreviation = NULL,
                           species = "leishmania major", savefile = NULL) {
  if (is.null(abbreviation) & is.null(species)) {
    stop("This requires either a species or 3 letter kegg id.")
  } else if (is.null(abbreviation)) {
    ## Then the species was provided.
    abbreviation <- get_kegg_orgn(species)
    message("The abbreviation detected was: ", abbreviation)
  }

  result <- NULL
  species <- gsub(pattern = " ", replacement = "_", x = as.character(species))
  savefile <- glue("kegg_{species}.rda.xz")
  kegg_data <- NULL
  if (file.exists(savefile)) {
    message("Reading from the savefile, delete ", savefile, " to regenerate.")
    result <- new.env()
    load(savefile, envir = result)
    result <- result[["result"]]
  } else {
    paths <- list()
    if (pathway == "all") {
      all_pathways <- unique(KEGGREST::keggLink("pathway", abbreviation))
      paths <- all_pathways
      paths <- gsub("path:", "", paths)
      ## all_modules = unique(KEGGREST::keggLink("module", abbreviation))
    } else if (class(pathway) == "list") {
      paths <- pathway
    } else {
      paths[1] <- pathway
    }
    total_genes <- 0
    result <- NULL
    for (count in seq_along(paths)) {
      path <- paths[count]
      path_name <- KEGGREST::keggGet(path)
      kegg_class <- path_name[[1]]$CLASS
      if (is.null(kegg_class)) {
        kegg_class <- ""
      }
      kegg_description <- path_name[[1]]$DESCRIPTION
      if (is.null(kegg_description)) {
        kegg_description <- ""
      }
      kegg_name <- path_name[[1]]$NAME
      kegg_name <- gsub("(.*) - .*", "\\1", kegg_name)
      kegg_name <- tolower(kegg_name)
      kegg_name <- gsub(" ", "_", kegg_name)
      ## RCurl is crap and fails sometimes for no apparent reason.
      kegg_geneids <- try(KEGGREST::keggLink(paste("path", path, sep = ":"))[, 2])
      kegg_geneids <- gsub(pattern = "^.*:", replacement = "", x = kegg_geneids)
      kegg_subst <- get_kegg_sub(abbreviation)
      tritryp_geneids <- kegg_geneids
      total_genes <- total_genes + length(kegg_geneids)
      patterns <- kegg_subst[["patterns"]]
      replaces <- kegg_subst[["replaces"]]
      message(count, "/", length(paths), ": Working on path: ",
              path, " which has ", length(kegg_geneids), " genes.")
      for (r in seq_along(kegg_subst[["patterns"]])) {
        tritryp_geneids <- gsub(pattern = patterns[r], replacement = replaces[r],
                                x = tritryp_geneids)
      }

      for (s in seq_along(tritryp_geneids)) {
        result <- rbind(result, data.frame(
                                    "GID" = tritryp_geneids[s],
                                    "KEGG_NAME" = kegg_name,
                                    "KEGG_CLASS" = kegg_class,
                                    "KEGG_DESCRIPTION" = kegg_description,
                                    "KEGG" = path))
      }
    } ## End iterating over pathways
    message("Total genes observed: ", total_genes)
    save(result, file = savefile)
  }
  return(result)
}

#' Provide a set of simple substitutions to convert geneIDs from KEGG->TriTryDB
#'
#' This function should provide 2 character lists which, when applied
#' sequentially, will result in a hopefully coherent set of mapped gene IDs
#' matching the TriTypDB/KEGG specifications.
#'
#' @param species 3 letter abbreviation for a given kegg type
#' @return 2 character lists containing the patterns and replace arguments for
#'  gsub(), order matters!
#' @seealso [KEGGREST]
#' @export
get_kegg_sub <- function(species = "lma") {
  patterns <- c()
  replaces <- c()
  if (species == "lma") {
    patterns <- c("LmjF",  "LMJF", "_")
    replaces <- c("LMJF_", "LmjF", ".")
  } else if (species == "tcr") {
    patterns <- c("TcCLB.")
    replaces <- c("")
  }

  ret <- list(
      "patterns" = patterns,
      "replaces" = replaces)
  return(ret)
}

#' Use gostats() against kegg pathways.
#'
#' This sets up a GSEABase analysis using KEGG pathways rather than gene
#' ontologies. Does this even work?  I don't think I have ever tested it yet.
#' oh, it sort of does, maybe if I export it I will rembmer it.
#'
#' @param organism The organism used to make the KEGG frame, human readable no taxonomic.
#' @param pathdb Name of the pathway database for this organism.
#' @param godb Name of the ontology database for this organism.
#' @return Results from hyperGTest using the KEGG pathways.
#' @seealso [AnnotationDbi] [GSEABase] [Category]
#' @export
gostats_kegg <- function(organism = "Homo sapiens",
                         pathdb = "org.Hs.egPATH", godb = "org.Hs.egGO") {
  org <- get0(pathdb)
  org_go <- get0(godb)
  frame <- AnnotationDbi::toTable(org)
  keggframedata <- data.frame(frame[["path_id"]], frame[["gene_id"]])
  keggFrame <- AnnotationDbi::KEGGFrame(keggframedata, organism = organism)
  gsc <- GSEABase::GeneSetCollection(keggFrame, setType = GSEABase::KEGGCollection())
  universe <- AnnotationDbi::Lkeys(org_go)
  genes <- universe[1:500]
  kparams <- Category::GSEAKEGGHyperGParams(name = "My Custom GSEA based annot Params",
                                            geneSetCollection = gsc,
                                            geneIds = genes,
                                            universeGeneIds = universe,
                                            pvalueCutoff = 0.05,
                                            testDirection = "over")
  kOver <- Category::hyperGTest(kparams)
  return(kOver)
}

#' Search KEGG identifiers for a given species name.
#'
#' KEGG identifiers do not always make sense.  For example, how am I supposed to
#' remember that Leishmania major is lmj?  This takes in a human readable string
#' and finds the KEGG identifiers that match it.
#'
#' @param species Search string (Something like 'Homo sapiens').
#' @param short Only pull the orgid?
#' @return Data frame of possible KEGG identifier codes, genome ID numbers,
#'  species, and phylogenetic classifications.
#' @seealso [RCurl]
#' @examples
#' \dontrun{
#'  fun = get_kegg_orgn('Canis')
#'  ## >     Tid     orgid      species                   phylogeny
#'  ## >  17 T01007   cfa Canis familiaris (dog) Eukaryotes;Animals;Vertebrates;Mammals
#' }
#' @export
get_kegg_orgn <- function(species = "Leishmania", short = TRUE) {
  all_organisms <- RCurl::getURL("https://rest.kegg.jp/list/organism")
  org_tsv <- textConnection(all_organisms)
  all <- read.table(org_tsv, sep = "\t", quote = "", fill = TRUE)
  close(org_tsv)
  colnames(all) <- c("Tid", "orgid", "species", "phylogeny")
  ## Look for the search string in the species column
  candidates <- all[grepl(species, all[["species"]]), ]
  if (isTRUE(short)) {
    candidates <- as.character(candidates[["orgid"]])
  }
  return(candidates)
}

#' Extract the percent differentially expressed genes for all KEGG pathways.
#'
#' KEGGgraph provides some interesting functionality for mapping KEGGids and
#' examining the pieces. This attempts to use that in order to evaluate how many
#' 'significant' genes are in a given pathway.
#'
#' @param all_ids Set of all gene IDs in a given analysis.
#' @param sig_ids Set of significant gene IDs.
#' @param organism KEGG organism identifier.
#' @param pathways What pathways to look at?
#' @param pathdir Directory into which to copy downloaded pathway files.
#' @param verbose Talky talky?
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return Dataframe including the filenames, percentages, nodes included, and
#'  differential nodes.
#' @seealso [KEGGgraph] [KEGGREST]
#' @export
pct_all_kegg <- function(all_ids, sig_ids, organism = "dme", pathways = "all",
                         pathdir = "kegg_pathways", verbose = FALSE, ...) {
  arglist <- list(...)
  if (!file.exists(pathdir)) {
    dir.create(pathdir)
  }
  paths <- NULL
  if (class(pathways) == "character") {
    if (length(pathways) > 1) {
      paths <- glue("path:{organism}{pathways}")
    } else if (pathways == "all") {
      all_pathways <- unique(KEGGREST::keggLink("pathway", organism))
      paths <- all_pathways
      paths <- gsub("path:", "", paths)
      ## all_modules = unique(KEGGREST::keggLink("module", species))
    } else {
      paths <- pathways
    }
  } else {
    paths <- pathways
  }

  pct_nodes <- list()
  pct_edges <- list()
  path_ids <- list()
  filenames <- list()
  path_names <- list()
  diff_nodes <- list()
  path_nodes <- list()
  diff_edges <- list()
  path_edges <- list()
  last_path <- length(paths)
  for (count in seq_along(paths)) {
    path <- paths[count]
    path_name <- try(KEGGREST::keggGet(path), silent = TRUE)
    if (class(path_name) == "try-error") {
      path_ids[count] <- NA
      path_names[count] <- NA
      filenames[count] <- NA
      pct_nodes[count] <- NA
      pct_edges[count] <- NA
      path_nodes[count] <- NA
      diff_nodes[count] <- NA
      path_edges[count] <- NA
      diff_edges[count] <- NA
    } else {
      path_name <- path_name[[1]]$NAME
      path_name <- gsub("(.*) - .*", "\\1", path_name)
      path_name <- tolower(path_name)
      path_name <- gsub(" ", "_", path_name)
      path_names[count] <- path_name
      message("Extracting data for ", path, ": ", path_name, ".")
      if (isTRUE(verbose)) {
        pct_diff <- sm(pct_kegg_diff(all_ids, sig_ids, pathway = path, organism = organism,
                                     pathdir = pathdir))
      } else {
        pct_diff <- sm(pct_kegg_diff(all_ids, sig_ids, pathway = path, organism = organism,
                                     pathdir = pathdir))
      }
      path_ids[count] <- path
      filenames[count] <- pct_diff[["filename"]]
      pct_nodes[count] <- pct_diff[["percent_nodes"]]
      pct_edges[count] <- pct_diff[["percent_edges"]]
      path_nodes[count] <- pct_diff[["all_nodes"]]
      diff_nodes[count] <- pct_diff[["diff_nodes"]]
      path_edges[count] <- pct_diff[["all_edges"]]
      diff_edges[count] <- pct_diff[["diff_edges"]]
      message(count, "/", last_path, ": The path: ", path_names[count],
              " was written to ", filenames[count], " and has ", pct_nodes[count],
              "% nodesdiff ", pct_edges[count], "% edgesdiff.")
    }
  } ## End of the for() loop.
  path_data <- data.frame()
  for (c in seq_along(path_ids)) {
    a_row <- list(
        "pathway" = path_ids[[c]],
        "path_name" = path_names[[c]],
        "filename" = filenames[[c]],
        "percent_nodes" = pct_nodes[[c]],
        "percent_edges" = pct_edges[[c]],
        "path_nodes" = path_nodes[[c]],
        "diff_nodes" = diff_nodes[[c]],
        "path_edges" = path_edges[[c]],
        "diff_edges" = diff_edges[[c]])
    ## Remove c(
    a_row <- gsub(pattern = "c\\(", replacement = "", x = a_row)
    ## Remove ),
    a_row <- gsub(pattern = "\\),", replacement = "", x = a_row)
    ## Remove "
    a_row <- gsub(pattern = "\"", replacement = "", x = a_row)
    ## Remove )
    a_row <- gsub(pattern = "\\)", replacement = "", x = a_row)
    ## Remove character(0
    a_row <- gsub(pattern = "character\\(0", replacement = "", x = a_row)
    ## Remove xxx: where xxx is the species abbreviation
    a_row <- gsub(pattern = glue("{organism}:"), replacement = "", x = a_row)
    ## Remove all commas
    a_row <- gsub(pattern = "\\,", replacement = "", x = a_row)
    ## Remove \n's
    a_row <- gsub(pattern = "\\\n", replacement = "", x = a_row)
    ## Get rid of excess space
    a_row <- gsub(pattern = " +", replacement = " ", x = a_row)
    for (k in seq_along(a_row)) {
      a_row[[k]] <- toString(a_row[[k]])
    }
    a_row <- as.list(a_row)
    names(a_row) <- c("pathway", "path_name", "filename", "percent_nodes",
                      "percent_edges", "path_nodes", "diff_nodes", "path_edges", "diff_edges")
    a_row <- as.data.frame(a_row)
    path_data <- rbind(path_data, a_row)
  }
  path_data[["pathway"]] <- as.character(path_data[["pathway"]])
  path_data[["path_name"]] <- as.character(path_data[["filename"]])
  path_data[["percent_nodes"]] <- as.numeric(as.character(path_data[["percent_nodes"]]))
  path_data[["percent_edges"]] <- as.numeric(as.character(path_data[["percent_edges"]]))
  path_data[["path_nodes"]] <- as.character(path_data[["path_nodes"]])
  path_data[["diff_nodes"]] <- as.character(path_data[["diff_nodes"]])
  path_data[["path_edges"]] <- as.character(path_data[["path_edges"]])
  path_data[["diff_edges"]] <- as.character(path_data[["diff_edges"]])

  return(path_data)
}

#' Extract the percent differentially expressed genes in a given KEGG pathway.
#'
#' KEGGgraph provides some interesting functionality for mapping KEGGids and
#' examining the pieces. This attempts to use that in order to evaluate how many
#' 'significant' genes are in a given pathway.
#'
#' @param all_ids Set of all gene IDs in a given analysis.
#' @param sig_ids Set of significant gene IDs.
#' @param pathway Numeric pathway identifier.
#' @param organism KEGG organism identifier.
#' @param pathdir Directory into which to copy downloaded pathway files.
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return Percent genes/pathway deemed significant.
#' @seealso [KEGGgraph] [KEGGREST]
#' @export
pct_kegg_diff <- function(all_ids, sig_ids, pathway = "00500",
                          organism = "dme", pathdir = "kegg_pathways", ...) {
  ## warning("This function may not work, it seems to be missing some genes.")
  arglist <- list(...)
  if (!file.exists(pathdir)) {
    dir.create(pathdir)
  }
  pathway <- gsub(pattern = organism, replacement = "", x = pathway)
  pathway <- gsub(pattern = "path:", replacement = "", x = pathway)
  filename <- file.path(pathdir, glue::glue("{organism}{pathway}.xml"))
  pathwayid <- glue::glue("{organism}{pathway}")
  retrieved <- NULL
  if (file.exists(filename)) {
    message("The file already exists, loading from it.")
    retrieved <- filename
  } else {
    retrieved <- try(myretrieveKGML(pathway, organism = organism,
                                    destfile = filename))
    if (class(retrieved) == "try-error") {
      retlist <- list(
          "pathway" = pathway,
          "filename" = "unavailable",
          "percent_nodes" = NA,
          "percent_edges" = NA,
          "all_nodes" = NULL,
          "diff_nodes" = NULL,
          "all_edges" = NULL,
          "diff_edges" = NULL)
      return(retlist)
    }
  }
  parse_result <- try(KEGGgraph::parseKGML2Graph(filename, expandGenes = TRUE))
  if (class(parse_result) == "try-error") {
    if (grepl(pattern = "Document is empty", x = parse_result[[1]])) {
      message("Deleting the empty file and trying again.")
      file.remove(filename)
      retrieved <- try(myretrieveKGML(pathway, organism = organism,
                                      destfile = filename,
                                      quiet = TRUE))
      parse_result <- try(KEGGgraph::parseKGML2Graph(filename, expandGenes = TRUE))
      if (class(parse_result) == "try-error") {
        retlist <- list(
            "pathway" = pathway,
            "filename" = "unavailable",
            "percent_nodes" = NA,
            "percent_edges" = NA,
            "all_nodes" = NULL,
            "diff_nodes" = NULL,
            "all_edges" = NULL,
            "diff_edges" = NULL)
        return(retlist)
      }
    } else if (grepl(pattern = "Start tag expected", x = parse_result[[1]])) {
      message("This pathway does not have a complete specification.")
      retlist <- list(
          "pathway" = pathway,
          "filename" = "unavailable",
          "percent_nodes" = NA,
          "percent_edges" = NA,
          "all_nodes" = NULL,
          "diff_nodes" = NULL,
          "all_edges" = NULL,
          "diff_edges" = NULL)
      return(retlist)
    }
  }

  all_keggids <- KEGGgraph::translateGeneID2KEGGID(all_ids, organism = organism)
  de_keggids <- KEGGgraph::translateGeneID2KEGGID(sig_ids, organism = organism)

  possible_nodes <- KEGGgraph::nodes(parse_result)
  node_is_differential <- possible_nodes %in% de_keggids
  found_nodes <- possible_nodes[node_is_differential]
  possible_edges <- KEGGgraph::edges(parse_result)
  edge_is_differential <- possible_edges %in% de_keggids
  found_edges <- possible_edges[edge_is_differential]

  pct_node_diff <- signif(mean(node_is_differential) * 100.0, 4)
  pct_edge_diff <- signif(mean(edge_is_differential) * 100.0, 4)

  path_data <- KEGGREST::keggGet(pathwayid)
  path_name <- path_data[[1]][["NAME"]]
  message(pct_node_diff, "% nodes differentially expressed in pathway ",
          pathway, ": '", path_name, "'.")
  message(pct_edge_diff, "% edges differentially expressed in pathway ",
          pathway, ": '", path_name, "'.")
  retlist <- list(
      "pathway" = pathway,
      "filename" = filename,
      "percent_nodes" = pct_node_diff,
      "percent_edges" = pct_edge_diff,
      "all_nodes" = toString(possible_nodes),
      "diff_nodes" = toString(found_nodes),
      "all_edges" = toString(possible_edges),
      "diff_edges" = toString(found_edges)
  )
  return(retlist)
}

#' A couple functions from KEGGgraph that have broken
#'
#' Some material in KEGGREST is borken.
#'
#' @param pathway The path to query.
#' @param organism Which organism to query?
#' @param destfile File to which to download.
#' @param silent Send stdout and stderr to dev null?
#' @param hostname Host to download from (this is what is broken.)
#' @param ...  Arglist!
#' @export
myretrieveKGML <- function(pathway, organism, destfile, silent = TRUE,
                           hostname = "http://www.kegg.jp", ...) {
  kgml <- mygetKGMLurl(pathwayid = pathway, organism = organism, hostname = hostname)
  referer <- glue::glue("{hostname}/kegg-bin/show_pathway?org_name={organism}&mapno=\\
                     {pathway}&mapscale=&show_description=hide")
  cmdline <- glue::glue("wget --header={shQuote('Accept: text/html')} \\
--user-agent={shQuote('Mozilla/5.0 (X11; Linux x86_64; \\
rv:45.0) Gecko/20100101 Firefox/45.0')} \\
--referer={shQuote(referer)} {shQuote(kgml)} \\
-O {shQuote(destfile)}")
  if (isTRUE(silent)) {
    cmdline <- glue::glue("{cmdline} 2>/dev/null 1>&2")
  }
  status <- system(cmdline)
  return(invisible(kgml))
}

mygetKGMLurl <- function(pathwayid, organism = "hsa", hostname = "http://www.kegg.jp", ...) {
  arglist <- list(...)
  if (!is.null(arglist[["hostname"]])) {
    hostname <- arglist[["hostname"]]
  }
  baseurl <- glue("{hostname}/kegg-bin/download?entry=%s%s&format=kgml")
  pathwayid <- gsub("path", "", pathwayid)
  pathwayid <- gsub(":", "", pathwayid)
  pco <- grepl("^[a-z][a-z][a-z]", pathwayid)
  org.len <- length(organism)
  if (org.len == 1 & length(pathwayid) != 1) {
    organisms <- rep(organism, length(pathwayid))
    organisms[pco] <- sapply(pathwayid[pco],
                             function(x) substr(x, 1L, 3L))
  } else if (org.len == length(pathwayid)) {
    organisms <- organism
  } else {
    stop("The length of 'organism' must be either one or the length of 'pathwayid'\n")
  }
  ids <- pathwayid
  ids[pco] <- sapply(pathwayid[pco],
                     function(x) substr(x, 4L, nchar(x)))
  urls <- sprintf(baseurl, organisms, ids)
  return(urls)
}

## EOF
