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
#' @examples
#'  uniprot_sc_downloaded <- download_uniprot_proteome(species="Saccharomyces cerevisiae S288c")
#'  uniprot_sc_downloaded$filename
#'  uniprot_sc_downloaded$species
#' @export
download_uniprot_proteome <- function(accession=NULL, species=NULL,
                                      taxonomy=NULL, all=FALSE, first=FALSE) {
  final_species <- ""
  if (!is.null(taxonomy)) {
    request_url <- glue::glue("https://www.uniprot.org/proteomes/?query=taxonomy%3A{xml2::url_escape(taxonomy)}")
    destination <- glue("{taxonomy}.txt.gz")
    ##request <- curl::curl(request_url)
    tt <- download.file(url=request_url, destfile=destination, method="wget", quiet=TRUE)
    result <- xml2::read_html(destination)
    result_html <- rvest::html_nodes(result, "tr")
    accessions_text <- rvest::html_attr(result_html, "id")
    ## The first two elements are headers
    accessions_text <- accessions_text[3:length(accessions_text)]
    accessions <- gsub(x=accessions_text, pattern="^(UP[0-9]+)(.*$)", replacement="\\1")
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
      destination <- glue("{tempfile()}.txt.gz")
      request_url <- glue("https://www.uniprot.org/proteomes/?query={xml2::url_escape(species)}")
      ##request <- curl::curl(request_url)
      tt <- download.file(url=request_url, destfile=destination, method="wget", quiet=TRUE)
      result <- xml2::read_html(destination)
      result_html <- rvest::html_nodes(result, "tr")
      accessions_text <- rvest::html_attr(result_html, "id")
      ## The first two elements are headers
      accessions_text <- accessions_text[3:length(accessions_text)]
      accessions <- gsub(x=accessions_text, pattern="^(UP[0-9]+)(.*$)", replacement="\\1")
      species_text <- rvest::html_nodes(result, "td") %>%
        rvest::html_nodes("span") %>%
        rvest::html_text()
      final_species <- species_text[species_text != ""]
      removed <- file.remove(destination)
      if (length(accessions) == 1) {
        accession <- accessions
      } else if (isTRUE(all)) {
        for (a in 1:length(accessions)) {
          name <- final_species[a]
          accession <- accessions[a]
          message("Downloading the proteome for ", name, ".")
          tmp <- download_uniprot_proteome(accession=accession)
          Sys.sleep(time=3)
        }
      } else if (isTRUE(first)) {
        accession <- accessions[1]
        name <- final_species[1]
        message("Downloading the proteome for ", name, ".")
        tmp <- download_uniprot_proteome(accession=accession)
      } else {
        message("Here are the species found, please choose one and try again.")
        for (a in 1:length(accessions)) {
          name <- final_species[a]
          accession <- accessions[a]
          message(a, ") ", accession, ": ", name)
        }
        message(toString(final_species))
        return(NULL)
      }
    }
  }
  request_url <- glue(
    "https://www.uniprot.org/uniprot/?query=proteome:\\
     {accession}&compress=yes&force=true&format=txt")
  destination <- glue("{accession}.txt.gz")
  ##tt <- curl::curl_fetch_disk(url=request_url, path=destination)
  tt <- download.file(url=request_url, destfile=destination, method="wget", quiet=TRUE)
  retlist <- list(
    "filename" = destination,
    "species" = final_species,
    "accession" = accession)
  return(retlist)
}

#' Read a uniprot text file and extract as much information from it as possible.
#'
#' I spent entirely too long fighting with Uniprot.ws, finally got mad and wrote this.
#'
#' @param file Uniprot file to read and parse
#' @param species Species name to download/load.
#' @param savefile Do a save?
#' @return Big dataframe of annotation data.
#' @examples
#'  sc_uniprot_annot <- load_uniprot_annotations(file=uniprot_sc_downloaded$filename)
#'  dim(sc_uniprot_annot)
#' @export
load_uniprot_annotations <- function(file=NULL, species=NULL, savefile=TRUE) {
  if (is.null(file) & is.null(species)) {
    stop("This requires either a filename or species name.")
  } else if (is.null(file)) {
    info <- download_uniprot_proteome(species=species)
    message("Downloaded proteome for: ", info[["species"]], " accession: ",
            info[["accession"]], " to file: ", info[["filename"]], ".")
    file <- info[["filename"]]
  }

  if (isTRUE(savefile)) {
    savefile <- "uniprot.rda"
  }
  ##if (!is.null(savefile) & savefile != NULL) {
  ##  if (file.exists(savefile)) {
  ##    uniprot_data <- new.env()
  ##    loaded <- load(savefile, envir=retlist)
  ##    uniprot_data <- uniprot_data[["uniprot_data"]]
  ##    return(retlist)
  ##  }
  ##}
  read_vec <- readr::read_lines(file)
  gene_num <- 0
  num_genes <- length(grep(pattern="^ID", x=read_vec))
  ## Vectors for those elements which will only have 1 answer
  many_ids <- list()
  id_types <- c(
    "primary_id", "amino_acids",
    "primary_accession", "uniprot_accessions",
    "recnames", "loci", "orfnames", "shortnames", "synonyms",
    ## The very many DR IDs
    "embl", "ccds", "pir", "refseq", "unigene", "proteinmodelportal", "smr",
    "intact", "string", "iptmnet", "phosphosite", "biomuta", "dmdm",
    "epd", "pax", "peptideatlas", "pride", "proteomicsdb", "dnasu",
    "ensembl", "ensbact", "geneid", "kegg", "tuberculist", "ucsc",
    "ctd", "eupathdb", "genecards", "hgnc", "hpa", "mim", "nextprot",
    "opentargets", "pharmgkb", "eggnog", "ko", "genetreee", "hogenom",
    "hovergen", "inparanoid", "oma", "orthodb", "phylomedb", "treefam",
    "genewiki", "genomernai", "pro", "proteomes", "bgee", "cleanex",
    "unipathway", "expressionatlas", "genevisible", "go", "cdd",
    "gene3d", "hamap", "interpro", "panther", "pfam", "pirsf", "prints",
    "supfam", "tigrfam",
    ## Final stuff
    "mw", "aa_length", "aa_sequence")
  uniprot_data <- data.frame(row.names=1:num_genes)
  for (id in id_types) {
    uniprot_data[[id]] <- ""
  }

  reading_sequence <- FALSE
  show_progress <- interactive() && is.null(getOption("knitr.in.progress"))
  if (isTRUE(show_progress)) {
    bar <- utils::txtProgressBar(style=3)
  }
  for (i in 1:length(read_vec)) {
    if (isTRUE(show_progress)) {
      pct_done <- i / length(read_vec)
      utils::setTxtProgressBar(bar, pct_done)
    }

    ## Start by skipping field types that we will never use.
    ## DT: History of the entry
    ## OS: Species, we kind of already know that.
    ## OC: Taxonomy of the species.
    ## RN: Looks like an arbitrary number
    ## RA: Authors
    ## RT: vector?
    ## RP: Data source?
    ## RX: Pubmed ID and DOI, this might be useful.
    ## CC: Reaction information
    ## PE: Evidence information
    ## KW: More reaction information
    ## FT: PDB information
    line <- read_vec[i]
    typestring <- substr(line, 1, 2)
    switchret <- switch(
      typestring,
      "DT|OS|OC|RN|RA|RT|RP|RX|CC|PE|KW|FT|OX" = {
        next
      },
      "ID" = {
        ## The master ID:
        ## Example: ID   3MGH_MYCTU              Reviewed;         203 AA.
        gene_num <- gene_num + 1
        material <- strsplit(x=line, split="\\s+")[[1]]
        gene_id <- material[2]
        uniprot_data[gene_num, "primary_id"] <- gene_id
        uniprot_data[gene_num, "amino_acids"] <- material[4]
        next
      },
      "AC" = {
        ## Now pull the primary uniprot accesstions
        ## Example: AC   P9WJP7; L0TAC1; O33190; P65412;
        tmp_ids <- gsub(pattern=";", replacement="", x=strsplit(x=line, split="\\s+")[[1]])
        uniprot_data[gene_num, "primary_accession"] <- tmp_ids[2]
        tmp_ids <- toString(tmp_ids[2:length(tmp_ids)])
        uniprot_data[gene_num, "uniprot_accessions"] <- tmp_ids
        next
      },
      "DE" = {
        ## Get the record names if available
        ## Example:
        ## DE RecName: Full=Putative 3-methyl DNA glycosylase {ECO:0000255|HAMAP-Rule:MF_00527};
        ##          DE            EC=3.2.2.- {ECO:0000255|HAMAP-Rule:MF_00527};
        if (grepl(pattern="DE\\s+RecName:", x=line)) {
          tmp_ids <- gsub(pattern="^.*Full=(.*?);.*$", replacement="\\1", x=line)
          uniprot_data[gene_num, "recnames"] <- tmp_ids
          next
        } else {
          next
        }
      },
      "GN" = {
        ## The GN field has a few interesting pieces of information and I think
        ## makes the primary link between uniprot and the IDs available at genbank,
        ## ensembl, microbesonline, etc. We may find one or more of the above fields
        ## in the GN, so I should take into account the various possible iterations.
        ## Example: GN   Name=pgl; Synonyms=devB; OrderedLocusNames=Rv1445c;
        ##          GN   ORFNames=MTCY493.09;
        if (grepl(pattern="^GN\\s+.*OrderedLocusNames=(.*?);.*$", x=line)) {
          ## message("Got a locusname on line ", i, " for gene number ", gene_num)
          ## i=565 is first interesting one.
          tmp_ids <- gsub(pattern="^GN\\s+.*OrderedLocusNames=(.*?);.*$", replacement="\\1", x=line)
          tmp_ids <- gsub(pattern="^(.*?),.*", replacement="\\1", x=tmp_ids)
          uniprot_data[gene_num, "loci"] <- gsub(
            pattern="^(.*?) \\{.*", replacement="\\1", x=tmp_ids)
          next
        } else if (grepl(pattern="^GN\\s+.*ORFNames=(.*?);.*$", x=line)) {
          uniprot_data[gene_num, "orfnames"] <- gsub(
            pattern="^GN\\s+.*ORFNames=(.*?);.*$", replacement="\\1", x=line)
          next
        } else if (grepl(pattern="^GN\\s+.*Name=(.*?);.*$", x=line)) {
          tmp_ids <- gsub(
            pattern="^GN\\s+.*Name=(.*?);.*$", replacement="\\1", x=line)
          uniprot_data[gene_num, "shortnames"] <- gsub(
            pattern="^(.*?) .*", replacement="\\1", x=tmp_ids)
          next
        } else if (grepl(pattern="^GN\\s+.*Synonyms=(.*?);.*$", x=line)) {
          uniprot_data[gene_num, "synonyms"] <- gsub(
            pattern="^GN\\s+.*OrderedLocusNames=(.*?);.*$", replacement="\\1", x=line)
          next
        } else {
          next
        }
      },
      "DR" = {
        ## The DR field contains mappings to many other databases
        ## Sadly, it too is quite a mess
        ## This stanza looks for EMBL IDs:
        ## Example: DR   EMBL; AL123456; CCP44204.1; -; Genomic_DNA.
        rest <- substr(x=line, start=6, stop=nchar(line))
        matches <- stringr::str_match(rest, "^(\\w+);\\s+(.*?)(\\.$|; \\-.$|\\. \\[.*\\]$)")
        intype <- matches[1, 2]
        information <- gsub(pattern=";", replacement="\\,", x=matches[1, 3])
        inswitchret <- switch(
          intype,
          "EMBL" = {
            uniprot_data[gene_num, "embl"] <- information
          },
          "CCDS" = {
            ## Consensus CDS protein set: https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi
            uniprot_data[gene_num, "ccds"] <- information
          },
          "PIR" = {
            ## The protein information resource: https://pir.georgetown.edu/
            uniprot_data[gene_num, "pir"] <- information
          },
          "RefSeq" = {
            ## RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
            if (uniprot_data[gene_num, "refseq"] == "") {
              uniprot_data[gene_num, "refseq"] <- information
            } else {
              uniprot_data[gene_num, "refseq"] <- toString(
                c(uniprot_data[gene_num, "refseq"], information))
            }
          },
          "UniGene" = {
            ## UniGene: https://www.ncbi.nlm.nih.gov/unigene
            uniprot_data[gene_num, "unigene"] <- information
          },
          "ProteinModelPortal" = {
            uniprot_data[gene_num, "proteinmodelportal"] <- information
          },
          "SMR" = {
            ## Small Multidrug Resistance proteins: This is actually a boolean identifying SMR proteins.
            uniprot_data[gene_num, "smr"] <- information
          },
          "IntAct" = {
            ## IntAct: Molecular Interaction Database: https://www.ebi.ac.uk/intact/
            uniprot_data[gene_num, "intact"] <- information
          },
          "STRING" = {
            ## STRING: The protein-protein interaction network database: https://string-db.org/
            uniprot_data[gene_num, "string"] <- information
          },
          "iPTMnet" = {
            ## iPTMnet: Protein post-translational modification database: https://research.bioinformatics.udel.edu/iptmnet/
            uniprot_data[gene_num, "iptmnet"] <- information
          },
          "PhosphoSitePlus" = {
            ## PhosphoSitePlus: Phosphorylation site database: https://www.phosphosite.org/homeAction.action
            uniprot_data[gene_num, "phosphosite"] <- information
          },
          "BioMuta" = {
            ## BioMuta: Single Nucleotide Variants in Cancer: https://hive.biochemistry.gwu.edu/biomuta
            uniprot_data[gene_num, "biomuta"] <- information
          },
          "DMDM" = {
            ## DMDM: Domain Mapping of Disease Mutations: http://bioinf.umbc.edu/dmdm/
            uniprot_data[gene_num, "dmdm"] <- information
          },
          "EPD" = {
            ## EPD: The Eukaryotic Promoter Database: https://epd.epfl.ch/index.php
            uniprot_data[gene_num, "epd"] <- information
          },
          "PaxDB" = {
            ## PaxDB: Protein Abundance Database: https://pax-db.org/
            uniprot_data[gene_num, "pax"] <- information
          },
          "PeptideAtlas" = {
            ## PeptideAtlas: Compendium of peptides in tandem mass spec datasets: http://www.peptideatlas.org/
            uniprot_data[gene_num, "peptideatlas"] <- information
          },
          "PRIDE" = {
            ## PRIDE: PRoteomics IDentifications database: https://www.ebi.ac.uk/pride/archive/
            uniprot_data[gene_num, "pride"] <- information
          },
          "ProteomicsDB" = {
            ## ProteomicsDB: https://www.proteomicsdb.org/
            uniprot_data[gene_num, "proteomicsdb"] <- information
          },
          "DNASU" = {
            ## DNASU: Plasmid Repository: https://dnasu.org/DNASU/Home.do
            uniprot_data[gene_num, "dnasu"] <- information
          },
          "Ensembl" = {
            ## Ensembl: https://useast.ensembl.org/index.html
            uniprot_data[gene_num, "ensembl"] <- information
          },
          "EnsemblBacteria" = {
            ## https://bacteria.ensembl.org/index.html
            uniprot_data[gene_num, "ensbact"] <- information
          },
          "GeneID" = {
            ## GeneID: http://genome.crg.es/software/geneid/
            uniprot_data[gene_num, "geneid"] <- information
          },
          "KEGG" = {
            ## Kyoto Encyclopedia of Genes and Genomes: https://www.genome.jp/kegg/
            uniprot_data[gene_num, "kegg"] <- information
          },
          "TubercuList" = {
            ## Tuberculist! http://genolist.pasteur.fr/TubercuList/
            uniprot_data[gene_num, "tuberculist"] <- information
          },
          "UCSC" = {
            uniprot_data[gene_num, "ucsc"] <- information
          },
          "CTD" = {
            uniprot_data[gene_num, "ctd"] <- information
          },
          "EuPathDB" = {
            ## EuPathDB: Eukaryotic Pathogen Database: https://eupathdb.org/eupathdb/
            uniprot_data[gene_num, "eupathdb"] <- information
          },
          "GeneCards" = {
            ## GeneCards: The Human Gene Database: https://www.genecards.org/
            uniprot_data[gene_num, "genecards"] <- information
          },
          "HGNC" = {
            ## HGNC: The HUGO Gene Nomenclature Commit: https://www.genenames.org/
            uniprot_data[gene_num, "hgnc"] <- information
          },
          "HPA" = {
            ## HPA: The Human Protein Atlas: https://www.proteinatlas.org/
            uniprot_data[gene_num, "hpa"] <- information
          },
          "MIM" = {
            ## MIM: Online Mendelian Inheritance in Man: https://www.omim.org/
            uniprot_data[gene_num, "mim"] <- information
          },
          "neXtProt" = {
            ## neXtProt: The human protein database: https://www.nextprot.org/
            uniprot_data[gene_num, "nextprot"] <- information
          },
          "OpenTargets" = {
            ## OpenTargets: Evaluate validity of therapeutic targets: https://www.opentargets.org/
            uniprot_data[gene_num, "opentargets"] <- information
          },
          "PharmGKB" = {
            ## PharmGKB: The Pharmacogenomics Knowledgebase: https://www.pharmgkb.org/
            uniprot_data[gene_num, "pharmgkb"] <- information
          },
          "eggNOG" = {
            ## eggNOG: Orthology predictions and function annotation: http://eggnogdb.embl.de/#/app/home
            if (uniprot_data[gene_num, "eggnog"] == "") {
              uniprot_data[gene_num, "eggnog"] <- information
            } else {
              uniprot_data[gene_num, "eggnog"] <- toString(
                c(uniprot_data[gene_num, "eggnog"], information))
            }
          },
          "KO" = {
            ## KO: KEGG Orthology database: https://www.genome.jp/kegg/ko.html
            uniprot_data[gene_num, "ko"] <- information
          },
          "GeneTree" = {
            ## GeneTree: Ensembl Genomes
            uniprot_data[gene_num, "genetree"] <- information
          },
          "HOGENOM" = {
            ## HOGENOM: Database of Complete Genome Homologous Gene Families:
            ## http://doua.prabi.fr/databases/hogenom/home.php?contents=query
            uniprot_data[gene_num, "hogenom"] <- information
          },
          "HOVERGEN" = {
            ## HOVERGEN: Homologous Vertebrate Genes Database:
            ## http://pbil.univ-lyon1.fr/databases/hovergen.php
            uniprot_data[gene_num, "hovergen"] <- information
          },
          "InParanoid" = {
            ## InParanoid: Ortholog groups with inparalogs: http://inparanoid.sbc.su.se/cgi-bin/index.cgi
            uniprot_data[gene_num, "inparanoid"] <- information
          },
          "OMA" = {
            ## OMA: Ortholog browser: https://omabrowser.org/oma/home/
            uniprot_data[gene_num, "oma"] <- information
          },
          "OrthoDB" = {
            ## OrthoDB: The heirarchical catalog of orthologs
            uniprot_data[gene_num, "orthodb"] <- information
          },
          "PhylomeDB" = {
            ## PhylomeDB: Repository of large scale phylogenetic information: http://phylomedb.org/
            uniprot_data[gene_num, "phylomedb"] <- information
          },
          "TreeFam" = {
            ## TreeFam: Phylogenetic Tree Database: http://www.treefam.org/
            uniprot_data[gene_num, "treefam"] <- information
          },
          "GeneWiki" = {
            ## GeneWiki: Wikipedia Gene Database: https://en.wikipedia.org/wiki/Gene_Wiki
            uniprot_data[gene_num, "genewiki"] <- information
          },
          "GenomeRNAi" = {
            ## GenomeRNAi: RNAi phenotypes and reagents: http://www.genomernai.org/
            uniprot_data[gene_num, "genomernai"] <- information
          },
          "PRO" = {
            ## PRO: Proteomics DB? https://www.proteomicsdb.org/
            uniprot_data[gene_num, "pro"] <- information
          },
          "Proteomes" = {
            ## Proteomes: uniprot proteomes database
            uniprot_data[gene_num, "proteomes"] <- information
          },
          "Bgee" = {
            ## Bgee: Gene Expression Data in Animals
            uniprot_data[gene_num, "bgee"] <- information
          },
          "CleanEx" = {
            ## CleanEx: Database of gene expression profiles: https://cleanex.epfl.ch//
            uniprot_data[gene_num, "cleanex"] <- information
          },
          "UniPathway" = {
            ## Uniprot pathways
            uniprot_data[gene_num, "unipathway"] <- information
          },
          "ExpressionAtlas" = {
            ## ExpressionAtlas: Gene expression results across species: https://www.ebi.ac.uk/gxa/home
            uniprot_data[gene_num, "expressionatlas"] <- information
          },
          "Genevisible" = {
            ## Genevisible: https://genevisible.com/search
            uniprot_data[gene_num, "genevisible"] <- information
          },
          "GO" = {
            ## Gene Ontology
            if (uniprot_data[gene_num, "go"] == "") {
              uniprot_data[gene_num, "go"] <- information
            } else {
              uniprot_data[gene_num, "go"] <- toString(
                c(uniprot_data[gene_num, "go"], information))
            }
          },
          "CDD" = {
            uniprot_data[gene_num, "cdd"] <- information
          },
          "Gene3D" = {
            uniprot_data[gene_num, "gene3d"] <- information
          },
          "HAMAP" = {
            uniprot_data[gene_num, "hamap"] <- information
          },
          "InterPro" = {
            if (uniprot_data[gene_num, "interpro"] == "") {
              uniprot_data[gene_num, "interpro"] <- information
            } else {
              uniprot_data[gene_num, "interpro"] <- toString(
                c(uniprot_data[gene_num, "interpro"], information))
            }
          },
          "PANTHER" = {
            uniprot_data[gene_num, "panther"] <- information
          },
          "Pfam" = {
            uniprot_data[gene_num, "pfam"] <- information
          },
          "PIRSF" = {
            uniprot_data[gene_num, "pirsf"] <- information
          },
          "PRINTS" = {
            uniprot_data[gene_num, "prints"] <- information
          },
          "SUPFAM" = {
            uniprot_data[gene_num, "supfam"] <- information
          },
          "TIGRFAMs" = {
            uniprot_data[gene_num, "tigrfam"] <- information
          })
      },
      "SQ" = {
        aa_seq <- ""
        reading_sequence <- TRUE
        mweight <- gsub(
          pattern="^SQ\\s+SEQUENCE\\s+\\d+\\s+AA;\\s+(\\d+)\\s+MW.*$", replacement="\\1", x=line)
        uniprot_data[gene_num, "mw"] <- mweight
        aa_length <- gsub(
          pattern="^SQ\\s+SEQUENCE\\s+(\\d+)\\s+AA;\\s+\\d+\\s+MW.*$", replacement="\\1", x=line)
        uniprot_data[gene_num, "aa_length"] <- aa_length
        if (isTRUE(reading_sequence)) {
          if (grepl(pattern="^\\s+", x=line)) {
            aa_line <- gsub(pattern="\\s", replacement="", x=line)
            aa_seq <- paste0(aa_seq, aa_line)
          }
        }
        if (grepl(pattern="^\\/\\/", x=line)) {
          uniprot_data[gene_num, "aa_sequence"] <- aa_seq
          reading_sequence <- FALSE
        }
      })
  } ## End of the for loop
  if (isTRUE(show_progress)) {
    close(bar)
  }
  message("Finished parsing, creating data frame.")

  if (!is.null(savefile)) {
    if (savefile != FALSE) {
      saved <- save(list="uniprot_data", file=savefile)
    }
  }
  return(uniprot_data)
}

#' Extract ontology information from a uniprot dataframe.
#'
#' @param input uniprot filename or dataframe.
#' @return Ontology dataframe
#' @examples
#'  sc_uniprot_go <- load_uniprot_go(sc_uniprot_annot)
#'  head(sc_uniprot_go)
#' @export
load_uniprot_go <- function(input) {
  if ("character" %in% class(input)) {
    input <- load_uniprot_annotations(file=input)
  } else if ("data.frame" %in% class(input)) {
    input <- as.data.frame(input)
  }

  kept <- input[, c("uniprot_accessions", "go", "aa_length")] %>%
    tidyr::separate_rows("uniprot_accessions")

  kept[["go"]] <- kept[["go"]] %>%
    stringr::str_extract_all(pattern="GO:\\d+")
  kept[["go"]] <- I(kept[["go"]])
  kept[["go"]] <- as.character(kept[["go"]])
  kept <- kept %>%
    tidyr::separate_rows("go", sep=",")
  kept[["go"]] <- gsub(pattern='"', replacement="", x=kept[["go"]])
  kept[["go"]] <- gsub(pattern=")", replacement="", x=kept[["go"]])
  kept[["go"]] <- gsub(pattern="c\\(", replacement="", x=kept[["go"]])
  kept[["go"]] <- gsub(pattern="\\s+", replacement="", x=kept[["go"]])
  colnames(kept) <- c("ID", "GO", "length")
  return(kept)
}

## EOF
