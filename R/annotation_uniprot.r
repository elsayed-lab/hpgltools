
load_uniprot_annotations <- function(file) {
  ## file <-  "uniprot_3AUP000001584.txt.gz"
  read_vec <- readr::read_lines(file, progress=TRUE)
  gene_num <- 0
  ## Vectors for those elements which will only have 1 answer
  gene_ids <- vector()
  amino_acids <- vector()
  recnames <- vector()
  loci <- vector()
  orfnames <- vector()
  shortnames <- vector()
  synonyms <- vector()
  ## Lists for those with the potential for multiple hits
  uniprot_accessions <- list()
  primary_accessions <- list()
  ## Lets do all the DR data in one shot
  embl_ids <- list()
  pir_ids <- list()
  refseq_ids <- list()
  proteinmodelportals_ids <- list()
  smr_ids <- list()
  string_ids <- list()
  pax_ids <- list()
  pride_ids <- list()
  ensbact_ids <- list()
  geneid_ids <- list()
  kegg_ids <- list()
  tuberculist_ids <- list()
  eggnog_ids <- list()
  ko_ids <- list()
  oma_ids <- list()
  phylome_ids <- list()
  proteome_ids <- list()
  go_ids <- list()
  cdd_ids <- list()
  gene3d_ids <- list()
  hamap_ids <- list()
  interpro_ids <- list()
  panther_ids <- list()
  pfam_ids <- list()
  supfam_ids <- list()
  tigrfam_ids <- list()
  message(paste0("Starting to iterate over ", length(read_vec), " lines."))
  for (i in 1:length(read_vec)) {
    line <- read_vec[i]
    if ((i %% 10000) == 0) {
      message(paste0("On line ", i))
    }
    ## Start by skipping field types that we will never use.
    if (grepl(pattern="^(DT|OS|OC|OX|RN|RA|RT|RP|RX|CC|PE|KW|FT|SQ)\\s+", x=line)) {
      next
    }
    ## The master ID:
    ## Example: ID   3MGH_MYCTU              Reviewed;         203 AA.
    if (grepl(pattern="^ID\\s+", x=line)) {
      gene_num <- gene_num + 1
      ## Initialize the ith element of our various data structures.
      primary_accessions[[gene_num]] <- ""
      amino_acids[gene_num] <- ""
      recnames[gene_num] <- ""
      loci[gene_num] <- ""
      orfnames[gene_num] <- ""
      shortnames[gene_num] <- ""
      synonyms[gene_num] <- ""
      uniprot_accessions[[gene_num]] <- ""
      ## The dr initializers go here
      embl_ids[[gene_num]] <- ""
      pir_ids[[gene_num]] <- ""
      refseq_ids[[gene_num]] <- ""
      proteinmodelportals_ids[[gene_num]] <- ""
      smr_ids[[gene_num]] <- ""
      string_ids[[gene_num]] <- ""
      pax_ids[[gene_num]] <- ""
      pride_ids[[gene_num]] <- ""
      ensbact_ids[[gene_num]] <- ""
      geneid_ids[[gene_num]] <- ""
      kegg_ids[[gene_num]] <- ""
      tuberculist_ids[[gene_num]] <- ""
      eggnog_ids[[gene_num]] <- ""
      ko_ids[[gene_num]] <- ""
      oma_ids[[gene_num]] <- ""
      phylome_ids[[gene_num]] <- ""
      proteome_ids[[gene_num]] <- ""
      go_ids[[gene_num]] <- ""
      cdd_ids[[gene_num]] <- ""
      gene3d_ids[[gene_num]] <- ""
      hamap_ids[[gene_num]] <- ""
      interpro_ids[[gene_num]] <- ""
      panther_ids[[gene_num]] <- ""
      pfam_ids[[gene_num]] <- ""
      supfam_ids[[gene_num]] <- ""
      tigrfam_ids[[gene_num]] <- ""
      ## Done initializing, now fill in the data.
      material <-  strsplit(x=line, split="\\s+")[[1]]
      gene_id <- material[2]
      gene_ids[gene_num] <- gene_id
      amino_acids[gene_num] <- material[4]
      next
    }
    ## Now pull the primary uniprot accesstions
    ## Example: AC   P9WJP7; L0TAC1; O33190; P65412;
    if (grepl(pattern="^AC\\s+", x=line)) {
      material <- gsub(pattern=";", replacement="", x=strsplit(x=line, split="\\s+")[[1]])
      primary_accessions[gene_num] <- material[2]
      possible_accessions <- material[2:length(material)]
      uniprot_accessions[[gene_num]] <- possible_accessions
      next
    }
    ## Get the record names if available
    ## Example: DE   RecName: Full=Putative 3-methyladenine DNA glycosylase {ECO:0000255|HAMAP-Rule:MF_00527};
    ##          DE            EC=3.2.2.- {ECO:0000255|HAMAP-Rule:MF_00527};
    if (grepl(pattern="DE\\s+RecName:", x=line)) {
      recnames[gene_num] <- gsub(pattern="^.*Full=(.*?);.*$", replacement="\\1", x=line)
      next
    }
    ## The GN field has a few interesting pieces of information and I think makes the primary link
    ## between uniprot and the IDs available at genbank, ensembl, microbesonline, etc.
    ## We may find one or more of the above fields in the GN, so I should take into account
    ## the various possible iterations.
    ## Example: GN   Name=pgl; Synonyms=devB; OrderedLocusNames=Rv1445c;
    ##          GN   ORFNames=MTCY493.09;
    if (grepl(pattern="^GN\\s+", x=line)) {
      pat <- "^GN\\s+.*OrderedLocusNames=(.*?);.*$"
      if (grepl(pattern=pat, x=line)) {
        ## message(paste0("Got a locusname on line ", i, " for gene number ", gene_num)) ## i=565 is first interesting one.
        loci[gene_num] <- gsub(pattern=pat, replacement="\\1", x=line)
        loci[gene_num] <- gsub(pattern="^(.*?),.*", replacement="\\1", x=loci[gene_num])
        loci[gene_num] <- gsub(pattern="^(.*?) \\{.*", replacement="\\1", x=loci[gene_num])
      }
      pat <- "^GN\\s+.*ORFNames=(.*?);.*$"
      if (grepl(pattern=pat, x=line)) {
        orfnames[gene_num] <- gsub(pattern=pat, replacement="\\1", x=line)
      }
      pat <- "^GN\\s+.*Name=(.*?);.*$"
      if (grepl(pattern=pat, x=line)) {
        shortnames[gene_num] <- gsub(pattern=pat, replacement="\\1", x=line)
        shortnames[gene_num] <- gsub(pattern="^(.*?) .*", replacement="\\1", x=shortnames[gene_num])
      }
      pat <- "^GN\\s+.*Synonyms=(.*?);.*$"
      if (grepl(pattern=pat, x=line)) {
        synonyms[gene_num] <- gsub(pattern=pat, replacement="\\1", x=line)
      }
      next
    }
    ## The DR field contains mappings to many other databases
    ## Sadly, it too is quite a mess
    ## This stanza looks for EMBL IDs:
    ## Example: DR   EMBL; AL123456; CCP44204.1; -; Genomic_DNA.
    pat_prefix <- "^DR\\s+"
    pat_suffix <- "; (.*?)( \\-.*)\\.$"
    pat <- paste0(pat_prefix, "EMBL", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      embl_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "PIR", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      pir_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "RefSeq", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      if (is.null(refseq_ids[gene_num][[1]])) {
        refseq_ids[[gene_num]] <- tmp_ids
      } else {
        refseq_ids[[gene_num]] <- c(refseq_ids[[gene_num]], tmp_ids)
      }
      refseq_ids[[gene_num]] <- unique(refseq_ids[[gene_num]])
      next
    }
    pat <- paste0(pat_prefix, "ProteinModelPortal", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      proteinmodelportals_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "SMR", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      smr_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "STRING", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      string_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "PaxDB", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      pax_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "PRIDE", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      pride_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "EnsemlBacteria", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      ensbact_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "GeneID", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      geneid_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "KEGG", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      kegg_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "TubercuList", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      tuberculist_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "eggNOG", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      if (is.null(eggnog_ids[gene_num][[1]])) {
        eggnog_ids[[gene_num]] <- tmp_ids
      } else {
        eggnog_ids[[gene_num]] <- c(eggnog_ids[[gene_num]], tmp_ids)
      }
      next
    }
    pat <- paste0(pat_prefix, "KO", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      ko_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "OMA", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      oma_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "PhylomeDB", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      phylome_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "Proteomes", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      proteome_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "GO", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      if (is.null(go_ids[gene_num][[1]])) {
        go_ids[[gene_num]] <- tmp_ids
      } else {
        go_ids[[gene_num]] <- c(go_ids[[gene_num]], tmp_ids)
      }
      go_ids[[gene_num]] <- unique(go_ids[[gene_num]])
      next
    }
    pat <- paste0(pat_prefix, "CDD", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      cdd_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "Gene3D", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      gene3d_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "HAMAP", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      hamap_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "InterPro", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      if (is.null(go_ids[gene_num][[1]])) {
        interpro_ids[[gene_num]] <- tmp_ids
      } else {
        interpro_ids[[gene_num]] <- c(interpro_ids[[gene_num]], tmp_ids)
      }
      interpro_ids[[gene_num]] <- unique(interpro_ids[[gene_num]])
      next
    }
    pat <- paste0(pat_prefix, "PANTHER", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      panther_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "Pfam", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      pfam_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "SUPFAM", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      supfam_ids[[gene_num]] <- tmp_ids
      next
    }
    pat <- paste0(pat_prefix, "TIGRFAMs", pat_suffix)
    if (grepl(pattern=pat, x=line)) {
      tmp_ids <- gsub(pattern=pat, replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      tigrfam_ids[[gene_num]] <- tmp_ids
      next
    }
  } ## End of the for loop
  message("Finished parsing, creating data frame.")
  final <- data.frame()
  final <- data.frame(row.names=gene_ids)
  final[["uniprot_primary"]] <- primary_accessions
  final[["amino_acids"]] <- amino_acids
  final[["recnames"]] <- recnames
  final[["loci"]] <- loci
  final[["orfnames"]] <- orfnames
  final[["shortnames"]] <- shortnames
  final[["synonyms"]] <- synonyms
  final[["uniprot_ids"]] <- uniprot_accessions
  final[["embl_ids"]] <- embl_ids
  final[["pir_ids"]] <- pir_ids
  final[["reseq_ids"]] <- refseq_ids
  final[["proteinmodelportals_ids"]] <- proteinmodelportals_ids
  final[["smr_ids"]] <- smr_ids
  final[["string_ids"]] <- string_ids
  final[["pax_ids"]] <- pax_ids
  final[["pride_ids"]] <- pride_ids
  final[["ensbact_ids"]] <- ensbact_ids
  final[["geneid_ids"]] <- geneid_ids
  final[["kegg_ids"]] <- kegg_ids
  final[["tuberculist_ids"]] <- tuberculist_ids
  final[["eggnog_ids"]] <- eggnog_ids
  final[["ko_ids"]] <- ko_ids
  final[["oma_ids"]] <- oma_ids
  final[["phylome_ids"]] <- phylome_ids
  final[["proteome_ids"]] <- proteome_ids
  final[["go_ids"]] <- go_ids
  final[["cdd_ids"]] <- cdd_ids
  final[["gene3d_ids"]] <- gene3d_ids
  final[["hamap_ids"]] <- hamap_ids
  final[["interpro_ids"]] <- interpro_ids
  final[["panther_ids"]] <- panther_ids
  final[["pfam_ids"]] <- pfam_ids
  final[["supfam_ids"]] <- supfam_ids
  final[["tigrfam_ids"]] <- tigrfam_ids
  return(final)
}
