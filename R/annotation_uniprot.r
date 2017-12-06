
parse_uniprot_txt_table <- function(file) {
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
  embl_ids <- list()
  refseq_ids <- list()
  pir_ids <- list()
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
      amino_acids[gene_num] <- ""
      recnames[gene_num] <- ""
      loci[gene_num] <- ""
      orfnames[gene_num] <- ""
      shortnames[gene_num] <- ""
      synonyms[gene_num] <- ""
      uniprot_accessions[[gene_num]] <- ""
      primary_accessions[[gene_num]] <- ""
      embl_ids[[gene_num]] <- ""
      refseq_ids[[gene_num]] <- ""
      pir_ids[[gene_num]] <- ""
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
    if (grepl(pattern="^DR\\s+EMBL;", x=line)) {
      tmp_ids <- gsub(pattern="^.*EMBL; (.*?) \\-;.*$", replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      embl_ids[[gene_num]] <- tmp_ids
      next
    }
    if (grepl(pattern="^DR\\s+RefSeq;", x=line)) {
      tmp_ids <- gsub(pattern="^.*RefSeq; (.*?)$", replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";|\\.$", replacement="", x=tmp_ids)
      if (is.null(refseq_ids[gene_num][[1]])) {
        refseq_ids[[gene_num]] <- tmp_ids
      } else {
        refseq_ids[[gene_num]] <- c(refseq_ids[[gene_num]], tmp_ids)
      }
      refseq_ids[[gene_num]] <- unique(refseq_ids[[gene_num]])
      next
    }
    if (grepl(pattern="^DR\\s+PIR;", x=line)) {
      tmp_ids <- gsub(pattern="^.*PIR; (.*?)\\.$", replacement="\\1", x=line)
      tmp_ids <- strsplit(x=tmp_ids, split="\\s+")[[1]]
      tmp_ids <- gsub(pattern=";", replacement="", x=tmp_ids)
      tmp_ids <- unique(tmp_ids)
      pir_ids[[gene_num]] <- tmp_ids
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
  final[["uniprot_acc"]] <- uniprot_accessions
  final[["embl_acc"]] <- embl_ids
  final[["refseq_acc"]] <- refseq_ids
  final[["pir_acc"]] <- pir_ids
  return(final)
}
