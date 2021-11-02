get_annots <- function(annots = "TMRC2") {
  if (annots == "TMRC2"){
  pan_db <- org.Lpanamensis.MHOMCOL81L13.v46.eg.db
  all_fields <- columns(pan_db)

  all_lp_annot <- sm(load_orgdb_annotations(
    pan_db,
    keytype = "gid",
    fields = c("annot_gene_entrez_id", "annot_gene_name",
               "annot_strand", "annot_chromosome", "annot_cds_length",
               "annot_gene_product")))$genes

  lp_go <- sm(load_orgdb_go(pan_db))
  lp_lengths <- all_lp_annot[, c("gid", "annot_cds_length")]
  colnames(lp_lengths)  <- c("ID", "length")
  all_lp_annot[["annot_gene_product"]] <- tolower(all_lp_annot[["annot_gene_product"]])
  orthos <- sm(EuPathDB::extract_eupath_orthologs(db = pan_db))

  hisat_annot <- all_lp_annot
  return(hisat_annot)
  } else if (annots == "TMRC3") {
    hs_annot <- sm(load_biomart_annotations(year="2020"))
    hs_annot <- hs_annot[["annotation"]]
    hs_annot[["transcript"]] <- paste0(rownames(hs_annot), ".", hs_annot[["version"]])
    rownames(hs_annot) <- make.names(hs_annot[["ensembl_gene_id"]], unique=TRUE)
    tx_gene_map <- hs_annot[, c("transcript", "ensembl_gene_id")]

    return(hs_annot)
  }

}


get_data <- function(data = "TMRC2") {
  if (data == "TMRC2"){
    load("/mnt/cbcb/fs00_reesyxan/cbcb-lab/nelsayed/scratch/atb/rnaseq/lpanamensis_tmrc_2019/rda/expt.rda")
    lp_expt <- expt %>%
      set_expt_conditions(fact = "zymodemecategorical") %>%
      subset_expt(nonzero = 8550) %>%
      subset_expt(coverage = 5000000) %>%
      semantic_expt_filter(semantic = c("amastin", "gp63", "leishmanolysin"),
                           semantic_column = "annot_gene_product")

    return(lp_expt)
  } else if (data == "TMRC3"){
    load("/mnt/cbcb/fs00_reesyxan/cbcb-lab/nelsayed/scratch/atb/rnaseq/lpanamensis_tmrc_2019/rda/hs_expt_all-v202110.rda")
    hs_expt <- expt %>%
      exclude_genes_expt(column="gene_biotype", method="keep",
                         pattern="protein_coding", meta_column="ncrna_lost")

    return(hs_expt)

  }
}


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
