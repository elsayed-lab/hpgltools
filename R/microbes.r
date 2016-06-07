#' Use the publicly available microbesonline mysql instance to get species ids.
#'
#' The microbesonline mysql instance is more complex than I like.  Their id system is reminiscent of
#' KEGG's and similarly annoying.  Though I haven't figured out how the tables interact, a query to
#' get ids is simple enough.
#'
#' @param name Text string containing some part of the species name of interest.
#' @param exact Use an exact species name?
#' @return Dataframe of ids and names.
#' @export
get_microbesonline_ids <- function(name, exact=FALSE) {
    requireNamespace("RMySQL")
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
    query <- "SELECT taxonomyId, shortName FROM Taxonomy WHERE shortName "
    if (isTRUE(exact)) {
        query <- paste0(query, "= '", name, "'")
    } else {
        query <- paste0(query, "like '%", name, "%'")
    }
    message(query)
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

#' Use the publicly available microbesonline mysql instance to get species name(s).
#'
#' The microbesonline mysql instance is more complex than I like.  Their id system is reminiscent of
#' KEGG's and similarly annoying.  Though I haven't figured out how the tables interact, a query to
#' get ids is simple enough.
#'
#' @param id Text string containing some part of the species name of interest.
#' @return Dataframe of ids and names.
#' @export
get_microbesonline_name <- function(id) {
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
    query <- paste0("SELECT shortName FROM Taxonomy WHERE taxonomyId = '", id, "'")
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    disconnect <- DBI::dbDisconnect(connection)
    print(result_df)
    return(result_df)
}

#' Skip the db and download all the text annotations for a given species.
#'
#' Like I said, the microbesonline mysqldb is rather more complex than I prefer.  This shortcuts
#' that process and just grabs a tsv copy of everything and loads it into a dataframe.
#'
#' @param ids List of ids to query.
#' @param species Species name(s) to use instead.
#' @return List of dataframes with the annotation information.
#' @export
get_microbesonline_annotation <- function(ids="160490", species=NULL) {
    retlist <- list()
    id_list <- list()
    if (is.null(ids) & is.null(species)) {
        stop("Either ids or species must be defined.")
    } else if (is.null(ids)) {
        for (spec in species) {
            id_df <- get_microbesonline_ids(spec)
            id_shortlist <- as.list(id_df[[1]])
            names(id_shortlist) <- id_df[[2]]
            id_list <- append(x=id_list, values=id_df[[1]])
        }
    } else {
        for (id in ids) {
            idl <- as.list(id)
            current_name <- get_microbesonline_name(id)
            names(idl) <- current_name[1,1]
            id_list <- append(x=id_list, values=idl)
        }
    }

    retlist <- list()
    print(id_list)
    for (t in 1:length(id_list)) {
        name <- names(id_list)[[t]]
        id <- id_list[[t]]
        url <- paste0("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", id, ";export=tab")
        string <- RCurl::getURL(url)
        con <- textConnection(string)
        data <- read.csv(con, sep="\t", header=TRUE)
        retlist[[name]] <- data
    }
    return(retlist)
}

mdesc_table <- function(table="Locus2GO") {
    db_driver <- DBI::dbDriver("MySQL")
    connection <- DBI::dbConnect(db_driver, user="guest", password="guest", host="pub.microbesonline.org", dbname="genomics")
    query <- paste0("DESCRIBE TABLE ", table)
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

get_loci_go <- function(taxonid="160490") {
    db_driver <- DBI::dbDriver("MySQL")
    ## cheese and crackers their database is entirely too complex and poorly documented.
    query <- paste0("SELECT L.locusId, G.goID, T.acc_synonym FROM genomics.Scaffold S, genomics.term_synonym T, genomics.Locus L, genomics.Locus2Go G where S.TaxonomyId = '",
                    taxonid, "' and S.isGenomic = 1 and S.scaffoldId = L.scaffoldId  and G.locusId = L.locusId and T.term_id = G.goID")
    result <- DBI::dbSendQuery(connection, query)
    result_df <- DBI::fetch(result, n=-1)
    disconnect <- DBI::dbDisconnect(connection)
    return(result_df)
}

## To get the set of GOids, I need to select goID, evidence from Locus2Go where locusID = stuff

##                Tables_in_genomics
## 1                           AASeq
## 2                             ACL
## 3                      Annotation
## 4                AnnotationDetail
## 5                        AutoAnno
## 6                             COG
## 7                        COGCount
## 8                          COGFun
## 9                         COGInfo
## 10                    COGrpsblast
## 11                          Carts
## 12                   CrisprFamily
## 13                    Description
## 14                     DomainInfo
## 15                         ECInfo
## 16              FasthmmFamily2HMM
## 17                 FasthmmRawHits
## 18                        GoCount
## 19                     GroupUsers
## 20                         Groups
## 21                        HitInfo
## 22                         IPR2Go
## 23                        IPRInfo
## 24                       InterPro
## 25                           Jobs
## 26                  KEGG2Taxonomy
## 27                   KEGGCompound
## 28                       KEGGConf
## 29                       KEGGInfo
## 30                        KEGGMap
## 31                          Locus
## 32                   Locus2Domain
## 33                       Locus2Ec
## 34                       Locus2Go
## 35                      Locus2Ipr
## 36                   Locus2Operon
## 37                      Locus2Pdb
## 38                      Locus2RTB
## 39              Locus2RTBArticles
## 40               Locus2RegPrecise
## 41                     Locus2Seed
## 42                   Locus2String
## 43                Locus2SwissProt
## 44               Locus2TigrFunCat
## 45                     Locus2Tree
## 46                     LocusCount
## 47                       LocusSeq
## 48                      LocusType
## 49                            MOG
## 50                   MOGComponent
## 51                      MOGMember
## 52              MOGNeighborScores
## 53                       OGMember
## 54                    ObjectDescr
## 55                         Operon
## 56                  OrthologGroup
## 57                     PdbEntries
## 58                        PdbReps
## 59                         PdbSeq
## 60                       PfamClan
## 61                       Position
## 62                 RegulonCluster
## 63                   RegulonLinks
## 64                       Scaffold
## 65           ScaffoldIsActiveFlag
## 66              ScaffoldPosChunks
## 67                    ScaffoldSeq
## 68                      SwissProt
## 69                SwissProt2Locus
## 70               SwissProt2Pubmed
## 71                        Synonym
## 72                    SynonymType
## 73                       TIGRInfo
## 74                      TIGRroles
## 75                        TaxName
## 76                        TaxNode
## 77                 TaxParentChild
## 78                       Taxonomy
## 79                  Taxonomy2Tree
## 80                  TranscriptSeq
## 81                           Tree
## 82                          Users
## 83                      assoc_rel
## 84                    association
## 85           association_property
## 86          association_qualifier
## 87  association_species_qualifier
## 88                             db
## 89                         dbxref
## 90                       evidence
## 91                evidence_dbxref
## 92                   gene_product
## 93          gene_product_ancestor
## 94             gene_product_count
## 95          gene_product_homology
## 96          gene_product_homolset
## 97          gene_product_property
## 98               gene_product_seq
## 99            gene_product_subset
## 100          gene_product_synonym
## 101                    graph_path
## 102               graph_path2term
## 103                      homolset
## 104                 instance_data
## 105                           seq
## 106                    seq_dbxref
## 107                  seq_property
## 108                  source_audit
## 109                       species
## 110                          term
## 111                     term2term
## 112            term2term_metadata
## 113                    term_audit
## 114                   term_dbxref
## 115               term_definition
## 116                 term_property
## 117                   term_subset
## 118                  term_synonym
