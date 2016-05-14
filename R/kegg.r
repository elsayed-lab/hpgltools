## Time-stamp: <Sat May 14 03:47:35 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Print some data onto KEGG pathways.
#'
#' KEGGREST and pathview provide neat functions for coloring molecular pathways with arbitrary data.
#' Unfortunately they are somewhat evil to use.  This attempts to alleviate that.
#'
#' @param path_data Some differentially expressed genes.
#' @param indir Directory into which the unmodified kegg images will be downloaded (or already exist).
#' @param outdir Directory which will contain the colored images.
#' @param pathway Perform the coloring for a specific pathway?
#' @param species Kegg identifier for the species of interest.
#' @param string_from Regex to help in renaming KEGG categories/gene names from one format to another.
#' @param string_to Regex to help in renaming KEGG categories/gene names from one format to another.
#' @param suffix Add a suffix to the completed, colored files.
#' @param second_from Sometimes just one regex is not enough!
#' @param second_to Sometimes just one regex is not enough!
#' @param filenames Name the final files by id or name?
#' @return A list of some information for every KEGG pathway downloaded/examined.  This information includes:
#'   a. The filename of the final image for each pathway.
#'   b. The number of genes which were found in each pathway image.
#'   c. The number of genes in the 'up' category
#'   d. The number of genes in the 'down' category
#' @seealso \pkg{Ramigo} \pkg{pathview}
#' @examples
#' \dontrun{
#'  thy_el_comp2_path = hpgl_pathview(thy_el_comp2_kegg, species="spz", indir="pathview_in",
#'                                    outdir="kegg_thy_el_comp2", string_from="_Spy",
#'                                    string_to="_Spy_", filenames="pathname")
#' }
#' @export
hpgl_pathview <- function(path_data, indir="pathview_in", outdir="pathview", pathway="all", species="lma", string_from="LmjF", string_to="LMJF", suffix="_colored", second_from=NULL, second_to=NULL, filenames="id") {
    ## Please note that the KGML parser fails if other XML parsers are loaded into R
    ## eh = new.env(hash=TRUE, size=NA)
    ## There is a weird namespace conflict when using pathview, so I will reload it here
    try(detach("package:Rgraphviz", unload=TRUE), silent=TRUE)
    try(detach("package:topGO", unload=TRUE), silent=TRUE)
    try(detach("package:pathview", unload=TRUE), silent=TRUE)
    try(detach("package:KEGGgraph", unload=TRUE), silent=TRUE)
    try(detach("package:RamiGO", unload=TRUE), silent=TRUE)
    try(detach("package:graph", unload=TRUE), silent=TRUE)
    require.auto("pathview")
    ## If a table from limma was passed to this, just assume that one wants logFC
    ## Similar stanzas should probably be added for deseq/edger
    ## This is added because pathview() only works with dataframes/lists with only numbers.
    ## So it is wise to pull only the column of numbers one cares about.
    if (!is.null(path_data$logFC)) {
        tmp_data <- as.vector(path_data[, "logFC"])
        names(tmp_data) <- rownames(path_data)
        path_data <- tmp_data
        rm(tmp_data)
    }
    tmp_names <- names(path_data)
    tmp_names <- gsub(string_from, string_to, tmp_names)
    if (!is.null(second_from)) {
        tmp_names <- gsub(second_from, second_to, tmp_names)
    }
    ## tmp_names = gsub("\\.","_", tmp_names)
    names(path_data) <- tmp_names
    rm(tmp_names)

    ## First check that the input pathview directory exists
    if (!file.exists(indir)) {
        dir.create(indir)
    }
    if (!file.exists(outdir)){
        dir.create(outdir)
    }
    paths <- list()
    if (pathway == "all") {
        all_pathways <- unique(KEGGREST::keggLink("pathway", species))
        paths <- all_pathways
        paths <- gsub("path:", "", paths)
        ## all_modules = unique(KEGGREST::keggLink("module", species))
    } else if (class(pathway) == "list") {
        paths <- pathway
    } else {
        paths[1] <- pathway
    }
    return_list <- list()
    for (count in 1:length(paths)) {
        path <- paths[count]
        path_name <- KEGGREST::keggGet(path)
        path_name <- path_name[[1]]$NAME
        path_name <- gsub("(.*) - .*", "\\1", path_name)
        path_name <- tolower(path_name)
        path_name <- gsub(" ", "_", path_name)
        ## RCurl is crap and fails sometimes for no apparent reason.
        gene_examples <- try(KEGGREST::keggLink(paste("path", path, sep=":"))[,2])
        ## limits=c(min(path_data, na.rm=TRUE), max(path_data, na.rm=TRUE))
        limit_test <- c(abs(min(path_data, na.rm=TRUE)), abs(max(path_data, na.rm=TRUE)))
        limit_min <- -1.0 * max(limit_test)
        limit_max <- max(limit_test)
        limits <- c(limit_min, limit_max)
        print(paste("Here are some path gene examples: ", gene_examples, sep=""))
        print(paste("Here are your genes: ", head(names(path_data))), sep="")
        pv <- try(pathview::pathview(gene.data=path_data, kegg.dir=indir, pathway.id=path,
                                     species=species, limit=list(gene=limits, cpd=limits),
                                     map.null=TRUE, gene.idtype="KEGG", out.suffix=suffix,
                                     split.group=TRUE, expand.node=TRUE, kegg.native=TRUE,
                                     map.symbol=TRUE, same.layer=FALSE, res=1200,
                                     new.signature=FALSE, cex=0.05, key.pos="topright"))
        if (class(pv) == "numeric") {
            colored_genes <- NULL
            newfile <- NULL
            up <- NULL
            down <- NULL
        } else {
            colored_genes <- dim(pv$plot.data.gene)[1]
            ## "lma04070._proeff.png"
            oldfile <- paste(path, ".", suffix, ".png", sep="")
            ## An if-statement to see if the user prefers pathnames by kegg ID or pathway name
            ## Dr. McIver wants path names...
            if (filenames == "id") {
                newfile <- paste(outdir,"/", path, suffix, ".png", sep="")
            } else {  ## If filenames is not 'id', put in the path name...
                newfile <- paste(outdir, "/", path_name, suffix, ".png", sep="")
            }
            rename_try <- try(file.rename(from=oldfile, to=newfile), silent=TRUE)
            if (class(rename_try)[1] == 'try-error') {
                warning("There was an error renaming a png file, likely because it didn't download properly.")
                warning("It is likely easiest to just delete the pathview input directory.")
            }
            data_low <- summary(path_data)[2]
            data_high <- summary(path_data)[3]
            numbers_in_plot <- as.numeric(pv$plot.data.gene$mol.data)
            up <- sum(numbers_in_plot > data_high, na.rm=TRUE)
            down <- sum(numbers_in_plot < data_low, na.rm=TRUE)
        }
        return_list[[path]]$file <- newfile
        return_list[[path]]$genes <- colored_genes
        return_list[[path]]$up <- up
        return_list[[path]]$down <- down
        message(paste0(count, "/", length(paths), ": Finished ", path_name))
    } ## End for loop

    retdf <- data.frame(rep(NA, length(names(return_list))))
    rownames(retdf) <- names(return_list)
    retdf$genes <- NA
    retdf$up <- NA
    retdf$down <- NA
    colnames(retdf) <- c("file","genes","up","down")
    for (path in names(return_list)) {
        if (is.null(return_list[[path]]$genes)) {
            retdf[path,]$genes <- 0
        } else {
            retdf[path,]$genes <- as.numeric(return_list[[path]]$genes)
        }
        retdf[path,]$file <- as.character(return_list[[path]]$file)
        retdf[path,]$up <- as.numeric(return_list[[path]]$up)
        retdf[path,]$down <- as.numeric(return_list[[path]]$down)
    }
    retdf <- retdf[with(retdf, order(up, down)), ]
    return(retdf)
}

#' Use gostats() against kegg pathways.
#'
#' This sets up a GSEABase analysis using KEGG pathways rather than gene ontologies.
#' Does this even work?  I don't think I have ever tested it yet.
#' oh, it sort of does, maybe if I export it I will rembmer it.
#'
#' @param organism The organism used to make the KEGG frame, human readable no taxonomic.
#' @param pathdb Name of the pathway database for this organism.
#' @param godb Name of the ontology database for this organism.
#' @return Results from hyperGTest using the KEGG pathways.
#' @export
gostats_kegg <- function(organism="Homo sapiens", pathdb="org.Hs.egPATH", godb="org.Hs.egGO") {
    org <- get0(pathdb)
    org_go <- get0(godb)
    frame <- AnnotationDbi::toTable(org)
    keggframedata <- data.frame(frame$path_id, frame$gene_id)
    keggFrame <- AnnotationDbi::KEGGFrame(keggframedata, organism=organism)
    gsc <- GSEABase::GeneSetCollection(keggFrame, setType=GSEABase::KEGGCollection())
    universe <- AnnotationDbi::Lkeys(org_go)
    genes <- universe[1:500]
    kparams <- Category::GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
                                              geneSetCollection=gsc,
                                              geneIds=genes,
                                              universeGeneIds=universe,
                                              pvalueCutoff=0.05,
                                              testDirection="over")
    kOver <- Category::hyperGTest(kparams)
    return(kOver)
}

#' Search KEGG identifiers for a given species name.
#'
#' KEGG identifiers do not always make sense.  For example, how am I supposed to remember that
#' Leishmania major is lmj?  This takes in a human readable string and finds the KEGG identifiers
#' that match it.
#'
#' @param species Search string (Something like 'Homo sapiens').
#' @param short Only pull the orgid?
#' @return Data frame of possible KEGG identifier codes, genome ID numbers, species, and
#'     phylogenetic classifications.
#' @seealso \pkg{RCurl}
#' @examples
#' \dontrun{
#'  fun = kegg_get_orgn('Canis')
#' ## >     Tid     orgid      species                   phylogeny
#' ## >  17 T01007   cfa Canis familiaris (dog) Eukaryotes;Animals;Vertebrates;Mammals
#' }
#' @export
kegg_get_orgn <- function(species="Leishmania", short=TRUE) {
    all_organisms <- RCurl::getURL("http://rest.kegg.jp/list/organism")
    org_tsv <- textConnection(all_organisms)
    all <- read.table(org_tsv, sep="\t", quote="", fill=TRUE)
    close(org_tsv)
    colnames(all) <- c("Tid","orgid","species","phylogeny")
    candidates <- all[grepl(species, all$species), ]  ## Look for the search string in the species column
    if (isTRUE(short)) {
        candidates <- as.character(candidates$orgid)
    }
    return(candidates)
}

#' Extract the percent differentially expressed genes for all KEGG pathways.
#'
#' KEGGgraph provides some interesting functionality for mapping KEGGids and examining the
#' pieces. This attempts to use that in order to evaluate how many 'significant' genes are in a
#' given pathway.
#'
#' @param all_ids Set of all gene IDs in a given analysis.
#' @param sig_ids Set of significant gene IDs.
#' @param organism KEGG organism identifier.
#' @param pathways What pathways to look at?
#' @param pathdir Directory into which to copy downloaded pathway files.
#' @return Dataframe including the filenames, percentages, nodes included, and differential nodes.
#' @seealso \pkg{KEGGgraph} \pkg{KEGGREST}
#' @export
pct_all_kegg <- function(all_ids, sig_ids, organism="dme", pathways="all", pathdir="kegg_pathways", ...) {
    arglist <- list(...)
    if (!file.exists(pathdir)) {
        dir.create(pathdir)
    }
    paths <- NULL
    if (class(pathways) == "character") {
        if (length(pathways) > 1) {
            paths <- paste0("path:", organism, pathways)
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

    percentages <- list()
    filenames <- list()
    path_names <- list()
    diff_nodes <- list()
    path_nodes <- list()
    for (count in 1:length(paths)) {
        path <- paths[count]
        path_name <- try(KEGGREST::keggGet(path), silent=TRUE)
        if (class(path_name) == "try-error") {
            path_names[count] <- NA
            filenames[count] <- NA
            percentages[count] <- NA
            path_nodes[count] <- NA
            diff_nodes[count] <- NA
        } else {
            path_name <- path_name[[1]]$NAME
            path_name <- gsub("(.*) - .*", "\\1", path_name)
            path_name <- tolower(path_name)
            path_name <- gsub(" ", "_", path_name)
            path_names[count] <- path_name
            message(paste0("Extracting data for ", path, ": ", path_name, "."))
            pct_diff <- pct_kegg_diff(all_ids, sig_ids, pathway=path, organism=organism,
                                                 pathdir=pathdir)
            filenames[count] <- pct_diff[["filename"]]
            percentages[count] <- pct_diff[["percent"]]
            path_nodes[count] <- pct_diff[["all_nodes"]]
            diff_nodes[count] <- pct_diff[["diff_nodes"]]
            message(paste0("The path: ", path_names[count], " was written to ", filenames[count], " and has ", percentages[count], "% diff."))
        }
    }
    path_data <- as.data.frame(cbind(pathways, path_names, filenames, percentages, path_nodes, diff_nodes))
    colnames(path_data) <- c("pathway","path_name","filename","percent","all_nodes","diff_nodes")
    path_data[["pathway"]] <- as.character(path_data[["pathway"]])
    path_data[["path_name"]] <- as.character(path_data[["path_name"]])
    path_data[["filename"]] <- as.character(path_data[["filename"]])
    path_data[["percent"]] <- as.numeric(path_data[["percent"]])
    path_data[["all_nodes"]] <- as.character(path_data[["all_nodes"]])
    path_data[["diff_nodes"]] <- as.character(path_data[["diff_nodes"]])
    return(path_data)
}

#' Extract the percent differentially expressed genes in a given KEGG pathway.
#'
#' KEGGgraph provides some interesting functionality for mapping KEGGids and examining the
#' pieces. This attempts to use that in order to evaluate how many 'significant' genes are in a
#' given pathway.
#'
#' @param all_ids Set of all gene IDs in a given analysis.
#' @param sig_ids Set of significant gene IDs.
#' @param pathway Numeric pathway identifier.
#' @param organism KEGG organism identifier.
#' @param pathdir Directory into which to copy downloaded pathway files.
#' @return Percent genes/pathway deemed significant.
#' @seealso \pkg{KEGGgraph} \pkg{KEGGREST}
#' @export
pct_kegg_diff <- function(all_ids, sig_ids, pathway="00500", organism="dme", pathdir="kegg_pathways", ...) {
    arglist <- list(...)
    if (!file.exists(pathdir)) {
        dir.create(pathdir)
    }
    pathway <- gsub(pattern=organism, replacement="", x=pathway)
    pathway <- gsub(pattern="path:", replacement="", x=pathway)
    filename <- paste0(pathdir, "/", organism, pathway, ".xml")
    pathwayid <- paste0(organism, pathway)
    retrieved <- NULL
    if (file.exists(filename)) {
        message("The file already exists, loading from it.")
        retrieved <- filename
    } else {
        log <- capture.output(type="message", { retrieved <- try(suppressMessages(KEGGgraph::retrieveKGML(pathwayid=pathway, organism=organism, destfile=filename, method="internal")), silent=TRUE); })
        if (class(retrieved) == "try-error") {
            retlist <- list(
                "pathway" = pathway,
                "filename" = unavailable,
                "percent" = NA,
                "all_nodes" = NULL,
                "diff_nodes" = NULL)
            return(retlist)
        }
    }
    parse_result <- KEGGgraph::parseKGML2Graph(filename, expandGenes=TRUE)

    all_keggids <- KEGGgraph::translateGeneID2KEGGID(all_ids, organism=organism)
    de_keggids <- KEGGgraph::translateGeneID2KEGGID(sig_ids, organism=organism)

    possible_nodes <- KEGGgraph::nodes(parse_result)
    is_differential <- possible_nodes %in% de_keggids
    found_nodes <- possible_nodes[is_differential]
    pct_diff <- signif(mean(is_differential) * 100.0, 4)
    path_data <- KEGGREST::keggGet(pathwayid)
    path_name <- path_data[[1]][["NAME"]]
    message(paste0(pct_diff, "% genes differentially expressed in pathway ", pathway, ": '", path_name, "'."))
    retlist <- list(
        "pathway" = pathway,
        "filename" = filename,
        "percent" = pct_diff,
        "all_nodes" = toString(possible_nodes),
        "diff_nodes" = toString(found_nodes))
    return(retlist)
}

play_kegggraph <- function() {
    library(KEGGgraph)
    map <- system.file("extdata/hsa04010.xml", package="KEGGgraph")
    mapkpathway <- parseKGML(map)
    mapkpathway
    map_graph <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
    map_nodes <- nodes(map_graph)
    map_edges <- KEGGgraph::edges(map_graph)
    node_data <- getKEGGnodeData(map_graph)
    library(Rgraphviz)
    mapped_graph <- subGraph(map_nodes, map_graph)
    mapped_graph
    outs <- sapply(edges(mapped_graph), length) > 0
    ins <- sapply(inEdges(mapped_graph), length) > 0
    ios <- outs | ins
    library("org.Hs.eg.db")
    io_gene_ids <- translateKEGGID2GeneID(names(ios))
    node_names <- sapply(mget(io_gene_ids, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
    names(node_names) <- names(ios)
    nattrs <- list()
    nattrs$fillcolor <- Rgraphviz::makeNodeAttrs(mapped_graph, "lightgrey", list(orang=names(ios)[ios]))
    nattrs$label <- node_names
    plot(mapped_graph, "neato", nodeAttrs=nattrs, attrs=list(node=list(fillcolor="lightgreen", width=0.75, shape="ellipse"), edge=list(arrowsize=0.7)))
}


## EOF
