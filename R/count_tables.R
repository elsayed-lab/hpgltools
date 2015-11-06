## Time-stamp: <Thu Nov  5 14:24:32 2015 Ashton Trey Belew (abelew@gmail.com)>

#' create_expt()  Wrap bioconductor's expressionset to include some other extraneous
#' information.  This simply calls create_experiment and then does
#' expt_subset for everything
#'
#' @param file  a comma separated file describing the samples with
#' information like condition,batch,count_filename,etc
#' @param color_hash default=NULL  a hash which describes how to color the samples,
#' it will generate its own colors using colorBrewer
#' @param suffix default='.count.gz'  when looking for the count tables in processed_data
#' look for this suffix on the end of the files.
#' @param header default=FALSE  Does the csv metadata file have a header?
#' @param genes default=NULL  annotation information describing the rows of the data set, usually
#' this comes from a call to import.gff()
#' @param by_type default=FALSE  when looking for count tables, are they organized by type?
#' @param by_sample default=FALSE  or by sample?  I do all mine by sample, but others do by type...
#' @param sep default=','  some people prefer their csv files as tab or semicolon separated.
#' @param include_type default='all'  I have usually assumed that all gff annotations should be used,
#' but that is not always true, this allows one to limit.
#' @param include_gff default=NULL  A gff file to help in sorting which features to keep
#' @param count_dataframe default=NULL  If one does not wish to read the count tables from processed_data/
#' they may instead be fed here
#' @param meta_dataframe default=NULL  an optional dataframe containing the metadata rather than a file
#' @param savefile default='expt'  an Rdata filename prefix for saving the data of the resulting expt.
#' @param low_files default=FALSE  whether or not to explicitly lowercase the filenames when searching in processed_data/
#' this is relevant because the ceph object storage by default lowercases filenames.
#'
#' It is worth noting that this function has a lot of logic used to
#' find the count tables in the local filesystem.  This logic has been
#' superceded by simply adding a field to the .csv file called
#' 'file'.  create_expt() will then just read that filename, it may be
#' a full pathname or local to the cwd of the project.
#'
#' @return  experiment an expressionset
#' @seealso \code{\link{pData}}, \code{\link{fData}},
#' \code{\link{exprs}}, \code{\link{hpgl_read_files}},
#' \code{\link{as.list.hash}}
#'
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Remember that this depends on an existing data structure of gene annotations.
create_expt = function(file=NULL, color_hash=NULL, suffix=".count.gz", header=FALSE, gene_info=NULL, by_type=FALSE, by_sample=FALSE, sep=",", include_type="all", include_gff=NULL, count_dataframe=NULL, meta_dataframe=NULL, savefile="expt", low_files=FALSE, ...) {
    if (is.null(meta_dataframe) & is.null(file)) {
        stop("This requires either a csv file or dataframe of metadata describing the samples.")
    } else if (is.null(file)) {
        tmp_definitions = meta_dataframe
    }  else {
        if (file_ext(file) == 'csv') {
            tmp_definitions = read.csv(file=file, comment.char="#", sep=sep)
        } else if (file_ext(file) == 'xls' | file_ext(file) == 'xlsx') {
            xls = loadWorkbook(file, create=FALSE)
            tmp_definitions = readWorksheet(xls, 1)
        } else {
            tmp_definitions = read.table(file=file)
        }
    }
    colnames(tmp_definitions) = tolower(colnames(tmp_definitions))
    tmp_definitions = subset(tmp_definitions, sample.id != "")
    condition_names = unique(tmp_definitions$condition)
    if (is.null(condition_names)) {
        warning("There is no 'condition' field in the definitions, this will make many analyses more difficult/impossible.")
    }
    if (is.null(color_hash)) {
        if (is.null(tmp_definitions$color)) {
            num_colors = length(condition_names)
            colors = suppressWarnings(colorRampPalette(brewer.pal(num_colors,"Dark2"))(num_colors))
            color_hash = hash(keys=as.character(condition_names), values=colors)
        } else {
            color_hash = hash(keys=as.character(tmp_definitions$sample.id), values=tmp_definitions$color)            
        }
    }
    ## Sometimes, R adds extra rows on the bottom of the data frame using this command.
    ## Thus the next line
    print("This function needs the conditions and batches to be an explicit column in the sample sheet.")
    expt_list = create_experiment(file=file, color_hash, suffix=suffix, header=header, gene_info=gene_info, by_type=by_type, by_sample=by_sample, count_dataframe=count_dataframe, meta_dataframe=meta_dataframe, sep=sep, low_files=low_files, include_type=include_type, include_gff=include_gff)
    expt = expt_list$expt
    def = expt_list$def
    new_expt = expt_subset(expt)
    new_expt$definitions = def

    if (!is.null(expt_list$annotation)) {
        new_expt$annotation = expt_list$annotation
        new_expt$gff_file = include_gff
        tmp_genes = new_expt$annotation
        tmp_genes = tmp_genes[tmp_genes$type == "gene",]
        rownames(tmp_genes) = make.names(tmp_genes$Name, unique=TRUE)
        tooltip_data = tmp_genes
        tooltip_data = tooltip_data[,c(11,12)]
        tooltip_data$tooltip = paste(tooltip_data$Name, tooltip_data$description, sep=": ")
        tooltip_data$tooltip = gsub("\\+", " ", tooltip_data$tooltip)
        rownames(tooltip_data) = tooltip_data$Name
        tooltip_data = tooltip_data[-1]
        tooltip_data = tooltip_data[-1]
        colnames(tooltip_data) = c("name.tooltip")
        new_expt$genes = genes
        new_expt$tooltip = tooltip_data
    }
    ## These entries in new_expt are intended to maintain a record of
    ## the transformation status of the data, thus if we now call
    ## normalize_expt() it should change these.
    ## Therefore, if we call a function like DESeq() which requires
    ## non-log2 counts, we can check these values and convert accordingly
    new_expt$filtered = FALSE
    new_expt$transform = "raw"
    new_expt$norm = "raw"
    new_expt$convert = "raw"
    new_expt$original_libsize = colSums(exprs(new_expt$expressionset))
    if (!is.null(savefile)) {
        save(list = c("new_expt"), file=paste(savefile, ".Rdata", sep=""))
    }
    return(new_expt)
}

#' create_experiment()  Wrap bioconductor's expressionset to include some other extraneous
#' information.
#'
#' @param file default=NULL  a comma separated file describing the samples with
#' information like condition,batch,count_filename,etc.
#' @param color_hash  a hash which describes how to color the samples
#' @param suffix default='.count.gz'  when looking for the count tables in processed_data
#' look for this suffix on the end of the files.
#' @param header default=FALSE  Does the csv metadata file have a header?
#' @param genes default=NULL  annotation information describing the rows of the data set, usually
#' this comes from a call to import.gff()
#' @param by_type default=FALSE  when looking for count tables, are they organized by type?
#' @param by_sample default=FALSE  or by sample?  I do all mine by sample, but others do by type...
#' @param sep default=','  some people prefer their csv files as tab or semicolon separated.
#' @param include_type default='all'  I have usually assumed that all gff annotations should be used,
#' but that is not always true, this allows one to limit.
#' @param include_gff default=NULL  A gff file to help in sorting which features to keep
#' @param count_dataframe default=NULL  If one does not wish to read the count tables from processed_data/
#' they may instead be fed here
#' @param meta_dataframe default=NULL  an optional dataframe containing the metadata rather than a file
#'
#' @return  experiment an expressionset
#' @seealso \code{\link{pData}}, \code{\link{fData}},
#' \code{\link{exprs}}, \code{\link{hpgl_read_files}},
#' \code{\link{as.list.hash}}
#'
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
create_experiment = function(file=NULL, color_hash, suffix=".count.gz", header=FALSE, gene_info=NULL, by_type=FALSE, by_sample=FALSE, include_type="all", include_gff=NULL, count_dataframe=NULL, meta_dataframe=NULL, sep=",", ...) {
    print("Please note that thus function assumes a specific set of columns in the sample sheet:")
    print("The most important ones are: Sample.ID, Stage, Type.")
    print("Other columns it will attempt to create by itself, but if")
    print("batch and condition are provided, that is a nice help.")
    if (is.null(meta_dataframe) & is.null(file)) {
        stop("This requires either a csv file or dataframe of metadata describing the samples.")
    } else if (is.null(file)) {
        sample_definitions = meta_dataframe
    } else {
        if (file_ext(file) == 'csv') {
            sample_definitions = read.csv(file=file, comment.char="#", sep=sep)
        } else if (file_ext(file) == 'xls' | file_ext(file) == 'xlsx') {
            xls = loadWorkbook(file, create=FALSE)
            sample_definitions = readWorksheet(xls, 1)
        } else {
            sample_definitions = read.table(file=file)
        }
    }
    colnames(sample_definitions) = tolower(colnames(sample_definitions))
    ##sample_definitions = sample_definitions[grepl('(^HPGL|^hpgl)', sample_definitions$sample.id, perl=TRUE),]
    if (is.null(sample_definitions$condition)) {
        sample_definitions$condition = tolower(paste(sample_definitions$type, sample_definitions$stage, sep="_"))
        sample_definitions$batch = gsub("\\s+|\\d+|\\*", "", sample_definitions$batch, perl=TRUE)
    }
    design_colors_list = as.list.hash(color_hash)
    sample_definitions$colors = as.list(design_colors_list[as.character(sample_definitions$condition)])
    sample_definitions = as.data.frame(sample_definitions)
    rownames(sample_definitions) = make.names(sample_definitions$sample.id, unique=TRUE)

    ## The logic here is that I want by_type to be the default, but only
    ## if no one chooses either.
    filenames = NULL
    found_counts = NULL
    ## This stanza allows one to have a field 'file' in the csv with filenames for the count tables.
    if (!is.null(sample_definitions$file)) {
        filenames = sample_definitions$file
        all_count_tables = hpgl_read_files(as.character(sample_definitions$sample.id),
            as.character(filenames), header=header, suffix=NULL)
        ## This stanza allows one to fill in the count tables with an external data frame.
    } else if (!is.null(count_dataframe)) {
        all_count_tables = count_dataframe
        colnames(all_count_tables) = rownames(sample_definitions)
        ## If neither of these cases is true, start looking for the files in the processed_data/ directory
    } else if (!isTRUE(by_type) & !isTRUE(by_sample) & is.null(filenames)) {
        ## If neither by_type or by_sample is set, look first by sample
        sample_definitions$counts = paste("processed_data/count_tables/", as.character(sample_definitions$sample.id), "/", as.character(sample_definitions$sample.id), suffix, sep="")
        found_counts = try(hpgl_read_files(as.character(sample_definitions$sample.id), as.character(sample_definitions$counts), header=header, suffix=suffix))
        if (class(found_counts) == 'try-error') {
            ## Then try by-type
            sample_definitions$counts = paste("processed_data/count_tables/", tolower(sample_definitions$type), "/", tolower(sample_definitions$stage), "/", sample_definitions$sample.id, suffix, sep="")
            found_counts = try(hpgl_read_files(as.character(sample_definitions$sample.id), as.character(sample_definitions$counts), header=header, suffix=suffix))
            if (class(found_counts) == 'try-error') {
                stop("Unable to find count tables, either by sample id nor by type")
            } else {
                all_count_tables = found_counts
            }
        } else {
            all_count_tables = found_counts
        }
    } else if (isTRUE(by_sample)) {
        sample_definitions$counts = paste("processed_data/count_tables/", as.character(sample_definitions$sample.id), "/", as.character(sample_definitions$sample.id), suffix, sep="")
        all_count_tables = try(hpgl_read_files(as.character(sample_definitions$sample.id), as.character(sample_definitions$counts), header=header, suffix=suffix))
    } else if (isTRUE(by_type)) {
        sample_definitions$counts = paste("processed_data/count_tables/", tolower(sample_definitions$type), "/", tolower(sample_definitions$stage), "/", sample_definitions$sample.id, suffix, sep="")
        all_count_tables = try(hpgl_read_files(as.character(sample_definitions$sample.id), as.character(sample_definitions$counts), header=header, suffix=suffix))
    } ## End checking by_type/by_samples
    all_count_matrix = as.data.frame(all_count_tables)

    rownames(all_count_matrix) = gsub("^exon:","", rownames(all_count_matrix))
    rownames(all_count_matrix) = make.names(gsub(":\\d+","", rownames(all_count_matrix)), unique=TRUE)
    if (is.null(gene_info)) {
        gene_info = data.frame(all_count_matrix)
    } else {
        if (is.null(gene_info$ID)) {
            gene_info$ID = rownames(gene_info)
        }
        gene_info = gene_info[gene_info$ID %in% rownames(all_count_matrix),]
        all_count_matrix = all_count_matrix[rownames(all_count_matrix) %in% genes$ID,]
    }

    ## Make sure that all columns have been filled in for every gene.
    complete_index = complete.cases(all_count_matrix)
    all_count_matrix = all_count_matrix[complete_index,]

    annotation = NULL
    if (!is.null(include_gff)) {
        if (include_type == "all") {
            print("Reading the annotation information, this may take a while.")
            annotation = BiocGenerics::as.data.frame(rtracklayer::import(include_gff, asRangedData=FALSE))
            print("Finished reading annotations, we should be done soon.")
        } else {
            print(paste("Excluding entries from the annotation which are not: ", include_type, sep=""))
            print("Reading the annotation information, this may take a while.")
            annotation = BiocGenerics::as.data.frame(rtracklayer::import(include_gff, asRangedData=FALSE))
            print("Finished reading annotations, we should be done soon.")
            keepers = annotation[annotation$type==include_type,]$gene_id
            index = row.names(all_count_matrix) %in% keepers
            all_count_matrix = all_count_matrix[index,]
            gene_info = gene_info[index,]
        }
    }

    if (is.null(sample_definitions$stage)) {
        sample_definitions$stage = "unknown"
    }
    if (is.null(sample_definitions$type)) {
        sample_definitions$type = "unknown"
    }
    if (is.null(sample_definitions$condition)) {
        sample_definitions$condition = "unknown"
    }
    if (is.null(sample_definitions$batch)) {
        sample_definitions$batch = "unknown"
    }
    if (is.null(sample_definitions$colors)) {
        sample_definitions$colors = "unknown"
    }
    if (is.null(sample_definitions$counts)) {
        sample_definitions$counts = "unknown"
    }
    if (is.null(sample_definitions$intercounts)) {
        sample_definitions$intercounts = "unknown"
    }
    metadata = new("AnnotatedDataFrame", data.frame(sample=as.character(sample_definitions$sample.id),
        stage=as.character(sample_definitions$stage),
        type=as.character(sample_definitions$type),
        condition=as.character(sample_definitions$condition),
        batch=as.character(sample_definitions$batch),
        color=as.character(sample_definitions$colors),
        counts=sample_definitions$counts,
        intercounts=sample_definitions$intercounts))
    sampleNames(metadata) = colnames(all_count_matrix)
    feature_data = new("AnnotatedDataFrame", gene_info)
    featureNames(feature_data) = rownames(all_count_matrix)
    experiment = new("ExpressionSet", exprs=all_count_matrix,
                     phenoData=metadata, featureData=feature_data)
    ret = list(expt=experiment, def=sample_definitions, annotation=annotation)
    return(ret)
}

#' expt_subset()  Extract a subset of samples following some rule(s) from an
#' experiment class
#'
#' @param expt  an expt which is a home-grown class containing an
#' expressionSet, design, colors, etc.
#' @param subset  a valid R expression which defines a subset of the
#' design to keep.
#' @param by_definitions default=TRUE  whether to use the definitions dataframe or design from pData(expressionset).
#' The definitions dataframe often has much more information.
#'
#' @return  metadata an expt class which contains the smaller set of
#' data
#' @seealso \code{\link{ExpressionSet}}, \code{\link{pData}},
#' \code{\link{exprs}}, and \code{\link{fData}}
#'
#' @export
#' @examples
#' ## smaller_expt = expt_subset(big_expt, "condition=='control'")
#' ## all_expt = expt_subset(expressionset, "")  ## extracts everything
expt_subset = function(expt, subset=NULL) {
    if (class(expt) == "ExpressionSet") {
        expressionset = expt
    } else if (class(expt) == "expt") {
        expressionset = expt$expressionset
    } else {
        stop("expt is neither an expt nor ExpressionSet")
    }
    if (is.null(expt$definitions)) {
        warning("There is no expt$definitions, using the expressionset.")
        initial_metadata = pData(expressionset)
    } else {
        initial_metadata = expt$definitions
    }
    if (is.null(subset)) {
        samples = initial_metadata
    } else {
        r_expression=paste("subset(initial_metadata,", subset, ")")
        samples = eval(parse(text=r_expression))
        ## design = data.frame(sample=samples$sample, condition=samples$condition, batch=samples$batch)
    }
    design = as.data.frame(samples)
    ## This is to get around stupidity with respect to needing all factors to be in a DESeqDataSet
    conditions = as.factor(as.character(design$condition))
    batches = as.factor(as.character(design$batch))
    design$condition = conditions
    design$batch = batches
    if (is.null(samples$sample.id)) {
        samplenames = as.character(samples$sample)
    } else {
        samplenames = as.character(samples$sample.id)
    }
    colors = as.character(samples$color)
    names = paste(conditions, batches, sep="-")
    subset_definitions = expt$definitions[rownames(expt$definitions) %in% samplenames, ]
    subset_libsize = expt$original_libsize[names(expt$original_libsize) %in% samplenames]
    expressionset = expressionset[, sampleNames(expressionset) %in% samplenames]
    columns = data.frame(sample=colnames(exprs(expressionset)))
    rownames(columns) = colnames(exprs(expressionset))
    metadata = list(initial_metadata=initial_metadata,
        original_expressionset=expressionset,
        expressionset=expressionset,
        samples=samples,
        design=design,
        definitions=subset_definitions,
        stages=initial_metadata$stage,
        types=initial_metadata$type,
        conditions=conditions,
        batches=batches,
        samplenames=samplenames,
        colors=colors,
        names=names,
        filtered=expt$filtered,
        transform=expt$transform,
        norm=expt$norm,
        convert=expt$convert,
        original_libsize=subset_libsize,
        columns=columns)
    class(metadata) = "expt"
    return(metadata)
}

#' hpgl_read_files()  Read a bunch of count tables and create a usable data frame from
#' them.
#'
#' @param ids  a list of experimental ids
#' @param files  a list of files to read
#' @param header default=FALSE  whether or not the count tables include a header row.
#' @param include_summary_rows default=FALSE  whether HTSeq summary rows should be included.
#' @param suffix default=NULL  an optional suffix to add to the filenames when reading them.
#'
#' It is worth noting that this function has some logic intended for the elsayed lab's data storage structure.
#' It shouldn't interfere with other usages, but it attempts to take into account different ways the data might be stored.
#'
#' @return  count_table a data frame of count tables
#' @seealso \code{\link{create_experiment}}
#'
#' @export
#' @examples
#' ## count_tables = hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
hpgl_read_files = function(ids, files, header=FALSE, include_summary_rows=FALSE, suffix=NULL, ...) {
    ## load first sample
    lower_filenames = files
    dirs = dirname(lower_filenames)
    low_files = tolower(basename(files))
    if (!is.null(suffix)) {
        low_hpgl = gsub(suffix, "", basename(files))
        low_hpgl = tolower(low_hpgl)
        low_hpgl = paste(low_hpgl, suffix, sep="")
    } else {
        low_hpgl = gsub("HPGL","hpgl", basename(files))
    }
    lower_filenames = paste(dirs, low_files, sep="/")
    lowhpgl_filenames = paste(dirs, low_hpgl, sep="/")
    if (file.exists(tolower(files[1]))) {
        files[1] = tolower(files[1])
    } else if (file.exists(lowhpgl_filenames[1])) {
        files[1] = lowhpgl_filenames[1]
    } else if (file.exists(lower_filenames[1])) {
        files[1] = lower_filenames[1]
    }
    ##count_table = read.table(files[1], header=header, ...)
    count_table = try(read.table(files[1], header=header))
    if (class(count_table)[1] == 'try-error') {
        stop(paste0("There was an error reading: ", files[1]))
    }
    print(paste0(files[1], " contains ", length(rownames(count_table)), " rows."))
    colnames(count_table) = c("ID", ids[1])
    ## iterate over and append remaining samples
    for (table in 2:length(files)) {
        if (file.exists(tolower(files[table]))) {
            files[table] = tolower(files[table])
        } else if (file.exists(lowhpgl_filenames[table])) {
            files[table] = lowhpgl_filenames[table]
        } else if (file.exists(lower_filenames[table])) {
            files[table] = lower_filenames[table]
        }
        tmp_count = try(read.table(files[table], header=header))
        if (class(tmp_count)[1] == 'try-error') {
            stop(paste0("There was an error reading: ", files[table]))
        }        
        colnames(tmp_count) = c("ID", ids[table])
        pre_merge = length(rownames(tmp_count))
        count_table = merge(count_table, tmp_count, by="ID")
        post_merge = length(rownames(count_table))
        print(paste0(files[table], " contains ", pre_merge, " rows and merges to ", post_merge, " rows."))
    }

    rm(tmp_count)
    ## set row and columns ids
    rownames(count_table) = count_table$ID
    count_table = count_table[-1]
    colnames(count_table) = ids

    ## remove summary fields added by HTSeq
    if (!include_summary_rows) {
        htseq_meta_rows = c('__no_feature', '__ambiguous', '__too_low_aQual',
                            '__not_aligned', '__alignment_not_unique')
        count_table = count_table[!rownames(count_table) %in% htseq_meta_rows,]
    }
    return(count_table)
}


#' concatenate_runs()  Sum the reads/gene for multiple sequencing runs of a single condition/batch
#'
#' @param expt  an experiment class containing the requisite metadata and count tables
#' @param column default='replicate'  a column of the design matrix used to specify which samples are replicates
#'
#' @return the input expt with the new design matrix, batches, conditions, colors, and count tables.
#' @export
#' @examples
#' ## compressed = concatenate_runs(expt)
concatenate_runs = function(expt, column='replicate') {
    data = exprs(expt$expressionset)
    design = expt$definitions
    replicates = levels(as.factor(design[,column]))
    final_expt = expt
    final_data = NULL
    final_design = NULL
    column_names = list()
    colors = list()
    conditions = list()
    batches = list()
    names = list()
    for (rep in replicates) {
        expression = paste0(column, "=='", rep, "'")
        tmp_expt =  expt_subset(expt, expression)
        tmp_data =  rowSums(exprs(tmp_expt$expressionset))
        tmp_design = tmp_expt$design[1,]
        final_data = cbind(final_data, tmp_data)
        final_design = rbind(final_design, tmp_design)
        column_names[[rep]] = as.character(tmp_design$sample.id)
        colors[[rep]] = as.character(tmp_design$color)
        batches[[rep]] = as.character(tmp_design$batch)
        conditions[[rep]] = as.character(tmp_design$condition)
        names[[rep]] = paste(conditions[[rep]], batches[[rep]], sep='-')
        colnames(final_data) = column_names
    }
    final_expt$design = final_design
    metadata = new("AnnotatedDataFrame", final_design)
    sampleNames(metadata) = colnames(final_data)
    feature_data = new("AnnotatedDataFrame", fData(expt$expressionset))
    featureNames(feature_data) = rownames(final_data)
    experiment = new("ExpressionSet", exprs=final_data,
                     phenoData=metadata, featureData=feature_data)
    final_expt$expressionset = experiment
    final_expt$original_expressionset = experiment
    final_expt$samples = final_design
    final_expt$colors = as.character(colors)
    final_expt$batches = as.character(batches)
    final_expt$conditions = as.character(conditions)
    final_expt$names = as.character(names)
    return(final_expt)
}

## EOF
