## Make sure I didn't introduce any stupid syntax errors.
expect_error(expect_error(library(hpgltools)))
## Make sure its friends load.
expect_error(expect_error(autoloads_all()))
## Load pasilla
expect_error(expect_error(require.auto("pasilla")))

## Load the pasilla data set
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single-end","single-end","paired-end",
        "paired-end","single-end","paired-end","paired-end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

## Make sure it is still possible to create an expt
expect_error(expect_error(pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)))

cbcb_data = counts
hpgl_data = exprs(pasilla_expt$expressionset)
expect_equal(cbcb_data, hpgl_data)

## Check that normalization tools work similarly
cbcb_quantile = cbcbSEQ::qNorm(cbcb_data)
hpgl_quantile_data = hpgltools::hpgl_norm(expt=pasilla_expt, transform="raw", norm="quant", convert="raw", filter_low=FALSE)
hpgl_quantile = hpgl_quantile_data$counts
expect_equal(cbcb_quantile, hpgl_quantile)

## log2/cpm that
cbcb_l2cpm_data = log2CPM(cbcb_quantile)
cbcb_l2cpm = cbcb_l2cpm_data$y
hpgl_l2cpm_data = hpgl_norm(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", filter_low=FALSE)
hpgl_l2cpm = hpgl_l2cpm_data$counts
expect_equal(cbcb_l2cpm, hpgl_l2cpm)
hpgl_l2cpm_expt = normalize_expt(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", filter_low=FALSE)
hpgl_l2cpm = exprs(hpgl_l2cpm_expt$expressionset)
expect_equal(cbcb_l2cpm, hpgl_l2cpm)

## Check that PCA invocations are similar
cbcb_svd = cbcbSEQ::makeSVD(cbcb_l2cpm)
hpgl_pca_info = hpgl_pca(expt=hpgl_l2cpm_expt)
hpgl_svd = hpgl_pca_info$pca
expect_equal(cbcb_svd$v, hpgl_svd$v)
expect_equal(cbcb_svd$d, hpgl_svd$d)

cbcb_res = cbcbSEQ::pcRes(cbcb_svd$v, cbcb_svd$d, design$condition, design$libType)
hpgl_res = hpgl_pca_info$res
expect_equal(cbcb_res, hpgl_res)

## invocation of batch correction
cbcb_batch = cbcbSEQ::combatMod(cbcb_l2cpm, batch=design$libType, mod=design$condition, noScale=TRUE)
hpgl_batch_expt = normalize_expt(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", batch="combatmod", filter_low=FALSE)
hpgl_batch = exprs(hpgl_batch_expt$expressionset)
expect_equal(cbcb_batch, hpgl_batch)

## voom invocation
cbcb_libsize = cbcb_l2cpm_data$lib.size
hpgl_libsize = hpgl_l2cpm_data$lib.size
expect_equal(cbcb_libsize, hpgl_libsize)
