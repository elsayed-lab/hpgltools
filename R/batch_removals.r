## Time-stamp: <Thu Mar 10 17:58:28 2016 Ashton Trey Belew (abelew@gmail.com)>

## Going to try and recapitulate the analyses found at:
## https://github.com/jtleek/svaseq/blob/master/recount.Rmd
## and use the results to attempt to more completely understand batch effects in our data

compare_batch_estimates <- function(start) {
    start_norm <- normalize_expt(start, convert="cpm")
    start_l2 <- normalize_expt(start, convert="cpm", transform="log2", filter_low=TRUE)
    start_low <- normalize_expt(start, filter_low=TRUE)
    data <- as.data.frame(Biobase::exprs(start_low$expressionset))
    mtrx <- as.matrix(data)
    design <- start$design
    conditions <- as.factor(design[, "condition"])
    batches <- as.factor(design[, "batch"])
    conditional_model <- model.matrix(~ conditions, data=data)
    null_model <- conditional_model[, 1]

    surrogate_estimate <- sva::num.sv(dat=mtrx, mod=conditional_model)
    control_likelihoods <- sva::empirical.controls(dat=mtrx, mod=conditional_model, mod0=null_model, n.sv=surrogate_estimate)
    supervised_sva_batch <- sva::svaseq(mtrx, conditional_model, null_model, controls=control_likelihoods)
    supervised_sva_batch_result <- supervised_sva_batch$sv

    ## Use unsupervised sva to estimate batch effect
    ## The data must be a matrix
    unsupervised_sva_batch <- sva::svaseq(mtrx, conditional_model, null_model)
    unsupervised_sva_batch_result <- unsupervised_sva_batch$sv

    ## Use a PCA batch test with rowMeans
    l2_data <- Biobase::exprs(start_l2$expressionset)
    pca_batch <- corpcor::fast.svd(l2_data - rowMeans(l2_data))$v[, 1]

    ## RUVSeq with the controls from sva
    ruv_batch_controlled <- RUVSeq::RUVg(mtrx, cIdx=as.logical(control_likelihoods), k=1)$W

    ## Use RUVSeq and residuals
    conditions <- as.factor(design[, "condition"])
    conditional_model <- model.matrix(~ conditions)
    ruv_input <- edgeR::DGEList(counts=data, group=conditions)
    ruv_input_norm <- edgeR::calcNormFactors(ruv_input, method="upperquartile")
    ruv_input_glm <- edgeR::estimateGLMCommonDisp(ruv_input_norm, ruv_design)
    ruv_input_tag <- edgeR::estimateGLMTagwiseDisp(ruv_input_glm, ruv_design)
    ruv_fit <- edgeR::glmFit(ruv_input_tag, ruv_design)
    ruv_res <- residuals(ruv_fit, type="deviance")
    ruv_normalized <- betweenLaneNormalization(mtrx, which="upper")  ## This also gets mad if you pass it a df and not matrix
    controls = rep(TRUE, dim(data)[1])
    ruv_batch = RUVSeq::RUVr(ruv_normalized, controls, k=1, ruv_res)$W

    ## Use RUVSeq with empirical controls
    ## The previous instance of ruv_input should work here, and the ruv_input_norm
    ## Ditto for _glm and _tag, and indeed ruv_fit
    ## Thus repeat the first 7 lines of the previous RUVSeq before anything changes.
    ruv_lrt <- glmLRT(ruv_fit, coef=2)
    ruv_controls = rank(ruv_lrt$table$LR) <= 400  ## what is going on here?!
    ruv_emp_batch <- RUVSeq::RUVg(mtrx, ruv_controls, k=1)$W

    plot(as.numeric(batches), pch=19, main="Batch assignment")
    plot(unsupervised_sva_batch_result, pch=19, main="SVA")
    plot(pca_batch, pch=19, main="PCA")
    plot(ruv_batch, pch=19, main="RUV residuals")
    plot(ruv_emp_batch, pch=19, main="RUV empirical")
    plot(unsupervised_sva_batch_result ~ batches, pch=19, main="unsupervised sva")
    plot(pca_batch ~ batches, pch=19, main="pca")
    plot(ruv_batch ~ batches, main="residual ruv")
    plot(ruv_emp_batch ~ batches, main="empirical controls ruv")

    estimate_batches <- cbind(conditions, batches, supervised_sva_batch_result, unsupervised_sva_batch_result,
                              pca_batch, ruv_batch_controlled, ruv_batch, ruv_emp_batch)
    colnames(estimate_batches) <- c("condition", "batch", "svasup", "svaun", "pca", "ruvsup", "ruvres", "ruvemp")

    tropical <- RSkittleBrewer::RSkittleBrewer('tropical')
    ## Since these can only go from 0->1 I am going to double and subtract 1 to get the full range
    ## no I guess that is a dumb thing to do
    ##corr = (cor(estimate_batches) * 2) - 1
    corr = cor(estimate_batches)
    corr
    cols = colorRampPalette(c(tropical[2], "white", tropical[1]))
    par(mar=c(5,5,5,5))
    ## If I read this plot correctly, the ruv residual method sucks for this data.
    corrplot::corrplot(corr, method="ellipse", type="lower", col=cols(100), tl.pos="d")
    ## A recap, each of these 6 things returned are additional factors which may be added to a model matrix passed to limma/edgeR/DESeq2
    ## in order to lower the effect of batch in the data.

    
##    dge <- DGEList(counts=dat0)
##    dge <- calcNormFactors(dge)
##    catplots = tstats = vector("list",6)
##    adj = c("+ study", "+ batch_unsup_sva",
##            "+ batch_ruv_res", "+ batch_ruv_emp",
##            "+ batch_pca","")
##    for(i in 1:6) {
##        design = model.matrix(as.formula(paste0("~ group",adj[i])))
##        v <- voom(dge,design,plot=FALSE)
##        fit <- lmFit(v,design)
##        fit <- eBayes(fit)
##        tstats[[i]] = abs(fit$t[,2])
##        names(tstats[[i]]) = as.character(1:dim(dat0)[1])
##        catplots[[i]] = CATplot(-rank(tstats[[i]]),-rank(tstats[[1]]),maxrank=1000,make.plot=F)
##    }
##
##    plot(catplots[[2]],ylim=c(0,1),col=trop[2],lwd=3,type="l",ylab="Concordance Between using study and different methods",xlab="Rank")
##    lines(catplots[[3]],col=trop[1],lwd=3,lty=2)
##    lines(catplots[[4]],col=trop[1],lwd=3)
##    lines(catplots[[5]],col=trop[3],lwd=3,lty=3)
##    lines(catplots[[6]],col=trop[4],lwd=3)
##    legend(200,0.5,legend=c("Unsup. svaseq","RUV Res.", "RUV Emp.","PCA","No adjustment"),col=trop[c(2,1,1,3,4)],lty=c(1,2,1,3,1),lwd=3)
}
