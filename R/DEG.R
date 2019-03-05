#' Call DEGenes
#'
#' This function loads the data as a dataframe, and method as a string.
#' It assumes that each line contains gene expression profile of one single
#' cell, and each column contains the one single gene expression profile in
#' different cells. The dataframe should also contain the cell type information
#' with column name 'cell_type', as well as group information as 'compare_group'
#' Batch information as 'batch' is optional. If included, users may want to use
#' the raw count data for later analysis. Differential expressed genes will be
#' called within each cell type by the method users select. For bulk RNAseq,
#' we provide edgeR, DESeq2. And for scRNA-seq, popular methods in packages
#' scde, monocle, DEsingle and MAST are available.
#'
#' @param data Input raw or normalized count data with column 'cell_type'
#' and 'compare_group'
#' @param method Method used to call DEGenes. Available options are:
##' \itemize{
##'  \item{Wilcox}: Wilcoxon rank sum test
##'  \item{DESeq2}: Negative binomial model based differential analysis
##'  (Love et al, Genome Biology, 2014)
##'  \item{SCDE}: Bayesian approach to single-cell differential
##'  expression analysis (Kharchenko et al, Nature Method, 2014)
##'  \item{monocle}: Census based differential analysis (Qiu et al,
##'  Nature Methods, 2017)
##'  \item{edgeR}: Negative binomial distributions, including empirical
##'  Bayes estimation, exact tests, generalized linear models and
##'  quasi-likelihood tests based differential analysis (McCarthy et al,
##'  Nucleic Acids Research, 2012)
##'  \item{DESingle}: Zero-Inflated Negative Binomial model to estimate
##'  the proportion of real and dropout zeros and to define and detect
##'  the 3 types of DE genes (Miao et al, Bioinformatics, 2018)
##'  \item{MAST}: GLM-framework that treates cellular detection rate as a
##'  covariate (Finak et al, Genome Biology, 2015)
##'  }
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param q_cut Cut-off for q value
#' @param add Whether add genes that are not differentially expressed 
#' but highly expressed for finding the significant pairs later
#' @param top Same as in function rawParse
#' @param stats Same as in function rawParse
#' @return A matrix of the differential expressed genes
#' @importFrom utils installed.packages
#' @import dplyr
#' @import DESeq2
#' @import scde
#' @import monocle
#' @import edgeR
#' @import DEsingle
#' @import MAST
#' @import Biobase
#' @export
DEG<-function(data,method, min_gene_expressed=0, min_valid_cells=0,contrast=NULL,q_cut=0.05,add=TRUE,top=50,stats='mean',...){
  if(method %in% c('SCDE','monocle','DESingle','MAST') && dim(data)[1]>=400){
    print('Warning: It may take a long time. You can go and brew a cup of coffee...')
  }
  if(length(unique(data$cell_type))!=1){
    stop('Error: please compare data with sinlge cell type')
  }
  sub_data<-subset(data,select=-cell_type)
  combination<-combn(unique(sub_data$compare_group),2)
  res=NULL
  if(method=='Wilcox'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,WilcoxTest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='DESeq2'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,DESeq2Test(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='SCDE'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,SCDETest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='monocle'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,MonocleTest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='edgeR'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,edgeRTest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='DESingle'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,DESingleTest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else if(method=='MAST'){
    for(i in ncol(combination)){
      if(is.null(contrast)){
        contrast=c(combination[,i])
      }
      sub_data<-sub_data[sub_data$compare_group %in% combination[,i],]
      res<-rbind(res,MASTTest(sub_data,min_gene_expressed, min_valid_cells, contrast=contrast,...))
    }
  }else{
    stop('Error: method currently not available')
  }
  cell_type<-unique(data$cell_type)
  res<-data.frame(res,cell_type,stringsAsFactors = FALSE)
  res <- res %>% filter(q.value<q_cut)
  if(add){
    parsedData<-rawParse(data %>% select(-compare_group),top=top,stats=stats)%>% select(c(gene,cell_type))%>% mutate(logFC=0.0001,p.value=NA,q.value=NA)
    parsedData<-parsedData %>% anti_join(res,by=c('gene'='gene'))
    res<-rbind(res,parsedData)
  }    
  return(res)
}

#' Differential expression using wilcox
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a Wilcoxon Rank Sum test
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param data_type Type of data. Available options are:
##' \itemize{
##'  \item{'raw data'}: Raw count data without any pre-processing
##'  \item{'log count'}: Normalized and log-transformed data
##' }
#' @param verbose Whether show the progress of computing
#'
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#'
#' @export
WilcoxTest<-function(sub_data,min_gene_expressed, min_valid_cells,
                     contrast=unique(sub_data$compare_group),
                     datatype='raw count', verbose=0){
  counts<-t(subset(sub_data,select=-compare_group))
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  groups<-as.factor(sub_data[valid_cells,'compare_group'])
  counts<-counts[expressed_genes,valid_cells]
  mysapply <- if (verbose) {pbsapply} else {sapply}
  p_val <- mysapply(
    X = 1:nrow(counts),
    FUN = function(i) {
      return(wilcox.test(counts[i,]~groups)$p.value)
    }
  )
  if(datatype=='raw count'){
    logFC<-log(rowSums(counts[,groups==contrast[1]])/ncol(counts[,groups==contrast[1]])/(rowSums(counts[,groups==contrast[2]])/ncol(counts[,groups==contrast[2]])),2)
  }else if(datatype=='log count'){
    logFC<-rowSums(counts[,groups==contrast[1]])/ncol(counts[,groups==contrast[1]])-rowSums(counts[,groups==contrast[2]])/ncol(counts[,groups==contrast[2]])
  }else{
    stop('Error: invalid data type')
  }
  res<-data.frame(rownames(counts),logFC,p_val,stringsAsFactors = FALSE)
  colnames(res)<-c('gene','logFC','p.value')
  res <- res %>% mutate(q.value=p.adjust(res$p.value, method = "BH"))
  return(res)
}

#' Differential expression using DESeq2
#'
#' Identifies differentially expressed genes between two groups of cells using
#' DESeq2
#'
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param test either "Wald" or "LRT", which will then use either
#' Wald significance tests (defined by \code{\link{nbinomWaldTest}}),
#' or the likelihood ratio test on the difference in deviance between a
#' full and reduced model formula (defined by \code{\link{nbinomLRT}})
#' @param fitType either "parametric", "local", or "mean"
#' for the type of fitting of dispersions to the mean intensity.
#' See \code{\link{estimateDispersions}} for description.
#' @param sfType either "ratio", "poscounts", or "iterate"
#' for teh type of size factor estimation. See
#' \code{\link{estimateSizeFactors}} for description.
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients
#' See \code{\link{nbinomWaldTest}} for description of the calculation
#' of the beta prior. In versions \code{>=1.16}, the default is set
#' to \code{FALSE}, and shrunken LFCs are obtained afterwards using
#' \code{\link{lfcShrink}}.
#' @param quiet whether to print messages at each step
#' @param modelMatrixType either "standard" or "expanded", which describe
#' how the model matrix, X of the GLM formula is formed.
#' "standard" is as created by \code{model.matrix} using the
#' design formula. "expanded" includes an indicator variable for each
#' level of factors in addition to an intercept. for more information
#' see the Description of \code{\link{nbinomWaldTest}}.
#' betaPrior must be set to TRUE in order for expanded model matrices
#' to be fit.
#' @param minReplicatesForReplace the minimum number of replicates required
#' in order to use \code{\link{replaceOutliers}} on a
#' sample. If there are samples with so many replicates, the model will
#' be refit after these replacing outliers, flagged by Cook's distance.
#' Set to \code{Inf} in order to never replace outliers.
#' @param useT logical, passed to \code{\link{nbinomWaldTest}}, default is FALSE,
#' where Wald statistics are assumed to follow a standard Normal
#' @param minmu lower bound on the estimated count for fitting gene-wise dispersion
#' and for use with \code{nbinomWaldTest} and \code{nbinomLRT}
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see next argument \code{BPPARAM}.
#' A note on running in parallel using \code{BiocParallel}: it may be
#' advantageous to remove large, unneeded objects from your current
#' R environment before calling \code{DESeq},
#' as it is possible that R's internal garbage collection
#' will copy these files while running on worker nodes.
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#' @import DESeq2
#'
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' This test does not support pre-processed genes. To use this method, please
#' install DESeq2, using the instructions at
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @importFrom utils installed.packages
#'
#' @export
DESeq2Test<-function(sub_data, min_gene_expressed, min_valid_cells, contrast=unique(sub_data$compare_group), test='Wald',
                     fitType='parametric',sfType='ratio',betaPrior=FALSE,quiet=FALSE,modelMatrixType='standard',
                     minReplicatesForReplace=7,useT=FALSE,minmu=0.5,parallel=FALSE,BPPARAM=bpparam()){
  counts<-t(sub_data[,-ncol(sub_data)])
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  counts<-counts[expressed_genes,valid_cells]
  coldata<-as.data.frame(sub_data[valid_cells,'compare_group'])
  rownames(coldata)<-rownames(sub_data[valid_cells,])
  colnames(coldata)<-'compare_group'
  dds<-DESeqDataSetFromMatrix(countData = counts,colData = coldata, design = ~ compare_group)
  dds$condition <- factor(dds$compare_group)
  dds<-DESeq(dds,test=test,
             fitType=fitType,
             sfType=sfType,
             betaPrior=betaPrior,
             quiet=quiet,
             minReplicatesForReplace=minReplicatesForReplace, modelMatrixType=modelMatrixType,
             useT=useT, minmu=minmu,
             parallel=parallel, BPPARAM=BPPARAM)
  res <- results(dds,contrast=c('compare_group',contrast))
  res<-data.frame(rownames(counts),res$log2FoldChange,res$pvalue,res$padj,stringsAsFactors = FALSE)
  colnames(res)<-c('gene','logFC','p.value','q.value')
  return(res)
}

#' Differential expression using scde
#'
#' Identifies differentially expressed genes between two groups of cells using
#' scde
#'
#' @references "Bayesian approach to single-cell differential expression
#' analysis" (Kharchenko PV, Silberstein L, Scadden DT, Nature Methods,
#' doi:10.1038/nmeth.2967)
#' https://github.com/hms-dbmi/scde
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param batch Different batch identifier
#' @param @param n.randomizations number of bootstrap randomizations to be performed
#' @param n.cores number of cores to utilize
#' @param batch.models (optional) separate models for the batch data (if generated
#' using batch-specific group argument). Normally the same models are used.
#' @param return.posteriors whether joint posterior matrices should be returned
#' @param verbose integer verbose level (1 for verbose)
#' @import scde
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' This test does not support pre-processed genes. To use this method, please
#' install scde, using the instructions at
#' http://hms-dbmi.github.io/scde/tutorials.html
#'
#' @importFrom utils installed.packages
#'
#' @export
SCDETest<-function(sub_data,min_gene_expressed,min_valid_cells,contrast=unique(sub_data$compare_group),batch=NULL,
                   n.randomizations=150,n.cores=10,batch.models=models,return.posteriors=FALSE,verbose=1){
  if(is.null(batch)){
    batch<-rep(NaN,nrow(sub_data))
    sub_data<-data.frame(sub_data,batch)
    batch<-NULL
  }
  counts<-t(subset(sub_data,select=-c(compare_group,batch)))
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  counts<-counts[expressed_genes,valid_cells]

  groups<-as.factor(sub_data[valid_cells,'compare_group'])
  o.ifm <- scde.error.models(counts = counts, groups = groups, n.cores = 2, threshold.segmentation = TRUE,
                             min.size.entries = 100,save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde.expression.prior(models=o.ifm,counts=counts)
  if(is.null(batch)){
    ediff <- scde.expression.difference(o.ifm, counts, o.prior, groups=groups,
                                        n.randomizations=n.randomizations, n.cores=n.cores, verbose=verbose)
  }else{
    batch<-as.factor(sub_data[valid_cells,'batch'])
    ediff_batch<- scde.expression.difference(o.ifm, counts, o.prior, groups = groups, batch = batch,
                                             n.randomizations = n.randomizations, n.cores = n.cores,
                                             return.posteriors=return.posteriors,verbose = verbose)
    ediff<- ediff_batch$batch.adjusted
  }
  p.value <- 2*pnorm(abs(ediff$Z),lower.tail=FALSE) # 2-tailed p-value
  q.value <- 2*pnorm(abs(ediff$cZ),lower.tail=FALSE) # Adjusted to control for FDR
  res<-data.frame(rownames(counts),ediff$mle,p.value,q.value,stringsAsFactors = FALSE)
  colnames(res) <- c('gene','logFC','p.value','q.value')
  if(all(levels(groups)!=contrast)){
    res$logFC<--res$logFC
  }
  return(res)
}

#' Differential expression using monocle
#'
#' Identifies differentially expressed genes between two groups of cells using
#' monocle
#'
#' @references Qiu X, Hill A, Packer J, Lin D, Ma Y, Trapnell C (2017).
#' “Single-cell mRNA quantification and differential analysis with Census.”
#' Nature Methods.
#' https://github.com/cole-trapnell-lab/monocle-release
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param batch Different batch identifier
#' @param cores The number of cores to be used while testing each gene
#' for differential expression.
#' @import monocle
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' This test does not support pre-processed genes. To use this method, please
#' install monocle, using the instructions at
#' https://bioconductor.org/packages/release/bioc/html/monocle.html
#'
#' @importFrom utils installed.packages
#'
#' @export
MonocleTest<-function(sub_data,min_gene_expressed,min_valid_cells,contrast=unique(sub_data$compare_group),
                      batch=NULL,cores=4){
  if(is.null(batch)){
    batch<-rep(NaN,nrow(sub_data))
    sub_data<-data.frame(sub_data,batch)
    batch<-NULL
  }
  counts<-t(subset(sub_data,select=-c(compare_group,batch)))
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  if(is.null(batch)){
    pd<-sub_data[,'compare_group',drop=FALSE] %>% mutate(num_genes_expressed=colSums(counts!=0))
  }else{
    pd<-sub_data[,c('compare_group','batch')] %>% mutate(num_genes_expressed=colSums(counts!=0))
  }
  rownames(pd)<-colnames(counts)
  pd<-new('AnnotatedDataFrame',pd)
  fd<-as.data.frame(rownames(counts))
  colnames(fd)<-'gene_short_name'
  fd<-fd %>% mutate(num_cells_expressed=rowSums(counts!=0))
  rownames(fd)<-rownames(counts)
  fd<-new('AnnotatedDataFrame',fd)
  data <- newCellDataSet(as.matrix(counts), phenoData = pd, expressionFamily=negbinomial.size(),featureData=fd)
  expressed_genes <- row.names(subset(fData(data),num_cells_expressed >= min_gene_expressed))
  valid_cells <- row.names(subset(pData(data),num_genes_expressed >= min_valid_cells))
  data <- data[expressed_genes,valid_cells]
  data <- estimateSizeFactors(data)
  data<- estimateDispersions(data)
  gene<-rownames(counts)
  if(is.null(batch)){
    diff_test_res <- differentialGeneTest(data,fullModelFormulaStr = "~compare_group")
    res<-diff_test_res[,c('pval','qval')]
    colnames(res)<-c('p.value','q.value')
    norm_data<- exprs(data)/pData(data)[,'Size_Factor']
    logFC<-log(rowSums(norm_data[,data$compare_group==contrast[1]])/rowSums(norm_data[,data$compare_group==contrast[2]]),2)
    res<-data.frame(logFC,res,stringsAsFactors = FALSE)
  }else{
    data <- reduceDimension(data,residualModelFormulaStr = "~batch",
                            verbose = TRUE)
    diff_test_res <- differentialGeneTest(data,fullModelFormulaStr = "~compare_group")
    res<-diff_test_res[,c('pval','qval')]
    colnames(res)<-c('p.value','q.value')
    norm_data<- exprs(data)/pData(data)[,'Size_Factor']
    logFC<-log(rowSums(norm_data[,data$compare_group==contrast[1]])/rowSums(norm_data[,data$compare_group==contrast[2]]),2)
    res<-data.frame(logFC,res,stringsAsFactors = FALSE)
  }
  res<-data.frame(gene,res,stringsAsFactors = FALSE)
  return(res)
}

#' Differential expression using edgeR
#'
#' Identifies differentially expressed genes between two groups of cells using
#' edgeR
#'
#' @references McCarthy, J. D, Chen, Yunshun, Smyth, K. G (2012). “Differential
#' expression analysis of multifactor RNA-Seq experiments with respect to
#' biological variation.” Nucleic Acids Research, 40(10), 4288-4297.
#' @references Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor
#' package for differential expression analysis of digital gene expression data.”
#' Bioinformatics, 26(1), 139-140.
#' https://github.com/cole-trapnell-lab/monocle-release
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param calcNormMethod normalization method to be used
#' @param trend.method method for estimating dispersion trend. Possible values
#' are "none", "movingave", "loess" and "locfit" (default).
#' @param tagwise logical, should the tagwise dispersions be estimated
#' @param robust logical, should the estimation of prior.df be robustified
#' against outliers
#' @import edgeR
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' This test does not support pre-processed genes. To use this method, please
#' install edgeR, using the instructions at
#' http://bioconductor.org/packages/release/bioc/html/edgeR.html
#'
#' @importFrom utils installed.packages
#'
#' @export
edgeRTest<-function(sub_data,min_gene_expressed,min_valid_cells,contrast=unique(sub_data$compare_group),
                    calcNormMethod='TMM',trend.method='locfit',tagwise=TRUE,robust=FALSE){
  counts<-t(subset(sub_data,select=-compare_group))
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  counts<-counts[expressed_genes,valid_cells]
  groups<-as.factor(sub_data[valid_cells,'compare_group'])
  dgList <- DGEList(counts=counts, genes=rownames(counts),group=groups)
  dgList <- calcNormFactors(dgList, method=calcNormMethod)
  designMat <- model.matrix(~groups)
  dgList <- estimateDisp(dgList, design=designMat,trend.method=trend.method,tagwise=tagwise,robust=robust)
  et <- exactTest(dgList,pair=contrast)
  res<-et$table[,c('logFC','PValue')]
  colnames(res)<-c('logFC','p.value')
  res <- res %>% mutate(q.value=p.adjust(res$p.value, method = "BH"))
  gene<-rownames(counts)
  res<-data.frame(gene,res,stringsAsFactors = FALSE)
  return(res)
}

#' Differential expression using DEsingle
#'
#' Identifies differentially expressed genes between two groups of cells using
#' DEsingle
#'
#' @references Zhun Miao, Ke Deng, Xiaowo Wang, Xuegong Zhang (2018). DEsingle
#' for detecting three types of differential expression in single-cell RNA-seq
#' data. Bioinformatics, bty332. 10.1093/bioinformatics/bty332.
#'
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param parallel If FALSE (default), no parallel computation is used;
#' if TRUE, parallel computation using \code{BiocParallel}, with argument
#' \code{BPPARAM}.
#' @param BPPARAM An optional parameter object passed internally to
#' \code{\link{bplapply}} when \code{parallel=TRUE}. If not specified,
#' \code{\link{bpparam}()} (default) will be used.
#' @import DEsingle
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' This test does not support pre-processed genes. To use this method, please
#' install DEsingle, using the instructions at
#' https://github.com/miaozhun/DEsingle
#'
#' @importFrom utils installed.packages
#'
#' @export
DESingleTest<-function(sub_data,min_gene_expressed,min_valid_cells,contrast=unique(sub_data$compare_group),
                       parallel=FALSE,BPPARAM=bpparam()){
  counts<-t(subset(sub_data,select=-compare_group))
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  counts<-counts[expressed_genes,valid_cells]
  groups<-as.factor(sub_data[valid_cells,'compare_group'])
  results <- DEsingle(counts = counts, group = groups,parallel=parallel,BPPARAM=BPPARAM)
  res<-results[,c('foldChange','pvalue','pvalue.adj.FDR')]
  colnames(res)<-c('foldChange','p.value','q.value')
  res<-res %>% mutate(logFC=log(foldChange,2)) %>% select(-foldChange)
  if(all(levels(groups)!=contrast)){
    res$logFC<--res$logFC
  }
  gene<-rownames(counts)
  res<-data.frame(gene,res,stringsAsFactors = FALSE)
  return(res)
}

#' Differential expression using MAST
#'
#' Identifies differentially expressed genes between two groups of cells using
#' MAST
#'
#' @references  MAST: a flexible statistical framework for assessing transcriptional
#' changes and characterizing heterogeneity in single-cell RNA sequencing
#' data G Finak, A McDavid, M Yajima, J Deng, V Gersuk, AK Shalek, CK Slichter
#' et al Genome biology 16 (1), 278
#'
#' @param sub_data Count data removed cell_type and selected certain two
#' compare_group
#' @param min_gene_expressed Genes expressed in minimum number of cells
#' @param min_valid_cells Minimum number of genes detected in the cell
#' @param contrast String vector specifying the contrast to be
#' tested against the log2-fold-change threshold
#' @param method Character vector, either ’glm’, ’glmer’ or ’bayesglm’
#' @param Silence Common problems with fitting some genes
#' @param check_logged Set FALSE to override sanity checks that try to
#' ensure that the default assay is log-transformed and has at least one
#' exact zero
#' @import MAST
#' @return A matrix of differentially expressed genes and related statistics.
#'
#' @details
#' To use this method, please install MAST, using the instructions at
#' https://github.com/RGLab/MAST
#'
#' @importFrom utils installed.packages
#'
#' @export
MASTTest<-function(sub_data,min_gene_expressed,min_valid_cells,contrast=unique(sub_data$compare_group),
                   method='glm',silent=FALSE,check_logged=TRUE){
  counts<-t(subset(sub_data,select=-compare_group))
  groups<-as.factor(sub_data[valid_cells,'compare_group'])
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
  expressed_genes <- rownames(subset(counts,rowSums(counts) >= min_gene_expressed))
  valid_cells<-(colSums(counts)>=min_valid_cells)
  counts<-counts[expressed_genes,valid_cells]
  counts<-log(counts+1,2)
  cdat<-as.data.frame(data.frame(rownames(sub_data[valid_cells,]),sub_data[valid_cells,'compare_group']),stringsAsFactors = FALSE)
  rownames(cdat)<-rownames(sub_data)
  colnames(cdat)<-c('wellKey','compare_group')
  fdat<-as.data.frame(rownames(counts),stringsAsFactors = FALSE)
  colnames(fdat)<-'primerid'
  sca <- FromMatrix(counts, cdat, fdat,check_logged=check_logged)
  cond<-as.factor(cdat$compare_group)
  zlmCond <- zlm(~compare_group, sca, method=method, silent=silent)
  summaryCond <- summary(zlmCond, logFC=TRUE, doLRT=TRUE)
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[component=='logFC', .(primerid, coef)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  res<-data.frame(fcHurdle[,c('coef','Pr(>Chisq)','fdr')],stringsAsFactors = FALSE)
  colnames(res)<-c('logFC','p.value','q.value')
  rownames(res)<-fcHurdle$primerid
  if(all(levels(groups)!=contrast)){
    res$logFC<--res$logFC
  }
  gene<-rownames(counts)
  res<-data.frame(gene,res,stringsAsFactors = FALSE)
  return(res)
}
