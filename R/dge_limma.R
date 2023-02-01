#' dge_limma
#'
#' @description differential gene expression analysis using limma algorithm.
#' @param expressions \code{data.frame()}. Gene expressions table.
#' @param is_rawcount \code{logical()}. Whether expressions is raw count or not. Default set to FALSE.
#' @param is_logged \code{logical()}. Whether provided expressions was log2 transformed or not. Default set to TRUE.
#' @param normalize \code{logical()}. Normalize provided expressions by library size. Default set to FALSE.
#' @param sample_frequency_threshold \code{numeric()}. Genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed. Default 0.5.
#' @param clinic_info \code{data.frame()}. Clinical information table.
#' @param ID_col \code{character()}. Column name for Sample ID (or Patient ID), should be consistent with column names of expressions.
#' @param group_col \code{character()}. Column name for group factor that differential gene expression will be conducted on.
#' @param covariate_col \code{character()}. column name for covariate factor that effects should be removed from model.
#' @param block_col \code{character()}. Column name for block factor if test is conducted on block model, such as paired test.
#' @param contrasts \code{vector(mode="character")} Specific contrasts if preferred, elements should be exactly same as group factor.
#' @param method \code{character(1)}. Method to conduct test. One of:
#' \itemize{
#'  \item "limma_trend" for non raw count expressions
#'  \item "limma_voom" for raw count expressions.
#' }
#' Default set as "limma_trend".
#'
#' @import DESeq2 statmod
#' @return \code{list()}, contains expressions, method, design, contrasts, test, and statistics of limma test.
#' @export
#'
#' @examples
#' data(woodman)
#' clean_log2_protein_expressions=log2_expressions
#' results<-dge_limma(
#'   clean_log2_protein_expressions,
#'   clinical_info = clean_RNA_sample_info,
#'   ID_col = "Sample_ID",
#'   group_col = "APOLLO_TIMPOINT",
#'   contrasts = c("TP2-Baseline"),
#'   method ="limma_trend")
#'
#'
dge_limma<-function(expressions,is_rawcount=FALSE,is_logged=T,normalize=FALSE,sample_frequency_threshold=0.5,clinic_info,ID_col,group_col,covariate_col,block_col,contrasts,method=c("limma_trend","limma_voom")){
  require(limma)
  require(edgeR)
  require(DESeq2)
  stopifnot("ID column was not found in clinic_info"=ID_col %in% colnames(clinic_info))
  stopifnot("group column was not found in clinic_info"=group_col %in% colnames(clinic_info))
  clinic_info[[group_col]]<-factor(clinic_info[[group_col]])
  if(!missing(covariate_col)) {
    stopifnot("covariate column was not found in clinic_info"=covariate_col %in% colnames(clinic_info))
    clinic_info[[covariate_col]]<-factor(clinic_info[[covariate_col]])
    if(any(levels(clinic_info[[covariate_col]]) %in% levels(clinic_info[[group_col]]))){clinic_info[[covariate_col]]<-factor(paste(clinic_info[[covariate_col]],"_",sep=""))}
    formula<-paste(paste("~",covariate_col,sep="+"),group_col,sep="+")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c("Base",paste(levels(clinic_info[[covariate_col]])[-1],levels(clinic_info[[covariate_col]])[1],sep="-"),paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    if(!missing(contrasts)){
      contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
      colnames(contrast.matrix)<-contrasts
      rownames(contrast.matrix)<-colnames(design)
      for(contrast in contrasts){
        if(contrast %in% rownames(contrast.matrix)){
          contrast.matrix[contrast,contrast]<-1
        } else {
          contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
          contrast.matrix[contrast_,contrast]<-c(1,-1)
        }
      }

    }
  } else {
    formula<-paste("~0+",group_col,sep="")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-levels(clinic_info[[group_col]])
    contrast.matrix<-makeContrasts(contrasts=contrasts,levels = design)
  }
  if(!missing(block_col)) {
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    clinic_info[[block_col]]<-factor(clinic_info[[block_col]])
  }
  stopifnot("expression colnames wass not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  method=match.arg(method)
  if(is_rawcount){
    dge <- DGEList(counts=expressions)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    if(method=="limma_trend"){
      logCPM <- cpm(dge, log=TRUE, prior.count=3)
      is_rawcount<-FALSE
      is_logged<-TRUE
      expressions<-logCPM
    } else if(method=="limma_voom"){
      v <- voom(dge, design, plot=TRUE, normalize="quantile")
      expressions<-v$E
      if(!missing(block_col)){
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        v <- voom(v, design, plot = TRUE, block = clinic_info[[block_col]], correlation = cor$consensus)
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        fit <- lmFit(v, design, block = clinic_info[[block_col]], correlation = cor$consensus)
      } else {
        fit <- lmFit(v, design)
      }
    }
  }
  if(!is_rawcount) {
    if(!is_logged){
      expressions<-log2(expressions+1)
    }
    if(normalize){
      library_size<-apply(expressions,2,sum)
      expressions<-t(t(expressions)/library_size)*mean(library_size)
    }
    expressions<-expressions[(rowSums(expressions==0))<(sample_frequency_threshold*ncol(expressions)),]
    cat("non raw_count data will be analyzed with limma\n")
    if(!missing(block_col)){
      corfit <- duplicateCorrelation(expressions,design,block=clinic_info[[block_col]])
      fit <- lmFit(expressions,design,block=clinic_info[[block_col]],correlation=corfit$consensus)
    } else {
      fit <- lmFit(expressions,design)
    }
  }
  fit<-contrasts.fit(fit,contrast.matrix)
  fit<-eBayes(fit,robust = T)
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    statistics=topTable(fit,coef=contrast,number=nrow(expressions))
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_statistics<-topTable(fit,number = nrow(expressions))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:length(contrasts))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method=method,design=design,contrast.matrix=contrast.matrix,fit=fit,statistics=contrast_statistics))
}
