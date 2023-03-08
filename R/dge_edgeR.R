
#' dge_edgeR
#'
#' @description differential gene expression analysis for raw count expressions using edgeR algorithm
#' @param expressions \code{data.frame()}. Gene expressions table.
#' @param sample_frequency_threshold \code{numeric()}. Genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed. Default 0.5.
#' @param clinic_info \code{data.frame()}. Clinical information table.
#' @param ID_col \code{character()}. Column name for Sample ID (or Patient ID), should be consistent with column names of expressions.
#' @param group_col \code{character()}. Column name for group factor that differential gene expression will be conducted on.
#' @param covariate_col \code{character()}. column name for covariate factor that effects should be removed from model.
#' @param block_col \code{character()}. Column name for block factor if test is conducted on block model, such as paired test.
#' @param contrasts \code{vector(mode="character")} Specific contrasts if preferred, elements should be exactly same as group factor.
#'
#' @return \code{list()}, contains expressions, method, design, contrasts, test, and statistics of edgeR test.
#' @export
#'
#' @examples
#' \dontrun{
#' data(woodman)
#' protein_expressions=2^(log2_expressions)
#' results<-dge_edgeR(
#'   protein_expressions,
#'    clinical_info = clean_RNA_sample_info,
#'    ID_col = "Sample_ID",
#'    group_col = "APOLLO_TIMPOINT",
#'    contrasts = c("TP2-Baseline"))
#' }
#'

dge_edgeR<-function(expressions,sample_frequency_threshold=0.5,clinic_info,ID_col,group_col,covariate_col,block_col,contrasts){
  require(limma)
  require(edgeR)
  require(DESeq2)
  stopifnot("expression colnames wass not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  if(!missing(block_col)){
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    expressions_<-as.data.frame(t(expressions))
    expressions_[["Sample_ID"]]<-rownames(expressions_)
    expressions_[["block_col"]]<-clinic_info[[block_col]]
    sample_block_indices<-expressions_ %>% group_by(block_col) %>% group_indices()
    sample_id_df<-data.frame(Sample_ID=rownames(expressions_),block_col=clinic_info[[block_col]],sample_block_indices=sample_block_indices)
    unique_sample_id_df<-sample_id_df[(1:nrow(expressions_))[!duplicated(sample_id_df$sample_block_indices)],]
    expressions_<-expressions_ %>% group_by(block_col) %>% summarise_all(mean)
    expressions_<-as.data.frame(expressions_)
    rownames(expressions_)<-unique_sample_id_df$Sample_ID[match(expressions_$block_col,unique_sample_id_df$block_col)]
    expressions_<-floor(expressions_[,-1])
    expressions<-data.frame(t(expressions_))
    clinic_info<-clinic_info[match(unique_sample_id_df$Sample_ID,clinic_info[[ID_col]]),]
  }
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
  dge <- DGEList(counts=expressions)
  keep <- filterByExpr(dge,design)
  keep<-keep & (rowSums(dge$counts<=1)<(sample_frequency_threshold*ncol(dge$counts)))
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  dge<-estimateCommonDisp(dge,design,robust=T)
  fit<-glmQLFit(dge,design)

  expressions<-cpm(dge, log=TRUE, prior.count=3)
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    test<-glmQLFTest(fit,contrast=contrast.matrix[,contrast])
    statistics=as.data.frame(topTags(test,n=nrow(expression)))
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_test<-glmQLFTest(fit,contrast=contrast.matrix)
    F_statistics<-as.data.frame(topTags(F_test,n = nrow(expressions)))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:(length(contrasts)+1))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method="edgeR",design=design,contrast.matrix=contrast.matrix,fit=fit,statistics=contrast_statistics))
}
