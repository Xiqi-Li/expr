
#' prepare_unsupervised_data
#'
#' @description extract qualified subset of data with most variance for unsupervised analysis
#' @param expressions \code{data.frame()}. Gene expressions table, optimally clean and filtered by prepare_clean_RNA_sample_info_and_protein_expressions.
#' @param method \code{character()}. Methods for calculation of variance, and etc, and later used by filtering. CV for coefficient of variance, DQ for quantile of both mean and variance (dual quantile), GUMBLE for gumbel distribution of mad. Default set as "MAD".
#' @param mad_top_n \code{numeric()}. Number of genes to be kept with top variance, if it is -1 keep all genes. Default set as -1.
#' @param mad_top_quantile \code{numeric()}. Percentage of genes to be kept with top mad values if mad_top_n is -1. Default set as 0.75.
#' @param cv_top_n \code{numeric()}. Number of genes to be kept with top cv values, if it is -1 keep all genes. Default set as -1,
#' @param cv_top_mean_quantile \code{numeric()}. Percentage of genes to be kept with top mad values if cv_top_n is -1. Default set as 0.5.
#' @param dq_top_mean_quantile \code{numeric()}. Genes to be kept with mean value quantile above defined dq_top_mean_quantile. Default set as 0.5.
#' @param dq_top_var_quantile \code{numeric()}. Genes to be kept with variance value quantile above defined dq_top_var_quantile. Default set as 0.5.
#' @param gumbel_p_cutoff \code{numeric()}. Genes with p value of MAD gumbel distribution less than gumbel_p_cutoff will be kept. Default set as 0.1.
#' @param remove_outlier \code{logical(1)}. Whether to remove outlier value before calculation of variance, and etc. Default set as FALSE.
#' @import goeveg ordinal
#' @return a \code{data.frame()} of a subset of expressions.
#' @export
#' @examples
#' \dontrun{
#' results<-prepare_unsupervised_data(log2_expression,method="MAD",mad_top_n=1000,remove_outlier=F)
#' }

prepare_unsupervised_data<-function(expressions,method=c("MAD","CV","DQ","GUMBEL"),mad_top_n=-1,mad_top_quantile=0.75,cv_top_n=-1,cv_top_mean_quantile=0.5,dq_top_mean_quantile=0.5,dq_top_var_quantile=0.5,gumbel_p_cutoff=0.1,remove_outlier=F){
  library(goeveg)
  library(ordinal)
  removeoutlier<-function(data){
    quartiles<-quantile(data,probs=c(0.25,0.75),na.rm=T)
    IQR<-IQR(data)
    Lower <- quartiles[1] - 1.5*IQR
    Upper <- quartiles[2] + 1.5*IQR
    data_wo_outlier <- subset(data, data > Lower & data < Upper)
    return(data_wo_outlier)
  }

  method=method[1] #XL1
  #MAD
  if(method=='MAD'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=mad(dat,na.rm=T);return(mad)})
    if(mad_top_n==-1){
      results=expressions[mads>quantile(mads,mad_top_quantile),]
    } else{
      results=expressions[order(mads,decreasing = T)[1:mad_top_n],]
    }
  }
  #GUMBEL
  if(method=='GUMBEL'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=mad(dat,na.rm=T);return(mad)})
    gumbel_mean=mean(mads,na.rm=T)
    gumbel_var<-var(mads,na.rm=T)
    threshold=qgumbel(1-gumbel_p_cutoff,gumbel_mean,gumbel_var)
    results=expressions[mads>=threshold,]
  }
  #CV
  if(method=="CV"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    #      rowMeans(expressions,na.rm=T)
    top_mean_filter=quantile(means,cv_top_mean_quantile,na.rm=T)
    expressions_=expressions[means>=top_mean_filter,]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    cvs=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);cv=cv(dat,na.rm=T);return(cv)})
    if(cv_top_n==-1){
      results=expressions_
    } else {
      results=expressions_[order(cvs,decreasing = T)[1:min(cv_top_n,nrow(expressions_))],]
    }
  }
  #DQ
  if(method=="DQ"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    top_mean_filter=quantile(means,dq_top_mean_quantile,na.rm=T)
    expressions_=expressions[(means>=top_mean_filter & (!is.na(means))),]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    vars=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);var=stats::var(dat,na.rm=T);return(var)})
    top_var_filter=quantile(vars,dq_top_var_quantile)
    results<-expressions_[(!is.na(vars)) & vars>=top_var_filter,]
  }
  results<-results[rowSums(is.na(as.matrix(results)))<(0.25*ncol(results)),]
  return(results)
}
