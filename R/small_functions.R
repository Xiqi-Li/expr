#' Merging multiple data frames
#' @description aggregate all data frames in an ordered list by columns with same headers to avoid duplicated headers
#' @param a list of data frames
#'
#' @return new data frame
#' @export
#'
table_org=function(a=list()){Reduce(function(x,y) full_join(x,y,by=intersect(colnames(x), colnames(y))), a)}

#' reassign NA
#'
#' @param x vector
#' @param y a character or a vector to replace NA in x.
#'
#' @return new vector.
#' @export
#'
reassignNA=function(x,y){
  x[is.na(x)]=y
  return(x)}

#' get heatmap column annotation for complexheatmap
#'
#' @param sampleAttr Sample meta data table.
#' @param data_mx Omic data with features as row and samples as columns.
#' @param nGroupMax Maximum number of unique value of an attribute for this attribute to be included.
#' @param essentialOnly Only essential metadata will be selected for annotation
#' @param tracks Names of attributes selected to be shown as tracks. Will be override if essentialOnly is TRUE. Attributes still obeys nGroupMax rule.
#'
#' @import ComplexHeatmap dplyr
#' @importFrom cols4all c4a
#' @return heatmap annotation object
#' @export
#'
getHeatMapAnnotation=function(sampleAttr,data_mx=NULL,nGroupMax=10,essentialOnly=F,tracks=vector(mode = "character")){
  if(!is.null(data_mx)){
    sampleAttr$colSum=colSums(data_mx)
    sampleAttr$max=apply(data_mx,2,max)
    sampleAttr$median=apply(data_mx,2,median)
    sampleAttr$mad=apply(data_mx,2,mad)
  }
  if(length(tracks)>0){
    sampleAttr=sampleAttr[,tracks]
  }
  feature.names=apply(sampleAttr, c(2), function(x) length(table(x)))
  feature.names=feature.names[sapply(sampleAttr[names(feature.names)],is.numeric)|(feature.names<=nGroupMax&feature.names>1)]
  feature.names=names(feature.names)
  feature.names.surr=tolower(setNames(gsub("[-_ \t\n\r\v\f]","",feature.names),feature.names))
  if(essentialOnly){
    ind=grep("survival|age$|stage|diagn|event|progression|pfs|time|point|protocol|gender|RECIST|response|vital",feature.names.surr,ignore.case = T)
    ind=c(ind,which(feature.names.surr=="os"))
    feature.names=feature.names[ind]
    feature.names.surr=feature.names.surr[ind]
  }
  anno.args=col=list()
  fixed.names=c("overallsurvival|^OS","vitalstatus","progressionfree|^PFS","eventforprogression")
  fixed.name=list()
  for(i in 1:length(fixed.names)){
    fixed.name[[i]]=feature.names[grep(fixed.names[i],feature.names.surr,ignore.case = T)]
  }
  colind=setNames(c("red","green"),c("1","0"))
  if(length(fixed.name[[1]])>0){
    if(length(fixed.name[[2]])>0){
      anno.args2=list(gp=gpar(col=colind[as.character(sampleAttr[[fixed.name[[2]]]])]%>%reassignNA("lightgrey")))
    }else{anno.args2=NULL}
    anno.args[["OS"]]=do.call(
      anno_barplot,
      c(list(x=sampleAttr[[fixed.name[[1]]]],baseline=min(sampleAttr[[fixed.name[[1]]]],na.rm = T)), # when NA exists, baseline needs to be specified.
        anno.args2))
  }
  if(length(fixed.name[[3]])>0){
    if(length(fixed.name[[4]])>0){
      anno.args2=list(gp=gpar(col=colind[as.character(sampleAttr[[fixed.name[[4]]]])]%>%reassignNA("lightgrey")))
    }else{anno.args2=NULL}
    anno.args[["PFS"]]=do.call(
      anno_barplot,
      c(list(x=sampleAttr[[fixed.name[[3]]]],baseline=min(sampleAttr[[fixed.name[[3]]]],na.rm = T)),# when NA exists, baseline needs to be specified.
        anno.args2))
  }
  feature.names=feature.names[!feature.names %in% unlist(fixed.name)]
  if(length(feature.names)>0){
    for (feature.name in feature.names[sapply(sampleAttr[feature.names],is.numeric)]){
      baseline=ifelse(
        length(unique(sampleAttr[[feature.name]]))>.5*nrow(sampleAttr) &
          all(sampleAttr[[feature.name]]>=0,na.rm = T),
        min(sampleAttr[[feature.name]],na.rm = T),0)
      anno.args[[feature.name]]=anno_barplot(sampleAttr[[feature.name]],baseline = baseline)
      feature.names=feature.names[!feature.names%in%feature.name]
    }
    if(length(feature.names)>0){
      for (feature.name in feature.names){
        sampleAttr[[feature.name]][sampleAttr[[feature.name]]==""]=NA
        tmp=factor(sampleAttr[[feature.name]])
        anno.args[[feature.name]]=tmp
        col[[feature.name]]=setNames(cols4all::c4a("dark24",nlevels(tmp)),levels(tmp))
      }
      anno.args[["col"]]=col
      anno.args[["na_col"]] = "lightgrey"
    }
  }
    heatmapAnnotation=do.call(HeatmapAnnotation,anno.args)
    return(heatmapAnnotation)
}

#' MRNsurr De-identify patient ID
#' @description randomly generate surrogate patient ID to deidentify patient information
#' @param x vector of IDs
#' @param seed seed for randomization. Default is NULL.
#'
#' @return vector of surrogate IDs
#' @export
#'
MRNsurr=function(x,seed=NULL){
  x=factor(x)
  if(is.null(seed)){seed=round(runif(1, 1, 1000))}
  set.seed(seed)
  levels(x)=paste("pt", sample(1:nlevels(x)), sep = "")
  return(as.vector(x))
}

#' change column names
#'
#' @param x - Data matrix
#' @param ind - vector of index of colnames needed to be changed
#' @param newNames - new column names
#' @return data matrix with changed column names
#' @export
changeColNames=function(x,ind,newNames){
  colnames(x)[ind]=newNames
  return(x)
}

#' getFill
#'
#' @param x The data matrix or data frame
#' @param byColumn Boolean: whether the function done by column
#'
#' @return fill rate
#' @export
#'
getFill=function(x,byColumn=F){
  if(byColumn){
    fill=apply(x, 2, function(x) sum(is.na(x)))/nrow(x)
    names(fill)=colnames(x)
  }else{
    fill=apply(x, 1, function(x) sum(is.na(x)))/ncol(x)
    names(fill)=rownames(x)
  }

  return(1-fill)
}


#' 2-sample t-test with only group statistics
#'
#' @param m1 the sample means of group 1
#' @param m2 the sample means of group 2
#' @param s1 the sample standard deviations of group 1
#' @param s2 the sample standard deviations of group 2
#' @param n1 sample size of group 1
#' @param n2 sample size of group 2
#' @param m0 the null value for the difference in means to be tested for. Default is 0.
#' @param equal.variance whether or not to assume equal variance. Default is FALSE.
#'
#' @return a named vector consisting "Difference of means", "Std Error", "t", "p-value".
#' @export
#'
t_test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE )
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat)
}
