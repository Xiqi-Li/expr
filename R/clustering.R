#' estimate best number of clusters
#'
#' @description estimate the best number of cluster for clustering, using 30 indices from NbClust::NbClust
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data (or matrix).
#' @param scale \code{character()}. Direction for scale data. Default ="row".
#' @param data_preTreat \code{logical()}. If TRUE, data will be PCA transformed; Components whose cumulative sum of variance reach variance_explained_cutoff will be kept. Default FALSE.
#' @param removeVar \code{numeric()}. Remove this percentage of variables based on low variance. Default=0.2.
#' @param variance_explained_cutoff \code{numeric()}. If data_preTreat is TRUE, components in PCA transformed data matrix whose cumulative sum of variance should reach this value. default=0.8.
#' @param method \code{character()}.The clustering method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans". Default = "complete".
#' @param min.nc \code{numeric()}. Minimal number of clusters, between 1 and (number of objects - 1), default=2.
#' @param max.nc \code{numeric()}. Maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
#' @param ... parameters to pass on to NbClust.
#' @return \code{list()} containing \code{Best.nc()} (statistics for various indices) and \code{Best.NumberofCluster()} (best number of clusters).
#' @export
#'
#' @examples
#' result<-estimate_bestNumberofClusters(data=data,scale='none',data_preTreat = F,max.nc = 8)
estimate_bestNumberofClusters<-function(
    data,
    scale=c("row","column","none"),
    data_preTreat=T,
    removeVar=0.2,
    variance_explained_cutoff=0.8,
    method="complete",
    min.nc=2,
    max.nc=15
    ){
  if(ncol(data)<nrow(data)){cat("sample space is less than feature space, reduction of feature space using PCA would be recommended!\n\n")}
  scale=match.arg(scale)
  if(scale=="none"){
    data<-t(data)
  }
  if(scale=="row"){
    data<-scale(t(data))
  }
  if(scale=="column"){
    data=scale(data)
  }

  if(data_preTreat){
    pca_res<-PCAtools::pca(t(data),removeVar = removeVar)
    n_components<-min(sum(cumsum(pca_res$variance)<variance_explained_cutoff)+1,nrow(data))
    data<-scale(pca_res$rotated[,pca_res$components[1:n_components]])
    res<-NbClust(data,  min.nc=min.nc, max.nc=max.nc, method = method, index = "all",...)
    bestNumberofCluster<-as.integer(names(which.max(table(as.numeric(res$Best.nc["Number_clusters",])))))
  }else{
    index_for_NbClust<-c("kl","ch","hartigan",'ccc',"scott","marriot","trcovw","tracew","friedman","rubin",
                         "cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball","ptbiserial","gap",
                         "frey","mcclain","gamma","gplus","tau","dunn","hubert","sdindex",'dindex',"sdbw")
  res<-lapply(index_for_NbClust,function(ind){
      res<-tryCatch(
        {NbClust(data,method=method,index=ind,min.nc = min.nc,max.nc = max.nc)$Best.nc},
        error=function(e){return(c("Number_clusters"=NA,"Value_Index"=NA))})
      # if(is.null(res)){return(c("Number_clusters"=NA,"Value_Index"=NA))};
      return(res)})
    names(res)<-index_for_NbClust
    res<-as.data.frame(res)
    bestNumberofCluster<-as.integer(names(which.max(table(as.numeric(res["Number_clusters",])))))}
  return(list(Best.nc=res,Best.NumberofCluster=bestNumberofCluster))
}

#' map_cluster
#' @description map cluster_id from query clusters to cluster_id from target clusters
#' @param query_clusters integer vector.
#' @param target_clusters integer vector.
#' @param bestnumberofclusters integer vector, containing mapped cluster_id from target clusters
#'
#' @return integer vector, containing mapped cluster_id from target clusters
#'
map_clusters<-function(query_clusters,target_clusters,bestnumberofclusters){
  stopifnot("query clusters should be as same length as target clusters!"=length(query_clusters)==length(target_clusters))
  cluster_table<-table(query_clusters,target_clusters)
  cluster_table<-cluster_table[order(apply(cluster_table,1,max),decreasing=T),]
  assigned_target_cluster_ids<-NULL
  extra_k=1
  for(i in rownames(cluster_table)){
    if(length(assigned_target_cluster_ids)>=length(unique(target_clusters))){
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,as.character(length(unique(target_clusters))+extra_k))
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
      extra_k=extra_k+1
    }
    cluster_table<-cluster_table[,setdiff(colnames(cluster_table),assigned_target_cluster_ids),drop=F]
    if(max(cluster_table[i,])!=0){
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,colnames(cluster_table)[which.max(cluster_table[i,])])
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
    } else {
      j=bestnumberofclusters+1
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,as.character(j))
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
      j=j+1
    }
  }
  return(as.integer(assigned_target_cluster_ids[order(names(assigned_target_cluster_ids))]))
}

# unified cluster call (kmeans,pam,hclust,fuzzy,mclust,apclust,dbscan,mclclust,specc,kkmean,skmeans,nmf,som), not export for user
kmeansCluster<-function(data,k,...){return(stats::kmeans(t(data),centers=k,...)$cluster)}
pamCluster<-function(data,k,...){return(cluster::pam(t(data),k = k,...)$clustering)}
hclustCluster<-function(data,k,...){
  hclust.out <- stats::hclust(dist(t(data)))
  hclust_clusters<-stats::cutree(hclust.out,k=k,...)
  return(hclust_clusters)
}
fuzzyCluster<-function(data,k,...){return(cluster::fanny(t(data),k=k,...)$clustering)}
mclustCluster<-function(data,k,...){return(mclust::Mclust(t(data),G=k)$classification)}
apclustCluster<-function(data,k,...){
  apclust_clusters<-apcluster::apclusterK(s=negDistMat(r=2),x=t(data),K=k)
  apclust_clusters<-rep(1:k,times=sapply(apclust_clusters@clusters,length))[order(unlist(apclust_clusters@clusters))]
  names(apclust_clusters)<-colnames(data)
  return(apclust_clusters)
}
hdbscanCluster<-function(data,k,...){return(dbscan::hdbscan(t(data),minPts = k)$cluster+1)}
mclCluster<-function(data,k,...){
  cor_mat<-cor(data)
  mcl_clusters<-MCL::mcl(cor_mat,addLoops = T,ESM=T,allow1=T)$Cluster
  names(mcl_clusters)<-colnames(data)
  return(mcl_clusters)
}

speccCluster<-function(data,k,...){
  specc_clusters<-kernlab::specc(t(data),centers=k,...)@.Data
  names(specc_clusters)<-colnames(data)
  return(specc_clusters)
}

kkmeansCluster<-function(data,k,...){
  kkmeans_cluster<-kernlab::kkmeans(t(data),centers=k,...)@.Data
  names(kkmeans_cluster)<-colnames(data)
  return(kkmeans_cluster)
}

skmeansCluster<-function(data,k,...){
  skmeans_cluster<-skmeans::skmeans(x=t(data),k=k,...)$cluster
  names(skmeans_cluster)<-colnames(data)
  return(skmeans_cluster)
}

nmfCluster<-function(data,k,...){
  data<-(data-min(data,na.rm=T))/(max(data,na.rm=T)-min(data,na.rm=T))
  fit = NMF::nmf(data, rank = k, ...)
  nmf_cluster<-apply(fit@fit@H, 2, which.max)
  names(nmf_cluster)<-colnames(data)
  return(nmf_cluster)
}

somCluster<-function(data,k,...){
  kr = floor(sqrt(ncol(data)))
  somfit = kohonen::som(t(data), grid = somgrid(kr, kr, "hexagonal"), ...)
  m = somfit$codes[[1]]
  m = m[seq_len(nrow(m)) %in% somfit$unit.classif, ]
  cl = cutree(hclust(dist(m)), k)
  group = numeric(ncol(data))
  for(cl_unique in unique(cl)) {
    ind = as.numeric(gsub("V", "", names(cl)[which(cl == cl_unique)]))
    l = somfit$unit.classif %in% ind
    group[l] = cl_unique
  }
  names(group)<-colnames(data)
  return(group)
}

#' multiCluster
#'
#' @description estimate clustering using multiple clustering algorithms, such as "kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","SKmeans","NMF","SOM"
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data (or matrix).
#' @param scale \code{character()}. Direction for scale data. Default ="row".
#' @param bestnumberofclusters \code{integer(1)}. Optimal number of clusters, if missing, calculated by estimate_bestNumberofClusters.
#' @param data_preTreat \code{logical()}. If TRUE, data will be PCA transformed; Components whose cumulative sum of variance reach variance_explained_cutoff will be kept. Default FALSE.
#' @param removeVar \code{numeric()}. Remove this % of variables based on low variance. Default=0.2.
#' @param variance_explained_cutoff \code{numeric()}. If data_preTreat is TRUE, components in PCA transformed data matrix whose cumulative sum of variance should reach this value. default=0.8.
#' @param nbclust_method 	the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#' @param min.nc minimal number of clusters, between 1 and (number of objects - 1)
#' @param max.nc maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
#' @param cluster_methods \code{vector(mode="character")}. Methods used for clustering and evaluation. default=c("kmeans","pam","hclust","fuzzy","mclust","apcluster","hdbscan","MCL","specc","kkmeans","SKmeans","NMF","SOM").
#' @param base_cluster_method \code{character()}. if missing, the most correlated method will be used.
#' @importFrom stats kmeans hclust cutree
#' @importFrom cluster pam fanny
#' @importFrom mclust Mclust mclustBIC adjustedRandIndex
#' @importFrom apcluster apclusterK negDistMat
#' @importFrom dbscan hdbscan
#' @importFrom MCL mcl
#' @importFrom kernlab specc kkmeans
#' @importFrom skmeans skmeans
#' @import NMF
#' @importFrom kohonen som somgrid
#' @importFrom ggcorrplot ggcorrplot
#'
#' @return list, containing best partition of sample, and co-clusters derived from multiple clustering algorithms.

#' @export
#'
#' @examples multiCluster(data=data,scale='none',data_preTreat = T)
multiCluster<-function(
    data,
    scale=c("row","column","none"),
    bestnumberofclusters,
    data_preTreat=F,
    removeVar=0.2,
    variance_explained_cutoff=0.8,
    nbclust_method="complete",
    min.nc=2,
    max.nc=15,
    cluster_methods=c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som"), #
    base_cluster_method){

  scale=match.arg(scale)
  Best.partition=NULL
  result=list()
  if(missing(bestnumberofclusters)){
    bnc<-estimate_bestNumberofClusters(data,scale=scale,data_preTreat=data_preTreat,removeVar = removeVar, variance_explained_cutoff = variance_explained_cutoff,method = nbclust_method,min.nc=min.nc,max.nc=max.nc)
    bestnumberofclusters<-bnc$Best.NumberofCluster
    Best.partition=bnc$Best.nc$Best.partition
    cat("The best number of clusters is ",bestnumberofclusters,"\n\n")
    result[["Best.partition"]]=Best.partition
  }
  data_o<-data
  if(scale=="row"){
    data<-t(scale(t(data)))
  }
  if(scale=="column"){
    data=t(scale(data))
  }

  kmeans_clusters=pam_clusters=hclust_clusters=
    fuzzy_clusters=mclust_clusters=apcluster_clusters=
    hdbscan_clusters=mcl_clusters=specc_clusters=
    kkmeans_clusters= skmeans_clusters=nmf_clusters=
    som_clusters=NULL


  if("kmeans" %in% cluster_methods){
    kmeans_clusters<-kmeansCluster(data,k = bestnumberofclusters)
  }
  if("pam" %in% cluster_methods){
    pam_clusters<-pamCluster(data,k = bestnumberofclusters)
  }
  if("hclust" %in% cluster_methods){
    hclust_clusters<-hclustCluster(data,k=bestnumberofclusters)
  }
  if("fuzzy" %in% cluster_methods){
    fuzzy_clusters<-fuzzyCluster(data,k=bestnumberofclusters)
  }
  if("mclust" %in% cluster_methods){
    mclust_clusters<-mclustCluster(data,k=bestnumberofclusters)
  }

  if("apclust" %in% cluster_methods){
    apclust_clusters<-apclustCluster(data,k=bestnumberofclusters)
  }
  if("hdbscan" %in% cluster_methods){
    hdbscan_clusters<-hdbscanCluster(data,k = bestnumberofclusters)
  }
  if("MCL" %in% cluster_methods){
    mcl_clusters<-mclCluster(data,k=bestnumberofclusters)
  }
  if("specc" %in% cluster_methods){
    specc_clusters<-speccCluster(data,k=bestnumberofclusters)
  }
  if("kkmeans" %in% cluster_methods){
    kkmeans_clusters<-kkmeansCluster(data,k=bestnumberofclusters)
  }
  if("skmeans" %in% cluster_methods){
    skmeans_clusters<-skmeansCluster(data,k=bestnumberofclusters)
  }
  if("nmf" %in% cluster_methods){
    require(NMF)
    nmf_clusters<-nmfCluster(data_o,k=bestnumberofclusters)
  }
  if("som" %in% cluster_methods){

    som_clusters<-somCluster(data,k=bestnumberofclusters)
  }
  res<-rbind(kmeans_clusters,pam_clusters,hclust_clusters,fuzzy_clusters,mclust_clusters,apclust_clusters,hdbscan_clusters,mcl_clusters,specc_clusters,kkmeans_clusters,skmeans_clusters,nmf_clusters,som_clusters)
  rownames(res)<-cluster_methods
  res=res[apply(res,1,function(x) length(unique(x)))>1,] # XL fix error when methods output same cluster number for all samples (cannot set k for hdbscan) https://hdbscan.readthedocs.io/en/latest/faq.html#q-hdbscan-is-failing-to-separate-the-clusters-i-think-it-should

  if(missing(base_cluster_method)){ # XL base_cluster_method=names(which.max(colSums(cor(t(res)))))
    pairs=data.frame(utils::combn(rownames(res),2))
    rand.sim=data.frame(matrix(nrow = nrow(res),ncol = nrow(res), dimnames = list(rownames(res),rownames(res))))
    for (x in pairs){rand.sim[x[1],x[2]]=rand.sim[x[2],x[1]]=mclust::adjustedRandIndex(res[x[1],],res[x[2],])}
    base_cluster_method=names(which.max(rowSums(rand.sim,na.rm = T)))
    p=ggcorrplot::ggcorrplot(rand.sim,hc.order = TRUE,type = "lower",lab = TRUE,method = c("circle"),title = "Similarity matrix (ARI)")
    # The adjusted Rand Index (ARI) should be interpreted as follows: ARI >= 0.90 excellent recovery; 0.80 =< ARI < 0.90 good recovery; 0.65 =< ARI < 0.80 moderate recovery; ARI < 0.65 poor recovery
    result[["rand.plot"]]=p
    result[[" base_method"]]=base_cluster_method
  }
  res<-as.data.frame(t(as.data.frame(lapply(as.data.frame(t(res)),function(col){map_clusters(col,res[base_cluster_method,],bestnumberofclusters = bestnumberofclusters)[col]}))))
  res_<-as.data.frame(matrix(paste("Cluster_",as.matrix(res),sep=""),nrow=nrow(res)))
  rownames(res_)<-rownames(res)
  colnames(res_)<-colnames(data)
  result[["CoClusters"]]=res_
  return(result)
}

#' consensusCluster
#'
#' @description consensus clustering using bootstrap or(and) noise addition of expression
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data.
#' @param scale \code{character()}. Default ="row", Direction for scale data.
#' @param method \code{character()}. default="bootstrap",  method used to produce subset expression for clustering
#' @param subFeatureSize \code{numeric()}. Default=0.8. The subset fraction of features for bootstrap.
#' @param subSampleSize \code{numeric()}. Default=1. The subset fraction of samples for bootstrap.
#' @param noise \code{numeric()}. Default=1, noise with variance of noise defined unit will be added onto expression
#' @param cutFUN \code{function()}. A vector of clustering functions. Options include:
#' \itemize{
#'  \item "Methods in ClassDiscovery package" - "cutHclust", "cutKmeans", "cutPam", and "cutRepeatedKmeans"
#'  \item "Methods provided by expr" - kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster
#' }
#' @param nTimes \code{integer()}. Default=100. Times of iteration of Clustering
#' @param clusters \code{vector(mode = "integer")}. Default=2, predefined number of clusters.
#' @param verbose \code{logical()}. default \code{FALSE}, if \code{TRUE}, print detailed process information.
#' @param ... params passed on to cutFUN.
#'
#' @return \code{data.frame()}, containing possibility of co-clustering.
#' @export
#'
#' @examples
#' results<-consensusCluster(data,method='combine',clusters=4)
#'
consensusCluster<-function(data,scale=c("row","column","none"),method=c("bootstrap","perturb","combine"),subFeatureSize=0.8,subSampleSize=1,noise=1,cutFUN,nTimes=100,clusters=2,verbose=F,...){
  #	A function that, given a data matrix, returns a vector of cluster assignments. Examples of functions with this behavior are cutHclust, cutKmeans, cutPam, and cutRepeatedKmeans, or kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster
  scale=match.arg(scale)
  if(scale=="row"){
    data<-t(scale(t(data)))
  }
  if(scale=="column"){
    data=t(scale(data))
  }
  method=match.arg(method)
  if(missing(cutFUN)){
    cutFUN=c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster)
  }else{
    tmp=setNames(
      c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster),
      c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som"))
    cutFUN=tmp[cutFUN]
    }
  N <- ncol(data)
  subFeatureSize <- as.integer(nrow(data)*subFeatureSize)
  subSampleSize <- as.integer(ncol(data)*subSampleSize)
  #  stableSampling<-
  stableMatch <- matrix(0, nrow = N, ncol = N,)
  rownames(stableMatch)<-colnames(data)
  colnames(stableMatch)<-colnames(data)
  for (i1 in 1:nTimes) {
    if(method=="bootstrap"){
      tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F)]
    }
    if(method=="perturb"){
      tempData <- data + matrix(rnorm(N * nrow(data), 0, noise),ncol = N)
    }
    if(method=="combine"){
      tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F) ]
      tempData <- tempData + matrix(rnorm(ncol(tempData) * nrow(tempData), 0, noise),ncol = ncol(tempData))
    }

    if (verbose) {
      cat(paste("[", i1, "] ", nrow(tempData), " ", sep = ""))
      if (i1%%10 == 0)
        cat("\n")
    }
    for(k in clusters){
      for(cutfun in cutFUN){
        tempCut <- cutfun(tempData, k=k,...)
        tempMatch <- matrix(0, nrow = N, ncol = N)
        rownames(tempMatch)<-colnames(data)
        colnames(tempMatch)<-colnames(data)

        for (i2 in 1:k) {
          for(samp1 in names(tempCut)){
            for(samp2 in names(tempCut)){
              tempMatch[samp1, samp2] <- ifelse(tempCut[samp1]==tempCut[samp2],1,0)
            }
          }
        }
        stableMatch <- stableMatch + tempMatch
      }
    }
  }
  if (verbose)
    cat("\n")
  result=stableMatch/stableMatch[row(stableMatch)==col(stableMatch)]
  return(result)
}



