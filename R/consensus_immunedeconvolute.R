
# @param expressions, data.frame, no default,gene expression with column as samples, row as gene names , must be provided

# @param method, character vector, default =c("abis","bindea","cibersort","consensustme","danaher","davoli","dcq","deconseq,"epic","mcpcounter","quantiseq","timer","xcell"),various algorithms of deconvolution

#' consensus_immunedeconvolute
#' @description Immune deconvolution of gene expression using multiple algorithms and immune signatures, and them extract consensus immune cellular fraction
#'
#' @param expressions \code{data.frame()}. Gene expression with column as samples, row as gene names.
#' @param methods \code{vector(mode="character")}. default =c("abis","bindea","cibersort","consensustme","danaher","davoli","dcq","deconseq,"epic","mcpcounter","quantiseq","timer","xcell"),various algorithms of deconvolution
#' @param celltype_mapping \code{data.frame()}. A mapping table from method specific cell type onto consensus cell type.
#' @param bindea_reference \code{data.frame()}. Immune signature for various immune cell type, if bindea is in the method list, bindea_reference must be provided
#' @param cibersort_reference \code{data.frame()}. Immune signature for various immune cell type, if cibersort is in the method list, cibersort_reference must be provided.
#' @param consensustme_indication \code{character()}. Parameters passed to consensustme. A character vector with one indication per sample for TIMER. Indication - one of TCGA cancer types, such as "uvm", "hnsc","meso", and etc. if consensustme is in the method list, consensustme_indication must be provided.
#' @param danaher_reference \code{data.frame()}. Immune signature for various immune cell type, if danaher is in the method list, danaher_reference must be provided.
#' @param davoli_reference \code{data.frame()}. Immune signature for various immune cell type, if davoli is in the method list, davoli_reference must be provided.
#' @param deconseq_reference \code{data.frame()}. Immune signature for various immune cell type, if deconseq is in the method list, deconseq_reference must be provided
#' @param timer_indication \code{vector(mode="character")}. A character vector with one indication per sample for TIMER. Indication - one of the TCGA cancer type, such as "UVM", "HNSC", and etc. if timer is in the method list, timer_indication must be provided.
#' @param method_frequency_cutoff \code{numeric()}. Default 2, cell types were found in a number of methods not less than the cutoff.
#' @param backround_noise \code{numeric()}. Default 0.00001, cell fraction less than background_noise was regarded as no expression.
#'
#' @import ADAPTS
#' @importFrom corto ssgsea
#'
#' @return \code{list()}, containing method_frequency table, all deconvolution results from individual methods, and final consensus deconvolution results
#' @export
#'
#' @examples
#' background_noise=0.00001
#' celltype_mapping<-xlsx2df(xlsx_file = "/home/harryjerry/Desktop/MyRScripts/Immunecell_signature.xlsx",sheet = "celltype_mapping",endCol = 3,header = T)
#' bindea_reference<-xlsx2df(xlsx_file = "/home/harryjerry/Desktop/MyRScripts/Immunecell_signature.xlsx",sheet = "Bindea",header = T)
#' cibersort_reference<-read.table("/home/harryjerry/Desktop/cibersort/LM22.txt",header=T,stringsAsFactors = F,check.names = F,sep="\t",row.names = 1)
#' danaher_reference<-xlsx2df(xlsx_file = "/home/harryjerry/Desktop/MyRScripts/Immunecell_signature.xlsx",sheet = "Danaher",header = T)
#' davoli_reference<-xlsx2df(xlsx_file = "/home/harryjerry/Desktop/MyRScripts/Immunecell_signature.xlsx",sheet = "Davoli",header = T)
#' res<-consensus_immunedeconvolute(expressions,celltype_mapping = celltype_mapping,bindea_reference = bindea_reference,cibersort_reference = cibersort_reference,consensustme_indication = "meso",danaher_reference = danaher_reference,davoli_reference = davoli_reference,timer_indication = "meso",backround_noise = background_noise)
#
consensus_immunedeconvolute<-function(expressions,methods=c("abis","bindea","cibersort","consensustme","danaher","davoli","dcq","deconseq","epic","mcpcounter","quantiseq","timer","xcell"),celltype_mapping,bindea_reference,cibersort_reference,consensustme_indication,danaher_reference,davoli_reference,deconseq_reference,timer_indication,method_frequency_cutoff=2,backround_noise=0.00001){
  results<-list()
  immune_abis<-NULL
  immune_bindea<-NULL
  immune_cibersort<-NULL
  immune_consensustme<-NULL
  immune_danaher<-NULL
  immune_davoli<-NULL
  immune_dcq<-NULL
  immune_deconseq<-NULL
  immune_epic<-NULL
  immune_mcpcounter<-NULL
  immune_quantiseq<-NULL
  immune_timer<-NULL
  immune_xcell<-NULL
  if("abis" %in% methods) {
    cat("running abis\n")
    immune_abis<-immunedeconv::deconvolute(expressions,method="abis",arrays=FALSE)
    immune_abis[,-1][immune_abis[,-1]<0]<-0
    immune_abis[,-1]<-t(t(immune_abis[,-1])/colSums(immune_abis[,-1]))
    assertthat::are_equal(colnames(expressions),colnames(immune_abis)[-1])
    results[["abis"]]<-immune_abis
  }

  if("bindea" %in% methods){
    cat("running bindea\n")
    bindea_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Bindea" & (!is.na(celltype_mapping$cell_type)),]
    bindea<-lapply(unique(bindea_reference[["Cell_Type"]]),function(ct){return(bindea_reference[["Symbol"]][bindea_reference[["Cell_Type"]]==ct])})
    names(bindea)<-unique(bindea_reference[["Cell_Type"]])
    immune_bindea<-corto::ssgsea(inmat=expressions,groups=bindea,minsize = 0)
    immune_bindea<-immune_bindea[bindea_mapping$method_cell_type,]
    rownames(immune_bindea)<-bindea_mapping$cell_type
    colnames(immune_bindea)<-colnames(expressions)
    immune_bindea[immune_bindea<0]<-0
    immune_bindea<-t(t(immune_bindea)/colSums(immune_bindea))
    immune_bindea<-data.frame(cell_type=rownames(immune_bindea),immune_bindea,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_bindea)[-1])
    results[["bindea"]]<-immune_bindea
  }

  if("cibersort" %in% methods){
    cat("running cibersort\n")
    cibersort_mapping<-celltype_mapping[celltype_mapping$method_dataset=="cibersort" & (!is.na(celltype_mapping$cell_type)),]
    immune_cibersort<-CIBERSORT(X=expressions,Y=cibersort_reference,perm=10)
    immune_cibersort<-immune_cibersort[cibersort_mapping$method_cell_type,]
    rownames(immune_cibersort)<-cibersort_mapping$cell_type
    immune_cibersort<-immune_cibersort[,colnames(expressions)]
    immune_cibersort<-data.frame(cell_type=rownames(immune_cibersort),immune_cibersort,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_cibersort)[-1])
    results[["cibersort"]]<-immune_cibersort
  }

  if("consensustme" %in% methods){
    cat("running consensustme\n")
    consensustme_mapping<-celltype_mapping[celltype_mapping$method_dataset=="consensus_tme",]
    immune_consensustme<-immunedeconv::deconvolute_consensus_tme(expressions,indications=rep(consensustme_indication,ncol(expressions)),method="ssgsea")
    consensustme_mapping<-consensustme_mapping[match(rownames(immune_consensustme),consensustme_mapping$method_cell_type),]
    rownames(immune_consensustme)<-consensustme_mapping$cell_type
    immune_consensustme[immune_consensustme<0]<-0
    immune_consensustme<-t(t(immune_consensustme)/colSums(immune_consensustme))
    immune_consensustme<-data.frame(cell_type=rownames(immune_consensustme),immune_consensustme,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_consensustme)[-1])
    results[["consensustme"]]<-immune_consensustme
  }

  if("danaher" %in% methods){
    cat("running danaher\n")
    danaher_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Danaher" & (!is.na(celltype_mapping$cell_type)),]
    danaher<-lapply(unique(danaher_reference[["Cell_Type"]]),function(ct){return(danaher_reference[["Symbol"]][danaher_reference[["Cell_Type"]]==ct])})
    names(danaher)<-unique(danaher_reference[["Cell_Type"]])
    immune_danaher<-corto::ssgsea(inmat=expressions,groups=danaher,minsize = 0)
    danaher_mapping<-danaher_mapping[danaher_mapping$method_cell_type %in% rownames(immune_danaher),]
    immune_danaher<-immune_danaher[danaher_mapping$method_cell_type,]
    rownames(immune_danaher)<-danaher_mapping$cell_type
    colnames(immune_danaher)<-colnames(expressions)
    immune_danaher[immune_danaher<0]<-0
    immune_danaher<-t(t(immune_danaher)/colSums(immune_danaher))
    immune_danaher<-data.frame(cell_type=rownames(immune_danaher),immune_danaher,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_danaher)[-1])
    results[["danaher"]]<-immune_danaher
  }

  if("davoli" %in% methods){
    cat("running davoli\n")
    davoli_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Davoli" & (!is.na(celltype_mapping$cell_type)),]
    davoli<-lapply(unique(davoli_reference[["Cell_Type"]]),function(ct){return(davoli_reference[["Symbol"]][davoli_reference[["Cell_Type"]]==ct])})
    names(davoli)<-unique(davoli_reference[["Cell_Type"]])
    immune_davoli<-corto::ssgsea(inmat=expressions,groups=davoli,minsize = 0)
    davoli_mapping<-davoli_mapping[davoli_mapping$method_cell_type %in% rownames(immune_davoli),]
    immune_davoli<-immune_davoli[davoli_mapping$method_cell_type,]
    rownames(immune_davoli)<-davoli_mapping$cell_type
    colnames(immune_davoli)<-colnames(expressions)
    immune_davoli[immune_davoli<0]<-0
    immune_davoli<-t(t(immune_davoli)/colSums(immune_davoli))
    immune_davoli<-data.frame(cell_type=rownames(immune_davoli),immune_davoli,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_davoli)[-1])
    results[["davoli"]]<-immune_davoli
  }

  if("dcq" %in% methods){
    cat("running dcq\n")
    dcq_mapping=celltype_mapping[celltype_mapping$method_dataset=="dcq",]
    immune_dcq<-ADAPTS::estCellPercent.DCQ(refExpr=ADAPTS::LM22,geneExpr=expressions)
    rownames(immune_dcq)<-trimws(gsub("\\.+"," ",rownames(immune_dcq)))
    dcq_mapping<-dcq_mapping[na.omit(match(rownames(immune_dcq),dcq_mapping$method_cell_type)),]
    immune_dcq<-immune_dcq[dcq_mapping$method_cell_type,]
    rownames(immune_dcq)<-dcq_mapping$cell_type
    immune_dcq<-data.frame(cell_type=rownames(immune_dcq),immune_dcq/100,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_dcq)[-1])
    results[["dcq"]]<-immune_dcq
  }
  if("deconseq" %in% methods){
    cat("running deconseq\n")
    deconseq_mapping<-celltype_mapping[celltype_mapping$method_dataset=="deconseq" & (!is.na(celltype_mapping$cell_type)),]
    immune_deconseq<-DeconRNASeq::DeconRNASeq(as.data.frame(expressions),deconseq_reference)
    immune_deconseq<-immune_deconseq[deconseq_mapping$method_cell_type,]
    rownames(immune_deconseq)<-deconseq_mapping$cell_type
    immune_deconseq<-immune_deconseq[,colnames(expressions)]
    immune_deconseq<-data.frame(cell_type=rownames(immune_deconseq),immune_deconseq,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_deconseq)[-1])
    results[["deconseq"]]<-immune_deconseq

  }
  if("epic" %in% methods){
    cat("running epic\n")
    immune_epic<-immunedeconv::deconvolute(expressions,method="epic",tumor=T,scale_mrna = T)
    assertthat::are_equal(colnames(expressions),colnames(immune_epic)[-1])
    results[["epic"]]<-immune_epic
  }

  if("mcpcounter" %in% methods){
    cat("running mcpcounter\n")
    mcp_counter_mapping<-celltype_mapping[celltype_mapping$method_dataset=="mcp_counter" & (!is.na(celltype_mapping$cell_type)),]
    immune_mcpcounter<-immunedeconv::deconvolute_mcp_counter(expressions)
    immune_mcpcounter<-t(t(immune_mcpcounter)/colSums(immune_mcpcounter))
    mcp_counter_mapping<-mcp_counter_mapping[match(rownames(immune_mcpcounter),mcp_counter_mapping$method_cell_type),]
    rownames(immune_mcpcounter)<-mcp_counter_mapping$cell_type
    immune_mcpcounter<-data.frame(cell_type=rownames(immune_mcpcounter),immune_mcpcounter,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_mcpcounter)[-1])
    results[["mcpcounter"]]<-immune_mcpcounter
  }

  if("quantiseq" %in% methods){
    cat("running quantiseq\n")
    immune_quantiseq<-immunedeconv::deconvolute(expressions,method='quantiseq',tumor=T,arrays = F,scale_mrna = T)
    assertthat::are_equal(colnames(expressions),colnames(immune_quantiseq)[-1])
    results[["quantiseq"]]<-immune_quantiseq
  }

  if("timer" %in% methods){
    cat("running timer\n")
    immune_timer<-immunedeconv::deconvolute(expressions,method="timer",indications = rep(timer_indication,ncol(expressions)))
    assertthat::are_equal(colnames(expressions),colnames(immune_timer)[-1])
    results[["timer"]]<-immune_timer
  }

  if("xcell" %in% methods){
    cat("running xcell\n")
    immune_xcell<-immunedeconv::deconvolute(expressions,method="xcell",tumor=T,arrays=F)
    assertthat::are_equal(colnames(expressions),colnames(immune_xcell)[-1])
    results[["xcell"]]<-immune_xcell
  }

  immune_consensus<-rbind(immune_abis,immune_bindea,immune_cibersort,immune_consensustme,immune_danaher,immune_davoli,immune_dcq,immune_epic,immune_mcpcounter,immune_quantiseq,immune_timer,immune_xcell)
  method_frequency<-table(immune_consensus$cell_type)
  results[["method_frequency"]]<-method_frequency
  consensus_cell_type<-setdiff(names(method_frequency[method_frequency>=method_frequency_cutoff]),"uncharacterized cell")
  immune_consensus<-immune_consensus[immune_consensus$cell_type %in% consensus_cell_type,]
  immune_consensus[,-1]<-immune_consensus[,-1]+background_noise
  immune_consensus<-immune_consensus %>% group_by(cell_type) %>% summarise_all(function(x){exp(mean(log(x),na.rm=T))})
  immune_consensus<-immune_consensus[rowSums(immune_consensus[,-1]<background_noise)<(ncol(expressions)/3),]
  results[["consensus"]]<-immune_consensus

  return(results)
}

## tiny modified cibersort code, not exposed to user
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#'
CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#'
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#'
CIBERSORT <- function(X, Y, perm=0, QN=TRUE){

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}
