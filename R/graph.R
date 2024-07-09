getPCorMxNetwork=function(data_mx,plotIg=F,method=c("ppcor","corpcor"),seed=123,pth=0.05){
  method=method[1]
  if(method=="ppcor"){
    pcorr=ppcor::pcor(t(data_mx))
    pcor_c=pcorr$estimate
    pcor_p=pcorr$p.value
    pcor_c.ig=pcor_c
    pcor_c.ig[pcor_p>pth/(ncol(pcor_p)^2)]=0
    diag(pcor_c.ig) = 0
  }else if(method=="corpcor"){
    pcor_c=corpcor::cor2pcor(cor(t(data_mx)))
    dimnames(pcor_c)=list(rownames(data_mx),rownames(data_mx))
    pcor_p=cor_pmat(t(data_mx))
    pcor_c.ig=pcor_c
    pcor_c.ig[pcor_p>pth/(ncol(pcor_p)^2)]=0
    diag(pcor_c.ig) = 0
  }
  set.seed(seed)
  ig.0=ig = graph.adjacency(pcor_c.ig, mode="undirected", weighted=TRUE, add.colnames='name')
  E(ig)$lty=ifelse(E(ig)$weight<=0,3,1)
  E(ig)$weight=abs(E(ig)$weight)
  V(ig)
  set.seed(123)
  if(plotIg){
    print(plot(ig,layout=layout_with_fr(ig),
               edge.width=E(ig)$weight*10,
               vertex.color="Grey"))
  }
  return(ig.0)
}

cor_pmat=function (x, ...) {
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

ParCorRidge=function(data,lambdaMax,fold){
  require(GGMridge)
  data=zscoreData(data,data)
  n=ncol(data)
  p=nrow(data)
  lambda.array <- seq(from = 0.1, to = lambdaMax, length = 0.2) * (n-1)
  pcut.array <- seq(from = 0.01, to = 0.05, by = 0.01)
  tpe <- lambda.pcut.cv(x = t(data),
                        lambda = lambda.array,
                        pcut = pcut.array,
                        fold = fold)
  w.mintpe <- which(tpe == min(tpe), arr.ind = TRUE)
  lambda <- lambda.array[w.mintpe[1L]]
  alpha <- pcut.array[w.mintpe[2L]]
  w.upper <- which(upper.tri(diag(p)))
  partial <- solve(lambda * diag(p) + cor(t(data)))
  partial <- (-scaledMat(x = partial))[w.upper]
  ###############################
  # get p-values from empirical
  # null distribution
  ###############################
  efron.fit <- getEfronp(z = transFisher(x = partial),
                         bins = 50L,
                         maxQ = 13)
  ###############################
  # estimate the edge set of
  # partial correlation graph
  ###############################
  th <- alpha
  w.array <- which(upper.tri(diag(p)),arr.ind=TRUE)
  # wsig <- which(p.adjust(efron.fit$correctp, method="BH") < th )
  wsig <- which(efron.fit$correctp < th )
  E <- w.array[wsig,]
  dim(E)
  ###############################
  # structured estimation
  ###############################
  fit <- structuredEstimate(x = t(data), E = E)
  th.partial <- fit$R
  dimnames(th.partial)=list(rownames(data),rownames(data))
  diag(th.partial)=0
  return(th.partial)
}
# @examples
# pcor=ParCorRidge(data=data_mx,lambdaMax=30,fold=10)


`ridge.net` <-
  function(X,lambda=NULL,plot.it=FALSE,scale=TRUE,k=10,verbose=FALSE){
    if (is.null(lambda)==TRUE){
      ss<-seq(-10,-1,length=1000)
      ss<-10^ss
      n<-nrow(X)
      nn<-n- floor(n/k)
      lambda<-ss*nn*ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    X <- scale(X,scale=scale)  # data needs to be centered and standardized
    B<- matrix(0, nrow=p, ncol=p) # matrix of regression coefficients
    lambda.opt<-rep(0,p)
    cat(paste("Performing local ridge regressions\n"))
    cat(paste("Vertex no "))
    for (i in 1:p) ## visit all nodes
    {
      if ((i/10)==floor(i/10)) {cat(paste(i,"..."))}
      noti <- (1:p)[-i]
      yi <- X[ , i]       # response
      Xi <- X[ , noti]    # predicted by all other nodes with i missing
      r.cv= ridge.cv(Xi,yi,lambda=lambda,scale=scale,plot.it=plot.it,k=k)
      B[i,-i]=r.cv$coefficients
      lambda.opt[i]=r.cv$lambda.opt
    }

    pcor<- Beta2parcor(B,verbose=verbose)  # compute matrix of partial correlations
    return(list(pcor=pcor,lambda.opt=lambda.opt))
  }

`Beta2parcor` <-
  function(Beta,verbose=FALSE){
    Dummy=Beta*t(Beta)
    if (verbose==TRUE){
      cat("\nNumber of pairwise regression coefficients with conflicting signs:", (sum((Dummy) < 0))/2, "\n")
      cat("Number of partial correlation coefficients greater than 1 in absolute value:",(sum((Dummy) >1))/2 ,
          "\n\n")
    }
    Dummy[Dummy<0]=0 # if a product is <0 the partial correlation coefficient is set to 0
    Dummy[Dummy>1]=1 # partial correlation coefficients should be in the range of [-1,1]
    P=sign(Beta)*sqrt(Dummy)
    diag(P)=rep(1,ncol(P))
    return(P)
  }

ridge.cv<-
  function(X,y,lambda=NULL,scale=TRUE,k=10,plot.it=FALSE){
    if (is.vector(X)==TRUE){
      X<-matrix(X,ncol=1)
    }
    if (is.null(lambda)==TRUE){
      ss<-seq(-10,-1,length=1000)
      ss<-10^ss
      n<-nrow(X)
      nn<-n- floor(n/k)
      lambda<-ss*nn*ncol(X)
    }
    cv<-rep(0,length(lambda))
    n<-nrow(X)
    all.folds <- split(sample(1:n), rep(1:k,length=n))
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      Xtrain=X[-omit,,drop=FALSE]
      ytrain=y[-omit]
      Xtest=X[omit,,drop=FALSE]
      ytest=y[omit]
      if ((is.vector(X)==TRUE)| (ncol(X)==1)){
        xtrain<-as.vector(Xtrain)
        coef.ll<-lm.ridge.univariate(xtrain,ytrain,lambda=lambda,scale=scale)
      }else{
        ll<-lm.ridge(ytrain~Xtrain,scale=scale,lambda=lambda)
        coef.ll<-coef(ll)
      }
      res<-matrix(length(ytest),length(lambda))
      pred<-t(matrix(coef.ll[,1],nrow=length(lambda),ncol=length(ytest))) + Xtest%*%t(coef.ll[,-1])
      res<-pred-matrix(ytest,nrow=length(ytest),ncol=length(lambda))
      cv<-cv+apply(res^2,2,sum)

    }
    cv<-cv/n
    lambda.opt<-lambda[which.min(cv)]
    if (plot.it==TRUE){
      plot(lambda,cv,type="l")
    }
    if ((is.vector(X)==TRUE)| (ncol(X)==1)){
      x<-as.vector(X)
      coefficients<- as.vector(lm.ridge.univariate(x,y,scale=scale,lambda=lambda.opt))

    }else{
      rr<-lm.ridge(y~X,scale=scale,lambda=lambda.opt)
      coefficients<-coef(rr)
    }
    intercept<-coefficients[1]
    coefficients<-coefficients[-1]
    return(list(intercept=intercept,coefficients=coefficients,lambda.opt=lambda.opt))
  }

lm.ridge.univariate<-function(x,y,lambda=0,scale=TRUE){
  r<-length(lambda)
  xs<-scale(x,scale=scale)
  ys<-scale(y,scale=scale)
  sx<-1
  sy<-1
  if (scale==TRUE){
    sx<-sd(x)
    sy<-sd(y)
  }
  b<-sum(xs*ys)/(sum(xs^2) +lambda)
  b<-b*sy/sx
  inter<-mean(y) - b*mean(x)
  coefficients<-cbind(inter,b)
  return(coefficients)
}

plotUmap=function(data_mx,sampleAttr=data.frame(),label,nGroupMax=10,axistextsize=10,dotsize=2,legendtextsize=9,textlabelsize=0.5,n_neighbors_rate=0.25){
  if(nrow(data_mx)>10000){
    message("*** selecting top 1000 variable features ***")
    data_mx=data_mx[
      names(head(sort(apply(data_mx[sample(1:nrow(data_mx),nrow(data_mx)/10),], 1, var),decreasing=T),1000)),]
  }
  n_neighbors=round(n_neighbors_rate*as.numeric(ncol(data_mx)))
  min_dist=0.1
  config=umap.defaults
  config$n_neighbors=n_neighbors
  config$min_dist=min_dist

  message("***UMAP scoring ***")
  UMAP=umap(t(as.matrix(data_mx)),config = config)
  scores = data.frame(UMAP$layout)

  message("***UMAP wrapper function***")
  umap.plots=list()
  sampleAttr$colSum=colSums(data_mx)
  sampleAttr$max=apply(data_mx,2,max)
  sampleAttr$median=apply(data_mx,2,median)
  feature.names=apply(sampleAttr, c(2), function(x) length(table(x)))
  feature.names=feature.names[sapply(sampleAttr[names(feature.names)],class)== "numeric"|(feature.names<=nGroupMax&feature.names>1)]
  feature.names=names(feature.names)
  for (feature.name in feature.names){

    if(class(sampleAttr[,feature.name])== "numeric") {
      feature = sampleAttr[,feature.name]
    }else{
      feature = as.factor(sampleAttr[,feature.name])
    }
    scores$feature=feature
    umap.plot=
      ggplot(data = scores, aes(x = X1, y = X2, label=sampleAttr[,label])) +
      geom_point(aes(colour = feature), size = dotsize) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize,colour = "black"),
            axis.text.x = element_text(size = axistextsize,colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize),
            legend.title = element_text(size = legendtextsize),
            legend.text = element_text(size = legendtextsize)) +
      labs(colour = feature.name) +
      geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
    if(class(sampleAttr[,feature.name])== "numeric") {
      umap.plot=umap.plot +
        scale_colour_gradient(low = "blue", high = "red")
    }
    umap.plots[[feature.name]]=umap.plot
  }
  return(umap.plots)
}
