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

#' switchName
#'
#' @param v character vector
#' @param old elements to be switched in vector v
#' @param new elements to substitute in old in oen to one order
#'
#' @return new vector
switchName=function(v,old,new){
  for(i in 1:length(old)){
    v=replace(v,v==old[i],new[i])
  }
  return(v)
}

#' run pipeline setup
#'
#' @export
#'
#' @examples
#' pplSetUp()
pplSetUp=function(){
  # shiny::runApp(list.dirs(system.file("settingsUI",package = "expr")),launch.browser = TRUE)
  source(system.file("settingsUI/app.R",package = "expr"))
  app=shinyApp(ui, server)
  runApp(app,launch.browser = T)
}

#' run pipeline setup
#'
#' @export
#'
#' @examples
#' runExprPPL()
runExprPPL=function(){
  source(system.file("runPPL.R",package = "expr"))
}
