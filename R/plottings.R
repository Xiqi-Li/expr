#' plot box plots
#'
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param sampleAttr - sample attributes; meta data
#' @param MOI - molecule of interest. A subset of data_mx colnames.
#' @param header4x - column name of sampelAttr used for x-axis categories.
#' @param header4color - column name of sampelAttr used for color filling
#' @param header4ID -column name of sampelAttr indicating sample IDs.
#' @param palette color palette
#' @importFrom ggpubr ggboxplot
#' @export
#'
plotBoxPlot=function(data_mx,sampleAttr,MOI,header4x,header4color=NULL,header4ID,palette= "Dark2"){
  df=as.data.frame(t(data_mx[MOI,]))
  df[,header4ID]=rownames(df)
  df[,header4x]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4x]
  if(!is.null(header4color)){
    df[,header4color]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4color]
    df=df[!is.na(df[,header4x])&!is.na(df[,header4color]),]
    p=ggpubr::ggboxplot(df, x = header4x ,
                        y = MOI,
                        combine = TRUE,
                        # merge = "flip",
                        ylab = "Value",
                        add = "jitter",                               # Add jittered points
                        add.params = list(size = 0.5, jitter = 0.1),
                        color = header4color ,
                        palette = "Dark2")
  }else{
    df=df[!is.na(df[,header4x]),]
    p=ggpubr::ggboxplot(df, x = header4x ,
                        y = MOI,
                        combine = TRUE,
                        # merge = "flip",
                        ylab = "Value",
                        add = "jitter",                               # Add jittered points
                        add.params = list(size = 0.5, jitter = 0.1),
                        palette = "Dark2")
  }
  return(p)
}

