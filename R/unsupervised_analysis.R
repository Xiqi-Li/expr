
#' unsupervised analysis
#'
#' @description MAP, TSNE, PCA, MDS and unsupervised clustering of expressions
#' @param dat \code{data.frame()}. Gene expressions table.
#' @param labels  \code{vector(mode = "character")}. Sample labels corresponding to all samples (column names) of \code{dat}.
#' @param cols \code{list()}. Colors to control presentation of labels.
#' @param scale \code{character(1)}. Direction to scale \code{dat}. Use "row" to scale by row; "column" to scale by column; "none" for no scaling. Default set as "row".
#' @param run_umap \code{logical(1)}. Whether to run UMAP analysis. Default set as TRUE.
#' @param umap_params \code{list()}. Params passed to UMAP. Default is NULL.
#' @param umap_config \code{list()}. Extra config param passed to UMAP. Default is umap.defaults.
#' @param run_tsne \code{logical(1)}. Whether to run TSNE analysis. Default set as TRUE.
#' @param tsne_params \code{list()}. Params passed to TSNE. Default is NULL.
#' @param run_pca \code{logical(1)}. Whether to run PCA analysis. Default set as TRUE.
#' @param pca_params \code{list()}. Params passed to PCA. Default is NULL.
#' @param run_mds \code{logical(1)}. Whether to run MDS analysis. Default set as TRUE.
#' @param mds_params \code{list()}. Params passed to MDS. Default is NULL.
#' @param run_heatmap \code{logical(1)}. Whether to run clustering analysis with heatmap. Default set as TRUE.
#' @param heatmap_params list, default NULL, params passed to Heatmap
#'
#' @importFrom plotly plot_ly add_markers
#' @import ggplot2 gridExtra ggplotify RColorBrewer umap Rtsne PCAtools plotly pheatmap ComplexHeatmap
#' @return \code{list()} containing results from UMPA, TSNE, PCA, MDS, Cluster(Heatmap) analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' results<-unsupervised_analysis(unsupervised_data,labels=clean_RNA_sample_info$PROTOCOL)
#' # export figure
#' tiff(filename ="unsupervised_cluster.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
#' grid.arrange(arrangeGrob(
#'   results$plots$umap,
#'   results$plots$tsne,
#'   results$plots$pca$`2D`,
#'   results$plots$mds,
#'   results$plots$heatmap,
#'   layout_matrix = matrix(c(1:4,rep(5,4)),ncol = 2,byrow = T)),
#'   cowplot::get_legend(results$plots$umap+theme(legend.position = "top")),
#'   nrow = 2, heights = c(12, 1))
#' dev.off()
#' }
#'
unsupervised_analysis<-function(dat,labels,cols,scale=c("row","column","none"),run_umap=T,umap_params=NULL,umap_config=umap.defaults,run_tsne=T,tsne_params=NULL,run_pca=T,pca_params=NULL,run_mds=T,mds_params=NULL,run_heatmap=T,heatmap_params=NULL){
  library(ggplot2)
  library(gridExtra)
  library(ggplotify)
  library(RColorBrewer)

  scale=scale[1]
  if(scale=='row'){
    dat<-na.omit(t(scale(t(dat))))
  }
  if(scale=="column"){
    dat<-scale(dat)
  }
  if(missing(cols)){
    cols<-palette.colors(palette = "Paired")
    if(!is.null(labels)){
      cols<-sample(cols,length(unique(labels)))
      names(cols)<-unique(labels)
    }
  }

  plots=list()
  ress=list()
  if(run_umap){
    library(umap)
    umap_res<-do.call(umap::umap,c(list(d=t(dat),config = umap_config),umap_params))
    layout<-data.frame(umap_res$layout)
    colnames(layout)<-c("L1","L2")
    umap_plot<-ggplot(data=layout,mapping = aes(x=L1,y=L2))+geom_point()+labs(title="umap_plot")
    if(!is.null(labels)){
      layout<-cbind(layout,Label=labels)
      umap_plot<-ggplot(data=layout, aes(x=L1,y=L2,color=Label))+
        geom_point()+labs(title="umap_plot")+
        scale_color_manual(values = cols)+theme(legend.position = "none")
    }
    plots[["umap"]]=umap_plot
    ress[["umap"]]=umap_res
  }

  if(run_tsne){
    library(Rtsne)
    dat_tsne<-unique(t(dat))
    dat_tsne<-as.matrix(dat_tsne)
    tsne_res<-do.call(Rtsne,c(list(X=dat_tsne, perplexity = min(floor((nrow(dat_tsne) - 1) / 3),100)),tsne_params))
    Y=data.frame(tsne_res$Y)
    colnames(Y)<-c("Y1","Y2")
    tsne_plot<-ggplot(data=Y,mapping = aes(x=Y1,y=Y2))+geom_point()+labs(title="tsne_plot")
    if(!is.null(labels)){
      Y<-cbind(Y,Label=labels)
      tsne_plot<-ggplot(data=Y,mapping = aes(x=Y1,y=Y2,color=Label))+geom_point()+
        labs(title="tsne_plot")+scale_color_manual(values = cols)+theme(legend.position = "none")
    }
    plots[["tsne"]]=tsne_plot
    ress[["tsne"]]=tsne_res
  }

  if(run_pca){
    # library(gg3D) # XL
    library(PCAtools)
    require(plotly)
    #devtools::install_github("Bioconductor/MatrixGenerics")
    pca_res<-do.call(PCAtools::pca,c(list(mat=t(dat)),pca_params)) #XL
    loadings<-data.frame(pca_res$loadings)
    colnames(loadings)<-paste("L",1:ncol(loadings),sep="")
    pca_plot<-ggplot(data=loadings,mapping=aes(x=L1,y=L2))+geom_point()+labs(title="pca_plot")
    if(!is.null(labels)){
      loadings<-cbind(loadings,Label=labels)
      pca_plot<-ggplot(data=loadings,mapping = aes(x=L1,y=L2,color=Label))+geom_point()+
        labs(title="pca_plot")+scale_color_manual(values = cols)+theme(legend.position = "none")
      # pca_3dplot<-ggplot(data=loadings,mapping = aes(x=L1,y=L2,z=L3,col=Label))+geom_point()+axes_3D(theta=30,phi=20)+stat_3D(theta=30,phi=20)+labs(title="pca_plot")+scale_color_manual(values = cols,guide=F)
      pca_3dplot=plotly::plot_ly(loadings,x=~L1, y=~L2, z=~L3, color = ~Label, colors=cols[unique(labels)],
                                 text = ~paste('Sample:', rownames(loadings)),
                                 width = 500, height = 500) %>% plotly::add_markers() #XL
    }
    plots[["pca"]]=list(`2D`=pca_plot,`3D`=pca_3dplot)
    ress[["pca"]]=pca_res
  }


  if(run_mds){
    distances<-stats::dist(t(dat))
    fit_mds<-do.call(cmdscale,c(list(d=distances,eig = T,k=2),mds_params))
    points<-data.frame(fit_mds$points)
    colnames(points)<-paste("P",1:ncol(points),sep="")
    mds_plot<-ggplot(data=points,mapping=aes(x=P1,y=P2))+geom_point()+labs(title="mds_plot")
    if(!is.null(labels)){
      points<-cbind(points,Label=labels)
      mds_plot<-ggplot(data=points,mapping = aes(x=P1,y=P2,color=Label))+geom_point()+
        labs(title="mds_plot")+scale_color_manual(values = cols)+theme(legend.position = "none")
    }
    plots[["mds"]]=mds_plot
    ress[["mds"]]=fit_mds
  }

  if(run_heatmap){
    library(pheatmap)
    library(ComplexHeatmap)
    mads<-apply(dat,1,mad,na.rm=T)
    heatmap_dat<-dat[names(sort(mads,decreasing = T)),]
    heatmap_plot<-as.ggplot(do.call(Heatmap,c(list(matrix=heatmap_dat,show_row_names=F,show_column_names=F),heatmap_params)))+labs(title = "heatmap")
    if(!is.null(labels)){
      top_anno<-HeatmapAnnotation(df=data.frame(label=labels),col=list(label=cols))
      heatmap_plot<-do.call(
        Heatmap,
        c(list(matrix=heatmap_dat,
               name="expressions",
               show_row_names=F,
               show_column_names=F,
               top_annotation=top_anno,
               column_title = "Hierarchical Clustering",
               column_title_gp = gpar(fontface = "bold")),heatmap_params))
    }
    plots[["heatmap"]]=heatmap_plot
  }

  return(list(plots=plots,ress=ress))
}
