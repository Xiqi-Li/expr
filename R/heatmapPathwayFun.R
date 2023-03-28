
#' getGseaRes
#' for heatmapPathway app
#' @param input shiny input
#'
#' @return
#' @export
#'
getGseaRes=function(input){
  phenotype=input$phenotype
  gseaRes=gseaRess[[phenotype]]
  gseaRes=gseaRes[gseaRes[[input$pthwySig]]<0.05,]
  return(gseaRes)
}

#' getPthwy
#' for heatmapPathway app
#' @param gseaRes global variable gseaRes
#'
#' @return
#' @export
#'
getPthwy=function(gseaRes){
  pthwy=pathway[gseaRes$pathway]
  pthwy=lapply(pthwy, function(x) x[x %in% rownames(log2_expressions)])
  return(pthwy)
}

#' getGenes
#' for heatmapPathway app
#'
#' @param pthwy global variable pthwy (res of getPthwy)
#' @param input
#'
#' @return
#' @export
#'
getGenes=function(pthwy,input){
  genes=pthwy[[input$pathwyName]]
  genes=genes[genes%in%rownames(log2_expressions)]
  return(genes)
}

#' drawHeatmap
#' for heatmapPathway app
#'
#' @param gseaRes global variable gseaRes
#' @param input shiny inpuy
#' @param z_data_mx z-scored data matrix
#' @param split_group global variable split_group
#' @param dge.p global variable
#' @param DEGp_colFun global variable
#' @param dge.logFC global variable
#' @param logFC_colFun global variable
#' @param dge.logFC.breaks global variable
#' @param genes global variable
#'
#' @return
#' @export
#'
drawHeatmap=function(input,z_data_mx,split_group,dge.p,DEGp_colFun,dge.logFC,logFC_colFun,dge.logFC.breaks,gseaRes,genes){

  heatmap_mRNA<-Heatmap(
    z_data_mx,
    name="level\nz_score",
    top_annotation = getHeatMapAnnotation(
      sampleAttr = clean_RNA_sample_info,
      # data_mx = log2_expressions,
      tracks=input$topTracks[input$topTracks%in%colnames(clean_RNA_sample_info)]
    ),
    left_annotation = rowAnnotation(
      pval = anno_simple(-log10(dge.p),col=DEGp_colFun,pch=ifelse(dge.p<0.05,"*",NA)),
      logFC = anno_simple(
        dge.logFC,col = logFC_colFun,
        pch=ifelse(genes %in% gseaRes$leadingEdge[[which(gseaRes$pathway==input$pathwyName)]],20,NA),
        pt_size = unit(0.1,"npc"),
        pt_gp = gpar(fontface=2))),
    show_row_names = F,column_dend_reorder = T,show_column_names = F,
    clustering_distance_rows =function(x,y) 1 - cor(x, y),
    column_split = split_group)
  # column_names_gp = gpar(fontsize=1))

  lgd_pvalue = Legend(title = input$dgeSig, col_fun = DEGp_colFun, at = -log10(c(min(dge.p),0.05,0.15,1)),
                      labels = c("0", "0.05", "0.15", 1))
  lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
  lgd_logFC = Legend(title = "logFC", col_fun = logFC_colFun, at = dge.logFC.breaks,
                     labels = format(dge.logFC.breaks,digits=2))
  lgd_leadingEdge = Legend(pch = 1, type = "points", labels = "leading edge")
  ptwHeatmap=draw(heatmap_mRNA, annotation_legend_list = list(lgd_pvalue, lgd_sig,lgd_logFC,lgd_leadingEdge))
  return(ptwHeatmap)
}

#' getVariables
#'
#' @param input app input.
#'
#' @return
#' @import igraph
#' @export
#'
getVariables=function(input){
  gseaRes=getGseaRes(input)
  pthwy=getPthwy(gseaRes)
  genes=getGenes(pthwy,input)
  split_group=split_groups[[input$phenotype]]

  # dge
  dge=DGE[[input$phenotype]][genes,]
  dge.p=setNames(dge[,grep(paste(":",input$dgeSig,sep=""),colnames(dge))],rownames(dge))
  dge.logFC=setNames(dge[,grep(":logFC",colnames(dge))],rownames(dge))
  dge.logFC.breaks=c(quantile(dge.logFC[dge.logFC<=0],probs = c(0,0.5)),0,quantile(dge.logFC[dge.logFC>0],probs = c(0.5,1)))

  # color
  DEGp_colFun=circlize::colorRamp2(-log10(c(min(dge.p),0.05,0.15,1)),colors = c("#BB172A","#EE3A48","#FFC9CB","white"))
  logFC_colFun=circlize::colorRamp2(dge.logFC.breaks,colors = c("#559633","#ADDC7C","white","#FD8A7A","#C82912"))

  # partial correlation graph
  data_mx=log2_expressions[genes,]
  ig=expr:::getPCorMxNetwork(data_mx,plotIg=F,method=c("corpcor"),seed=123)

  ig=delete_vertices(ig,which(degree(ig)==0))
  E(ig)$lty=ifelse(E(ig)$weight<=0,3,1)
  E(ig)$weight=abs(E(ig)$weight)
  V(ig)$color=logFC_colFun(dge.logFC[V(ig)$name])
  V(ig)$shape=ifelse(V(ig)$name %in% gseaRes$leadingEdge[[which(gseaRes$pathway==input$pathwyName)]],"circle","square")
  V(ig)$size=-log10(dge.p[V(ig)$name])
  V(ig)$label.cex=0.7
  set.seed(123)
  layout=layout_with_fr(ig)

  # plot heatmap
  z_data_mx=zscoreData(data_mx,data_mx)
  z_data_mx=data_mx[rownames(data_mx),]
  split_group=split_group[colnames(z_data_mx)]
  ptwHeatmap=drawHeatmap(input,z_data_mx,split_group,dge.p,DEGp_colFun,dge.logFC,logFC_colFun,dge.logFC.breaks,gseaRes,genes)

  return(list(ig=ig,layout=layout,ptwHeatmap=ptwHeatmap,ig=ig,dge=dge,genes=genes,dge.logFC=dge.logFC,logFC_colFun=logFC_colFun))
}

#' brush_action
#'
#' @param df selected area
#' @param output app output
#'
#' @return
#' @export
#'
brush_action = function(df, output) {
  row_index = unique(unlist(df$row_index))
  output[["info"]]= renderDT({
    if(!is.null(df)) {
      dge[match(genes[row_index],dge$Gene_Symbol),grep("logFC|P.Value|Gene_Info|band",colnames(dge),value = T)]%>%
        datatable(options = list(
          scrollX=T,scrollY="800px",scrollCollapse=T,
          autoWidth = TRUE,
          columnDefs = list(list(width = '400px', targets = 2)),
          initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'font-size': '12px'});","}"))) %>%
        formatStyle(columns = colnames(.$x$data), `font-size` = '12px')
    }
  })

  output$corPlot=renderUI({
    htmltools::plotTag(
      {
        # plot graph
        if(is.null(df)) {
          plot(ig,layout=layout,
               edge.width=E(ig)$weight*3,
               vertex.frame.color=logFC_colFun(dge.logFC[V(ig)$name]),
               #vertex.frame.color=ifelse("c8c500")
               vertex.size=V(ig)$size*2+1,
               vertex.label.font=1,
               vertex.label.cex = 1,
               vertex.label.dist = 0.4,
               vertex.shape=V(ig)$shape,
               vertex.label.color = "black",
               height = 900,width = 900)
        }else{
          plot(ig,layout=layout,
               edge.width=E(ig)$weight*3,
               vertex.frame.color=ifelse(V(ig)$name %in% genes[row_index],"black",NA),
               vertex.size=V(ig)$size*2+1,
               vertex.label.font=1,
               vertex.label.cex = 1,
               vertex.label.dist = 0.4,
               vertex.shape=V(ig)$shape,
               vertex.label.color = "black")
        }
      },
      "corPlot", width = 900, height = 900, pixelratio=3)
  })
}
