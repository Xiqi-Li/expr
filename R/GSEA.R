# pathway_gene_table
# produce a table containing log2 value of fold change of gene for each pathway, internally used by gsea
# not exposed to users
pathway_gene_table<-function(pathways,stats){
  stats<-data.frame(Gene_ID=names(stats),logFC=stats)
  genes_ordered_by_freq<-names(genes_freq<-sort(table(Reduce(c,pathways))[intersect(stats[["Gene_ID"]],unique(Reduce(c,pathways)))],decreasing = T))
  results<-data.frame(Gene=genes_ordered_by_freq)
  for(pathway in names(pathways)){
    genes_of_pathway<-intersect(genes_ordered_by_freq,pathways[[pathway]])
    results[[pathway]]=NA
    results[[pathway]][match(genes_of_pathway,genes_ordered_by_freq)]<-stats[["logFC"]][match(genes_of_pathway,stats[["Gene_ID"]])]
  }
  return(results)
}

# fgseatableplot
# produce a gsea tableplot,including"Pathway", "Gene ranks", "NES", "pval", "padj" columns, internally used by gsea
# not exposed to users
# replacement for fgsea::plotGseaTable
fgseatableplot<-function (pathways, stats, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 3, 1.2, 1.2), render = TRUE)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  pathways <- pathways[sapply(pathways, length) > 0]
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(pn, just = "right", x = unit(0.95, "npc")),
         ggplot() +
           geom_segment(aes(x = p, xend = p, y = 0, yend = statsAdj[p]), size = 0.2) +
           scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         ggplot()+geom_rect(aes(xmin=0,xmax=fgseaRes$NES[fgseaRes$pathway==pn],ymin=-0.6,ymax=0.6),fill=ifelse(sign(fgseaRes$NES[fgseaRes$pathway==pn])==1,"green","red")) +
           scale_x_continuous(limits = c(min(fgseaRes$NES)*1.1, max(max(fgseaRes$NES)*1.1,0)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         #         textGrob(sprintf("%.2f", annotation$NES)),
         textGrob(sprintf("%.1e", annotation$pval)),
         textGrob(sprintf("%.1e", annotation$padj)))
  })
  rankPlot_segment <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  rankPlot_rect <- ggplot() + geom_blank() + scale_x_continuous(limits = c(min(fgseaRes$NES)*1.1, max(fgseaRes$NES)*1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(), rankPlot_segment, rankPlot_rect, nullGrob(), nullGrob()))

  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))

  p <- arrangeGrob(grobs = grobs[grobsToDraw], ncol = sum(colwidths != 0), widths = colwidths[colwidths != 0])
  if (render) {
    grid.draw(p)
  }
  else {
    p
  }
}

# @param enrichment_padj_cutoff, numeric, default=0.05, geneset will be regarded as significant if padj value less or equal to enrichment_padj_cutoff
# @param pathway_gene_table_file, character, no default, file name for saving gene logFC information for each significant pathway
# @param fgseatableplot_file, character, default="fgseatableplot.tiff", file name for saving gsea table plot
# @param enrich_pathways, character, default="up2down2", if TRUE, information about how many up-regulated(or down-regulated) geneset will be plot individually

#' run gsea
#' @description gene set enrichment analysis using fast gsea.
#' @param statistics \code{data.frame()} or \code{matrix()}. No default, statistics from comparisions such as dge_limma,must contain logFold and p.value, must be provided.
#' @param output_dir \code{character()}. Path to the directory where outputs should be saved. Default is current working directory.
#' @param logFC_col \code{character()}. Column for log FoldChange. Default="logFC".
#' @param pval_col \code{character()}. Column for P value. Default="P.Value".
#' @param pathways \code{list()}. No default, Gene sets for gsea, refer to fgsea::gmtPathways.
#' @param enrichment_padj_cutoff \code{numeric()}. Default=0.05. Gene set will be regarded as significant if padj value is less or equal to enrichment_padj_cutoff.
#' @param pathway_gene_table_file \code{character()}. Name of the file where gene logFC information for each significant pathway be saved.
#' @param fgseatableplot_file \code{character()}. default="fgseatableplot.tiff", file name for saving gsea table plot
#' @param enrich_pathways \code{logic()}. Default="up2down2", if TRUE, information about how many up-regulated(or down-regulated) geneset will be plot individually
#' @param enrichment_map \code{logic()}. Default TRUE, if TRUE, perform pathway enrichment.
#' @param pval_cutoff \code{numeric()}. Default=0.05. Genes with p.value less or equal to pval_cutoff are considered significantly changed gene and fed to enrichment_map
#' @param FC_cutoff \code{numeric()}. Default=1.5. Genes with logFC more than or equal to FC_cutoff will be considered as up regulated, and ones less than or equal to negative FC_cutoff as down regulated.
#' @param kappa_cutoff \code{numeric(}). Default=0.3. Pathways with kappa values more than or equal to kappa_cutoff will be considered as close pathways, otherwise, they are separated
#' @param ...  Params passed to fgsea
#'
#' @return \code{list()}, containing gsea table/plot, pathway gene table/plot and enrich map table/plot.
#' @export
#'
#' @examples
#' results<-gsea(statistics,pathways=hallmark_pathway)
#'
gsea<-function(statistics,output_dir=".",logFC_col="logFC",pval_col="P.Value",pathways,enrichment_padj_cutoff=0.05,pathway_gene_table_file,fgseatableplot_file="fgseatableplot.tiff",enrich_pathways="up2down2",enrichment_map=T,pval_cutoff=0.05,FC_cutoff=1.5,kappa_cutoff=0.3,...){
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
  if(class(statistics)=="data.frame"){
    stats<-statistics[[logFC_col]]
    names(stats)<-rownames(statistics)
  }
  fgseaRes <- fgsea(pathways = pathways, stats = stats,minSize=15,maxSize=500,nperm=10000,...)
  fgseaRes<-fgseaRes[order(fgseaRes$NES,decreasing=T)]
  sig_pathways<-pathways[fgseaRes$pathway[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj))]]
  if(length(sig_pathways)==0) {
    cat("no significant signaling pathway was found.")
    return(NULL)
  }
  if(sum(fgseaRes$padj<=enrichment_padj_cutoff  & (!is.na(fgseaRes$padj)))>30){
    sig_pathways<-pathways[fgseaRes$pathway[order(fgseaRes$padj)][1:30]]
  }
  pathway_gene_table=pathway_gene_table(sig_pathways,stats)
  if(missing(pathway_gene_table_file)){
    pathway_gene_table_file="pathway_gene_table.csv"
  }
  write.csv(pathway_gene_table,file.path(output_dir,pathway_gene_table_file))

  tiff(filename = file.path(output_dir,fgseatableplot_file),units="in", width=15, height=5, res=300, compression = 'lzw')
  fgseatableplot(pathways = sig_pathways,stats = stats,fgseaRes = fgseaRes)
  dev.off()

  if(!is.na(enrich_pathways)){
    up_enrich_pathways<-NULL
    down_enrich_pathways<-NULL
    enrich_pathways_<-NULL
    if(grepl("up",enrich_pathways)){
      n_up_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=up)[0-9]+",enrich_pathways,perl=T))))
      up_enrich_pathways<-head(fgseaRes$pathway,n_up_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(grepl("down",enrich_pathways)){
      n_down_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=down)[0-9]+",enrich_pathways,perl=T))))
      down_enrich_pathways<-tail(fgseaRes$pathway,n_down_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(is.null(enrich_pathways_)){enrich_pathways_<-enrich_pathways}
  }
  enrichpathwaysplot<-list()
  for(enrich_pathway in enrich_pathways_){
    enrich_filename= paste(enrich_pathway,".tiff",sep="")
    tiff(filename = file.path(output_dir,enrich_filename),units="in",width=6,height=4,res=300,compression='lzw')
    plot_enrichment<-plotEnrichment(pathways[[enrich_pathway]],stats)+labs(title=enrich_pathway)
    print(plot_enrichment)
    dev.off()
    enrichpathwaysplot[[enrich_pathway]]<-paste(enrich_pathway,".tiff",sep="")
  }
  enrichment_map_data<-NULL
  if(enrichment_map){
    sig_genes<-rownames(statistics)[statistics[[pval_col]]<=pval_cutoff]
    enrichment_map_data<-fgseaRes[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj)),]
    enrichment_map_data<-data.frame(ID=paste("pathway",1:nrow(enrichment_map_data),sep="_"),
                                    Term_Description=enrichment_map_data[["pathway"]],
                                    lowest_p=enrichment_map_data[["padj"]])
    changed_genes<-sapply(enrichment_map_data[["Term_Description"]],function(path){up_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]>=FC_cutoff]);
    up_regulated<-up_regulated[up_regulated %in% sig_genes];
    down_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]<=(-FC_cutoff)]);
    down_regulated<-down_regulated[down_regulated %in% sig_genes];
    return(c("Up_regulated"=paste(up_regulated,collapse = ", "),"Down_regulated"=paste(down_regulated,collapse = ", ")))
    })
    changed_genes<-as.data.frame(t(changed_genes))
    assertthat::are_equal(enrichment_map_data$Term_Description,rownames(changed_genes))
    enrichment_map_data<-cbind(enrichment_map_data,changed_genes)
    enrichment_map_data<-enrichment_map_data[(enrichment_map_data[["Up_regulated"]]!="") | (enrichment_map_data[["Down_regulated"]]!=""),]
    if(nrow(enrichment_map_data)>=2){
      tiff(filename = file.path(output_dir,"enrichment_map.tiff"),units="in", width=8, height=8, res=300, compression = 'lzw')
      enrichment_map_data<-cluster_enriched_terms(enrichment_map_data,use_description = T,plot_clusters_graph = T,plot_dend=T,kappa_threshold=kappa_cutoff)
      dev.off()
    } else{
      cat("rows of enrichment_map_data is less than 2. therefore, no enriched terms plot will be produced.")
    }
  }

  return(list(fgseatable=fgseaRes,sig_pathways=sig_pathways,pathway_gene_table=pathway_gene_table,enrichment_map_data=enrichment_map_data,pathway_gene_table_file=pathway_gene_table_file,fgseatableplot_file=fgseatableplot_file,enrichpathwaysplot=enrichpathwaysplot))
}

