## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy=TRUE,
  tidy.opts=list(arrow=TRUE, indent=2,width.cutoff=60),
  message = FALSE,
  warning = FALSE
)

require(expr)
require(dplyr)
require(htmltools)

## ----load-and-clean-data------------------------------------------------------
expressions_file = system.file("extdata/expressions.csv",package = "exprClean")
RNA_sample_info_file = system.file("extdata/RNA_sample_info.csv",package = "exprClean")
protein_coding_ensemble2symbol_file = system.file("extdata/protein_coding_ensemble2symbol.csv",package = "exprClean")

RNA_sample_info<-read.csv(RNA_sample_info_file,header=T,stringsAsFactors = F,check.names = F)
expressions<-read.csv(expressions_file,header=T,stringsAsFactors = F,check.names = F)
protein_coding_ensemble2symbol<-read.csv(protein_coding_ensemble2symbol_file,header=T,stringsAsFactors = F,check.names = F)

# set param
paramInput=data.frame(
  sample_id_col = "Sample_ID",gene_id_col = "gene_id",
  mark_library_size = T,adjusted_by_library_size = T,tolerant_library_size_factor = 1.1,
  ens_id_col = "Gene_ID",symbol_col = "Gene_Symbol",
  mark_purity = T,remove_low_purity_sample = T,low_purity_threshold = 0.3,
  mark_lowexpressgene_pct = T,lowexpression_threshold = 1,remove_sample_with_intense_lowexpressgene = T,outlier_lowexpressgene_pct_factor = 1.5,
  remove_lowexpressgene = T,sample_frequency_threshold = 0.5)

# clean
tmp=expression(prepare_clean_RNA_sample_info_and_protein_expressions(
  RNA_sample_info=RNA_sample_info,sample_id_col = paramInput[["sample_id_col"]],expressions = expressions,gene_id_col = paramInput[["gene_id_col"]],
  mark_library_size = paramInput[["mark_library_size"]],adjusted_by_library_size = paramInput[["adjusted_by_library_size"]],tolerant_library_size_factor = paramInput[["tolerant_library_size_factor"]],
  protein_coding_ensemble2symbol_table = protein_coding_ensemble2symbol,ens_id_col = paramInput[["ens_id_col"]],symbol_col = paramInput[["symbol_col"]],
  mark_purity = paramInput[["mark_purity"]],remove_low_purity_sample = paramInput[["remove_low_purity_sample"]],low_purity_threshold = paramInput[["low_purity_threshold"]],
  mark_lowexpressgene_pct = paramInput[["mark_lowexpressgene_pct"]],lowexpression_threshold = paramInput[["lowexpression_threshold"]],remove_sample_with_intense_lowexpressgene = paramInput[["remove_sample_with_intense_lowexpressgene"]],outlier_lowexpressgene_pct_factor = paramInput[["outlier_lowexpressgene_pct_factor"]],
  remove_lowexpressgene = paramInput[["remove_lowexpressgene"]],sample_frequency_threshold = paramInput[["sample_frequency_threshold"]]))

consoleOutput=capture.output(results<-eval(tmp),type = c("output"))

# data loading summary
SumInfo=t(data.frame(
  paramInput[,grep("factor|threshold",colnames(paramInput))],
  `Number_of_Cleaned_Samples`=dim(results$clean_RNA_sample_info)[1],
  `Number_of_Cleaned_Genes`=dim(results$clean_protein_expressions)[1],
  `Number_of_Deleted_Samples`=dim(results$marked_ori_RNA_sample_info)[1]-dim(results$clean_RNA_sample_info)[1],
  Notes=grep("too many low expression genes",consoleOutput,value = T)))
t1=knitr::kable(SumInfo,caption = "Data loading summary")%>%kableExtra::kable_styling(font_size = "200%")


# cleaned RNA sample info
t2=format(results$clean_RNA_sample_info,nsmall=2) %>%
  DT::datatable(options = list(scrollX=T,scrollY=T),
                caption = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left; color:black; font-size:200% ;',
                  'Sample Info (cleaned)') )

# log2 transformation
clean_RNA_sample_info<-results[["clean_RNA_sample_info"]]
clean_protein_expressions<-results[["clean_protein_expressions"]] # it is not log2 transformed yet
log2_clean_protein_expressions<-log2(clean_protein_expressions+1)



## ----load-and-clean-data-res, echo=FALSE, results = 'asis'--------------------
cat(paste0("#### ","Show data loading summary."), sep = "\n")
t1

cat(paste0("\n#### ","Show cleaned RNA sample info."), sep = "\n")
t2

## ----batch-effect-visualization-----------------------------------------------
rm(list = setdiff(ls(),c(lsf.str(), grep("clean",ls(),value = T))))

# set param
paramInput=data.frame(method="MAD",mad_top_n = 1000,remove_outlier = F)

# prep data
unsupervised_data<-prepare_unsupervised_data(
  expressions=log2_clean_protein_expressions,
  method=paramInput[["method"]],
  mad_top_n = paramInput[["mad_top_n"]],
  remove_outlier = paramInput[["remove_outlier"]])
assertthat::are_equal(colnames(unsupervised_data),clean_RNA_sample_info[["Sample_ID"]])

# run analysis
unsupervised_results=unsupervised_analysis(unsupervised_data,labels = clean_RNA_sample_info$PROTOCOL,cols = c("2014-0938"="green","2016-0861"="blue","LAB08-0380"="red"))

# generate plots
p1=unsupervised_results$plots$pca$`3D`%>%layout(title = "PCA", autosize = F)

p2=ggplotly(
  unsupervised_results$plots$umap +
    theme_bw() + theme(legend.position='none') +
    geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
  tooltip = c("x","y","Label","text")) %>%
    layout(title = "UMAP")

p3=ggplotly(
  unsupervised_results$plots$tsne +
    theme_bw() + theme(legend.position='none') +
    geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
  tooltip = c("x","y","Label","text")) %>%
  layout(title = "tSNE")

p4=ggplotly(
  unsupervised_results$plots$mds +
    theme_bw()+ theme(legend.position='none') +
    geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
  tooltip = c("x","y","Label","text"))%>%
  layout(title = "MDS")

p5=unsupervised_results$plots$heatmap

## ----batch-effect-visualization-plots,class.source = 'fold-hide'--------------
a=data.frame(x=c(0.2,0.8,0.2,0.8),
           y=c(1.0,1.0,0.45,0.45),
           text=c("PCA","UMAP","tSNE","MDS"),
           xref = "paper",
           yref = "paper",
           xanchor = "center",
           yanchor = "bottom",
           showarrow = FALSE 
           )
annotations = apply(a,1, function(x) lapply(as.list(x),type.convert, as.is=TRUE))

# grid.draw(cowplot::get_legend(unsupervised_results$plots$umap+theme(legend.position = "top")))
p=subplot(list(style(p1,showlegend = F),p2,p3,p4),nrows = 2,margin = 0.1) %>%
  layout(title = "Unsupervised clustering before correction",
         annotations = annotations,
         scene = list(domain=list(x=c(0,0.5),y=c(0.5,1))))

browsable(tagList(tags$div(
  style = 'width:100%;height=12;display:block;float:left;',
  p)))

## ----batch-effect-visualization-plots2,results='asis',fig.width=6,fig.height=3,out.width="100%",dpi=130,class.source = 'fold-hide'----
cat("<br>")
p5

## ----batch-effect-determination, results = 'asis'-----------------------------
kmean.df=(unsupervised_results$ress$umap$layout)
# calculate Silhouette score for k from 2 to 2 times number of batches and determine the optimal k for clustering
avg_sil <- function(k) {
  km.res <- kmeans(kmean.df, centers = k, nstart = 25)
  ss <- cluster::silhouette(km.res$cluster, dist(kmean.df))
  mean(ss[, 3])}
k.values = 2:(2*(length(unique(clean_RNA_sample_info$PROTOCOL))))
avg_sil_values <- sapply(k.values, avg_sil)
km.res <- kmeans(kmean.df, centers = k.values[which.max(avg_sil_values)], nstart = 25)

# test if good clustering
good=all(c(km.res$betweenss/km.res$totss>0.5,
      max(avg_sil_values)>0.1))

# test if k-means clustering match with batch definition
clusters=lapply(unique(km.res$cluster), function(x) names(km.res$cluster)[km.res$cluster==x])
batches=lapply(unique(clean_RNA_sample_info$PROTOCOL),
               function(x) clean_RNA_sample_info$Sample_ID[clean_RNA_sample_info$PROTOCOL==x])
hits=c()
i=0
for (ncluster in 1:length(clusters)){
  for (nbatches in 1:length(batches)){
    hit=length(setdiff(clusters[[ncluster]],batches[[nbatches]]))==0 &
      length(setdiff(batches[[nbatches]],clusters[[ncluster]]))==0
    hits=c(hits,hit)
  }
}
matched=any(hits)==T

batch_effect_exist=all(matched,good)
cat(sprintf("Clustering good?: **%s**<br>
            Clustering matched?: **%s**<br>
            Batch effect exists: **%s**",good,matched,batch_effect_exist))


## ----batch-effect-correction--------------------------------------------------
unsupervised_resultss=list()
if(batch_effect_exist){
  #limma
  limma_rbe_log2_protein_expressions<-limma::removeBatchEffect(
    log2_clean_protein_expressions,batch=clean_RNA_sample_info$PROTOCOL)
  uns_limma_rbe_log2_protein_expressions<-prepare_unsupervised_data(limma_rbe_log2_protein_expressions,method="MAD",mad_top_n=1000)
  unsupervised_resultss[["limma"]]=unsupervised_analysis(uns_limma_rbe_log2_protein_expressions,labels=clean_RNA_sample_info$PROTOCOL,cols = c("2014-0938"="green","2016-0861"="blue","LAB08-0380"="red"))
  #combat
  combat_rbe_log2_protein_expressions<-sva::ComBat(log2_clean_protein_expressions,batch=clean_RNA_sample_info$PROTOCOL)
  uns_combat_rbe_log2_protein_expressions<-prepare_unsupervised_data(combat_rbe_log2_protein_expressions,method="MAD",mad_top_n=1000)
  unsupervised_resultss[["combat"]]=unsupervised_analysis(uns_combat_rbe_log2_protein_expressions,labels=clean_RNA_sample_info$PROTOCOL,cols = c("2014-0938"="green","2016-0861"="blue","LAB08-0380"="red"))


# generate plots
p1=p2=p3=p4=p5=list()
for (method in c("limma","combat")){
  unsupervised_results=unsupervised_resultss[[method]]
  p1[[method]]=unsupervised_results$plots$pca$`3D`%>%layout(title = sprintf("PCA %s",method), autosize = F)
  
  p2[[method]]=ggplotly(
  unsupervised_results$plots$umap +
    theme_bw() + theme(legend.position='none') +
    geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
  tooltip = c("x","y","Label","text")) %>%
    layout(title = sprintf("UMAP %s",method))

  p3[[method]]=ggplotly(
    unsupervised_results$plots$tsne +
      theme_bw() + theme(legend.position='none') +
      geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
    tooltip = c("x","y","Label","text")) %>%
    layout(title = sprintf("tSNE %s",method))
  
  p4[[method]]=ggplotly(
    unsupervised_results$plots$mds +
      theme_bw()+ theme(legend.position='none') +
      geom_point(aes(text=sprintf("sample: %s", colnames(unsupervised_data)))),
    tooltip = c("x","y","Label","text"))%>%
    layout(title = sprintf("MDS %s",method))
  
  p5[[method]]=unsupervised_results$plots$heatmap
}
}

## ----limma-plots,class.source = 'fold-hide'-----------------------------------
if (batch_effect_exist) {
a=data.frame(x=c(0.2,0.8,0.2,0.8),
           y=c(1.0,1.0,0.45,0.45),
           text=c("PCA","UMAP","tSNE","MDS"),
           xref = "paper",
           yref = "paper",
           xanchor = "center",
           yanchor = "bottom",
           showarrow = FALSE 
           )
annotations = apply(a,1, function(x) lapply(as.list(x),type.convert, as.is=TRUE))

pp1=p1[["limma"]]
pp2=p2[["limma"]]
pp3=p3[["limma"]]
pp4=p4[["limma"]]
p=subplot(list(style(pp1,showlegend = F),pp2,pp3,pp4),nrows = 2,margin = 0.1) %>%
  layout(title = "Unsupervised clustering after correction",
         annotations = annotations,
         scene = list(domain=list(x=c(0,0.5),y=c(0.5,1))))

browsable(tagList(tags$div(
  style = 'width:100%;height=12;display:block;float:left;',
  p)))
}

## ----limma-plots2,results='asis',fig.width=6,fig.height=3,out.width="100%",dpi=130,class.source = 'fold-hide'----
if (batch_effect_exist) {
cat("<br>")
p5[["limma"]]
}

## ----combat-plots,class.source = 'fold-hide'----------------------------------
if (batch_effect_exist) {
a=data.frame(x=c(0.2,0.8,0.2,0.8),
           y=c(1.0,1.0,0.45,0.45),
           text=c("PCA","UMAP","tSNE","MDS"),
           xref = "paper",
           yref = "paper",
           xanchor = "center",
           yanchor = "bottom",
           showarrow = FALSE 
           )
annotations = apply(a,1, function(x) lapply(as.list(x),type.convert, as.is=TRUE))

pp1=p1[["combat"]]
pp2=p2[["combat"]]
pp3=p3[["combat"]]
pp4=p4[["combat"]]
p=subplot(list(style(pp1,showlegend = F),pp2,pp3,pp4),nrows = 2,margin = 0.1) %>%
  layout(title = "Unsupervised clustering after correction",
         annotations = annotations,
         scene = list(domain=list(x=c(0,0.5),y=c(0.5,1))))

browsable(tagList(tags$div(
  style = 'width:100%;height=12;display:block;float:left;',
  p)))
}

## ----combat-plots2,fig.width=6,fig.height=3,out.width="100%",dpi=130,,class.source = 'fold-hide'----
if (batch_effect_exist) {
cat("<br>")
p5[["combat"]]
}

## ----downstream,out.width="100%",fig.width=7,fig.height=4,dpi=130-------------
if (batch_effect_exist){log2_expressions=limma_rbe_log2_protein_expressions}else{expressions=log2_clean_protein_expressions}

# clear environment
rm(list = setdiff(ls(),c(grep("log2_expressions|info",ls(),value = T))))
unsupervised_expressions<-prepare_unsupervised_data(log2_expressions,method="MAD",mad_top_n=1000)
Heatmap(t(scale(t(unsupervised_expressions))),show_row_names = F,name="expression\nz_score")

