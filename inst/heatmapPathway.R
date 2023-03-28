args = commandArgs(trailingOnly=TRUE)
TypeOfComparison=args[1]

require(ComplexHeatmap)
require(expr)
require(igraph)
require(dplyr)
require(DT)

#load data
if (TypeOfComparison=="delta"){
  load("delta_log2.RData")
  log2_expressions=delta_log2_expressions
  clean_RNA_sample_info=paired_pt_info
  clean_RNA_sample_info$Sample_ID=rownames(clean_RNA_sample_info)
}else{
  load("clean_batchexamined_logRNA.RData")
  if (correctionMetod=="limma"){
    log2_expressions=limma_rbe_log2_protein_expressions
    cat("<br>Using limma batch-corrected expression matrix...<br>")
  }else if(correctionMetod=="combat"){
    log2_expressions=combat_rbe_log2_protein_expressions
    cat("<br>Using ComBat batch-corrected expression matrix...<br>")
  }else{
    log2_expressions=log2_clean_protein_expressions
    cat("<br>Using uncorrected expression matrix...<br>")
  }
}



rm(list = setdiff(ls(),c(lsf.str(),"clean_RNA_sample_info","log2_expressions","TypeOfComparison")))


load(sprintf("gsea_%s.RData",TypeOfComparison))

# if baseline analysis then select only baseline samples
clean_RNA_sample_info=clean_RNA_sample_info[clean_RNA_sample_info$Sample_ID%in%names(split_groups[[1]]),]
log2_expressions=log2_expressions[,clean_RNA_sample_info$Sample_ID]

tracks.d=c("OS","vital_status","PFS","event_for_progression","best_response_perc","time_point","gender","best_response","ORGAN")
tracks.d=tracks.d[tracks.d %in% colnames(clean_RNA_sample_info)]

# pathway gene annotation
load(system.file("extdata/ens2symbol_anno.RData",package = "expr"))

require(shiny)
require(shinyWidgets)
require(InteractiveComplexHeatmap)
ui = fluidPage(
  h3("Zoomed-in view of significant pathway (GSEA)"),
  br(),
  dropdownButton(
    tags$h3("Select pathway."),
    selectInput(inputId = 'phenotype',
                label = 'Type of comparison',
                choices = names(DGE),selected=names(DGE)[1]),
    selectInput(inputId = 'pthwySig',
                label = 'P-value type to determine GSEA significance',
                choices = c("pval","padj"),selected=c("pval","padj")[2]),
    selectInput(inputId = 'dgeSig',
                label = 'P-value type to determine DGE significance',
                choices = c("P.Value","adj.P.Val"),selected=c("P.Value","adj.P.Val")[1]),
    selectInput(inputId = 'pathwyName',
                label = 'Pathway name',
                choices = c(NA)),
    selectInput(inputId = "topTracks",
                label = "Heatmap tracks",
                choices = colnames(clean_RNA_sample_info),
                selected = tracks.d,
                multiple = TRUE),
    circle = TRUE, status = "danger",
    icon = icon("gear"), width = "300px",
    tooltip = tooltipOptions(title = "Select pathway.")
  ),
  tags$h5("Interactive heatmap visualization of gene expression levels."),
  textOutput("text"),
  InteractiveComplexHeatmapOutput(),
  tags$hr(style = "border-top: 1px solid #ccc;"),
  tags$h5("Differential expression of slected genes"),
  DTOutput("info"),
  tags$hr(style = "border-top: 1px solid #ccc;"),
  tags$h5("Ridge regularized partial correlation network of slected genes"),
  tags$p("Nodes in larger size indicate higher DGE significance; \n
         Circle nodes are leading edges of the pathway and square nodes are not; \n
         Colors of the nodes are equavilant to colors indicated in heatmap logFC legend; \n
         Framed nodes are selected in the heatmap." ),
  uiOutput("corPlot",width = "100%")
)

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

server = function(input, output, session) {

  # update pathway names choices
  toListen=reactive({list(input$phenotype,input$pthwySig)})
  observeEvent(toListen(), {
    pthwy=getPthwy(getGseaRes(input))
    updateSelectInput(session, "pathwyName", choices = names(pthwy),selected =names(pthwy[1]))
  })

  observe({
    if(input$pathwyName=="NA"|is.na(input$pathwyName)){
      output$text=renderText("meawmeaw")
    }else{
      tmp=getVariables(input)
      for (a in names(tmp)){assign(a,tmp[[a]],envir = .GlobalEnv)}
      output$text=renderText(ls())
      makeInteractiveComplexHeatmap(input, output, session, ptwHeatmap,brush_action = brush_action)
    }
  })

  session$onSessionEnded(function() {
    stopApp()
  })
}

shinyApp(ui, server)

