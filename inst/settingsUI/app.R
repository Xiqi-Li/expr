library(shiny)
require(DT)
require(plotly)
require(shinyWidgets)
require(umap)

# Define UI
ui <- navbarPage("expr-SetUp: inputs and parameters",
    sidebarLayout(
      sidebarPanel(
        fileInput("data_file", "Upload RData File",accept = c(".RData")),
        prettyCheckbox("logged","log-transformed data?", value = TRUE,
                       status = "danger",shape = "curve",outline = TRUE),

        tags$hr(style = "border-top: 1px solid #ccc;"),
        tags$h5("Select the metadata headers corresponding to below phenotypes. If not available, select NA.", style="color:#505050"),
        pickerInput("Sample_ID", "Sample ID", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("MRN", "Patient ID/MRN", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("headerRm", "Columns that should not be shown or used", choices = c(NA), multiple = TRUE,
                    options = list(`live-search` = TRUE)),
        pickerInput("OS", "Overall Survival", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("vital_status", "Vital Status", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("PFS", "Progression-free survival", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("event_for_progression", "Event for Progression", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("time_point", "Collection time point", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("best_response", "Best_Response_(RECIST)", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("best_response_perc", "Best_Response_(%)", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        pickerInput("batch", "Used for batch correction - reference umaps", choices = c(NA),
                    options = list(`live-search` = TRUE)),
        switchInput("forceCorrection","Force batch correction", value = FALSE, size = "mini"),

        tags$hr(style = "border-top: 1px solid #ccc;"),
        tags$h5("Set parameters for dataset cleanup.", style="color:#505050"),
        tags$p("Add library size column to metadata?"),
        switchInput(inputId = "mark_library_size", value = TRUE, size = "mini"),
        sliderInput("tolerant_library_size_factor","The ratio of maximal to minimal library size that riggers library size normalization.",
                    min=1,max=2,value = 1.1,step = 0.05),
        tags$p("Add tumor purity column to metadata?"),
        switchInput(inputId = "mark_purity", value = TRUE, size = "mini"),
        sliderInput("low_purity_threshold","Exclude samples with tumor purity lower than: ",
                    min=0,max=1,value = 0.3,step = 0.05),
        tags$p("Add sparsely expressed gene percentage column to metadata?"),
        switchInput(inputId = "mark_lowexpressgene_pct", value = TRUE, size = "mini"),
        sliderInput("lowexpression_threshold","A gene is defined as sparsely expressed if below: (default - 1/4 quantile, max - max count)",
                    min=0,max=50,value = 1,step = 0.05),
        tags$p("Remove sample with high percentage of sparsely expressed gene?"),
        switchInput(inputId = "remove_sample_with_intense_lowexpressgene",value = TRUE,size = "mini"),
        sliderInput("outlier_lowexpressgene_pct_factor","Exclude samples with low expressed gene percentage over cohort average times:",
                    min=1,max=3,value = 1.5,step = 0.05),
        tags$p("Remove low expressed genes?"),
        switchInput(inputId = "remove_lowexpressgene",value = TRUE,size = "mini"),
        sliderInput("sample_frequency_threshold","Remove genes expressed sparsely in more than x (rate) of all samples:",
                    min=0,max=1,value = 0.5,step = 0.05),

        tags$hr(style = "border-top: 1px solid #ccc;"),
        actionButton("save_settings", "Save Settings")
        ),
      mainPanel(
        tags$h3("UMAP embedding of Samples: Select batch factor."),
        textOutput('test'),
        br(),
        dropdownButton(
          tags$h3("Sample Metadata"),
          selectInput(inputId = 'feature',
                      label = 'Color guide',
                      choices = NA),
          sliderInput("nGroupMax",
                      "Maximum number of categories for categorical meta data",
                      value = 10,min = 2, max = 20),
          circle = TRUE, status = "danger",
          icon = icon("gear"), width = "300px",

          tooltip = tooltipOptions(title = "Select color guide!")
        ),
        plotlyOutput('umap',height = "400px"),
        tags$hr(style = "border-top: 1px solid #ccc;"),
        tags$h3("Sample Metadata table"),
        DTOutput('sampleAttrDT'),
        )
      )
    )


# Define server
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  #define functions
  getUmapScore=function(data_mx){
    if(nrow(data_mx)>10000){
      message("*** selecting top 1000 variable features ***")
      data_mx=data_mx[
        names(head(sort(apply(data_mx[sample(1:nrow(data_mx),nrow(data_mx)/10),], 1, var),decreasing=T),1000)),]
    }
    n_neighbors=round(0.25*as.numeric(ncol(data_mx)))
    min_dist=0.1
    config=umap.defaults
    config$n_neighbors=n_neighbors
    config$min_dist=min_dist

    message("***UMAP scoring ***")
    UMAP=umap(t(as.matrix(data_mx)),config = config)
    scores = data.frame(UMAP$layout)
    return(scores)
  }

  getFeatureNames=function(sampleAttr,nGroupMax=10){
    feature.names=apply(sampleAttr, c(2), function(x) length(table(x)))
    feature.names=feature.names[sapply(sampleAttr[names(feature.names)],class)== "numeric"|(feature.names<=nGroupMax&feature.names>1)]
    feature.names=names(feature.names)
    return(feature.names)
  }

  plotlyUmaps=function(scores,sampleAttr,feature.name, label,axistextsize=10,dotsize=2,legendtextsize=9,textlabelsize=0.5){
    message("***UMAP wrapper function***")
    # umap.plots=list()

    # for (feature.name in feature.names){

    if(class(sampleAttr[,feature.name])== "numeric") {
      feature = sampleAttr[,feature.name]
    }else{
      feature = as.factor(sampleAttr[,feature.name])
    }
    scores$feature=feature
    umap.plot=
      ggplot(data = scores, aes(x = X1, y = X2, label=sampleAttr[,label])) +
      geom_point(aes(colour = feature,text=sampleAttr[,label]), size = dotsize) +
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
    umap.plots=ggplotly(umap.plot,tooltip = c("x","y","Label","text"))#[[feature.name]]
    # }
    return(umap.plots)
  }

  getSampleAttr.tmp=function(sampleAttr,data_mx){
    cbind(sampleAttr,
          colSum=colSums(data_mx),
          max=apply(data_mx,2,max),
          median=apply(data_mx,2,median))
  }

  headers=c("Sample_ID","MRN","headerRm","batch",
            "OS","vital_status","PFS","event_for_progression",
            "time_point",
            "best_response","best_response_perc")

  cleanUpSettings=c("mark_library_size","tolerant_library_size_factor",
                    "mark_purity","low_purity_threshold",
                    "mark_lowexpressgene_pct","lowexpression_threshold",
                    "remove_sample_with_intense_lowexpressgene","outlier_lowexpressgene_pct_factor",
                    "remove_lowexpressgene","sample_frequency_threshold")
  # Load data
  observeEvent(input$data_file, {
    load(input$data_file$datapath)
    # define class of ID attributes
    coerce_ind=names(which(sapply(sampleAttr[grep("ID|MRN",colnames(sampleAttr),ignore.case = T)],is.numeric)))
    if(length(coerce_ind)>0){
      for(ind in coerce_ind){
        sampleAttr[,coerce_ind]=as.character(sampleAttr[,coerce_ind])
      }}

    choices=colnames(sampleAttr)

    observeEvent(input$logged,{
      if(input$logged){
        expressions=as.data.frame(2^gene_expressions-1)
      }else{
        expressions=as.data.frame(gene_expressions)
      }

    #   Generate UMAP distribution
      scores = getUmapScore(data_mx = log2(expressions+1))

      #Update cleanup settings options
      # tmp=na.omit(unlist(expressions))
      # updateSliderInput(session, "lowexpression_threshold",max = quantile(tmp,prob=0.75))

      # Update select input options
      lapply(headers, function(x) updatePickerInput(session, x, choices = c(choices,"NA"),selected ="NA"))


    # update metadata
    sampleAttr_updated <- reactive({
      sampleAttr_copy <- sampleAttr
      for(header in headers[headers!="headerRm"]){
        eval(parse(
          text=sprintf("colnames(sampleAttr_copy)[colnames(sampleAttr_copy) == input$%s] <- '%s'",header,header)
          ))
      }
      if(length(na.omit(input$headerRm))>0) {sampleAttr_copy=sampleAttr_copy[,!colnames(sampleAttr_copy) %in% input$headerRm]}
      sampleAttr_copy
    })

    # DT updated metadata
    output$sampleAttrDT=renderDT({
      format(sampleAttr_updated(),nsmall=2) %>%
        DT::datatable(options = list(scrollX=T,scrollY="750px",scrollCollapse=T),rownames = F)
    })

    #update plot inputs
    toListenHeader=eval(parse(
      text = sprintf("reactive({list(%s)})",paste(paste("input$",headers,sep=""),collapse = ","))
      ))

    observeEvent(toListenHeader(),{
    updateSelectInput(session,"feature",
                      choices = getFeatureNames(sampleAttr_updated(),nGroupMax = input$nGroupMax),
                      selected = getFeatureNames(sampleAttr_updated())[1])
      })

    # Generate plot
    toListenUmap=eval(parse(
      text = sprintf("reactive({list(%s)})",
                     paste(paste("input$",c("datafile","feature","Sample_ID","nGroupMax"),sep=""),collapse = ","))
    ))

    observeEvent(toListenUmap(),{
      output$umap=renderPlotly({
        plotlyUmaps(
          scores=scores,
          sampleAttr = getSampleAttr.tmp(sampleAttr_updated(),data_mx = expressions),
          feature.name = input$feature, label = "Sample_ID")
        })
      })

    # Save data when settings are saved
    cleanUpSet=eval(parse(
      text = sprintf("reactive({list(%s)})",paste(paste(cleanUpSettings,"=input$",cleanUpSettings,sep=""),collapse = ","))
    ))
    metaHeaders=eval(parse(
      text = sprintf("reactive({list(%s)})",paste(paste(headers,"=input$",headers,sep=""),collapse = ","))
    ))

    # output$test=renderPrint({metaHeaders()})
    observeEvent(input$save_settings, {
      assign("sampleAttr", sampleAttr_updated(), envir = .GlobalEnv)
      assign("expressions", expressions, envir = .GlobalEnv)
      assign("settings",cleanUpSet(), envir = .GlobalEnv)
      assign("metaHeaders",metaHeaders(), envir = .GlobalEnv)
      sampleAttr=sampleAttr_updated()
      settings=cleanUpSet()
      metaHeaders=metaHeaders()
      forceCorrection=input$forceCorrection
      save(sampleAttr,expressions,settings,metaHeaders,forceCorrection,file="exprSetup.RData")
      stopApp()
    })
    # session$onSessionEnded(function() {
      # assign("sampleAttr", sampleAttr_updated(), envir = .GlobalEnv)
      # assign("expressions", expressions, envir = .GlobalEnv)
      # assign("settings",cleanUpSet(), envir = .GlobalEnv)
      # assign("metaHeaders",metaHeaders(), envir = .GlobalEnv)
      # save(sampleAttr,expressions,settings,metaHeaders,file="exprSetup.RData")
      # stopApp()
    # })

    })
    })
}

# Run the app
shinyApp(ui, server)
