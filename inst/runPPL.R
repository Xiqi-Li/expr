require(expr)
require(dplyr)
rm(list=ls())
pplSetUp()
# save.image("exprSetup.RData")


# Run Clean Up
load("exprSetup.RData")
file.copy(system.file("exprCleanUp.Rmd",package = "expr"), "exprCleanUp.Rmd",overwrite=T)
rmarkdown::render(
  "exprCleanUp.Rmd", params = c(settings,metaHeaders[names(metaHeaders)!="headerRm"],list(forceCorrection=forceCorrection)),
  output_file = paste0("exprCleanUp", format(Sys.time(), "%B%d%Y"),".html",collapse = "")
)
file.remove("exprCleanUp.Rmd")

# Run Main Module
load("exprSetup.RData")
load("clean_batchexamined_logRNA.RData")
if(metaHeaders$time_point!="NA"&!is.na(metaHeaders$time_point)){
  titleStr="of baseline samples"
}else{
  titleStr=""
}
OSthresh=tryCatch({format(quantile(sampleAttr[,"OS"],na.rm=T,probs=c(0,0.5,1)),digits=2)},
  error=function(e){return(c(NA,NA,NA))})
PFSthresh=tryCatch({format(quantile(sampleAttr[,"PFS"],na.rm=T,probs=c(0,0.5,1)),digits=2)},
                  error=function(e){return(c(NA,NA,NA))})
tracks.d=c("OS","vital_status","PFS","event_for_progression","best_response_perc","time_point","gender","best_response")[
  c("OS","vital_status","PFS","event_for_progression","best_response_perc","time_point","gender","best_response") %in% colnames(clean_RNA_sample_info)]
file.copy(system.file("exprMain.Rmd",package = "expr"), "exprMain.Rmd",overwrite=T)
rmarkdown::render(
  "exprMain.Rmd", params = "ask",
  output_file = paste0("exprMain", format(Sys.time(), "%B%d%Y"),".html",collapse = "")
)
file.remove("exprMain.Rmd")

# Run Timeline Module
rm(list=ls())
metaHeaders=R.utils::loadToEnv("exprSetup.RData")[["metaHeaders"]]
if(metaHeaders$time_point!="NA"&!is.na(metaHeaders$time_point)){
  load("clean_batchexamined_logRNA.RData")
  paired_RNA_sample_info=clean_RNA_sample_info[clean_RNA_sample_info$MRNsurr%in%names(which(table(clean_RNA_sample_info$MRNsurr)>1)),] %>% arrange(time_point)%>%arrange(order(gtools::mixedorder(MRNsurr)))
  paired_pt_info=reshape2::dcast(paired_RNA_sample_info[,c("MRNsurr","time_point","Sample_ID")],MRNsurr~time_point,value.var="Sample_ID",fun.aggregate = function(x) paste(unique(x),collapse = ";")) %>% tibble::column_to_rownames("MRNsurr")
  num.ind=colnames(paired_RNA_sample_info)[sapply(paired_RNA_sample_info, is.numeric)]
  for (header in colnames(paired_RNA_sample_info)[!grepl("Sample_ID|MRN|full|name",colnames(paired_RNA_sample_info),ignore.case = T)]){
    tmp=reshape2::dcast(paired_RNA_sample_info[,c("MRNsurr","time_point",header)],MRNsurr~time_point,value.var=header,fun.aggregate = function(x) paste(unique(x),collapse = ";")) %>% tibble::column_to_rownames("MRNsurr")
    if(all(apply(tmp, 1, function(x) length(unique(na.omit(x[x!=""]))))==1)){
      paired_pt_info[[header]]=apply(tmp,1,function(x) unique(na.omit(x[x!=""])))
    }
  }
  num.ind=num.ind[num.ind%in%colnames(paired_pt_info)]
  paired_pt_info[,num.ind]=apply(paired_pt_info[,num.ind],2,as.numeric)
  paired_pt_info[paired_pt_info==""]=NA
  paired_pt_info=paired_pt_info[gtools::mixedorder(rownames(paired_pt_info)),]
  paired_pt_info=paired_pt_info[complete.cases(paired_pt_info[,c(1,2)]),]
  OSthresh=tryCatch({format(quantile(paired_pt_info[,"OS"],na.rm=T,probs=c(0,0.5,1)),digits=2)},
                    error=function(e){return(c(NA,NA,NA))})
  PFSthresh=tryCatch({format(quantile(paired_pt_info[,"PFS"],na.rm=T,probs=c(0,0.5,1)),digits=2)},
                     error=function(e){return(c(NA,NA,NA))})
  tracks.d=c("OS","vital_status","PFS","event_for_progression","best_response_perc","time_point","gender","best_response")[c("OS","vital_status","PFS","event_for_progression","best_response_perc","time_point","gender","best_response") %in% colnames(paired_pt_info)]
  file.copy(system.file("exprLongitudinal.Rmd",package = "expr"), "exprLongitudinal.Rmd",overwrite=T)
  rmarkdown::render(
    "exprLongitudinal.Rmd", params = "ask",
    output_file = paste0("exprLongitudinal", format(Sys.time(), "%B%d%Y"),".html",collapse = "")
  )
  file.remove("exprLongitudinal.Rmd")
}

