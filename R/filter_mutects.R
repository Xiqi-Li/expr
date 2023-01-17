
#' filter_mutects
#' @description filter mutects using various filter sets
#'
#' @param mutects \code{data.frame()} containing SNV annotation.
#' @param traced_columns \code{vector(mode="character")}. Used for removing duplicated SNV annotations.
#'  Default set as c("Sample_ID","gene","type","chr","start","end").
#' @param low_purity_filter \code{logical()}. Remove low purity sample or keep. Default is TRUE.
#' @param low_purity_samples \code{vector(mode="character")}. A vector of samples with low purity.
#' @param exonutronly_filter \code{logical()}. Keep only SNVs located in exon, UTR or not. Default is TRUE.
#' @param keep \code{vector(mode="character")}. A vector of SNV loci (column func.knowngene) to keep.
#'  Default c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3").
#' @param off_target_filter \code{logical()}. Remove off-target SNVs or not. Default is TRUE.
#' @param intervals \code{GRanges object}. Providing the targeted regions information.
#' @param germline_filter \code{logical()}. Remove germline SNVs or not. Default is TRUE.
#' @param cutoff_n_vaf \code{numeric()}. Cutoff for normal variant allele fraction. Default 0.01.
#' @param treat_known_cancer_gene_specially \code{logical()}. Split known cancer genes and treat differently. Default set to TRUE.
#' @param known_cancer_genes \code{vector(mode="character")}. A vector of cancer gene names.
#' @param special_mutect_genes \code{vector(mode="character")}. A vector of genes whose SNVs should be exempt from filter sets.
#' @param mapping_quality_filter \code{logical()}. Remove low mapping quality SNVs or not. Default TRUE.
#' @param cutoff_mq \code{numeric()}. Cutoff for maximal mapping quality of reads. Default 30,
#' @param statistical_filter \code{logical()}. Remove SNVs not passing proportional test or not. Default TRUE.
#' @param cutoff_prop_test_p \code{numeric()}. Cutoff for proportional test pvalue. Default 0.01.
#' @param tumor_f_filter \code{logical()}. Remove SNVs with low tumor fraction or not. Default TRUE.
#' @param cutoff_tumor_f \code{numeric()}. Cutoff for tumor fraction. Default 0.0.
#' @param tumor_coverage_filter \code{logical()}. Remove SNVs with low tumor coverage or not. Default TRUE.
#' @param cutoff_t_coverage \code{numeric()}. Cutoff for tumor coverage (column t_alt+ column t_ref). Default 10.
#' @param normal_coverage_filter \code{logical()}. Remove SNVs with low normal coverage or not. Default TRUE.
#' @param cutoff_n_coverage \code{numeric()}. Cutoff for normal coverage (column n_alt+ column n_ref). Default 10.
#' @param lodt_filter \code{logical()}. Remove SNVs with low log likelihood of tumor event or not. Default TRUE.
#' @param cutoff_lodt \code{numeric()}. Cutoff for number of log likelihood of tumor event. Default 6.3.
#' @param consistent_mutect_statistic_filter \code{logical()}. Remove SNVs failed to pass consistency statistical test (binom test). Default set to TRUE.
#' @param mutect_for_consistent_statistic \code{data.frame()} object containing gene, number of SNVs, number of samples (columns: gene,n_pindel,n_sample), used for binom test with targeted SNVs
#' @param cutoff_binom_pval \code{numeric()}. Cutoff for binom test p-value. Default 0.05.
#' @param description_file \code{character()}. File name for deposit of information of each operation (filter sets). Default "description.txt".
#' @param save_traced_mutectsave \code{logical()}. Save traced_mutect_file or not. Default set to TRUE.
#' @param traced_mutect_file \code{character()}. File name for deposit results after each operation (filter sets) if save_traced_mutect is TRUE.
#'  Default "traced_pindel.txt".

#'
#' @return \code{data.frame()}. Filtered SNVs.
#' @export
#'
filter_mutects<-function(mutects, traced_columns=c("Sample_ID","gene","type","chr","start","end"),
                         low_purity_filter=T,low_purity_samples,
                         exonutronly_filter=T,keep=c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),
                         off_target_filter=T,intervals=intervals,
                         germline_filter=T,cutoff_n_vaf=0.01,
                         treat_known_cancer_gene_specially=T,known_cancer_genes=c("TP53","GNA11","GNAQ"),special_mutect_genes=c(),
                         mapping_quality_filter=T,cutoff_mq=30,
                         statistical_filter=T,cutoff_prop_test_p=0.01,
                         tumor_f_filter=T,cutoff_tumor_f=0,
                         tumor_coverage_filter=T,cutoff_t_coverage=10,
                         normal_coverage_filter=T,cutoff_n_coverage=10,
                         lodt_filter=T,cutoff_lodt=6.3,
                         consistent_mutect_statistic_filter=T,mutect_for_consistent_statistic=NA,cutoff_binom_pval=0.05,
                         description_file="description.txt",
                         save_traced_mutect=TRUE,traced_mutect_file="traced_mutect.txt"){
  require(plyranges)
  options(warn=-1)
  mutects<-mutects[!duplicated(do.call(paste,mutects[,traced_columns])),]
  fileConn<-description_file
  write(paste(paste("Starting from ",nrow(mutects),sep="")," mutects.",sep=""),fileConn,append=F)
  mutects_traced<-list()
  traced_names<-c()
  mutects_traced[["original_mutects"]]<-mutects[,traced_columns]
  traced_names<-c(traced_names,"original_mutects")

  #low purity filter
  if(low_purity_filter){
    mutects<-mutects[!mutects[["Sample_ID"]] %in% low_purity_samples,]
    write("After filter out low purity (0.5);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_low_purity"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_low_purity")
  }

  #exon UTRs only
  if(exonutronly_filter){
    mutects<-mutects[mutects$func.knowngene %in% c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),]
    write("After filter out non exonic or UTRs", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_nonexonutr"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #offtarget filter
  if(off_target_filter){
    library(rtracklayer)
    library(plyranges)
    mutects_gr<-df2granges(mutects,genome = "hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,seqnames_col = "chr",start_col = "start",end_col = "end",meta_cols = colnames(mutects))
    seqlevels(intervals)<-seqlevels(mutects_gr)
    seqinfo(intervals)<-seqinfo(mutects_gr)
    targets<-filter_by_overlaps(mutects_gr,intervals)
    mutects<-as.data.frame(mcols(targets))
    write("After filter out off-targets mutects;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_off_target"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_off_target")
  }

  #germline cutoff 0.01
  if(germline_filter){
    mutects<-mutects[(mutects$n_alt_count/(mutects$n_alt_count+mutects$n_ref_count))<=cutoff_n_vaf,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out germline mutects with cutoff_n_vaf ",cutoff_n_vaf,sep=" "),fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_germline"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_germline")
  }

  #treat known cancer gene specially
  if(treat_known_cancer_gene_specially){
    if(length(known_cancer_genes)>0){
      known_cancer_gene_mutects<-mutects[mutects$gene %in% setdiff(known_cancer_genes,special_mutect_genes),]
    }
    if(length(special_mutect_genes)>0){
      special_cancer_gene_mutects<-mutects[mutects$gene %in% special_mutect_genes,,drop=F]
    }

    mutects<-mutects[!(mutects$gene %in% known_cancer_genes),]
    write("After split known cancer gene (cosmic or cancer_genetic);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," non known cancer gene mutects.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(known_cancer_gene_mutects),sep="")," known cancer gene mutects.",sep=""), fileConn,append = T)
    mutects_traced[["non_known_cancer_gene_mutect"]]<-mutects[,traced_columns]
    mutects_traced[["known_cancer_gene_mutect"]]<-known_cancer_gene_mutects[,traced_columns]
    traced_names<-c(traced_names,"non_known_cancer_gene_mutect")
    traced_names<-c(traced_names,"known_cancer_gene_mutect")
    if(length(special_mutect_genes)>0){
      write(paste(paste("There are still ",nrow(special_cancer_gene_mutects),sep="")," special cancer gene mutects.",sep=""), fileConn,append = T)
      mutects_traced[["special_cancer_gene_mutect"]]<-special_cancer_gene_mutects[,traced_columns]
      traced_names<-c(traced_names,"special_cancer_gene_mutect")
    }
  }

  #mapping quality filter
  if(mapping_quality_filter){
    mutects<-mutects[(mutects$t_ref_max_mapq>=cutoff_mq & mutects$t_alt_max_mapq>=cutoff_mq),]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out mutects with low average mapping quality with cutoff_mq ",cutoff_mq,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_mapping_quality"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_mapping_quality")
  }

  #statistical filter
  #browser()
  if(statistical_filter){
    pval<-apply(mutects[,c("t_ref_count","t_alt_count",'n_ref_count',"n_alt_count")],1,function(row){prop.test(x=c(as.integer(row["t_alt_count"]),as.integer(row["n_alt_count"])),n=c(as.integer(row["t_alt_count"])+as.integer(row["t_ref_count"])+1,as.integer(row["n_alt_count"])+as.integer(row["n_ref_count"])+1))$p.value})
    mutects<-mutects[pval<=cutoff_prop_test_p,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out prop test failed mutects with cutoff_prop_test_p ",cutoff_prop_test_p,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_prop_statistic"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_prop_statistic")
  }

  # tumor_f filter
  if(tumor_f_filter){
    mutects<-mutects[mutects$tumor_f>=cutoff_tumor_f,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low tumor fraction mutects with  cutoff_tumor_f ",cutoff_tumor_f,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_tumor_f"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_f")
  }

  #tumor depth coverage filter
  if(tumor_coverage_filter){
    mutects<-mutects[(mutects$t_alt_count+mutects$t_ref_count)>=cutoff_t_coverage,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low tumor coverage mutects with cutoff_t_coverage ",cutoff_t_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_tumor_coverage"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_coverage")
  }

  #normal depth coverage filter
  if(normal_coverage_filter){
    mutects<-mutects[(mutects$n_alt_count+mutects$n_ref_count)>=cutoff_n_coverage,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low normal coverage mutects with cutoff_n_coverage ",cutoff_n_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_normal_coverage"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_normal_coverage")
  }

  #lodt filter
  if(lodt_filter){
    mutects<-mutects[mutects$t_lod_fstar>=cutoff_lodt,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low t_lod_fstar mutects with cutoff_lodt ",cutoff_lodt,sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_lodt"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_lodt")
  }


  #combine known_cancer gene mutects
  if(treat_known_cancer_gene_specially & length(known_cancer_genes)>0){
    mutects<-rbind(mutects,known_cancer_gene_mutects)
    write("After combine known cancer gene mutects;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["combined_known_mutects"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"combined_known_mutects")
  }

  #browser()
  if(consistent_mutect_statistic_filter){
    n_sample_subject<-length(unique(mutects[["Sample_ID"]]))
    mutects_subject<-mutects %>% group_by(gene) %>% summarise(n_mutect_subject=n())
    mutects_subject[["n_sample_subject"]]<-n_sample_subject

    mutects_for_statistic<-merge(mutects_subject,mutect_for_consistent_statistic,by=c("gene"),all.x=T)
    mutects_for_statistic[["n_sample"]][is.na(mutects_for_statistic[["n_sample"]])]=max(mutect_for_consistent_statistic$n_sample)
    mutects_for_statistic[["n_mutect"]][is.na(mutects_for_statistic[["n_mutect"]])]<-0
    mutects_for_statistic<-set_column_as_rownames(mutects_for_statistic,"gene")
    # mutects_for_statistic_<-mutects_for_statistic[,c("n_mutect_subject","n_sample_subject","n_mutect","n_sample")]
    mutects_for_statistic$n_mutect_subject<-apply(mutects_for_statistic[,c("n_mutect_subject","n_sample_subject")],1,min)
    mutects_for_statistic$n_mutect<-apply(mutects_for_statistic[,c("n_mutect","n_sample")],1,min)
    pvals<-apply(mutects_for_statistic,1,function(r){test<-binom.test(r["n_mutect_subject"],r["n_sample_subject"],p=r["n_mutect"]/r["n_sample"]);return(test$p.value)})
    adjpvals<-p.adjust(pvals,method="BH")
    keep_genes<-names(adjpvals)[adjpvals>cutoff_binom_pval]
    mutects<-mutects[mutects$gene %in% keep_genes,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out non consitent mutect with tcga mutect with cutoff_binom_pval ",cutoff_binom_pval,sep=" "), fileConn,append = T)
    write(paste(paste("There are still totally ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",sum(mutects$gene %in% known_cancer_gene_mutects$gene),sep="")," known cancer gene mutects.",sep=""), fileConn,append = T)
    write(paste(paste("There are ",nrow(known_cancer_gene_mutects)-sum(mutects$gene %in% known_cancer_genes),sep="")," known cancer gene mutects removed.",sep=""), fileConn,append = T)
    mutects_traced[["filter_nonconsistent_mutect_statistic"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_nonconsistent_mutect_statistic")
  }

  #combine special cancer gene mutects
  if(treat_known_cancer_gene_specially & length(special_mutect_genes)>0){
    mutects<-rbind(mutects,special_cancer_gene_mutects)
    write("After combine special cancer gene mutects;", fileConn,append = T)
    write(paste(paste("There are still totally",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["combined_final_mutects"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"combined_final_mutects")
  }

  mutect_order=NA
  for(traced_name in rev(traced_names)){
    if(traced_name %in% c("special_cancer_gene_mutect","known_cancer_gene_mutect")){
      next
    } else{
      mutects_temp=mutects_traced[[traced_name]]
      rownames(mutects_temp)<-do.call(paste,mutects_temp)
      mutects_temp<-mutects_temp[c(intersect(mutect_order,rownames(mutects_temp)),setdiff(rownames(mutects_temp),mutect_order)),]
      mutect_order=rownames(mutects_temp)
      rownames(mutects_temp)<-NULL
      mutects_traced[[traced_name]]=mutects_temp
    }
  }

  if(save_traced_mutect){
    for(traced_name in traced_names){
      mutects_temp<-mutects_traced[[traced_name]]
      colnames(mutects_temp)<-paste(traced_name,colnames(mutects_temp),sep=":")
      write.table(mutects_temp,file = traced_mutect_file,append = T,row.names = F,sep="\t")
    }
  }
  return(mutects)
}
