
#' filter_pindels
#' @description filter pindels using various filter sets
#'
#' @param pindels \code{data.frame()}. Pindels containing insertion and deletion annotation.
#' @param traced_columns \code{vector(mode="character")}. Used for removing duplicated pindel annotations.
#'  Default set as c("Sample_ID","gene","type","chr","start","end").
#' @param low_purity_filter \code{logical()}. Remove low purity sample or keep. Default is TRUE.
#' @param low_purity_samples  \code{vector(mode="character")}. A vector of samples with low purity.
#' @param exonutronly_filter \code{logical()}. keep only pindels located in exon, UTR or not. Default is TRUE.
#' @param keep \code{vector(mode="character")}. Pindel locus (column func.knowngene) to keep. Default set as c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),
#' @param remove_td_filter \code{logical()}. Remove tandem duplicated or not. Default is TRUE.
#' @param td_coverage_filter \code{logical()}. Remove low coverage tandom duplicate. Default is TRUE.
#' @param cutoff_td_coverage \code{numeric()}. Cutoff for ratio of tandem duplicate to average genome coverage. Default 1.5.
#' @param off_target_filter \code{logical()}. Remove off-target pindels or not. Default is TRUE.
#' @param intervals \code{GRanges object} Providing the targeted regions information.
#' @param germline_filter \code{logical()}. Remove germline pindels or not. Default is TRUE.
#' @param cutoff_n_vaf \code{numeric()}. Cutoff for normal variant allele fraction. Default 0.01.
#' @param treat_known_cancer_gene_specially \code{logical()}. Split known cancer genes and treat differently. Default set to TRUE.
#' @param known_cancer_genes \code{vector(mode="character")}. A vector of cancer gene names.
#' @param special_pindel_genes \code{vector(mode="character")}. A vector of genes whose pindel no through filter sets.
#' @param mapping_quality_filter \code{logical()}. Remove low mapping quality pindels or not. Default set to TRUE.
#' @param cutoff_mq \code{numeric()}. Cutoff for averaged mapping quality of anchor reads. Default 25.
#' @param statistical_filter \code{logical()}. Remove pindels not passing proportional test or not. Default set to TRUE.
#' @param cutoff_prop_test_p \code{numeric()}. Cutoff for proportional test p-value. Default 0.01.
#' @param t_vaf_filter \code{logical()}. Remove pindels with low tumor variant allele fraction or not. Default set to TRUE.
#' @param cutoff_t_vaf \code{numeric()}. Cutoff for tumor variant allele fraction. Default 0.0.
#' @param tumor_f_filter \code{logical()}. Remove pindels with low tumor fraction or not. Default set to TRUE.
#' @param cutoff_tumor_f \code{numeric()}. Cutoff for tumor fraction. Default 0.0.
#' @param deletion_length_filter \code{logical()}. Remove or not remove very long deletions.
#' @param cutoff_d_length \code{numeric()}. Cutoff for deletion length. Default 1000.
#' @param insertion_length_filter \code{logical()}. Remove or not remove very short insertions. Default set to TRUE.
#' @param cutoff_i_length \code{numeric()}. Cutoff for insertion length. Default set to 6.
#' @param tumor_coverage_filter \code{logical()}. Remove pindels with low tumor coverage or not. Default set to TRUE.
#' @param cutoff_t_coverage \code{numeric()}. Cutoff for tumor coverage (column t_alt+ column t_ref). Default 10.
#' @param normal_coverage_filter \code{logical()}. Remove pindels with low normal coverage or not. Default set to TRUE.
#' @param cutoff_n_coverage \code{numeric()}. Cutoff for normal coverage (column n_alt+ column n_ref). Default set to 10.
#' @param support_filter \code{logical()}. Remove pindels with low support reads or not. Default set to TRUE.
#' @param cutoff_support \code{numeric()}. Cutoff for number of support reads. Default 3.
#' @param homopolymer_filter \code{logical()}. Remove pindels with too many homopolymer or not. Default set to TRUE.
#' @param maxn_homopolymer \code{numeric()}. Allowed maximal number of homopolymer. Default 6.
#' @param repeat_pindel_filter \code{logical()}.  Remove repeated pindels or not. Default set to TRUE.
#' @param cutoff_repeat \code{numeric()}. Cutoff for repeating frequency. Default 2.
#' @param consistent_pindel_statistic_filter \code{logical()}. Remove pindels failed to pass consistency statistical test (binom test). Default set to TRUE.
#' @param pindel_for_consistent_statistic \code{data.frame()} containing gene, number of pindels, number of samples (columns: gene,n_pindel,n_sample), used for binom test with targeted pindels
#' @param cutoff_binom_pval \code{numeric()}. Cutoff for binom test pvalue. Default 0.05.
#' @param description_file \code{character()}. File name for deposit of information of each operation (filter sets). Default description.txt.
#' @param save_traced_pindel \code{logical()}. Save traced_pindel_file. Default set as TRUE.
#' @param traced_pindel_file \code{character()}. File name of saved results after each operation (filter sets) if save_traced_pindel=TRUE, Default set as traced_pindel.txt.
#' @import rtracklayer plyranges
#'
#' @return \code{data.frame()}. Filtered pindels.
#' @export
#'
#' @examples
filter_pindels<-function(pindels,traced_columns=c("Sample_ID","gene","type","chr","start","end"),
                         low_purity_filter=T,low_purity_samples,
                         exonutronly_filter=T,keep=c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),
                         remove_td_filter=T,
                         td_coverage_filter=T,cutoff_td_coverage=1.5,
                         off_target_filter=T,intervals=intervals,
                         germline_filter=T,cutoff_n_vaf=0.01,
                         treat_known_cancer_gene_specially=T,known_cancer_genes=c("TP53","GNA11","GNAQ"),special_pindel_genes=c("BAP1","PRKDC"),
                         mapping_quality_filter=T,cutoff_mq=25,
                         statistical_filter=T,cutoff_prop_test_p=0.01,
                         t_vaf_filter=T,cutoff_t_vaf=0.0,
                         tumor_f_filter=T,cutoff_tumor_f=0,
                         deletion_length_filter=T,cutoff_d_length=1000,
                         insertion_length_filter=T,cutoff_i_length=6,
                         tumor_coverage_filter=T,cutoff_t_coverage=10,
                         normal_coverage_filter=T,cutoff_n_coverage=10,
                         support_filter=T,cutoff_support=3,
                         homopolymer_filter=T,maxn_homopolymer=6,
                         repeat_pindel_filter=T,cutoff_repeat=2,
                         consistent_pindel_statistic_filter=T,pindel_for_consistent_statistic=NA,cutoff_binom_pval=0.05,
                         description_file="description.txt",
                         save_traced_pindel=T,traced_pindel_file="traced_pindel.txt"){
  #filter derived from SVFilter,gatk,and github.com/genome/pindel
  options(warn=-1)
  require(plyranges)
  pindels<-pindels[!duplicated(do.call(paste,pindels[,traced_columns])),]
  fileConn<-description_file
  write(paste(paste("Starting from ",nrow(pindels),sep="")," indels.",sep=""),fileConn,append=F)
  genome_coverage_tf<- pindels %>% group_by(Sample_ID) %>% summarise(genome_coverage=sum(t_ref_count+t_alt_count)/n())
  genome_coverage<-genome_coverage_tf[["genome_coverage"]]
  names(genome_coverage)<-genome_coverage_tf[["Sample_ID"]]
  pindels_traced<-list()
  traced_names<-c()
  pindels_traced[["original_pindels"]]<-pindels[,traced_columns]
  traced_names<-c(traced_names,"original_pindels")

  #low purity filter
  if(low_purity_filter){
    pindels<-pindels[!pindels[["Sample_ID"]] %in% low_purity_samples,]
    write("After filter out low purity (0.5);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_low_purity"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_low_purity")
  }

  #exon UTRs only
  if(exonutronly_filter){
    pindels<-pindels[pindels$func.knowngene %in% c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),]
    write("After filter out non exonic or UTRs", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_nonexonutr"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #remove TD
  if(remove_td_filter){
    pindels<-pindels[pindels$type!="TD",]
    write("After filter out tandem duplicate indels;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_TD"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  } else if(td_coverage_filter){
    #sequencing relative depth (ratio between the average of sequencing depth in the duplicated region and that over entire genome) filter >=1.5
    genome_coverage<-genome_coverage[pindels[['Sample_ID']]]
    pindels<-pindels[!(pindels$type=="TD" & ((pindels$t_ref_count+pindels$t_alt_count)>=(cutoff_td_coverage*genome_coverage))),]
    write(paste("After filter out low coverage tandem duplicate indels with cutoff_td_coverage ",cutoff_td_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_TD_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #offtarget filter
  if(off_target_filter){
    library(rtracklayer)
    library(plyranges)
    pindels_gr<-df2granges(pindels,genome = "hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,seqnames_col = "chr",start_col = "start",end_col = "end",meta_cols = colnames(pindels))
    seqlevels(intervals)<-seqlevels(pindels_gr)
    seqinfo(intervals)<-seqinfo(pindels_gr)
    targets<-filter_by_overlaps(pindels_gr,intervals)
    pindels<-as.data.frame(mcols(targets))
    write("After filter out off-targets indel;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_off_target"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_off_target")
  }

  #germline cutoff 0.01
  if(germline_filter){
    pindels<-pindels[pindels$n_vaf<=cutoff_n_vaf,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out germline indel with cutoff_n_vaf ",cutoff_n_vaf,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_germline"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_germline")
  }

  #treat known cancer gene specially
  if(treat_known_cancer_gene_specially){
    if(length(known_cancer_genes)>0){
      known_cancer_gene_pindels<-pindels[pindels$gene %in% setdiff(known_cancer_genes,special_pindel_genes),]
    }
    if(length(special_pindel_genes)>0){
      special_cancer_gene_pindels<-pindels[pindels$gene %in% special_pindel_genes,,drop=F]
    }
    pindels<-pindels[!(pindels$gene %in% known_cancer_genes),]
    write("After split known cancer gene (cosmic or cancer_genetic);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," non known cancer gene pindels.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(known_cancer_gene_pindels),sep="")," known cancer gene pindels.",sep=""), fileConn,append = T)
    pindels_traced[["non_known_cancer_gene_pindel"]]<-pindels[,traced_columns]
    pindels_traced[["known_cancer_gene_pindel"]]<-known_cancer_gene_pindels[,traced_columns]
    traced_names<-c(traced_names,"non_known_cancer_gene_pindel")
    traced_names<-c(traced_names,"known_cancer_gene_pindel")
    if(length(special_pindel_genes)>0){
      write(paste(paste("There are still ",nrow(special_cancer_gene_pindels),sep="")," special cancer gene pindels.",sep=""), fileConn,append = T)
      pindels_traced[["special_cancer_gene_pindel"]]<-special_cancer_gene_pindels[,traced_columns]
      traced_names<-c(traced_names,"special_cancer_gene_pindel")
    }
  }

  #mapping quality filter
  if(mapping_quality_filter){
    pindels<-pindels[(pindels$mq/pindels$support)>=cutoff_mq,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low mapping quality indels with cutoff_mq ",cutoff_mq,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_mapping_quality"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_mapping_quality")
  }

  #statistical filter
  if(statistical_filter){
    pval<-apply(pindels[,c("t_ref_count","t_alt_count",'n_ref_count',"n_alt_count")],1,function(row){prop.test(x=c(as.integer(row[2]),as.integer(row[4])),n=c(as.integer(row[1])+as.integer(row[2])+1,as.integer(row[3])+as.integer(row[4])+1))$p.value})
    pindels<-pindels[pval<=cutoff_prop_test_p,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out prop test failed indels with cufoff_prop_test_p ",cutoff_prop_test_p,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_prop_statistic"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_prop_statistic")
  }

  # tumor_f filter
  if(tumor_f_filter){
    pindels<-pindels[pindels$tumor_f>=cutoff_tumor_f,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out lwo tumor fraction indels with cutoff_tumor_f ",cutoff_tumor_f,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_tumor_f"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_f")
  }

  # t_vaf filter
  if(t_vaf_filter){
    pindels<-pindels[pindels$t_vaf>=cutoff_t_vaf,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low t_vaf indels with cutoff_t_vaf ",cutoff_t_vaf,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_t_vaf"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_t_vaf")
  }

  #deletion length filter
  if(deletion_length_filter){
    pindels<-pindels[pindels$length<=cutoff_d_length,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out  too long deletions with cutoff_d_length",cutoff_d_length,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_deletion_length"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_deletion_length")
  }

  #insertion length filter
  if(insertion_length_filter){
    pindels<-pindels[pindels$insertlen<=cutoff_i_length,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out too long insertion with cutoff_i_length ",cutoff_i_length,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_insertion_length"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_insertion_length")
  }

  #tumor depth coverage filter
  if(tumor_coverage_filter){
    pindels<-pindels[(pindels$t_alt_count+pindels$t_ref_count)>=cutoff_t_coverage,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low tumor coverage indels with cutoff_t_coverage ",cutoff_t_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_tumor_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_coverage")
  }

  #normal depth coverage filter
  if(normal_coverage_filter){
    pindels<-pindels[(pindels$n_alt_count+pindels$n_ref_count)>=cutoff_n_coverage,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low normal coverage indels with cutoff_n_coverage ",cutoff_n_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_normal_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_normal_coverage")
  }

  #support filter
  if(support_filter){
    pindels<-pindels[pindels$support>=cutoff_support,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low support indels with cutoff_support",cutoff_support,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_low_support"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_low_support")
  }

  #homopolymer filter
  if(homopolymer_filter){
    #homopolymer(many copies of a single repeating unit) filer maxlength=6
    homopolymer_pat<-paste(c(paste("A{",maxn_homopolymer,",}",sep=""),paste("T{",maxn_homopolymer,",}",sep=""),paste("C{",maxn_homopolymer,",}",sep=""),paste("G{",maxn_homopolymer,",}",sep="")),collapse="|")
    pindels<-pindels[!grepl(homopolymer_pat,pindels$refs,perl = T) | grepl(homopolymer_pat,pindels$sams,perl = T),]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out too much homopolymer indels with cutoff_maxn_homopolymer ",maxn_homopolymer,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_homopolymer"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_homopolymer")
  }

  #combine known_cancer gene pindels
  if(treat_known_cancer_gene_specially & length(known_cancer_genes)>0){
    pindels<-rbind(pindels,known_cancer_gene_pindels)
    write("After combine known cancer gene pindels;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["combined_known_pindels"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"combined_known_pindels")
  }


  #reapeat pindel filter
  if(repeat_pindel_filter){
    pindels_<-pindels
    pindels_[["Pindel_ID"]]<-do.call(paste,pindels_[,c("gene","chr","start")])
    repeat_pindel_IDs<-names(table(pindels_[["Pindel_ID"]]))[table(pindels_[["Pindel_ID"]])>=cutoff_repeat]
    pindels_<-pindels_[!(pindels_[["Pindel_ID"]] %in% repeat_pindel_IDs),]
    pindels<-pindels_[,-ncol(pindels_)]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out high repeated indels with cutoff_repeat ",cutoff_repeat,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",sum(pindels$gene %in% known_cancer_gene_pindels$gene),sep="")," known cancer gene indels.",sep=""), fileConn,append = T)
    write(paste(paste("There are ",nrow(known_cancer_gene_pindels)-sum(pindels$gene %in% known_cancer_genes),sep="")," known cancer gene indels removed.",sep=""), fileConn,append = T)
    pindels_traced[["filter_repeat_pindel"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_repeat_pindel")
  }

  #browser()
  if(consistent_pindel_statistic_filter){
    n_sample_subject<-length(unique(pindels[["Sample_ID"]]))
    pindels_subject<-pindels %>% group_by(gene) %>% summarise(n_pindel_subject=n())
    pindels_subject[["n_sample_subject"]]<-n_sample_subject
    pindels_for_statistic<-merge(pindels_subject,pindel_for_consistent_statistic,by="gene",all.x=T)
    pindels_for_statistic[["n_sample"]][is.na(pindels_for_statistic[["n_sample"]])]=max(pindel_for_consistent_statistic$n_sample)
    pindels_for_statistic[["n_pindel"]][is.na(pindels_for_statistic[["n_pindel"]])]<-0
    pindels_for_statistic<-set_column_as_rownames(pindels_for_statistic,"gene")
    pindels_for_statistic<-pindels_for_statistic[pindels_for_statistic[["n_pindel_subject"]]<=pindels_for_statistic[["n_sample_subject"]],]
    pindels_for_statistic$n_pindel_subject<-apply(pindels_for_statistic[,c("n_pindel_subject","n_sample_subject")],1,min)
    pindels_for_statistic$n_pindel<-apply(pindels_for_statistic[,c("n_pindel","n_sample")],1,min)
    pvals<-apply(pindels_for_statistic,1,function(r){test<-binom.test(r["n_pindel_subject"],r["n_sample_subject"],p=r["n_pindel"]/r["n_sample"]);return(test$p.value)})
    adjpvals<-p.adjust(pvals,method="BH")
    keep_genes<-names(adjpvals)[adjpvals>cutoff_binom_pval]
    pindels<-pindels[pindels[["gene"]] %in% keep_genes,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out non consitent pindel with tcga pindel with cutoff_binom_pval ",cutoff_binom_pval,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",sum(pindels$gene %in% known_cancer_gene_pindels$gene),sep="")," known cancer gene indels.",sep=""), fileConn,append = T)
    write(paste(paste("There are ",nrow(known_cancer_gene_pindels)-sum(pindels$gene %in% known_cancer_genes),sep="")," known cancer gene indels removed.",sep=""), fileConn,append = T)
    pindels_traced[["filter_nonconsistent_pindel_statistic"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonconsistent_pindel_statistic")
  }

  #combine special cancer gene pindels
  if(treat_known_cancer_gene_specially & length(special_pindel_genes)>0){
    pindels<-rbind(pindels,special_cancer_gene_pindels)
    write("After combine special cancer gene pindels;", fileConn,append = T)
    write(paste(paste("There are still totally ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["combined_final_pindels"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"combined_final_pindels")
  }

  pindel_order=NA
  for(traced_name in rev(traced_names)){
    if(traced_name %in% c("special_cancer_gene_pindel","known_cancer_gene_pindel")){
      next
    } else{
      pindels_temp=pindels_traced[[traced_name]]
      rownames(pindels_temp)<-do.call(paste,pindels_temp)
      pindels_temp<-pindels_temp[c(intersect(pindel_order,rownames(pindels_temp)),setdiff(rownames(pindels_temp),pindel_order)),]
      pindel_order=rownames(pindels_temp)
      rownames(pindels_temp)<-NULL
      pindels_traced[[traced_name]]=pindels_temp
    }
  }
  if(save_traced_pindel){
    for(traced_name in traced_names){
      pindels_temp<-pindels_traced[[traced_name]]
      colnames(pindels_temp)<-paste(traced_name,colnames(pindels_temp),sep=":")
      write.table(pindels_temp,file = traced_pindel_file,append = T,row.names = F,sep="\t")
    }
  }
  return(pindels)
}

