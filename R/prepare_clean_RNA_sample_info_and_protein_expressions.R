# prepare_clean_RNA_sample_info_and_protein_expressions
# mark RNA sample information table with is_sequenced information, sample library size,  sample tumorpurity,sample low expressed gene percentage, and extract protein coding genes, and finally remove low express gene
# and prepare clean RNA sample information table and consistent protein coding mRNA expressions

# contact: hzhu2@mdanderson.org

# version:0.1.0

# OS:ubuntu 16.04
# R:4.1.2

#' prepare_clean_RNA_sample_info_and_protein_expressions
#' @description Mark RNA sample information table with is_sequenced information, sample library size,  sample tumorpurity,sample low expressed gene percentage, and extract protein coding genes, and finally remove low express gene and prepare clean RNA sample information table and consistent protein coding mRNA expressions.
#' @param RNA_sample_info \code{data.frame()}. RNA sample information table, must be provided and contain sample id column
#' @param sample_id_col \code{character()}. default ="Sample_ID",column name for sample_id column in RNA_sample info
#' @param expressions \code{data.frame()}. no default, gene expressions table,must be provided
#' @param gene_id_col \code{character()}. default="rowname" if rownames of expressions are gene id, otherwise, provide a column name containing gene id
#' @param mark_library_size \code{logical(1)}. default TRUE, mark RNA sample information table with sample library size
#' @param adjusted_by_library_size \code{logical(1)}. default TRUE, normalize expressions with corresponding sample library size
#' @param tolerant_library_size_factor  \code{numeric()}. default=1.1, if maximal sample library size is over tolerant_library_size_factor multiple minimal sample library size, trigger normalization of expression by sample library size, otherwise keep unchanged
#' @param protein_coding_ensemble2symbol_table \code{data.frame()}.  Protein coding gene table containing ensemble gene id and hugo symbol, used for extraction of protein coding expressions and convert ensemble gene id to hugo symbol
#' @param ens_id_col \code{character()}. default ="Gene_ID",column name for ensemble gene id column in protein coding ensemble2symbol table
#' @param symbol_col \code{character()}. default ="Gene_Symbol",column name for hugo symbol column in protein coding ensemble2symbol table
#' @param mark_purity \code{logical(1)}. default TRUE, mark RNA sample information table with sample purity from estimate algorithm
#' @param remove_low_purity_sample \code{logical(1)}. default TRUE, remove sample with low purity than low_purity_threshold

#' @param low_purity_threshold \code{numeric()}. default=0.3, if remove_low_purity_sample as TRUE, sample with purity lower than low_purity_threshold  will be excluded
#' @param mark_lowexpressgene_pct \code{logical(1)}.default TRUE, mark RNA sample information table with sample low expressed gene percentage
#' @param lowexpression_threshold \code{numeric()}. default=1, define whether gene is low expressed based on the threshold
#' @param remove_sample_with_intense_lowexpressgene \code{logical(1)}. default TRUE, remove sample with too high percentage of low expressed gene
#' @param outlier_lowexpressgene_pct_factor \code{numeric()}. default=1.5, if remove_sample_with_intense_lowexpressgene is TRUE, sample with low expressed gene percentage over the product of outlier_lowexpressgene_pct_factor and average low expressed gene percentage across all samples will be excluded
#' @param remove_lowexpressgene \code{logical(1)}. default TRUE, remove low expressed genes
#' @param sample_frequency_threshold \code{numeric()}. default=0.5, if remove_lowexpressgene is TRUE, genes with low expression occurred in at least sample_frequency_threshold fration of all samples will be removed

#'
#' @return \code{list()}, contains marked original RNA sample information table, filtered clean RNA sample information table and filtered clean and consistent protein expressions.

#' @export
#'
#' @examples
#' results<-prepare_clean_RNA_sample_info_and_protein_expressions(RNA_sample_info=RNA_sample_info,sample_id_col = "Sample_ID",expressions = expressions,gene_id_col = "gene_id",                                     #
#'                                                                mark_library_size = T,adjusted_by_library_size = T,tolerant_library_size_factor = 1.1,                                                             #
#'                                                                protein_coding_ensemble2symbol_table = protein_coding_ensemble2symbol,ens_id_col = "Gene_ID",symbol_col = "Gene_Symbol",                           #
#'                                                                mark_purity = T,remove_low_purity_sample = T,low_purity_threshold = 0.3,                                                                           #
#'                                                                mark_lowexpressgene_pct = T,lowexpression_threshold = 1,remove_sample_with_intense_lowexpressgene = T,outlier_lowexpressgene_pct_factor = 1.5,     #
#'                                                                remove_lowexpressgene = T,sample_frequency_threshold = 0.5)
#'
prepare_clean_RNA_sample_info_and_protein_expressions<-function(RNA_sample_info,sample_id_col=c("Sample_ID","Tumor_Sample_Barcode","sample_id","id"),expressions,gene_id_col=c("rowname","gene_id","Gene_ID","ID"),
                                                                mark_library_size=T,adjusted_by_library_size=T,tolerant_library_size_factor=1.1,
                                                                protein_coding_ensemble2symbol_table,ens_id_col=c("Gene_ID","gene_id","id"),symbol_col=c("Gene_Symbol","Gene_Name","gene_symbol","gene_name"),
                                                                mark_purity=T,remove_low_purity_sample=T,low_purity_threshold=0.3,
                                                                mark_lowexpressgene_pct=T,lowexpression_threshold=1,remove_sample_with_intense_lowexpressgene=T,outlier_lowexpressgene_pct_factor=1.5,
                                                                remove_lowexpressgene=T,sample_frequency_threshold=0.5
){
  # sample_id_col=match.arg(sample_id_col) #XL
  # gene_id_col<-match.arg(gene_id_col)#XL
  # ens_id_col<-match.arg(ens_id_col)#XL
  # symbol_col<-match.arg(symbol_col)#XL
  #XL
  sample_id_col=sample_id_col[1]
  gene_id_col<gene_id_col[1]
  ens_id_col<-ens_id_col[1]
  symbol_col<-symbol_col[1]

  sequenced_RNA_samples<-colnames(expressions)[-match(gene_id_col,colnames(expressions))]
  stopifnot("some sequencing sample come without sample infomation!"=all( sequenced_RNA_samples %in% RNA_sample_info[[sample_id_col]]))
  if(length(RNA_sample_info[[sample_id_col]])>length(sequenced_RNA_samples)){
    cat(paste(setdiff(RNA_sample_info[[sample_id_col]],sequenced_RNA_samples),collapse = ", ")," are not sequenced!\n")
  }
  RNA_sample_info[["IS_Sequenced"]]<-ifelse(RNA_sample_info[[sample_id_col]] %in% sequenced_RNA_samples,"YES","NO")

  #Library size balance
  library_size<-colSums(expressions[,sequenced_RNA_samples],na.rm=T)
  if(mark_library_size){
    RNA_sample_info[["Library_Size"]]<-NA
    RNA_sample_info[["Library_Size"]][match(sequenced_RNA_samples,RNA_sample_info[[sample_id_col]])]<-library_size
  }

  if(adjusted_by_library_size){
    if(max(library_size,na.rm=T)>(tolerant_library_size_factor*(min(library_size,na.rm=T)))){
      expressions[,sequenced_RNA_samples]<-as.data.frame(t(t(expressions[,sequenced_RNA_samples])/library_size*mean(library_size,na.rm=T)))
    } else {
      cat("all sample libraries are quite stable, no need to adjust them!")
    }
  }

  #extract protein coding gene
  if(gene_id_col=="rowname"){
    if(any(grepl("ENSG",rownames(expressions)))){
      if(missing(protein_coding_ensemble2symbol_table)) {stop("Since expression presented with ensemble id,protein coding gene ensemble to symbol table is required!")}
      protein_expressions<-expressions[rownames(expressions) %in% protein_coding_ensemble2symbol_table[[ens_id_col]],]
      rownames(protein_expressions)<-protein_coding_ensemble2symbol_table[[symbol_col]][match(rownames(protein_expressions),protein_coding_ensemble2symbol_table[[ens_id_col]])]
    } else {protein_expressions<-expressions[rownames(expressions) %in% protein_coding_ensemble2symbol_table[[symbol_col]]]}
  } else {
    if(any(grepl("ENSG",expressions[[gene_id_col]]))){
      if(missing(protein_coding_ensemble2symbol_table)) {stop("Since expression presented with ensemble id,protein coding gene ensemble to symbol table is required!")}
      protein_expressions<-expressions[expressions[[gene_id_col]] %in% protein_coding_ensemble2symbol_table[[ens_id_col]],]
      rownames(protein_expressions)<-protein_coding_ensemble2symbol_table[[symbol_col]][match(protein_expressions[[gene_id_col]],protein_coding_ensemble2symbol_table[[ens_id_col]])]
      protein_expressions<-protein_expressions[,-match(gene_id_col,colnames(protein_expressions))]
    } else {
      protein_expressions<-expressions[expressions[[gene_id_col]] %in% protein_coding_ensemble2symbol_table[[symbol_col]]]
      rownames(protein_expressions)<-protein_expressions[[gene_id_col]]
      protein_expressions<-protein_expressions[,-match(gene_id_col,colnames(protein_expressions))]
    }
  }

  #sample tumor purity
  estimate_res<-as.data.frame(t(immunedeconv::deconvolute_estimate(protein_expressions)))
  estimate_res<-data.frame(Sample_ID=rownames(estimate_res),estimate_res)
  if(mark_purity){
    RNA_sample_info<-merge(RNA_sample_info,estimate_res,by.x=sample_id_col,by.y="Sample_ID",all.x=T)
  }

  #calculate percentage of lowexpressgene
  tatal_low_protein_expression<-sum(protein_expressions<lowexpression_threshold)
  average_lowexpressgene_pct<-tatal_low_protein_expression/(nrow(protein_expressions)*ncol(protein_expressions))*100
  cat("On average, every sample have ",average_lowexpressgene_pct,"% low express genes.\n\n")
  sample_lowexpressgene_pct<-colSums(protein_expressions<lowexpression_threshold)/nrow(protein_expressions)*100
  outlier_sample_lowexpressgene_pct<-sample_lowexpressgene_pct[sample_lowexpressgene_pct>average_lowexpressgene_pct*outlier_lowexpressgene_pct_factor]
  if(length(outlier_sample_lowexpressgene_pct)==0){
    cat("distribution of low expression gene for every sample is similar.")
  } else{
    cat("Sample",paste(names(outlier_sample_lowexpressgene_pct),collapse = ","),"has (have) too many low expression genes: ",paste(outlier_sample_lowexpressgene_pct,collapse=","),"% respectively.\n\n")
  }
  if(mark_lowexpressgene_pct){
    RNA_sample_info[["lowexpressgene_pct"]]<-NA
    RNA_sample_info[["lowexpressgene_pct"]][match(names(sample_lowexpressgene_pct),RNA_sample_info[[sample_id_col]])]<-sample_lowexpressgene_pct
  }

  #prepare clean RNA_sample_info and protein expressions
  clean_RNA_sample_info<-RNA_sample_info[RNA_sample_info[["IS_Sequenced"]]=="YES",]
  clean_protein_expressions<-protein_expressions[,clean_RNA_sample_info[[sample_id_col]]]
  if(remove_low_purity_sample){
    clean_RNA_sample_info<-clean_RNA_sample_info[clean_RNA_sample_info[["TumorPurity"]]>=low_purity_threshold,]
    clean_protein_expressions<-clean_protein_expressions[,clean_RNA_sample_info[[sample_id_col]]]
  }
  if(remove_sample_with_intense_lowexpressgene){
    clean_RNA_sample_info<-clean_RNA_sample_info[clean_RNA_sample_info[["lowexpressgene_pct"]]<outlier_lowexpressgene_pct_factor*average_lowexpressgene_pct,]
    clean_protein_expressions<-clean_protein_expressions[,clean_RNA_sample_info[[sample_id_col]]]
  }

  #remove low express gene if low express gene is present across defined pencentage of samples
  if(remove_lowexpressgene){
    clean_protein_expressions<-clean_protein_expressions[rowSums(clean_protein_expressions<lowexpression_threshold)<(sample_frequency_threshold*ncol(clean_protein_expressions)),]
  }

  assertthat::are_equal(colnames(clean_protein_expressions),clean_RNA_sample_info[[sample_id_col]])

  return(list(marked_ori_RNA_sample_info=RNA_sample_info,clean_RNA_sample_info=clean_RNA_sample_info,clean_protein_expressions=clean_protein_expressions))
}
