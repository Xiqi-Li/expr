context("Test functions to prepare clean RNA sample info and protein expressions")

expressions_file = system.file("extdata/expressions.csv",package = "exprClean")
RNA_sample_info_file = system.file("extdata/RNA_sample_info.csv",package = "exprClean")
protein_coding_ensemble2symbol_file = system.file("extdata/protein_coding_ensemble2symbol.csv",package = "exprClean")


RNA_sample_info<-read.csv(RNA_sample_info_file,header=T,stringsAsFactors = F,check.names = F)
expressions<-read.csv(expressions_file,header=T,stringsAsFactors = F,check.names = F)
protein_coding_ensemble2symbol<-read.csv(protein_coding_ensemble2symbol_file,header=T,stringsAsFactors = F,check.names = F)

test_that("prepare_clean_RNA_sample_info_and_protein_expressions can produce matching info and expression tables", {
  expr=expression(prepare_clean_RNA_sample_info_and_protein_expressions(RNA_sample_info=RNA_sample_info,sample_id_col = "Sample_ID",expressions = expressions,gene_id_col = "gene_id",
                                                                mark_library_size = T,adjusted_by_library_size = T,tolerant_library_size_factor = 1.1,
                                                                protein_coding_ensemble2symbol_table = protein_coding_ensemble2symbol,ens_id_col = "Gene_ID",symbol_col = "Gene_Symbol",
                                                                mark_purity = T,remove_low_purity_sample = T,low_purity_threshold = 0.3,
                                                                mark_lowexpressgene_pct = T,lowexpression_threshold = 1,remove_sample_with_intense_lowexpressgene = T,outlier_lowexpressgene_pct_factor = 1.5,
                                                                remove_lowexpressgene = T,sample_frequency_threshold = 0.5))
  expect_output(eval(expr),"all sample libraries are quite stable, no need to adjust them!" , fixed = TRUE)
  expect_output(eval(expr),"On average, every sample have  25.75076 % low express genes.")
  expect_output(eval(expr),"sample Cap2115-6-ID04 looks having too many low expression genes:  84.6735024284943 % respectively.")
  expect_output(eval(expr),"Merged dataset includes")
  results=eval(expr)
  expect_true(dim(results$clean_protein_expressions)[2]==dim(results$clean_RNA_sample_info)[1])
  expect_true(all(colnames(results$clean_protein_expressions)==results$clean_RNA_sample_info$Sample_ID))
})

test_that("prepare_clean_RNA_sample_info_and_protein_expressions can adjust expressions by library size", {
  expr=expression(prepare_clean_RNA_sample_info_and_protein_expressions(RNA_sample_info=RNA_sample_info,sample_id_col = "Sample_ID",expressions = expressions,gene_id_col = "gene_id",
                                                                        mark_library_size = T,adjusted_by_library_size = T,tolerant_library_size_factor = 1,
                                                                        protein_coding_ensemble2symbol_table = protein_coding_ensemble2symbol,ens_id_col = "Gene_ID",symbol_col = "Gene_Symbol",
                                                                        mark_purity = T,remove_low_purity_sample = T,low_purity_threshold = 0.3,
                                                                        mark_lowexpressgene_pct = T,lowexpression_threshold = 1,remove_sample_with_intense_lowexpressgene = T,outlier_lowexpressgene_pct_factor = 1.5,
                                                                        remove_lowexpressgene = T,sample_frequency_threshold = 0.5))

  results=eval(expr)
  expect_true(results$clean_protein_expressions[1,1]>40.49)
})

