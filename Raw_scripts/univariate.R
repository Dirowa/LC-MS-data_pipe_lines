#BiocManager::install("ropls")
#BiocManager::install('MultiDataSet')
#devtools::install_github("SciDoPhenIA/phenomis")




library( phenomis)
###############################
# variables #
##########################

input_folder <- 'F:/avans/stage MM/Sherloktest_data_2/peakpick_output/batch_correction_main_output'
output_folder <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output/batch_correction_main_output/"
final_output <-"F:/avans/stage MM/Sherloktest_data_2/significant_features"
  
dir.create(final_output, showWarnings = F)
setwd(output_folder)
#c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
#"limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter"
test <- "ttest"
factor_of_interest <- 'gender'
P_value_treshhold <- 0.20
max_features_output <- NA
graph_title <- NA
pre_fix_of_report <- ""


# create output name#
result <- paste0(test,'_',factor_of_interest)
#retrieve best info
path_Corrected_data_matrix <- paste0(input_folder,'/',"data_matrix.tsv")
path_sample_meta_data <- paste0(input_folder,'/',"sample_metadata.tsv")
path_variable_meta_data <- paste0(input_folder,'/',"variable_metadata.tsv")

# reading in data
set <- reading(NA,
                         files.ls = list(dataMatrix = file.path(path_Corrected_data_matrix),
                                         sampleMetadata = file.path(path_sample_meta_data),
                                         variableMetadata = file.path(path_variable_meta_data)))



#removing pooled and blanks


samples_Set <- set[, Biobase::pData(set)[, "sampleType"] != "QC"]
samples_Set <- samples_Set[, Biobase::pData(samples_Set)[, "sampleType"] != "blank"	]


univariate_set <- phenomis::hypotesting(
                                        samples_Set,
                                        test.c = test,
                                        factor_of_interest,
                                        adjust_thresh.n = P_value_treshhold,
                                        signif_maxprint.i = max_features_output,
                                        title.c = graph_title,
                                        prefix.c = pre_fix_of_report,
                                        figure.c = paste0(result,'.pdf'),
                                        report.c = paste0(result,'.txt'))


##############################################################
# rewriting the output files to be annotated by teh databases#
##############################################################

significant_features <- read.table(paste0(output_folder,result,'.txt'), sep = '')
variable_metadata <- read.table(path_variable_meta_data, sep = '')
rownames(variable_metadata) <- variable_metadata[,1]
variable_metadata <- variable_metadata[,-1]
colnames(variable_metadata)<- variable_metadata[1,]
variable_metadata <- variable_metadata[-1,]
variable_metadata <- select(variable_metadata, mzmed, rtmed,npeaks,blank,pool,sample)

#selecting which rows to get
significant_features <-cbind(variable_metadata[(rownames(significant_features)),],significant_features)


write.csv(significant_features, paste0(final_output,'/',result,".csv"))





