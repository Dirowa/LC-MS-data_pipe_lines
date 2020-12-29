#BiocManager::install("ropls")
#BiocManager::install('MultiDataSet')
#devtools::install_github("SciDoPhenIA/phenomis")
#BiocManager::install("biosigner")

########################################################
# ote to self incoperate list of variables of interests#
##########################################################

library( phenomis)
library(biosigner)
###############################
# variables #
##########################

input_folder <- 'F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default/batch_correction'
Corrected_data_matrix <- "Data_matrix_xcms_default.tsv_batchcorrected.tsv"
sample_meta_data <- "sample_meta_data_XCMS_default.tsv_batchcorrected.tsv"
variable_meta_data <- "Variable_metaData_XCMS_default.tsv_batchcorrected.tsv"
#c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
#"limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter"
test <- "ttest"

heat_map_statistics = c("pearson", "kendall", "spearman")[1]
variables_of_interest <- c("age",'bmi', 'gender', 'osmolality')
correcting_data_set_according_variable <- NULL
factor_of_interest <- 'gender'
P_value_treshhold <- 0.05
max_features_output <- NA
graph_title <- NA
pre_fix_of_report <- ""




output_folder  <- paste0(input_folder,'/',test,'_',factor_of_interest,'/')
dir.create(output_folder, showWarnings = T)
setwd(output_folder)


# create output name#
result <- paste0(test,'_',factor_of_interest)

# remove file to resolve bug
file.remove(paste0(output_folder,result,'.txt'))
#retrieve best info
path_Corrected_data_matrix <- paste0(input_folder,'/',Corrected_data_matrix)
path_sample_meta_data <- paste0(input_folder,'/',sample_meta_data)
path_variable_meta_data <- paste0(input_folder,'/',variable_meta_data)

# reading in data
set <- phenomis::reading(NA,
                         files.ls = list(dataMatrix = file.path(path_Corrected_data_matrix),
                                         sampleMetadata = file.path(path_sample_meta_data),
                                         variableMetadata = file.path(path_variable_meta_data)))


#correcting dataset according variable information
if(!(is.null(correcting_data_set_according_variable))){
  Biobase::exprs(set) <- sweep(Biobase::exprs(set),
                                         2,
                                         Biobase::pData(set)[, correcting_data_set_according_variable],
                                         "/")

}

set <- phenomis::inspecting(set)
samples_Set <- set[, Biobase::pData(set)[, "sampleType"] != "QC"]
samples_Set <- samples_Set[, Biobase::pData(samples_Set)[, "sampleType"] != "blank"	]
set <- phenomis::inspecting(samples_Set)

png(filename = paste0(output_folder,result,'overview.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
sacurine.eset <- phenomis::inspecting(samples_Set)

dev.off()



png(filename = paste0(output_folder,result,'signif_chart.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)



try(univariate_set <- phenomis::hypotesting(
                                        sacurine.eset,
                                        test.c = test,
                                        factor_of_interest,
                                        adjust_thresh.n = P_value_treshhold,
                                        signif_maxprint.i = max_features_output,
                                        title.c = graph_title,
                                        prefix.c = pre_fix_of_report,
                                        adjust.c = "none",
                                        report.c = paste0(result,'.txt'))
)


dev.off()
##########################
# Generating plots########
##########################


#update variables_of_interest

#variables_of_interest <- variables_of_interest[variables_of_interest != factor_of_interest]


png(filename = paste0(output_folder,result,'_Global_PCA.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
sacPca <- ropls::opls(univariate_set, info.txt = 'interactive')

dev.off()




for (i in 1:length(variables_of_interest)){
  
  png(filename = paste0(output_folder,result,'_',variables_of_interest[i],'.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  
  ropls::plot(sacPca,
              parAsColFcVn = Biobase::pData(univariate_set)[, variables_of_interest[i]],
              typeVc = "x-score",
              figure.c = paste0(result,'.pdf'),
              plotSubC = variables_of_interest[i])
  
  
  
  dev.off()
}

##### heatmap ##########
png(filename = paste0(output_folder,result,'_heatmap.png'),
    width = 1840, height = 1840, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)

sacurine.eset <- ropls::getEset(sacPca)
sacurine.eset <- phenomis::clustering(sacurine.eset, correl.c = heat_map_statistics,
                                      clusters.vi = c(5, 3))

dev.off()


## Supervised modeling
### (O)PLS(-DA) modeling
png(filename = paste0(output_folder,result,'_Testing.png'),
    width = 1840, height = 1840, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
sacPlsda <- ropls::opls(sacurine.eset, factor_of_interest)
sacurine.eset <- ropls::getEset(sacPlsda)

sacurine.biosign <- biosigner::biosign(sacurine.eset, factor_of_interest, seedI = 123)
sacurine.eset <- biosigner::getEset(sacurine.biosign)


phenomis::writing(sacurine.eset, dir.c = getwd(),overwrite.l = TRUE)
dev.off()

##############################################################
# rewriting the output files to be annotated by teh databases#
##############################################################

file.rename(paste0(output_folder,"/dataMatrix.tsv"), paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
file.rename(paste0(output_folder,"/sampleMetadata.tsv"), paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
file.rename(paste0(output_folder,"/variableMetadata.tsv"), paste0(output_folder,variable_meta_data,'_',test,"_",factor_of_interest,'.tsv'))


