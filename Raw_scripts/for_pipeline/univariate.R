#BiocManager::install("ropls")
#BiocManager::install('MultiDataSet')
#devtools::install_github("SciDoPhenIA/phenomis")
#BiocManager::install("biosigner")

########################################################
# ote to self incoperate list of variables of interests#
##########################################################

library(phenomis)
library(biosigner)
library(ropls)
library(dendextend)
library(data.table)
library(dplyr)
library(stringr)

###############################
# variables #
##########################

input_folder <- 'F:/avans/stage MM/xcms_pipeline/_XCMS_default/batch_correction/Noice_reducted/'
Corrected_data_matrix <- "XCMS_default_batchcorrected_noice_reduced_matrix.tsv"
sample_meta_data <- "XCMS_default_batchcorrected_noice_reduced_sample_metadata.tsv"
variable_meta_data <- "XCMS_default_batchcorrected_noice_reduced_variable_metadata.tsv"


test <- c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman","limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[1]


variables_of_interest <- "gender,age,bmi"
correcting_data_set_according_variable <- "NULL"
factor_of_interest <- 'gender'
second_factor_of_interest <- 'age'

P_value_treshhold <- 0.05
max_features_output <- "NA"
graph_title <- ''
pre_fix_of_report <- ""



pixelsize1 <- 20
pixelsize2 <- 12


######### filtering data ##########
# keeps data who is bigger then these numbers#
cutoff_hotelPval <-  0.001
cutoff_missPval <-  0.001 
cutoffDecipval <- 0.001


######## heatmap ###########
cluster_groups_of_samples<- 6
cluster_groups_of_features<- 6

heat_map_statistics <- c("pearson", "kendall", "spearman")[1]
heat_map_algorithm <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
                        "1-cor", "1-abs(cor)")[7]






###############################
# correcting some dataformats #
###############################
if (correcting_data_set_according_variable == "NULL"){
  correcting_data_set_according_variable <- NULL
}

if (max_features_output == "NA"){
    max_features_output <- NA
}


try(variables_of_interest <- strsplit(variables_of_interest,",")[[1]])

items_to_filter <- c('blank', 'QC')


##############################################
sampleType <- "sample_type"

output_folder  <- paste0(input_folder,'/',test,'_',factor_of_interest,'/')
dir.create(output_folder, showWarnings = T)
setwd(output_folder)

numerical_factor_of_interest <- second_factor_of_interest

# create output name#
result <- paste0(test,'_',factor_of_interest)

# remove file to resolve bug
try(file.remove(paste0(output_folder,result,'.txt')))
#retrieve best info
path_Corrected_data_matrix <- paste0(input_folder,'/',Corrected_data_matrix)
path_sample_meta_data <- paste0(input_folder,'/',sample_meta_data)
path_variable_meta_data <- paste0(input_folder,'/',variable_meta_data)




#############################################
# changing datasets for better graph display#
#############################################
sample_metadata <- (read.table(path_sample_meta_data, sep = '\t', header = TRUE, row.names = 1))
data_matrix <- (read.table(path_Corrected_data_matrix, sep = '\t', header = TRUE, row.names = 1))


#######################################################################
# edit datamatrix 1/5 of lowest feature in dataset will be replace NA#
#######################################################################

for (i in 1:length(rownames(data_matrix))){
  data1 <- data_matrix[i,]
  data <- data1[!is.na(data1)]
  data <- min(data)-log10(5)
  data <- data1[is.na(data1)] <- data
  data_matrix[i,] <- data1
}










factor <- sample_metadata[, factor_of_interest]
sample_names <- rownames(sample_metadata)
eset_names <- list()
for (i in 1:length(factor)){
  string <- paste0(sample_names[[i]],"-",factor[[i]])
  eset_names[i] <- string
  
}


colnames(data_matrix) <- eset_names
rownames(sample_metadata) <- eset_names

datamatrix_path <- paste0(output_folder,Corrected_data_matrix,'_tmp_',factor_of_interest,'.tsv')
write.table(data_matrix, datamatrix_path, sep ='\t')

sample_metadata_path <- paste0(output_folder,sample_meta_data,'_tmp_',factor_of_interest,'.tsv')
write.table(sample_metadata, sample_metadata_path, sep ='\t')


line <- readLines(datamatrix_path,)
line[1] <- paste0('""\t',line[1])
writeLines(line,datamatrix_path,)


line <- readLines(sample_metadata_path,)
line[1] <- paste0('""\t',line[1])
writeLines(line,sample_metadata_path, )

##################
# reading in data#
##################


set <- phenomis::reading(NA,
                         files.ls = list(dataMatrix = file.path(datamatrix_path),
                                         sampleMetadata = file.path(sample_metadata_path),
                                         variableMetadata = file.path(path_variable_meta_data)))


png(filename = paste0(output_folder,result,'overview.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)



sacurine.eset <- phenomis::inspecting(set,
                                      pool_as_pool1.l = FALSE,
                                      pool_cv.n = 0.3,
                                      span.n = 1,
                                      sample_intensity.c = c("median", "mean")[1],
                                      title.c = "Inital Overview",
                                      plot_dims.l = TRUE,
                                      figure.c = c("none", "interactive", "myfile.pdf")[2],
                                      
                                      report.c = c("none", "interactive", "myfile.txt")[3])

dev.off()

#######################
# Filtering of samples#
#######################


if (!is.null(items_to_filter)){
  for (i in 1:length(items_to_filter)){
    set <- set[, Biobase::pData(set)[, sampleType] != items_to_filter[i]]
    
  }
  png(filename = paste0(output_folder,result,'overview_after_sample_filtering.png'),
      width = 960, height = 960, units = "px", pointsize = pixelsize1,
      bg = "white",  res = NA)
  set <- sacurine.eset <- phenomis::inspecting(set,
                                               pool_as_pool1.l = FALSE,
                                               pool_cv.n = 0.3,
                                               span.n = 1,
                                               sample_intensity.c = c("median", "mean")[1],
                                               title.c = "Inital Overview",
                                               plot_dims.l = TRUE,
                                               figure.c = c("none", "interactive", "myfile.pdf")[2],
                                               
                                               report.c = c("none", "interactive", "myfile.txt")[3])
  dev.off()
}


####################################################
#correcting dataset according variable information##
####################################################
if(!(is.null(correcting_data_set_according_variable))){
  Biobase::exprs(set) <- sweep(Biobase::exprs(set),
                               2,
                               Biobase::pData(set)[, correcting_data_set_according_variable],
                               "/")
  set <- phenomis::inspecting(set)
}

#####################################################
######### Filtering according data statistics########
#####################################################

set <- set[, Biobase::pData(set)[, "hotel_pval"] >= cutoff_hotelPval &
             Biobase::pData(set)[, "miss_pval"] >= cutoff_missPval &
             Biobase::pData(set)[, "deci_pval"] >= cutoffDecipval]

phenomis::inspecting(set)


######################
# Univariate Testing #
######################

png(filename = paste0(output_folder,result,'signif_chart.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize2,
    bg = "white",  res = NA,
)


univariate_set <- phenomis::hypotesting(
  set,
  test.c = test,
  factor_of_interest,
  adjust_thresh.n = P_value_treshhold,
  signif_maxprint.i = max_features_output,
  title.c = 'Significant levels of univariate testing',
  prefix.c = pre_fix_of_report,
  adjust.c = "none",
  report.c = paste0(result,'.txt'),
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2])




dev.off()


################################################
##################### PCA ######################
################################################

png(filename = paste0(output_folder,result,'_Global_PCA.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)
sacPca <- ropls::opls(univariate_set, info.txt = 'interactive',
                      y = NULL,
                      predI = NA,
                      orthoI = 0,
                      algoC = c("default", "nipals", "svd")[1],
                      crossvalI = 1,
                      log10L = FALSE,
                      permI = 7,
                      scaleC = c("none", "center", "pareto", "standard")[4],
                      subset = NULL,
                      plotSubC = 'Variance between Principal Components',
                      fig.pdfC = c("none", "interactive", "myfile.pdf")[2],
                      info.txtC = c("none", "interactive", "myfile.txt")[2],
                      printL = TRUE,
                      plotL = TRUE,
                      .sinkC = NULL,)

dev.off()

#update variables_of_interest

variables_of_interest <- variables_of_interest[variables_of_interest != factor_of_interest]


for (i in 1:length(variables_of_interest)){
  print(i)
  png(filename = paste0(output_folder,result,'_',variables_of_interest[i],'.png'),
      width = 960, height = 960, units = "px", pointsize = pixelsize1,
      bg = "white",  res = NA,
  )
  
  ropls::plot(sacPca,
              parAsColFcVn = Biobase::pData(univariate_set)[, variables_of_interest[i]],
              typeVc = "x-score",
              
              figure.c = paste0(result,'.pdf'),
              plotSubC = variables_of_interest[i])
  
  
  
  dev.off()
}








########################
##### heatmap ##########
########################



png(filename = paste0(output_folder,result,'_heatmap.png'),
    width = 1840, height = 1840, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)

sacurine.eset <- ropls::getEset(sacPca)
sacurine.eset <- phenomis::clustering(sacurine.eset, correl.c = heat_map_statistics,
                                      cex.vn = c(1,1),
                                      clusters.vi = c(cluster_groups_of_samples, cluster_groups_of_features),
                                      dissym.c = heat_map_algorithm ,
                                      agglo.c = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
                                                  "median", "centroid")[2],
                                      title.c = paste0('Clustering of Features and sampels ',heat_map_statistics))




dev.off()

###########################################
########### Supervised modeling############
########### (O)PLS(-DA) modeling###########
###########################################
png(filename = paste0(output_folder,result,'_Testing.png'),
    width = 1840, height = 1840, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)

sacPlsda <- ropls::opls(sacurine.eset, predI = 5, factor_of_interest, seedI = 123)
sacurine.eset <- ropls::getEset(sacPlsda)




dev.off()


##############################################################
# rewriting the output for filtering#
##############################################################
######## removing tmp files ############


phenomis::writing(sacurine.eset, dir.c = getwd(),overwrite.l = TRUE)


if (!is.null(items_to_filter)){
  library(plyr)
  metadata <- read.table(path_sample_meta_data, sep = '\t', header = TRUE, row.names = 1)
  metadata1 <- read.table(paste0(output_folder,"/sampleMetadata.tsv"), sep = '\t', header = TRUE, row.names = 1)
  colnames1 <- setdiff(colnames(metadata1), colnames(metadata))
  
  metadata <- metadata[  metadata$sample_type %in% items_to_filter, ]
  
  
  for (i in 1:length(colnames1)){
    metadata[[colnames1[i]]] <- NA
  }
  
  
  metadata2 <- rbind(metadata1, metadata)
  write.table(metadata2, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  file.remove(paste0(output_folder,"/sampleMetadata.tsv"))
  
  sample_meta_data1 <- read.table(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), header = TRUE, row.names = 1)
  rownames(sample_meta_data1) <- gsub('(.*)-\\w+', '\\1', rownames(sample_meta_data1))
  write.table(sample_meta_data1, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  
  
  file.copy(path_Corrected_data_matrix,paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv') )
  file.copy(paste0(output_folder,"/dataMatrix.tsv"), paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  file.remove(paste0(output_folder,"/dataMatrix.tsv"))
  
  line <- readLines(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  
} else if(is.null(items_to_filter)) {
  file.rename(paste0(output_folder,"/sampleMetadata.tsv"), paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  
  sample_meta_data1 <- read.table(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), header = TRUE, row.names = 1)
  rownames(sample_meta_data1) <- gsub('(.*)-\\w+', '\\1', rownames(sample_meta_data1))
  write.table(sample_meta_data1, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  
  file.rename(paste0(output_folder,"/dataMatrix.tsv"), paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  
  ####################
  # correcting output#
  ####################
  
  line <- readLines(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line,paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  line <- readLines(paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line,paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  
}





file.rename(paste0(output_folder,"/variableMetadata.tsv"), paste0(output_folder,variable_meta_data,'_',test,"_",factor_of_interest,'.tsv'))


file.remove(datamatrix_path)
file.remove(sample_metadata_path)

