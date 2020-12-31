########################################################
# installing softawre#############
##########################################
#BiocManager::install("ProteoMM")
#library(ProteoMM)
#data("hs_peptides") # loads variable hs_peptides
#dim(hs_peptides)  # 695 x 13
#remotes::install_github("ricoderks/Rcpm")
#library(Rcpm)
#test_check("Rcpm")

#metanr_packages <- function(){
#  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn")
#  list_installed <- installed.packages()
#  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
#  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#    BiocManager::install(new_pkgs)
#    print(c(new_pkgs, " packages added..."))
#  }
#  
#  if((length(new_pkgs)<1)){
#    print("No new packages added...")
#  }
#}
#metanr_packages()




#install.packages("pacman")

library(pacman)
library("impute")
library("pcaMethods")
library("globaltest")
library("GlobalAncova")
library("Rgraphviz")
library("preprocessCore")
library("genefilter")
library("SSPA")
library("sva")
library("limma")
library("KEGGgraph")

library("BiocParallel")
library("MSnbase")
library("multtest")
library("RBGL")
library("fgsea")
library(MetaboAnalystR)
library(dplyr)


library(devtools)
library( phenomis)


# Step 2: Install MetaboAnalystR with documentation
#devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T, force = TRUE)

########################
# important variables #
#######################
# do not change#

input_folder <- 'F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default'


Data_matrix_xcms <- "Data_matrix_xcms_default.tsv"
samplemetadata <- "sample_meta_data_XCMS_default.tsv"
variable_metadata <- "Variable_metaData_XCMS_default.tsv"

final_output_folder <- "batch_correction"
# columns names in dataframe of metadata
sample_name <- 'Sample_name'
class <- 'sampleType'
batch <- 'batch'
injectionOrder <- 'injectionOrder'

#of class (what is the pooled sample code ( QC, pooled etc.))
QC <- 'pool'
sample_qc_part_name <- 'QC'


Filter_on_RSD <- TRUE
RSD_treshhold <- 25



# setting up some directories#
final_output_folder <- paste0(input_folder,'/',final_output_folder)
dir.create(final_output_folder, showWarnings = F)
setwd(final_output_folder)

variable_metadata_naam <- variable_metadata





#################################################################
# IMPORT DATA SET#
###############################################

path_to_Data_matrix_xcms <- paste0(input_folder,'/',Data_matrix_xcms )
path_to_samplemetadata <- paste0(input_folder,'/',samplemetadata )
Path_to_variable_metadata <- paste0(input_folder,'/',variable_metadata )
# the sample metadata 
# the sample metadata is required to have these column names:  sample_name, class, batch, injectionOrder
metadata <- read.csv(file = path_to_samplemetadata, header = FALSE, sep = '\t')

#standartd intensity list retrieved from XCMS 
intensities <- read.csv(file = path_to_Data_matrix_xcms, header = FALSE, sep = '\t')



####################
# correcting format#
####################
# problem with XCMS is that the headers is shifted. this has to be corrected 
# switching rows and columns of metadata
metadata1 <- as.data.frame(t(metadata))
metadata1['V1'] <- c("Sample_name", head(metadata1['V1'], dim(metadata1)[1] - 1)[[1]])
metadata <- as.data.frame(t(metadata1))
# same counts for the data_matrix
intensities <- as.data.frame(t(intensities))
intensities['V1'] <- c("Sample_name", head(intensities['V1'], dim(intensities)[1] - 1)[[1]])
intensities <- as.data.frame(t(intensities))
#editing the colnames
# add rowname variavle, sorting the dataframe by samplename
rownames(metadata1)<- metadata1[,1]
rownames(intensities) <- intensities[,1]
intensities <- intensities[,-1]
metadata1 <- metadata1[,-1]
# order the dataframes
intensities <- intensities[,order((intensities[1,]))]
metadata1 <- metadata1[,order((metadata1[1,]))]
#add colnames, the colnames is taken from another table but this is to prevent typos
#thereby the files are directly generated by XCMS so no error should be there
colnames(metadata1)<- intensities[1,]
colnames(intensities) <- intensities[1,]
#remove first kolumn data, no need for it
metadata1 <- metadata1[-1,]
intensities <- intensities[-1,]






###############
# COMBINING FOR OTHER BATCH CORRECTIONS#
###############\
identical(names(intensities), names( metadata1))
combined <- rbind(metadata1[1,], intensities)
#combined[is.na(combined)] <- 0
write.csv(combined, 'BatchCorrect_import_data_1')

#######################
# editing sample data#
#####################
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,]
metadata <- select(metadata, sample_name, class, batch, injectionOrder)
write.csv(metadata, 'batchCorrect_import_metadata')

rownames(metadata) <- metadata[,1]
metadata <- t(metadata)

names <- metadata[ , grepl( sample_qc_part_name , colnames( metadata ) ) ]


####################################
# Generate overview of current data#
####################################
rm(mSet)
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, 'BatchCorrect_import_data_1', "col", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
if (Filter_on_RSD == TRUE){
mSet<-FilterVariable(mSet, "iqr", "F", RSD_treshhold)
}
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ref = names[1,])

normalized_set <- mSet[["dataSet"]][["norm"]]
ordered_normalized_set1 <- normalized_set[order(row.names(normalized_set)), ]




#####################
# ordering datasets #
#####################
meta_data1 <- read.csv('batchCorrect_import_metadata')
new_normalized_set1 <- cbind(meta_data1[,2:5], ordered_normalized_set1);
new_normalized_set1 <- new_normalized_set1[,-1]

new_normalized_set1$sampleType <- as.character(new_normalized_set1$sampleType)
new_normalized_set1$sampleType[new_normalized_set1$sampleType == QC] <- "QC"

write.csv(new_normalized_set1, "new_normalized_set.csv")



#############################
# correcting for loes method#
#############################

metadata2 <- as.data.frame(t(metadata1))
metadata2$sampleType[metadata2$sampleType == QC] <- "QC"

row_itensities <- rownames(ordered_normalized_set1)
Feature_table <- as.data.frame(read.csv(Path_to_variable_metadata, sep = '\t'))

#### possible correcting for removed features#####

a <-rownames(Feature_table)
b <-colnames(ordered_normalized_set1)
diffrences <- setdiff(a,b)
diffrences <- diffrences[diffrences != ""]

#Remove the differences
for (i in 1:length(diffrences)){
  Feature_table<-Feature_table[!(rownames(Feature_table)==diffrences[[i]]),]
}

ordered_normalized_set1 <- t(ordered_normalized_set1)

loes_file_dataMatrix.tsv <- paste0(final_output_folder,'/loes_dataMatrix.tsv')
loes_file_sampleMetadata.tsv <- paste0(final_output_folder,'/loes_sampleMetadata.tsv')
loes_file_variableMetadata.tsv <- paste0(final_output_folder,'/loes_variableMetadata.tsv')



write.table(ordered_normalized_set1, loes_file_dataMatrix.tsv, sep ='\t')
write.table(metadata2, loes_file_sampleMetadata.tsv, sep ='\t')
write.table(Feature_table, loes_file_variableMetadata.tsv, sep ='\t')





# rewiting becouse it needs "" \t for, cant be done by playing with dataframes
line <- readLines(loes_file_dataMatrix.tsv,)
line[1] <- paste0('""\t',line[1])
writeLines(line,loes_file_dataMatrix.tsv )

line <- readLines(loes_file_sampleMetadata.tsv,)
line[1] <- paste0('""\t',line[1])
writeLines(line,loes_file_sampleMetadata.tsv )

line <- readLines(loes_file_variableMetadata.tsv,)
line[1] <- paste0('""\t',line[1])
writeLines(line,loes_file_variableMetadata.tsv )



################
# calculate PCA#
################
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)



#############################
# performing batchcorrection#
#############################
rm(mSet)
# there are more but those require intenalstandards
methods<- c("Loess", "Combat", "WaveICA","EigenMS","QC_RLSC","ANCOVA")
#methods<- c("Combat","EigenMS","QC_RLSC","ANCOVA")
counter <- 1
items <- list()


for (i in 1:length(methods))  {
  counter <- counter + 1
  if (methods[i] == "Loess"){
    eset <- reading(NA,
                    files.ls = list(dataMatrix = loes_file_dataMatrix.tsv,
                                    sampleMetadata = loes_file_sampleMetadata.tsv,
                                    variableMetadata = loes_file_variableMetadata.tsv))
    
    corrcted_set <- phenomis::correcting(
      eset,
      reference.c = "QC",
      col_batch.c = batch,
      col_injectionOrder.c = injectionOrder,
      col_sampleType.c = class,
      span.n = 1,
      title.c = NA,
      figure.c = c("none", "interactive", "Loess_correction.pdf")[3],
      report.c = c("none", "interactive", "Loess_correction.pdf")[3]
    )
    #corrcted_set <- inspecting(corrcted_set)
    phenomis::writing(corrcted_set, dir.c = final_output_folder, prefix.c = 'loessALL_BatchCorrected_data',
                      overwrite.l = TRUE)
    #filter for only QC sets

    
    file.remove(paste0(final_output_folder,'/',"loessALL_BatchCorrected_data_sampleMetadata.tsv"))
    file.remove(paste0(final_output_folder,'/',"loessALL_BatchCorrected_data_variableMetadata.tsv"))

    #items[counter] <- 'loessQC_BatchCorrected_data.tsv'
  }
  else {
    rm(mSet)
    mSet <- InitDataObjects("pktable", "utils", FALSE)
    mSet <- Read.BatchDataTB(mSet, "new_normalized_set.csv", "row")
    try(mSet <- PerformBatchCorrection(mSetObj = mSet, imgName = methods[[i]], Method = methods[[i]]), silent = F) 
    #try(info <- (mSet$dataSet$interbatch_dis))
    safe_name <- paste0("BatchCorrected_data",'_',methods[[i]],'.csv')
    (file.rename("MetaboAnalyst_batch_data.csv", safe_name))
    
  }
}



if (counter == 1){
  print('NO batch correction worked')
}



###############################################
# calculate average variance across QC samples#
###############################################
items <- list.files(path = final_output_folder, pattern = "BatchCorrected_data", all.files = TRUE,
                    full.names = FALSE, recursive = FALSE,
                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

items[[(length(items) + 1)]] <- "new_normalized_set.csv"
items
means <- list()
means2 <- list()

for (i in 1:length(items)){
  print(i)
  items <- Filter(Negate(is.null), items)
  try(rm(variance_dataframe))
  print(items[[i]])
  if (items[[i]] == "new_normalized_set.csv"){
    variance_dataframe <- read.csv(paste0(final_output_folder,'/','new_normalized_set.csv'), header = TRUE, sep = ',')
    variance_dataframe1 <- variance_dataframe[,-1:-4]
    means2[[i]] <- mean(sapply(variance_dataframe1,  var))
    
    variance_dataframe <- variance_dataframe[grep('QC', variance_dataframe$sampleType),]
    variance_dataframe <- variance_dataframe[,-1:-4]
    #variance_dataframe[is.na(variance_dataframe)] <- 0
    print(mean(sapply(variance_dataframe,  var)))
    means[[i]] <- mean(sapply(variance_dataframe,  var))
    rm(variance_dataframe)
    rm(variance_dataframe1)
  }
  
  
  else if (as.character(items[[i]]) == 'loessALL_BatchCorrected_data_dataMatrix.tsv'){
    
    variance_dataframe <- read.csv(paste0(final_output_folder,'/',items[[i]]), header = TRUE, row.names(1), sep = '\t')
    rownames(variance_dataframe) <- variance_dataframe[,1]
    variance_dataframe <- variance_dataframe[,-1]
    variance_dataframe1 <- as.data.frame(t(variance_dataframe))
    
    vars <- list()
    for (ii in 1: length(names( variance_dataframe ))){
      variance_dataframe_subset <- variance_dataframe1[,ii]
      variance_dataframe_subset <- variance_dataframe_subset[!is.na(variance_dataframe_subset)]
      vars[[ii]] <- var(variance_dataframe_subset)
      
    }
    
    var <- (unlist(vars))
    var <- var[!is.na(var)]
    print(mean(var))
    means2[[i]] <- mean(var)
    
    
    
    variance_dataframe <- variance_dataframe[ , grepl( sample_qc_part_name , names( variance_dataframe ) ) ]
    variance_dataframe <- as.data.frame(t(variance_dataframe))
    
    vars <- list()
    for (ii in 1: length(names( variance_dataframe ))){
      variance_dataframe_subset <- variance_dataframe[,ii]
      variance_dataframe_subset <- variance_dataframe_subset[!is.na(variance_dataframe_subset)]
      vars[[ii]] <- var(variance_dataframe_subset)
      
    }
    
    var <- (unlist(vars))
    var <- var[!is.na(var)]
    print(mean(var))
    means[[i]] <- mean(var)
    
    
    
    rm(variance_dataframe)
    rm(variance_dataframe1)
    rm(var)
    rm(vars)
    
  } 
  else {
    variance_dataframe <- read.csv(paste0(final_output_folder,'/',items[[i]]), header = TRUE, sep = ',')
    variance_dataframe1 <- variance_dataframe[,-1:-3]
    var <- sapply(variance_dataframe1,  var)
    means2[[i]] <- mean(var)
    
    
    variance_dataframe <- variance_dataframe[grep("QC", variance_dataframe$CLASS),]
    variance_dataframe <- variance_dataframe[,-1:-3]
    
    #variance_dataframe[is.na(variance_dataframe)] <- 0
    var <- sapply(variance_dataframe,  var)
    print(mean(var))
    means[[i]] <- mean(var)
    rm(variance_dataframe)
  }
}


## ordering output
means1 <- do.call(rbind, Map(data.frame,mean_var_total_log10=means2, mean_var_QC_log10=means, file=items))
means <- means1[order(means1$mean_var_QC_log10),]
print(means)
#write down output of benchmarking
write.csv(means, file = 'Benchmarking_batch_correction.csv')


#selecting Best dataset


if (sub('.*\\.', '', means[1,3]) == "tsv"){
  try(means <- gsub('QC', 'ALL',means[1,3]))
  best_batch_corrected <- (read.csv(means, sep = '\t'))
  best_batch_corrected <- as.data.frame((best_batch_corrected))
  rownames(best_batch_corrected) <- best_batch_corrected[,1]
  best_batch_corrected <- best_batch_corrected[,-1]
  
  write.table(best_batch_corrected,'Best_batchcorrected.tsv', sep = '\t')
  
  
} else if ((sub('.*\\.', '', means[1,3]) == "csv") ) {
  best_batch_corrected <- read.csv(means[1,3])
  best_batch_corrected <- as.data.frame((best_batch_corrected))
  rownames(best_batch_corrected) <- best_batch_corrected[,1]
  best_batch_corrected <- t(best_batch_corrected[,-1:-3])
  
  write.table(best_batch_corrected,'Best_batchcorrected.tsv', sep = '\t')
  
}




##################################################
# perform updated PCA with best batch corrections#
##################################################

rm(mSet)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "Best_batchcorrected.tsv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_Best_corrected_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca__Best_corrected_scree_0_", "png", 72, width=NA, 5)


##################################################
# Redisgn output as input for statistical testing#
##################################################
# samplemetada
file.copy(from=loes_file_sampleMetadata.tsv, to=paste0(final_output_folder,'/',samplemetadata,'_batchcorrected.tsv' )
          , 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)



#datamatrix and variable metadata needs to be corrected due filtering
variable_metadata <-  read.table(paste0(final_output_folder,'/loes_variableMetadata.tsv'), sep = '\t')

data_matrix <- read.table(paste0(final_output_folder,'/Best_batchcorrected.tsv'), sep = '\t')

#retrieve differences#
a <-variable_metadata[,1]
b <-rownames(data_matrix)
diffrences <- setdiff(a,b)
diffrences <- diffrences[diffrences != ""]

#Remove the differences
for (i in 1:length(diffrences)){
  variable_metadata<-variable_metadata[!(variable_metadata$V1==diffrences[[i]]),]
}

# make the files pretty enough for statistical analyis
colnames(variable_metadata) <- variable_metadata[1,]
variable_metadata <- variable_metadata[-1,]
TMP <- variable_metadata[,1]
variable_metadata <- variable_metadata[,-1]
variable_metadata <- mutate_all(variable_metadata, function(x) as.numeric(as.character(x)))
rownames(variable_metadata) <- TMP

#rownames(data_matrix) <- data_matrix[,1]
#data_matrix <- data_matrix[,-1]


output_datamatrix <- paste0(final_output_folder,'/',Data_matrix_xcms,'_batchcorrected.tsv')
output_variable_metadata <- paste0(final_output_folder,'/',variable_metadata_naam,'_batchcorrected.tsv')

write.table(data_matrix, output_datamatrix, sep ='\t')
write.table(variable_metadata, output_variable_metadata, sep ='\t')


# rewiting becouse it needs "" \t for, cant be done by playing with dataframes
line <- readLines(output_datamatrix,)
line[1] <- paste0('""\t',line[1])
writeLines(line,output_datamatrix,)


line <- readLines(output_variable_metadata,)
line[1] <- paste0('""\t',line[1])
writeLines(line,output_variable_metadata, )
