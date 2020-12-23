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




library(colSu)
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

path_to_Data_matrix_xcms <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output/peakpick_outputData_matrix_xcmsdefault.tsv"
path_to_samplemetadata <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output/peakpick_outputsample_meta_data_XCMSdefault.tsv"
Path_to_variable_metadata <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output/peakpick_outputVariable_metaDatadefault.tsv"

output <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output"
# columns names in dataframe of metadata
sample_name <- 'X'
class <- 'sampleType'
batch <- 'batch'
injectionOrder <- 'injectionOrder'

#of class (what is the pooled sample code ( QC, pooled etc.))
QC <- 'pool'

setwd(output)




#################################################################
# IMPORT DATA SET#
###############################################
# the sample metadata 
# the sample metadata is required to have these column names:  sample_name, class, batch, injectionOrder
metadata <- read.csv(file = path_to_samplemetadata, header = TRUE, sep = '\t')

#standartd intensity list retrieved from XCMS 
intensities <- read.csv(file = path_to_Data_matrix_xcms, header = FALSE, sep = '\t')

####################
# correcting format#
####################

# switching rows and columns of metadata
metadata1 <- as.data.frame(t(metadata))
# add rowname variavle, sorting the dataframe by samplename
colnames(metadata1)<- metadata1[1,]
metadata1 <- metadata1[,order((metadata1[1,]))]
metadata1 <- metadata1[-1,]

#####################
# intensities #
###################

colnames(intensities) <- intensities[1,]
intensities <- intensities[,order((intensities[1,]))]
intensities <- intensities[-1,]
rownames(intensities) <- intensities[,1]
intensities <- intensities[,-1]


###############
# combining#
###############\

combined <- rbind(metadata1[2,], intensities)
combined[is.na(combined)] <- 0
write.csv(combined, 'BatchCorrect_import_data_1')

#######################
# editing sample data#
#####################
metadata <- select(metadata, sample_name, class, batch, injectionOrder)
write.csv(metadata, 'batchCorrect_import_metadata')


####################################
# Generate overview of current data#
####################################

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, 'BatchCorrect_import_data_1', "col", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-FilterVariable(mSet, "iqr", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

normalized_set <- mSet[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]

#####################
# ordering datasets #
#####################
meta_data <- read.csv('batchCorrect_import_metadata')
new_normalized_set <- cbind(meta_data[,2:5], ordered_normalized_set);
new_normalized_set <- new_normalized_set[,-1]

new_normalized_set$sampleType <- as.character(new_normalized_set$sampleType)
new_normalized_set$sampleType[new_normalized_set$sampleType == QC] <- "QC"

write.csv(new_normalized_set, "new_normalized_set.csv")

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
#mSet <- InitDataObjects("pktable", "utils", FALSE)
#mSet <- Read.BatchDataTB(mSet, "new_normalized_set.csv", "row")
methods<- c("Loess", "Combat", "WaveICA","EigenMS","QC_RLSC","ANCOVA")
counter <- 0
items <- list()


for (i in 1:length(methods))  {
  counter <- counter + 1
  if (methods[i] == "Loess"){
    eset <- reading(NA,
                             files.ls = list(dataMatrix = file.path(path_to_Data_matrix_xcms),
                                             sampleMetadata = file.path(path_to_samplemetadata),
                                             variableMetadata = file.path(Path_to_variable_metadata)))
    
    corrcted_set <- correcting(
      eset,
      reference.c = QC,
      col_batch.c = batch,
      col_injectionOrder.c = injectionOrder,
      col_sampleType.c = class,
      span.n = 1,
      title.c = NA,
      figure.c = c("none", "interactive", "Loess_correction.pdf")[3],
      report.c = c("none", "interactive", "Loess_correction.pdf")[3]
    )
    corrcted_set <- inspecting(corrcted_set)
    corrcted_set <- phenomis::transforming(corrcted_set, method.c = "log10")
    phenomis::writing(corrcted_set, dir.c = output, prefix.c = 'loessALL',
                      overwrite.l = TRUE)
    #filter for only QC sets
    corrcted_set <- corrcted_set[, Biobase::pData(corrcted_set)[, "sampleType"] == QC]
    phenomis::writing(corrcted_set, dir.c = output, prefix.c = 'loessQC_BatchCorrected_data',
                      overwrite.l = TRUE)
    
    file.remove(paste0(output,'/',"loessQC_BatchCorrected_data_sampleMetadata.tsv"))
    file.remove(paste0(output,'/',"loessQC_BatchCorrected_data_variableMetadata.tsv"))
    
  }
  else {
  rm(mSet)
  mSet <- InitDataObjects("pktable", "utils", FALSE)
  mSet <- Read.BatchDataTB(mSet, "new_normalized_set.csv", "row")
  try(mSet <- PerformBatchCorrection(mSetObj = mSet, imgName = methods[[i]], Method = methods[[i]]), silent = T) 
  #try(info <- (mSet$dataSet$interbatch_dis))
  safe_name <- paste0("BatchCorrected_data",'_',methods[[i]],'.csv')
  (file.rename("MetaboAnalyst_batch_data.csv", safe_name))
  
  }
}



if (counter == 0){
  print('NO batch correction worked')
}



###############################################
# calculate average variance across QC samples#
###############################################
items <- list.files(path = output, pattern = "BatchCorrected_data", all.files = FALSE,
                    full.names = FALSE, recursive = FALSE,
                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

items[[(length(items) + 1)]] <- "new_normalized_set.csv"


items
means <- list()
for (i in 1:length(items)){
  items <- Filter(Negate(is.null), items)
  try(rm(variance_dataframe))
  print(items[[i]])
  if (items[[i]] == "new_normalized_set.csv"){
    variance_dataframe <- read.csv(paste0(output,'/','new_normalized_set.csv'), header = TRUE, sep = ',')
    variance_dataframe <- variance_dataframe[grep('QC', variance_dataframe$sampleType),]
    variance_dataframe <- variance_dataframe[,-1:-4]
    #variance_dataframe[is.na(variance_dataframe)] <- 0
    means[[i]] <- mean(rowVars(variance_dataframe))
    rm(variance_dataframe)

  }
 
  
  else if (as.character(items[[i]]) == 'loessQC_BatchCorrected_data_dataMatrix.tsv'){

    variance_dataframe <- read.csv(paste0(output,'/',"loessQC_BatchCorrected_data_dataMatrix.tsv"), header = FALSE, sep = '\t')
    variance_dataframe <- as.data.frame(t(variance_dataframe))
    variance_dataframe <- variance_dataframe[,-1]
    variance_dataframe <- variance_dataframe[-1,]
    variance_dataframe <- na.exclude(variance_dataframe)
    variance_dataframe <- sapply(variance_dataframe,  var)
    means[[i]] <- mean(variance_dataframe)
    rm(variance_dataframe)
  }
  else {
    variance_dataframe <- read.csv(paste0(output,'/',items[[i]]), header = TRUE, sep = ',')
    variance_dataframe <- variance_dataframe[grep("QC", variance_dataframe$CLASS),]
    variance_dataframe <- variance_dataframe[,-1:-3]
    #variance_dataframe[is.na(variance_dataframe)] <- 0
    variance_dataframe <- sapply(variance_dataframe,  var)
    means[[i]] <- mean(variance_dataframe)
    rm(variance_dataframe)
  }
}


## ordering output
means1 <- do.call(rbind, Map(data.frame, mean=means, file=items))
means <- means1[order(means1$mean),]
print(means)
#write down output of benchmarking
write.csv(means, file = 'Benchmarking_batch_correction.csv')


#selecting Best dataset


if (sub('.*\\.', '', means[1,2]) == "tsv"){
  try(means <- gsub('QC', 'ALL',means[1,2]))
  best_batch_corrected <- (read.csv(means, sep = '\t'))
  best_batch_corrected <- as.data.frame((best_batch_corrected))
  rownames(best_batch_corrected) <- best_batch_corrected[,1]
  best_batch_corrected <- best_batch_corrected[,-1]
  
  write.csv(best_batch_corrected,'Best_batchcorrected1.csv')
  
  
} else if ((sub('.*\\.', '', means[1,2]) == "csv") ) {
  best_batch_corrected <- read.csv(means[1,2])
  best_batch_corrected <- as.data.frame((best_batch_corrected))
  rownames(best_batch_corrected) <- best_batch_corrected[,1]
  best_batch_corrected <- t(best_batch_corrected[,-1:-3])
  
  write.csv(best_batch_corrected,'Best_batchcorrected.csv')

}


##################################################
# perform updated PCA with best batch corrections#
##################################################

rm(mSet)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "Best_batchcorrected.csv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_Best_corrected_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca__Best_corrected_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_Best_corrected__score2d_0_", "png", 72, width=NA,  1,2,0.95,0,1)


