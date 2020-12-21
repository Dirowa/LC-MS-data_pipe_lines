########################################################
# installing softawre#############
##########################################
BiocManager::install("ProteoMM")
library(ProteoMM)
data("hs_peptides") # loads variable hs_peptides
dim(hs_peptides)  # 695 x 13
remotes::install_github("ricoderks/Rcpm")
library(Rcpm)
test_check("Rcpm")

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()

install.packages("pacman")

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

install.packages("devtools")
library(devtools)

# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T, force = TRUE)

########################
# important variables #
#######################

path_to_metadata <- "F:/avans/stage MM/peakpicking 2/sample_meta_data_csv_XCMSdefault.txt"
path_to_feature_intensity_table <- "F:/avans/stage MM/peakpicking 2/feature_intensity_XCMSdefault.txt"

# columns names in dataframe of metadata
metadata <- 'metadata'
sample_name <- 'sample_name'
class <- 'class'
batch <- 'batch'
injectionOrder <- 'injectionOrder'

#of class (what is the pooled sample code ( QC, pooled etc.))
QC <- 'QC'





#################################################################
# IMPORT DATA SET#
###############################################
# the sample metadata 
# the sample metadata is required to have these column names:  sample_name, class, batch, injectionOrder
metadata <- read.csv(file = path_to_metadata, header = TRUE, sep = ',')

#standartd intensity list retrieved from XCMS 
intensities <- read.csv(file = path_to_feature_intensity_table, header = FALSE, sep = ',')

####################
# correcting format#
####################

### add sample_name to intensities#
intensities[1,1] <- 'sample_name'

# switching rows and columns of metadata
metadata1 <- as.data.frame(t(metadata))
# add rowname variavle, sorting the dataframe by samplename
metadata1 <- add_rownames(metadata1, var = "name")
metadata1 <- metadata1[-1,]
metadata2 <- metadata1[,-1]
metadata2 <- metadata2[,order((metadata2[1,]))]
rownames <- metadata1[,1]
metadata2 <- cbind(names = rownames, metadata2)
rm(metadata1)
colnames(metadata2) <- metadata2[1,]

#####################
# intensities #
###################

intensities1 <- intensities[,-1]
rownames(intensities1) <- intensities[,1]
intensities1 <- intensities1[,order((intensities1[1,]))]
intensities1 <- cbind(names = row.names(intensities1),intensities1)
intensities1[1,] <- gsub("\\..*","",intensities1[1,])
colnames(intensities1) <- intensities1[1,]
###############
# combining#
###############\

combined <- rbind(metadata2[2,], intensities1[-1,])
combined <- rbind(intensities1[1,],combined)
combined <- combined[-1,]
rownames(combined) <- combined[,1]
combined <- combined[,-1]

###########################
# change all NA data in 0##
###########################
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
new_normalized_set <- cbind(meta_data[,3:5], ordered_normalized_set);
write.csv(new_normalized_set, "new_normalized_set.csv")

################
# calculate PCA#
################
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
rm(mSet)



#############################
# performing batchcorrection#
#############################
mSet <- InitDataObjects("pktable", "utils", FALSE)
mSet <- Read.BatchDataTB(mSet, "new_normalized_set.csv", "row")
methods<- c("Combat", "WaveICA","EigenMS","QC_RLSC","ANCOVA","RUV_random","RUV_2","RUV_s","RUV_r","RUV_g","NOMIS","CCMN")
counter <- 0
items <- list()


# try function becouse not all datasets accepts the other functions (to many arguments for the LOG2 function)
for (i in 1:length(methods))  {
  try(mSet <- PerformBatchCorrection(mSetObj = mSet, imgName = methods[[i]], Method = methods[[i]]), silent = T)
  (info <- (mSet$dataSet$interbatch_dis))
  counter <- counter + 1
  safe_name <- paste0("MetaboAnalyst_batch_data.csv",'_',methods[[i]],'.csv')
  file.rename("MetaboAnalyst_batch_data.csv", safe_name)
  items[[counter]] <- safe_name
}

items[[(counter + 1 )]] <- "new_normalized_set.csv"
if (counter == 0){
  print('NO batch correction worked')
}



###############################################
# calculate average variance across QC samples#
###############################################

means <- list()
CLASS <- (toupper(class))
for (i in 1:length(items)){
  if (items[[i]] == 'new_normalized_set.csv'){
    variance_dataframe <- read.csv(items[[i]], header = TRUE, sep = ',')
    variance_dataframe <- variance_dataframe[grep(QC, variance_dataframe$class),]
    variance_dataframe <- variance_dataframe[,-1:-4]
    variance_dataframe <- sapply(variance_dataframe, var)
    means[[i]] <- mean(variance_dataframe)
  }
  else {
    print(items[[i]])
    variance_dataframe <- read.csv(items[[i]], header = TRUE, sep = ',')
    variance_dataframe <- variance_dataframe[grep(QC, variance_dataframe$CLASS),]
    variance_dataframe <- variance_dataframe[,-1:-3]
    variance_dataframe <- sapply(variance_dataframe, var)
    means[[i]] <- mean(variance_dataframe)
  }
}

means <- do.call(rbind, Map(data.frame, mean=means, file=items))
means <- means[order(means$mean),]
means
#write down output of benchmarking
write.csv(means, file = 'Benchmarking_batch_correction.csv')



#selecting Best dataset
best_batch_corrected <- read.csv(means[1,2])
write.csv(best_batch_corrected,'Best_batchcorrected.csv')



data_corrected_new <- best_batch_corrected[,-3]
write.csv(data_corrected_new,file = "MetaboAnalyst_batch_data_stats.csv",row.names =  F)

##################################################
# perform updated PCA with best batch corrections#
##################################################

rm(mSet)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "MetaboAnalyst_batch_data_stats.csv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_Best_corrected_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca__Best_corrected_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_Best_corrected__score2d_0_", "png", 72, width=NA,  1,2,0.95,0,0)

