########################################################
# installing softawre#############
##########################################


BiocManager::install("ProteoMM")
library(ProteoMM)

data("hs_peptides") # loads variable hs_peptides
dim(hs_peptides)  # 695 x 13


remotes::install_github("ricoderks/Rcpm")
library(Rcpm)
??Rcpm-package


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

#################################################################
# IMPORT DATA SET#
###############################################
# the sample metadata 
# the sample metadata is required to have these column names:  sample_name, class, batch, injectionOrder
metadata <- read.csv(file = "F:/avans/stage MM/peakpicking 2/sample_meta_data_csv_XCMSdefault.txt", header = TRUE, sep = ',')

#standartd intensity list retrieved from XCMS 
intensities <- read.csv(file = "F:/avans/stage MM/peakpicking 2/feature_intensity_XCMSdefault.txt", header = FALSE, sep = ',')

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



write.csv(combined, 'BatchCorrect_import_data_1')


#######################
# editing sample data#
######################


metadata <- select(metadata, sample_name, class, batch, injectionOrder)
write.csv(metadata, 'batchCorrect_import_metadata')


##################################################################
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, 'BatchCorrect_import_data_1', "col", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-FilterVariable(mSet, "iqr", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

normalized_set <- mSet[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]

###########
meta_data <- read.csv('batchCorrect_import_metadata')
new_normalized_set <- cbind(meta_data[,3:5], ordered_normalized_set);
write.csv(new_normalized_set,file = "new_normalized_set.csv")

mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
rm(mSet)


rm(mSet)
mSet <- InitDataObjects("pktable", "utils", FALSE)
mSet <- Read.BatchDataTB(mSet, "new_normalized_set.csv", "row")


methods<- c("Combat", "WaveICA","EigenMS","QC_RLSC","ANCOVA","RUV_random","RUV_2","RUV_s","RUV_r","RUV_g","NOMIS","CCMN")
counter <- 0

# try function becouse not all datasets accepts the other functions (to many arguments for the LOG2 function)
for (i in 1:length(methods))  {
  try(mSet <- PerformBatchCorrection(mSetObj = mSet, imgName = methods[[i]], Method = methods[[i]]), silent = T)
  info <- (mSet$dataSet$interbatch_dis)
  counter <- counter + 1
  safe_name <- paste0("MetaboAnalyst_batch_data.csv",'_',methods[[i]],'.csv')
  file.rename("MetaboAnalyst_batch_data.csv", safe_name)
  
}

if (counter == 0){
  print('NO batch correction worked')
}


#selecting Best dataset
best_batchcorrection <- paste0("MetaboAnalyst_batch_data.csv",'_',gsub("_.*","",(names(info)[which.min(info)])),'.csv')
data_corrected <- read.csv(paste0(best_batchcorrection))
data_corrected_new <- data_corrected[,-3]
write.csv(data_corrected_new,file = "MetaboAnalyst_batch_data_stats.csv",row.names =  F)


rm(mSet)
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "MetaboAnalyst_batch_data_stats.csv", "row", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20) # No need to normalize again
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA,  1,2,0.95,0,0)



