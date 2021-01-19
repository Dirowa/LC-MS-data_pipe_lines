#################
# install things#
#################
#install.packages("BiocManager")
#BiocManager::install(c("xcms","RColorBrewer","pander","pheatmap","SummarizedExperiment","magrittr","MSnbase"))
#install.packages('doParallel','stats')
#BiocManager::install("IPO")
#BiocManager::install("CAMERA")

###############
# Variables   #
###############

# path to data & metadata (tsv format)
path <- "F:/avans/stage MM/Sherloktest_data_2/"
output_folder <- 'F:/avans/stage MM/xcms_pipeline/'
unique_name <- 'default'
setwd(path)

sample_name <- 'sample_name'
sample_type <- 'sampleType'

sample_metadata <-'sample_metadata.tsv'
data_file_extention <- "*.mzXML"


data_file_extention <- paste0("*.",data_file_extention)
output_folder <- paste0(output_folder,'_','XCMS_',unique_name)
dir.create(output_folder, showWarnings = F)
################
# Library's####
##############
library(xcms)
library(faahKO)
library(pander)
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(stats)
library(CAMERA)
library(tibble)
library(dplyr)
library(data.table)


####################
# Setting variables#
####################


#peakpicking Centwave
cwp <- CentWaveParam(
  ppm = 25,
  peakwidth = c(20, 50),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1L,
  mzdiff = -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = numeric(),
  extendLengthMSW = FALSE
)

mpp <- MergeNeighboringPeaksParam(
  expandRt = 2,
  expandMz = 0,
  ppm = 10,
  minProp = 0.75
)

bwdpd <- 30
minFractiondpd <- 0.5
minSamplesdpd <- 1
binSizedpd <- 0.25
maxFeaturespdp <- 2000



obi <- ObiwarpParam(
  binSize = 1,
  centerSample = 3,
  response = 1L,
  distFun = "cor",
  gapInit = 0.5,
  gapExtend = 2.5,
  factorDiag = 2,
  factorGap = 1,
  localAlignment = FALSE,
  initPenalty = 0,
  subset = integer(),
  subsetAdjust = c("average", "previous")
)



fill <- FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)




#############################################
##  READING DATA IN & correctiong format ####
###########################################

# importing the files from the folder into R. all !!!! all files must be .mzML and one .tsv <-= sample meta data
data_files <- list.files(path = path, pattern = data_file_extention, full.names = TRUE, recursive = FALSE)
Sample_metadata <- list.files(path = path, pattern = "*.tsv", full.names = TRUE, recursive = FALSE)            
# create an table of the meta data
data_frame <- read.table(file = Sample_metadata, sep = '\t', header = TRUE)
# order the dataframe by name (so all pathfiles will be going to the right row)
data_frame <- data_frame[order(data_frame[sample_name]),]
# convert the table to an dataframe
data_frame <- as.data.frame.matrix(data_frame) 
#order the datafiles oin name
data_files <- sort(data_files)



data_files
data_frame[sample_type]


############### Analysing data   ##########################
raw_data <- readMSData(files = data_files, msLevel = 1, pdata =  new("NAnnotatedDataFrame", data_frame), mode = "onDisk")
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- refineChromPeaks(xdata, mpp)
xdata <- adjustRtime(xdata, param = obi)



pdp <- PeakDensityParam(sampleGroups = unlist(data_frame[sample_type]),
                        bwdpd,
                        minFractiondpd,
                        minSamplesdpd,
                        binSizedpd,
                        maxFeaturespdp)


xdata <- groupChromPeaks(xdata, param = pdp)
xdata <- fillChromPeaks(xdata, param = fill)



#############
# annotation#
#############
xs <- as(xdata, "xcmsSet")
an <- xsAnnotate(xs)

anF <- groupFWHM(an)
xsaC <- groupCorr(anF)
xsaFI <- findIsotopes(xsaC)
xsaFA <- findAdducts(xsaFI, polarity="negative")

peaklist <- getPeaklist(xsaFA)
peaklist <- peaklist[with(peaklist, order(rt, mz)),]

peaklist <- peaklist[, -grep("X", colnames(peaklist))]


##############################
#Quantifying output and data #
##############################
res <- quantify(xdata, value = "into")
#meta data samples




#feature list
safe_name <- paste0(output_folder,'/','Variable_metaData_XCMS_',unique_name,'.tsv')
feature_list<- as.data.frame(featureDefinitions(xdata))
rownames(peaklist) <- rownames(feature_list)

write.table(peaklist,safe_name, sep = '\t') 

line <- readLines(safe_name,)
line[1] <- paste0('""\t',line[1])
writeLines(line,safe_name,)



#intensity of found features
#head(assay(res))
safe_name <- paste0(output_folder,'/','Data_matrix_XCMS_',unique_name,'.tsv')
data_matrix <- featureValues(xdata, value = "into")

d <- data_matrix
features <- rownames(d)
rownames(d) <- NULL
data_matrix <- cbind(features,d)
for ( col in 1:ncol(data_matrix)){
  colnames(data_matrix)[col] <-  sub(data_file_extention, "", colnames(data_matrix)[col])
}

rownames(data_matrix) <- data_matrix[,1]
data_matrix <- data_matrix[,-1]


write.table(data_matrix, safe_name, sep = '\t')

line <- readLines(safe_name,)
line[1] <- paste0('""\t',line[1])
writeLines(line,safe_name,)


safe_name <- paste0(output_folder,'/','sample_meta_data_XCMS_',unique_name,'.tsv')
sample_metadata1 <- as.data.frame(colData(res))
rownames(sample_metadata1) <- sample_metadata1[,1]
sample_metadata1 <- sample_metadata1[,-1]
rownames(sample_metadata1) <- colnames(data_matrix)
write.table(sample_metadata1, safe_name, sep = '\t')



line <- readLines(safe_name,)
line[1] <- paste0('"',sample_name,'"\t',line[1])
writeLines(line,safe_name,)




