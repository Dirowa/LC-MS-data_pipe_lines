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
unique_name <- '!@#$%^&1&^%$#@!'
path <- "!@#$%^&2&^%$#@!"
output_folder <- '!@#$%^&3&^%$#@!'

try(setwd(path))
sample_name <- '!@#$%^&4&^%$#@!'
sample_type <- '!@#$%^&5&^%$#@!'
QC_in_sample_type <- '!@#$%^&6&^%$#@!'
blank_in_sampleType <- "!@#$%^&7&^%$#@!"
sample_in_sampleType <- "!@#$%^&8&^%$#@!"

injection_order <- '!@#$%^&9&^%$#@!'
batch_name <- '!@#$%^&10&^%$#@!'

sample_metadata <-'!@#$%^&11&^%$#@!'
data_file_extention <- "!@#$%^&12&^%$#@!"

polarity1 <- c('negative','positive')[!@#$%^&13&^%$#@!]



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
  ppm = !@#$%^&14&^%$#@!,
  peakwidth = c(!@#$%^&15&^%$#@!),
  snthresh = !@#$%^&16&^%$#@!,
  prefilter = c(3, 100),
  mzCenterFun = c("wMean","mean","apex","wMeanApex3")[!@#$%^&17&^%$#@!],
  integrate = !@#$%^&18&^%$#@!,
  mzdiff = !@#$%^&19&^%$#@!,
  fitgauss = !@#$%^&20&^%$#@!,
  noise = !@#$%^&21&^%$#@!,
  verboseColumns = !@#$%^&22&^%$#@!,
  roiList = list(),
  firstBaselineCheck = !@#$%^&23&^%$#@!,
  roiScales = numeric(),
  extendLengthMSW = !@#$%^&24&^%$#@!
)

mpp <- MergeNeighboringPeaksParam(
  expandRt = !@#$%^&25&^%$#@!,
  expandMz = !@#$%^&26&^%$#@!,
  ppm = !@#$%^&27&^%$#@!,
  minProp = !@#$%^&28&^%$#@!
)

bwdpd <- !@#$%^&29&^%$#@!
minFractiondpd <- !@#$%^&30&^%$#@!
minSamplesdpd <- !@#$%^&31&^%$#@!
binSizedpd <- !@#$%^&32&^%$#@!
maxFeaturespdp <- !@#$%^&33&^%$#@!



obi <- ObiwarpParam(
  binSize = !@#$%^&34&^%$#@!,
  centerSample = !@#$%^&35&^%$#@!,
  response = !@#$%^&36&^%$#@!,
  distFun = c("cor","cor_opt","cov","prd","euc")[!@#$%^&37&^%$#@!],
  gapInit = !@#$%^&38&^%$#@!,
  gapExtend = !@#$%^&39&^%$#@!,
  factorDiag = !@#$%^&40&^%$#@!,
  factorGap = !@#$%^&41&^%$#@!,
  localAlignment = !@#$%^&42&^%$#@!,
  initPenalty = !@#$%^&43&^%$#@!,
  subset = integer(),
  subsetAdjust = c("average", "previous")
)



fill <- FillChromPeaksParam(
  expandMz = !@#$%^&44&^%$#@!,
  expandRt = !@#$%^&45&^%$#@!,
  ppm = !@#$%^&46&^%$#@!,
  fixedMz = !@#$%^&47&^%$#@!,
  fixedRt = !@#$%^&48&^%$#@!
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
xsaFA <- findAdducts(xsaFI, polarity=polarity1)

peaklist <- getPeaklist(xsaFA)
peaklist <- peaklist[with(peaklist, order(rt, mz)),]

peaklist <- peaklist[, -grep("X", colnames(peaklist))]


##############################
#Quantifying output and data #
##############################
res <- quantify(xdata, value = "into")
#meta data samples

sample_name <- 'sample_name'


#feature list
safe_name <- paste0(output_folder,'/','Variable_metaData_XCMS_',unique_name,'.tsv')
feature_list<- as.data.frame(featureDefinitions(xdata))
rownames(peaklist) <- rownames(feature_list)


colnames(peaklist)[colnames(peaklist) == QC_in_sample_type] <- "QC"
colnames(peaklist)[colnames(peaklist) == blank_in_sampleType] <- "blank"
colnames(peaklist)[colnames(peaklist) == sample_in_sampleType] <- "sample"

write.table(peaklist,safe_name, sep = '\t') 

line <- readLines(safe_name,)
line[1] <- paste0('""\t',line[1])
writeLines(line,safe_name,)



#intensity of found features
#head(assay(res))
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







safe_name <- paste0(output_folder,'/','sample_meta_data_XCMS_',unique_name,'.tsv')
sample_metadata1 <- as.data.frame(colData(res))
rownames(sample_metadata1) <- sample_metadata1[,1]
sample_metadata1 <- sample_metadata1[,-1]
rownames(sample_metadata1) <- colnames(data_matrix)
sample_metadata1[sample_type][sample_metadata1[sample_type] == QC_in_sample_type] <- "QC"
sample_metadata1[sample_type][sample_metadata1[sample_type] == blank_in_sampleType] <- "blank"
sample_metadata1[sample_type][sample_metadata1[sample_type] == sample_in_sampleType] <- "sample"

colnames(sample_metadata1)[colnames(sample_metadata1) == sample_type] <- "sample_type"
colnames(sample_metadata1)[colnames(sample_metadata1) == batch_name] <- 'batch'
colnames(sample_metadata1)[colnames(sample_metadata1) == injection_order] <- "Injection_order"


write.table(sample_metadata1, safe_name, sep = '\t')



line <- readLines(safe_name,)
line[1] <- paste0('"',sample_name,'"\t',line[1])
writeLines(line,safe_name,)







safe_name <- paste0(output_folder,'/','Data_matrix_XCMS_',unique_name,'.tsv')
colnames(data_matrix) <- rownames(sample_metadata1)


write.table(data_matrix, safe_name, sep = '\t')

line <- readLines(safe_name,)
line[1] <- paste0('""\t',line[1])
writeLines(line,safe_name,)
