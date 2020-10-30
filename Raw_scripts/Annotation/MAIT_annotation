######################
# Camera Annotation###
######################

install.packages("BiocManager")
BiocManager::install(c("xcms","RColorBrewer","pander","pheatmap","SummarizedExperiment","magrittr","MSnbase"))
install.packages('doParallel','stats')
BiocManager::install("CAMERA")

######################
# important variables#
######################

# path to data
path <- "F:/avans/stage MM/test_data_sherlok/test_data/"

############
# library's#
############

################
# Library's####
##############
library(BiocStyle)
library(xcms)
library(pander)
library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(stats)
library(CAMERA)


##############
# Peakpicking#
##############


#############################################
##  READING DATA IN & correctiong format ####
###########################################

# importing the files from the folder into R. all !!!! all files must be .mzML and one .tsv <-= sample meta data
data_files <- list.files(path = path, pattern = "*.mzML", full.names = TRUE, recursive = TRUE)
Sample_metadata <- list.files(path = path, pattern = "*.tsv", full.names = TRUE, recursive = TRUE)            
# create an table of the meta data
data_frame <- read.table(file = Sample_metadata, sep = '\t', header = TRUE)
# order the dataframe by name (so all pathfiles will be going to the right row)
data_frame <- data_frame[order(data_frame$sample_name),]
# convert the table to an dataframe
data_frame <- as.data.frame.matrix(data_frame) 
#order the datafiles oin name
data_files <- sort(data_files)
#data_files
#data_frame



###################
# XCMS peakpicking#
###################

############### Reading data in   ##########################
raw_data <- readMSData(files = data_files, msLevel = 1, pdata =  new("NAnnotatedDataFrame", data_frame), mode = "onDisk")
head(rtime(raw_data))

#################################################
# first peak selection with m/z #################
################################################

#peakpicking Centwave
cwp <- CentWaveParam(
  ppm = 42.5,
  peakwidth = c(11, 57.2),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1L,
  mzdiff = 0.0045,
  fitgauss = FALSE,
  noise = 5000,
  verboseColumns = FALSE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = numeric(),
  extendLengthMSW = FALSE
)


xdata <- findChromPeaks(raw_data, param = cwp)

#head(chromPeaks(xdata))
##########################################
# chrom peaks adjustment          ########
#########################################


mpp <- MergeNeighboringPeaksParam(
  expandRt = 4,
  expandMz = 0,
  ppm = 10,
  minProp = 0.75
)

xdata <- refineChromPeaks(xdata, mpp)


##########################
#adjust r time###########
########################


obi <- ObiwarpParam(
  binSize = 0.682,
  centerSample = 3,
  response = 1L,
  distFun = "cor",
  gapInit = 0.544,
  gapExtend = 2.4,
  factorDiag = 2,
  factorGap = 1,
  localAlignment = FALSE,
  initPenalty = 0,
  subset = integer(),
  subsetAdjust = c("average", "previous")
)

xdata <- adjustRtime(xdata, param = obi)

######################
# group into features#
######################

pdp <- PeakDensityParam(sampleGroups = xdata$sampleType,
                        bw = 0.25,
                        minFraction = 1,
                        minSamples = 1,
                        binSize = 0.25)


xdata <- groupChromPeaks(xdata, param = pdp)

#####################
# fill group peaks##
###################
fill <- FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)
xdata <- fillChromPeaks(xdata, param = fill)

res <- quantify(xdata, value = "into")


#meta data samples
write.csv(colData(res), 'F:/avans/stage MM/XCMS/output/sample_meta_data_CSV_XCMS.txt')
#feature list
write.csv(rowData(res),'F:/avans/stage MM/XCMS/output/Feature_list_centwave_N5000.txt') 

#intensity of found features
#head(assay(res))
write.csv((featureValues(xdata, value = "into")), "F:/avans/stage MM/XCMS/output/feature_intensity_centwave_N5000.txt")

head(featureSummary(xdata, group = xdata$class))

########
# MAIT #
########


####################
# MAIT#############
###################

library(MAIT)


peaks <- (featureValues(xdata))
masses <- getPeaklist(xsaFA)$mz
rt <- getPeaklist(xsaFA)$rt/60


importMAIT <- MAITbuilder(data = peaks, masses = masses,
                          rt = rt, significantFeatures = TRUE,
                          spectraEstimation = TRUE,rtRange=0.2,
                          corThresh=0.7)

importMAIT <- Biotransformations(MAIT.object = importMAIT,
                                 adductAnnotation = TRUE,
                                 peakPrecision = 0.005, adductTable = NULL)

signTable <- sigPeaksTable(MAIT.object = importMAIT, printCSVfile = FALSE)



