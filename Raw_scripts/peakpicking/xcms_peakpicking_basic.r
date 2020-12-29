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
output_folder <- 'F:/avans/stage MM/Sherloktest_data_2/peakpick_output'
sample_metadata <-'sample_metadata.tsv'
unique_name <- 'default'
setwd(path)




output_folder <- paste0(output_folder,'_','XCMS_',unique_name)


dir.create(output_folder, showWarnings = F)
plot_allignment <- TRUE
################
# Library's####
##############
library(BiocStyle)
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

bwdpd <- 30
minFractiondpd <- 0.5
minSamplesdpd <- 1
binSizedpd <- 0.25
maxFeaturespdp <- 2000


mpp <- MergeNeighboringPeaksParam(
  expandRt = 4,
  expandMz = 0,
  ppm = 10,
  minProp = 0.75
)

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
data_files <- list.files(path = path, pattern = "*mzXML", full.names = TRUE, recursive = TRUE)
Sample_metadata <- list.files(path = path, pattern = "sample_metadata.tsv", full.names = TRUE, recursive = TRUE)            

# create an table of the meta data
data_frame <- read.table(file = Sample_metadata, sep = '\t', header = TRUE)

# order the dataframe by name (so all pathfiles will be going to the right row)
data_frame <- data_frame[order(data_frame$sample_name),]

# convert the table to an dataframe
data_frame <- as.data.frame.matrix(data_frame) 

#order the datafiles oin name
data_files <- sort(data_files)

data_files
data_frame


############### Analysing data   ##########################
raw_data <- readMSData(files = data_files, msLevel = 1, pdata =  new("NAnnotatedDataFrame", data_frame), mode = "onDisk")
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- refineChromPeaks(xdata, mpp)
xdata <- adjustRtime(xdata, param = obi)

if (plot_allignment == 'TRUE'){
  
  
    png(filename = paste0(path,'RT_adjustment.png'),
        width = 1840, height = 1840, units = "px", pointsize = 12,
        bg = "white",  res = NA,
    )

  bpis_adj <- chromatogram(xdata, aggregationFun = "max", include = "none")
  par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
  plot(bpis_adj, col = group_colors[bpis_adj$sampleType])
  ## Plot also the difference of adjusted to raw retention time.
  plotAdjustedRtime(xdata, col = group_colors[xdata$sampleType])
  
  
  dev.off()
}
pdp <- PeakDensityParam(sampleGroups = xdata$sampleType,
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
#xs <- as(xdata, "xcmsSet")
#an <- xsAnnotate(xs)

#anF <- groupFWHM(an)
#xsaC <- groupCorr(anF)
#xsaFI <- findIsotopes(xsaC)
#xsaFA <- findAdducts(xsaFI, polarity="positive")

#peaklist <- getPeaklist(xsaFA)
#peaklist <- peaklist[with(peaklist, order(rt, mz)),]

#peaklist1 <- select(peaklist, mz, rt, pcgroup, adduct, isotopes)

#isolist <- getIsotopeCluster(xsaFA)


#safe_name <- paste0(output_folder,'peaklistGrouped',unique_name,'.txt')
#write.csv(peaklist1, safe_name)
##############################
#Quantifying output and data #
##############################
res <- quantify(xdata, value = "into")
#meta data samples



safe_name <- paste0(output_folder,'/','sample_meta_data_XCMS_',unique_name,'.tsv')
sample_metadata1 <- as.data.frame(colData(res))
rownames(sample_metadata1) <- sample_metadata1[,1]
sample_metadata1 <- sample_metadata1[,-1]
write.table(sample_metadata1, safe_name, sep = '\t')

#feature list
safe_name <- paste0(output_folder,'/','Variable_metaData_XCMS_',unique_name,'.tsv')
feature_list<- as.data.frame(featureDefinitions(xdata))
feature_list <- feature_list[ , !(names(feature_list) %in% 'peakidx')]

write.table(feature_list,safe_name, sep = '\t') 



#intensity of found features
#head(assay(res))
safe_name <- paste0(output_folder,'/','Data_matrix_XCMS_',unique_name,'.tsv')
data_matrix <- featureValues(xdata, value = "into")

d <- data_matrix
features <- rownames(d)
rownames(d) <- NULL
data_matrix <- cbind(features,d)
for ( col in 1:ncol(data_matrix)){
  colnames(data_matrix)[col] <-  sub(".mzXML", "", colnames(data_matrix)[col])
}

rownames(data_matrix) <- data_matrix[,1]
data_matrix <- data_matrix[,-1]


write.table(data_matrix, safe_name, sep = '\t')







