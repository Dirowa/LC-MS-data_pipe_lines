
###############
# Variables   #
###############


#################
# install things#
#################
install.packages("BiocManager")
BiocManager::install(c("xcms","RColorBrewer","pander","pheatmap","SummarizedExperiment","magrittr","MSnbase"))
install.packages('doParallel','stats')
BiocManager::install("IPO")
###############
# Variables   #
###############

# path to data
path <- "F:/avans/stage MM/test_data_sherlok/test_data/"
output_folder <- 'F:/avans/stage MM/XCMS/output/'
name_safe <- 'XCMS_default'
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
library(IPO)
###########################
# activate speedy conzales#
###########################

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(3)))
} else {
  register(bpstart(SnowParam(3)))
}
register(SerialParam())



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

data_files
data_frame




############### Reading data in   ##########################
raw_data <- readMSData(files = data_files, msLevel = 1, pdata =  new("NAnnotatedDataFrame", data_frame), mode = "onDisk")


head(rtime(raw_data))


########################
# Retrieving mz value's#
########################

#mzs <- mz(raw_data)
#mzs_by_file <- split(mzs, f = fromFile(raw_data))
#length(mzs_by_file)

##########################
# plotting the basic data#
##########################

########################
# creating chromatogram#
########################
#bpis <- chromatogram(raw_data, aggregationFun = "max")
## Define colors for the two groups
#group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
#names(group_colors) <- c("sample", "QC")
## Plot all chromatograms.
#plot(bpis, col = group_colors[raw_data$class])

###################
# generate boxplot#
###################

#tc <- split(tic(raw_data), f = fromFile(raw_data))
#boxplot(tc, col = group_colors[raw_data$class],
#        ylab = "intensity", main = "Total ion current")

####################
# generate heatmaph#
#####################

#bpis_bin <- bin(bpis, binSize = 2)
## Calculate correlation on the log2 transformed base peak intensities
#cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
#colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
#ann <- data.frame(group = raw_data$class)
#rownames(ann) <- raw_data$sample_name
## Perform the cluster analysis
#pheatmap(cormat, annotation = ann,
#         annotation_color = list(group = group_colors))



########################
# generatiing ion graph#
########################

#raw_data %>%
#  filterRt(rt = rtr) %>%
#  filterMz(mz = mzr) %>%
#  plot(type = "XIC")


###########################################3#
#making a smaller graph for noice settings##
# peak picking for adjustment of parameters#  
#############################################

#rtr <- c(0, 1000)
#mzr <- c(500, 100000)
## extract the chromatogram
#chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
#plot(chr_raw, col = group_colors[chr_raw$class])


#centwave parameters
#cwp <-CentWaveParam(
#  ppm = 25,
#  peakwidth = c(20, 80),
#  snthresh = 10,
#  prefilter = c(3, 100),
#  mzCenterFun = "wMean",
#  integrate = 2L,
#  mzdiff = -0.001,
#  fitgauss = FALSE,
#  noise = 5000,
#  verboseColumns = FALSE,
#  roiList = list(),
#  firstBaselineCheck = TRUE,
#  roiScales = numeric(),
#  extendLengthMSW = FALSE
#)

#chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
#xchr <- findChromPeaks(chr_raw, param = cwp)
#sample_colors <- group_colors[xchr$class]
#plot(xchr, col = sample_colors,
#     peakBg = sample_colors[chromPeaks(xchr)[, "column"]])



#################################################
# first peak selection with m/z #################
################################################

#peakpicking
cwp <- CentWaveParam()
#  ppm = 25,
#  peakwidth = c(20, 80),
#  snthresh = 10,
#  prefilter = c(3, 100),
#  mzCenterFun = "wMean",
#  integrate = 2L,
#  mzdiff = -0.001,
#  fitgauss = FALSE,
#  noise = 5000,
#  verboseColumns = FALSE,
#  roiList = list(),
#  firstBaselineCheck = TRUE,
#  roiScales = numeric(),
#  extendLengthMSW = FALSE
#)


xdata <- findChromPeaks(raw_data, param = cwp)

#head(chromPeaks(xdata))
##########################################
# chrom peaks adjustment          ########
#########################################


mpp <- MergeNeighboringPeaksParam()
#  expandRt = 4,
#  expandMz = 0,
#  ppm = 10,
#  minProp = 0.75,
#  msLevel = 1
  
#)

xdata <- refineChromPeaks(xdata, mpp)



##########################
# Retrieving a summary####
#########################
#xdata <- xdata_pp

#summary_fun <- function(z)
#  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))

# <- lapply(split.data.frame(
#  chromPeaks(xdata), f = chromPeaks(xdata)[, "sample"]),
#  FUN = summary_fun)
#T <- do.call(rbind, T)
#rownames(T) <- basename(fileNames(xdata))
#pandoc.table(
#  T,
#  caption = paste0("Summary statistics on identified chromatographic",
#                   " peaks. Shown are number of identified peaks per",
#                   " sample and widths/duration of chromatographic ",
#                   "peaks."))


##########################
#adjust r time###########
########################
obi <- ObiwarpParam()
#  binSize = 1,
#  centerSample = integer(),
#  response = 100L,
#  distFun = "cor_opt",
#  gapInit = 0.3,
#  gapExtend = 2.4,
#  factorDiag = 2,
#  factorGap = 1,
#  localAlignment = FALSE,
#  initPenalty = 0,
#  subset = integer(),
#  subsetAdjust = c("average", "previous")
#)

xdata <- adjustRtime(xdata, param = obi)


## Get the base peak chromatograms.
#bpis_adj <- chromatogram(xdata, aggregationFun = "max", include = "none")
#par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
#plot(bpis_adj, col = group_colors[bpis_adj$class])
## Plot also the difference of adjusted to raw retention time.
#plotAdjustedRtime(xdata, col = group_colors[xdata$class])


# graph the adjustment
## Plot also the difference of adjusted to raw retention time.
#plotAdjustedRtime(xdata, col = group_colors[xdata$class])

## ----alignment-peak-groups-example-peak, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Example extracted ion chromatogram before (top) and after alignment (bottom)."----
#par(mfrow = c(2, 1))
## Plot the raw data
#plot(chr_raw, col = group_colors[chr_raw$class])

## Extract the chromatogram from the adjusted object
#chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
#plot(chr_adj, col = group_colors[chr_raw$class], peakType = "none")


#############################################################################
# This part is for when there is a big data, and will be done on QC samples#
###########################################################################


#drop the Adjustmennt
#xdata <- dropAdjustedRtime(xdata)

#define the layout
#xdata$class <- "sample"
#xdata$class[c(1, 2, 3)] <- "QC"

# set the parameters on what
#pdp_subs <- PeakDensityParam(sampleGroups = xdata$class, minFraction = 0.9)

#xdata <- groupChromPeaks(xdata, param = pdp_subs)

# group the parametrs
#pgp_subs <- PeakGroupsParam(minFraction = 0.85,
#subset = which(xdata$sample_type == "QC"),
#subsetAdjust = "average", span = 0.4)



#xdata <- adjustRtime(xdata, param = pgp_subs)

#make a nice plot
#clrs[xdata$sample_type == "QC"] <- c("#00ce0080")
#par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 0.5))
#plot(chromatogram(xdata, aggregationFun = "sum"),
#     col = clrs, peakType = "none")
#plotAdjustedRtime(xdata, col = clrs, peakGroupsPch = 1,
#peakGroupsCol = "#00ce0040")



######################
# group into features#
######################

pdp <- PeakDensityParam(sampleGroups = xdata$sampleType)
                        #,
                        #bw = 30,
                        #minFraction = 0.4,
                        #minSamples = 1,
                        #binSize = 0.25)


xdata <- groupChromPeaks(xdata, param = pdp)

#####################
# fill group peaks##
###################
fill <- FillChromPeaksParam()
#  expandMz = 0,
#  expandRt = 0,
#  ppm = 0,
#  fixedMz = 0,
#  fixedRt = 0
#)
xdata <- fillChromPeaks(xdata, param = fill)

res <- quantify(xdata, value = "into")
#meta data samples
colData(res)

#feature list
rowData(res)

#intensity of found features
head(assay(res))
head(featureValues(xdata, value = "into"))
head(featureSummary(xdata, group = xdata$class))


#meta data samples
filename = paste0(output_folder,'Sample_metadata_',name_safe,'.txt')
write.csv(colData(res), filename)

#feature list
filename = paste0(output_folder,'Feature_list_',name_safe,'.txt')
write.csv(rowData(res),filename) 

#intensity of found features
#head(assay(res))
filename = paste0(output_folder,'feature_intensity_',name_safe,'.txt')
write.csv((featureValues(xdata, value = "into")), "F:/avans/stage MM/XCMS/output/feature_intensity_centwave_N5000.txt")

head(featureSummary(xdata, group = xdata$class))
