library(xcms)
library(dplyr)



path_to_files <- "!@#$%^&1&^%$#@!"
data_file_extention <- "!@#$%^&2&^%$#@!"

path_to_variable_metadata <- "!@#$%^&3&^%$#@!"
path_to_output_folder <- '!@#$%^&4&^%$#@!'
unique_name <- '!@#$%^&5&^%$#@!'




####################

#peakpicking Centwave
cwp <- CentWaveParam(
  ppm = !@#$%^&6&^%$#@!,
  peakwidth = c(!@#$%^&7&^%$#@!),
  snthresh = !@#$%^&8&^%$#@!,
  prefilter = c(!@#$%^&9&^%$#@!),
  mzCenterFun = c("wMean","mean","apex","wMeanApex3")[!@#$%^&10&^%$#@!],,
  integrate = !@#$%^&11&^%$#@!,
  mzdiff = !@#$%^&12&^%$#@!,
  fitgauss = !@#$%^&13&^%$#@!,
  noise = !@#$%^&14&^%$#@!,
  verboseColumns = !@#$%^&15&^%$#@!,
  roiList = list(),
  firstBaselineCheck = !@#$%^&16&^%$#@!,
  roiScales = numeric(),
  extendLengthMSW = !@#$%^&17&^%$#@!
)


path_to_output_folder <- paste0(path_to_output_folder,'/Ms2_feature_list', unique_name, '.tsv')
data_file_extention <- paste0('*.',data_file_extention)
data_files <- list.files(path = path_to_files, pattern = data_file_extention, full.names = TRUE, recursive = FALSE)

swath_data <- readMSData(data_files, mode = "onDisk")


table(msLevel(swath_data))
head(fData(swath_data)[, c("isolationWindowTargetMZ",
                           "isolationWindowLowerOffset",
                           "isolationWindowUpperOffset",
                           "msLevel", "retentionTime")])


table(isolationWindowTargetMz(swath_data))


#################
# finding peaks #

swath_data <- findChromPeaks(swath_data, param = cwp)
cwp@snthresh <- cwp@snthresh/2
swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp)

chromPeakData(swath_data)

table(chromPeakData(swath_data)$isolationWindow)

###############################
# Reconstructing a MS2 spectra#
###############################
data <- read.table(path_to_variable_metadata, header = T, row.names(1), sep = '\t')

fenamiphos_mz <- data[,1]
Rt_ms1_feature <- data[,4]
feature_names <- rownames(data)
##########################
# find matching peaks MS1#
##########################
frame <- data.frame()
swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.9)
dataframe <- chromPeaks(swath_data)

for (ii in 1:length(fenamiphos_mz)){
  
  fenamiphos_ms1_peaka <- chromPeaks(swath_data, mz = fenamiphos_mz[ii], ppm = 5)
  
  ########## retrieving rownames, debuggging for if multiple match, chooses the best
  fenamiphos_ms1_peak1 <- as.data.frame(fenamiphos_ms1_peaka)
  fenamiphos_ms1_peak1$rt_dt <- abs(fenamiphos_ms1_peak1$rt - Rt_ms1_feature[ii])
  fenamiphos_ms1_peak1 <- fenamiphos_ms1_peak1[order(fenamiphos_ms1_peak1$rt_dt),]
  name <- rownames(fenamiphos_ms1_peak1[1,])
  idx_nr <- which(rownames(fenamiphos_ms1_peaka) == name) 
  peak_ID <- rownames(fenamiphos_ms1_peaka)[1]

    ###########################################################################
  fenamiphos_swath_spectrum <- swath_spectra[mcols(swath_spectra)$peak_id == peak_ID]
  pk_ids <- NULL
  try(pk_ids <- mcols(fenamiphos_swath_spectrum)$ms2_peak_id[[1]])
  #print(length(pk_ids))
  if (!(length(pk_ids) == 0)){
     Ms2_peaks <- (dataframe[rownames(dataframe) %in% pk_ids,]) 
    if (is.null(rownames(Ms2_peaks))){
      Ms2_peaks <- t(as.data.frame(Ms2_peaks))
      rownames(Ms2_peaks) <- fenamiphos_swath_spectrum@elementMetadata@listData[["ms2_peak_id"]]@listData[[1]]
      print(Ms2_peaks)
    }else{
      Ms2_peaks <- (as.data.frame(Ms2_peaks))
      
    }
     Ms2_peaks <- as.data.frame(Ms2_peaks)
     #Ms2_peaks$precursor_MZ <- fenamiphos_mz
     tmp  <- length(rownames(Ms2_peaks))
     mz_column <- rep(fenamiphos_mz[ii], tmp)
     rt_column <- rep(Rt_ms1_feature[ii], tmp)
     feature <- rep(feature_names[ii], tmp)
     
     Ms2_peaks$pre_cursor_mz <- mz_column
     Ms2_peaks$pre_cursor_rt <- rt_column
     Ms2_peaks$pre_cursor_feature_name <- feature
     
     print(Ms2_peaks)
     frame <- rbind(frame, Ms2_peaks)
     
  }
}

write.table(frame, file = path_to_output_folder, sep = '\t')




