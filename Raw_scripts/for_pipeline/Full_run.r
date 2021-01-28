####################
##### Full Run #####
####################

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


################
## Peakpicking #
################

  ###############
  # Variables   #
  ###############
  
  # path to data & metadata (tsv format)
  path <- "F:/avans/stage MM/Sherloktest_data_2/"
  output_folder <- 'F:/avans/stage MM/full_test_run/'
  unique_name <- 'default'
  setwd(path)
  sample_name <- 'sample_name'
  sample_type <- 'sampleType'
  QC_in_sample_type <- 'pool'
  blank_in_sampleType <- "blank"
  sample_in_sampleType <- "sample"
  
  injection_order <- 'injectionOrder'
  batch_name <- 'batch'
  
  sample_metadata <-'sample_metadata.tsv'
  data_file_extention <- "mzXML"
  
  polarity1 <- c('negative','positive')[2]
  
  data_file_extention <- paste0("*.",data_file_extention)
  dir.create(output_folder, showWarnings = F)
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
  Sample_metadata <- list.files(path = path, pattern = sample_metadata, full.names = TRUE, recursive = FALSE)            
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
  
  
  #########################################
  ############### Analysing data ########## 
  ########################################
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
  
  
  
  ################
  # annotation####
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
  
##############################
# unload packages #
######################
  
  detach_package(xcms)
  detach_package(faahKO)
  detach_package(pander)
  detach_package(xcms)
  detach_package(faahKO)
  detach_package(RColorBrewer)
  detach_package(pander)
  detach_package(magrittr)
  detach_package(pheatmap)
  detach_package(SummarizedExperiment)
  detach_package(stats)
  detach_package(CAMERA)
  detach_package(tibble)
  detach_package(dplyr)
  detach_package(data.table)

  rm(an, anF, cwp, d, data_frame, data_matrix, feature_list, fill, mpp, obi, pdp, peaklist, raw_data
     ,res, sample_metadata1, xdata, xs, xsaC, xsaFA, xsaFI, line, data_file_extention)
  
#############################
###### Batch_correction #####
  ###########################
  
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
  
  #######################
  # Remove uneeded items#
  #######################
  

  ########################
  # important variables #
  #######################
  # do not change#
  
  input_folder <- output_folder

  Data_matrix_xcms <- paste0('Data_matrix_XCMS_',unique_name,'.tsv')
  samplemetadata <- paste0('sample_meta_data_XCMS_',unique_name,'.tsv')
  variable_metadata <- paste0('Variable_metaData_XCMS_',unique_name,'.tsv')
  
  final_output_folder <- "batch_correction"
  # columns names in dataframe of metadata
  
  Best_peakpicked_based_on <- c("mean_var_QC_log10","Ratio_VAR_QC_log")[2]
  
  
  
  
  
  
  
  pixelsize1 <- 20
  pixelsize2 <- 12
  # setting up some directories#
  final_output_folder <- paste0(input_folder,'/',final_output_folder)
  dir.create(final_output_folder, showWarnings = F)
  setwd(final_output_folder)
  
  variable_metadata_naam <- variable_metadata
  
  
  
  ############# settins which should not be edited ################
  sample_name <- 'Sample_name'
  sample_type <- 'sample_type'
  batch <- 'batch'
  injectionOrder <- 'Injection_order'
  #################################################################
  # IMPORT DATA SET#
  ###############################################
  
  path_to_Data_matrix_xcms <- paste0(input_folder,'/',Data_matrix_xcms )
  path_to_samplemetadata <- paste0(input_folder,'/',samplemetadata )
  Path_to_variable_metadata <- paste0(input_folder,'/',variable_metadata )
  # the sample metadata 
  # the sample metadata is required to have these column names:  sample_name, sample_type, batch, injectionOrder
  metadata <- read.csv(file = path_to_samplemetadata, header = TRUE, sep = '\t', row.names = 1)
  #standartd intensity list retrieved from XCMS 
  intensities <- read.csv(file = path_to_Data_matrix_xcms, header = TRUE, sep = '\t', row.names = 1)
  variable_metadata_matrix <- read.csv(file = Path_to_variable_metadata, header = TRUE, sep = '\t', row.names = 1)
  
  
  #######################################
  # Editing datasets for batchcorrection#
  #######################################
  # switch dataframe
  metadata1 <- as.data.frame(t(metadata))
  
  # order to be sure
  metadata1<- metadata1[ , order(names(metadata1))]
  intensities<- intensities[ , order(names(intensities))]
  
  Index <- which(rownames(metadata1) == sample_type) 
  sample_types  <- metadata1[Index, ]  ## subsets dataframe
  
  colnames(sample_types) <- colnames(intensities)
  combined <- rbind(sample_types,intensities)
  
  write.csv(combined, 'BatchCorrect_import_data_1')
  ####################################
  # Generate overview of current data#
  ####################################
  try(rm(mSet))
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, 'BatchCorrect_import_data_1', "col", "disc")
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet)
  # here was RSD filtering before #
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ref = names[1,])
  
  normalized_set <- mSet[["dataSet"]][["norm"]]
  ordered_normalized_set1 <- normalized_set[order(row.names(normalized_set)), ]
  
  
  ########################################
  # creating dataset for batch correction#
  ########################################
  
  Samplemetadata1 <- select(metadata, sample_type, batch, injectionOrder)
  new_normalised_set <- cbind(Samplemetadata1,ordered_normalized_set1)
  
  
  #####################
  # ordering datasets #
  #####################
  
  write.csv(new_normalised_set, "new_normalized_set.csv")
  
  #############################
  # correcting for loes method#
  #############################
  #### possible correcting for removed features#####
  
  a <-rownames(variable_metadata_matrix)
  b <-colnames(new_normalised_set)
  diffrences <- setdiff(a,b)
  if (length(diffrences) != 0){
    diffrences <- diffrences[diffrences != ""]
    
    #Remove the differences
    for (i in 1:length(diffrences)){
      variable_metadata_matrix<-variable_metadata_matrix[!(rownames(variable_metadata_matrix)==diffrences[[i]]),]
    }}
  
  
  
  loes_file_dataMatrix.tsv <- paste0(final_output_folder,'/loes_dataMatrix.tsv')
  loes_file_sampleMetadata.tsv <- paste0(final_output_folder,'/loes_sampleMetadata.tsv')
  loes_file_variableMetadata.tsv <- paste0(final_output_folder,'/loes_variableMetadata.tsv')
  
  ordered_normalized_set1 <- t((ordered_normalized_set1))
  
  write.table(ordered_normalized_set1, loes_file_dataMatrix.tsv, sep ='\t')
  write.table(metadata, loes_file_sampleMetadata.tsv, sep ='\t')
  write.table(variable_metadata_matrix, loes_file_variableMetadata.tsv, sep ='\t')
  
  
  
  
  
  # rewiting becouse it needs "" \t for, cant be done by playing with dataframes
  line <- readLines(loes_file_dataMatrix.tsv,)
  line[1] <- paste0('""\t',line[1])
  writeLines(line,loes_file_dataMatrix.tsv )
  
  line <- readLines(loes_file_sampleMetadata.tsv,)
  line[1] <- paste0('"',sample_name,'"\t',line[1])
  #line[1] <- paste0('""\t',line[1])
  writeLines(line,loes_file_sampleMetadata.tsv )
  
  line <- readLines(loes_file_variableMetadata.tsv,)
  line[1] <- paste0('""\t',line[1])
  writeLines(line,loes_file_variableMetadata.tsv )
  
  
  ######################################
  # inspecting #
  #######################################
  
  
  set <- phenomis::reading(NA,
                           files.ls = list(dataMatrix = file.path(loes_file_dataMatrix.tsv),
                                           sampleMetadata = file.path(loes_file_sampleMetadata.tsv),
                                           variableMetadata = file.path(loes_file_variableMetadata.tsv)))
  
  
  png(filename = paste0(final_output_folder,'/inspecting dataset.png'),
      width = 960, height = 960, units = "px", pointsize = pixelsize1,
      bg = "white",  res = NA
  )
  
  
  
  sacurine.eset <- phenomis::inspecting(set,
                                        pool_as_pool1.l = FALSE,
                                        pool_cv.n = 0.3,
                                        span.n = 1,
                                        sample_intensity.c = c("median", "mean")[1],
                                        title.c = "Inital Overview",
                                        plot_dims.l = TRUE,
                                        figure.c = c("none", "interactive", "myfile.pdf")[2],
                                        
                                        report.c = c("none", "interactive", "myfile.txt")[3])
  
  
  dev.off()
  try(rm(set))
  try(rm(sacurine.eset))
  
  ################
  # calculate PCA#
  ################
  mSet <- PCA.Anal(mSet)
  mSet <- PlotPCAPairSummary(mSet, "pca_pair_before_correction", "png", 72, width=NA, 5)
  mSet <- PlotPCAScree(mSet, "pca_scree_before_correction_", "png", 72, width=NA, 5)
  mSet <- PlotPCA2DScore(mSet, "pca_score2d_before_correction_", "png", 72, width=NA, 1,2,0.95,0,0)
  
  
  
  #############################
  # performing batchcorrection#
  #############################
  try(rm(mSet))
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
      
      png(filename = paste0(final_output_folder,'/loes_batchcorrection.png'),
          width = 960, height = 960, units = "px", pointsize = pixelsize2,
          bg = "white",  res = NA,
      )
      
      
      
      
      corrcted_set <- phenomis::correcting(
        eset,
        reference.c = "QC",
        col_batch.c = batch,
        col_injectionOrder.c = injectionOrder,
        col_sampleType.c = sample_type,
        span.n = 1,
        title.c = 'Loess Batch Correction',
        figure.c = c("none", "interactive", "Loess_correction.pdf")[2],
        report.c = c("none", "interactive", "Loess_correction.pdf")[3]
      )
      
      
      dev.off()
      
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
      print(items[[i]])
      variance_dataframe <- read.csv(paste0(final_output_folder,'/','new_normalized_set.csv'), header = TRUE, sep = ',')
      variance_dataframe1 <- variance_dataframe[,-1:-4]
      print(mean(sapply(variance_dataframe1,  var)))
      means2[[i]] <- mean(sapply(variance_dataframe1,  var))
      
      variance_dataframe <- variance_dataframe[grep('QC', variance_dataframe[[sample_type]]),]
      variance_dataframe <- variance_dataframe[,-1:-4]
      #variance_dataframe[is.na(variance_dataframe)] <- 0
      print(mean(sapply(variance_dataframe,  var)))
      means[[i]] <- mean(sapply(variance_dataframe,  var))
      rm(variance_dataframe)
      rm(variance_dataframe1)
    }
    
    
    else if (as.character(items[[i]]) == 'loessALL_BatchCorrected_data_dataMatrix.tsv'){
      print(items[[i]])
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
      
      
      
      variance_dataframe <- variance_dataframe[ , grepl( 'QC' , names( variance_dataframe ) ) ]
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
      print(items[[i]])
      variance_dataframe1 <- variance_dataframe[,-1:-3]
      var <- sapply(variance_dataframe1,  var)
      print(mean(var))
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
  means <- do.call(rbind, Map(data.frame,mean_var_total_log10=means2, mean_var_QC_log10=means, file=items))
  means <- transform(means, Ratio_VAR_QC_log =  mean_var_QC_log10 / mean_var_total_log10)
  means <- means[order(means[Best_peakpicked_based_on]),]
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
  a <-rownames(variable_metadata)
  b <-rownames(data_matrix)
  diffrences <- setdiff(a,b)
  if (length(diffrences) != 0){
    diffrences <- diffrences[diffrences != ""]
    
    #Remove the differences
    for (i in 1:length(diffrences)){
      variable_metadata<-variable_metadata[!(variable_metadata$V1==diffrences[[i]]),]
    }}
  
  # make the files pretty enough for statistical analyis
  colnames(variable_metadata) <- variable_metadata[1,]
  variable_metadata <- variable_metadata[-1,]
  TMP_features_names <- variable_metadata[,1]
  variable_metadata <- variable_metadata[,-1]
  
  isotopes <- variable_metadata$isotopes
  adducts <- variable_metadata$adduct
  
  variable_metadata <- mutate_all(variable_metadata, function(x) as.numeric(as.character(x)))
  rownames(variable_metadata) <- TMP_features_names
  
  variable_metadata$isotopes <- isotopes
  variable_metadata$adduct <- adducts
  
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
  
  
################################
# cleaning up #################
############################
  
  detach_package(pacman)
  detach_package("impute")
  detach_package("pcaMethods")
  detach_package("globaltest")
  detach_package("GlobalAncova")
  detach_package("Rgraphviz")
  detach_package("preprocessCore")
  detach_package("genefilter")
  detach_package("SSPA")
  detach_package("sva")
  detach_package("limma")
  detach_package("KEGGgraph")
  
  detach_package("BiocParallel")
  detach_package("MSnbase")
  detach_package("multtest")
  detach_package("RBGL")
  detach_package("fgsea")
  detach_package(MetaboAnalystR)
  detach_package(dplyr)
  detach_package(devtools)
  detach_package( phenomis)
  
  
  rm (combined, mSet, new_normalised_set, data_matrix, intensities, best_batch_corrected, ordered_normalized_set1, corrcted_set, eset,
      normalized_set, variable_metadata, line, variable_metadata_matrix, a, diffrences,b,features , TMP_features_names, metadata, adducts, isotopes, metadata1
      ,data_files, sample_types, sample_metadata1, means, items, means2, Samplemetadata1)
  
  
################################
# pre_univariate filtering#
################################  #
  
  
  #####################
  #  filter the noise#
  ####################
  #install.packages("randomcoloR")
  library(dplyr)
  library( randomcoloR)
  library(reshape2)
  library(MetaboAnalystR)
  
  folder <- final_output_folder
  outputfolder_name <- "pre_univariate_filterd"
  
  
  
  variable_metadata <- paste0(variable_metadata_naam,'_batchcorrected.tsv')
  data_matrix <- paste0(Data_matrix_xcms,'_batchcorrected.tsv')
  sample_metadata <-paste0(samplemetadata,'_batchcorrected.tsv' )
  
  
  Filter_on_isotopes <- TRUE
  Filter_on_RSD <- TRUE
  noise_reduction <- TRUE
  
  ### isotope reduction #
  ppm_mz = 0.003
  ppm_rt = 0.003
  
  ### RSD filtering ###
  RSD_treshhold <- 25
  
  ### noise filtering ####
  min_amoutn_of_sample_hit <- 1
  
  
  
  #####################
  # creating the paths#
  #####################
  outputfolder <- (paste0(folder,'/',outputfolder_name,'/'))
  try(dir.create(outputfolder, showWarnings = T))
  
  
  
  
  Variable_metadata_path <- paste0(folder,'/',  variable_metadata)
  data_matrix_path <- paste0(folder, '/',data_matrix)
  sample_metadata_path <- paste0(folder, '/',sample_metadata)
  
  ####################
  # reading in data###
  ####################
  
  
  variable_metadata <- read.table(Variable_metadata_path, sep = '\t', header = TRUE, row.names = 1)
  data_matrix <- read.table(data_matrix_path, sep = '\t', header = TRUE, row.names = 1)
  sample_metadata <- (read.table(sample_metadata_path, sep = '\t', header = TRUE, row.names = 1))
  
  
  colnames(data_matrix) <- rownames(sample_metadata)
  
  
  
  
  ###########################
  #  ordering the dataframes#
  ###########################
  data_matrix <- data_matrix[ , order(names(data_matrix))]
  sample_metadata <- sample_metadata[order(row.names(sample_metadata)), ]
  
  
  if (isTRUE(Filter_on_isotopes)){
    ###################################
    # Filter on isotopes #
    ##################################
    
    isotope_table <- variable_metadata %>% select(mz, rt)
    isotope_table <- isotope_table[order(isotope_table$rt),]
    #screen for the sum of all samples intensities#
    sample_names <- rownames(sample_metadata[sample_metadata$sample_type == 'sample',])
    
    data_orginal <- 10**data_matrix[sample_names]
    isotope_table <-  cbind(isotope_table,sumINT=(rowMeans(data_orginal)))
    
    # order the isotope table from mz high to low to find possible isotopes#
    isotope_table <- isotope_table[order(isotope_table$mz),]
    
    # lets find isotopes
    isotope_table[,"isotopes"] <- NA
    
    #################################################
    # grouping the isotopes by MZ and Retention time#
    #################################################
    
    
    i <- 0
    while (i <= length(rownames(isotope_table))){
      i = i + 1
      if ( is.na(isotope_table[i,'isotopes'])){
        mz_value <- isotope_table$mz[i]
        rt_value <- isotope_table$rt[i]
        original_feat <- rownames(isotope_table[isotope_table$mz == mz_value,])
        
        mz_value1 <- mz_value - 1 - ppm_mz
        mz_value2 <- mz_value + 1 + ppm_mz
        rt_value1 <- rt_value - 1 - ppm_rt
        rt_value2 <- rt_value + 1 - ppm_rt
        
        ###### with the first Mz value is checked weither there are any MZ values near it on 1 +/- 1 and ppm
        possible_isotopes <- isotope_table[isotope_table$rt >= rt_value1 & isotope_table$rt <= rt_value2 ,]
        # this list getrs updated with near Rt time ######
        possible_istopes1 <- (possible_isotopes[possible_isotopes$mz >= mz_value1 & possible_isotopes$mz <= mz_value2 ,])
        
        #remove features who are already assigned a isotope group #
        possible_istopes1 <- possible_istopes1[is.na(possible_istopes1$isotopes), ]
        
        # assinging isotopes groups #
        tmp_features_to_assign <- rownames(possible_istopes1)
        for (ii in 1:length(tmp_features_to_assign)){
          rowN <- which(rownames(isotope_table) == tmp_features_to_assign[ii]) 
          isotope_table[rowN,'isotopes'] <- i
        }
        
        if (length(tmp_features_to_assign) != 1){
          
          mz_value1 <- mz_value - 2 - ppm_mz
          mz_value2 <- mz_value + 2 + ppm_mz
          possible_istopes1 <- rownames(possible_isotopes[possible_isotopes$mz >= mz_value1 & possible_isotopes$mz <= mz_value2 ,])
          
          for (ii in 1:length(possible_istopes1)){
            rowN <- which(rownames(isotope_table) == possible_istopes1[ii]) 
            isotope_table[rowN,'isotopes'] <- i
          }
          
          if (length(tmp_features_to_assign) != length(possible_istopes1)){
            
            mz_value1 <- mz_value - 3 - ppm_mz
            mz_value2 <- mz_value + 3 + ppm_mz
            possible_istopes1 <- rownames(possible_isotopes[possible_isotopes$mz >= mz_value1 & possible_isotopes$mz <= mz_value2 ,])
            
            for (ii in 1:length(possible_istopes1)){
              rowN <- which(rownames(isotope_table) == possible_istopes1[ii]) 
              isotope_table[rowN,'isotopes'] <- i
            } 
          }}}}
    
    ############################
    # checking the intensities #
    ############################
    variable_metadata <- variable_metadata[ order(row.names(variable_metadata)), ]
    isotope_table <- isotope_table[ order(row.names(isotope_table)), ]
    
    
    variable_metadata <- cbind(variable_metadata, isotope_group =isotope_table$isotopes)
    
    features_to_keep <- list()
    unique_groups <- unique(isotope_table$isotopes)
    
    
    # chose the monoisopic mass based on the highest intensities#
    i <- 0
    while (i < length(unique_groups)){
      i <- 1 + i
      group_tot_test <- unique_groups[[i]]
      totest <- isotope_table[isotope_table$isotopes ==  group_tot_test,]
      #print(totest)
      N <- rownames(totest)
      
      if (length(N) == 1){
        features_to_keep[i] <- N
        
      }else if (length(N) >= 1)  {
        N <- rownames(totest[totest$sumINT  == max(totest$sumINT ),])
        
      }else{
        break
      }
      #print(N)
      #print(group_tot_test)
      features_to_keep[i] <- N
    }  
    
    print(length(features_to_keep))
    ##################################################################
    
    
    
    features_to_keep <- unlist(features_to_keep)
    # update variable_metadata
    
    
    # updating variable metadata and datamatrix
    variable_metadata <- variable_metadata[(rownames(variable_metadata) %in% features_to_keep),]
    data_matrix <- data_matrix[ (rownames(data_matrix) %in% features_to_keep), ]
    
    
    
    rm(data_orginal)
    rm(isotope_table)
    rm(possible_isotopes)
  }
  if (isTRUE(Filter_on_RSD)){
    
    ####################################
    # filtering on RSD#
    ####################################
    metadata1 <- as.data.frame(t(sample_metadata))
    
    # order to be sure
    metadata1<- metadata1[ , order(names(metadata1))]
    intensities<- data_matrix[ , order(names(data_matrix))]
    intensities <- 10**intensities
    Index <- which(rownames(metadata1) == "sample_type") 
    sample_types  <- metadata1[Index, ]  ## subsets dataframe
    
    colnames(sample_types) <- colnames(intensities)
    combined <- rbind(sample_types,intensities)
    
    write.csv(combined, 'import for RSD filtering')
    
    
    try(rm(mSet))
    mSet<-InitDataObjects("pktable", "stat", FALSE)
    mSet<-Read.TextData(mSet, 'import for RSD filtering', "col", "disc")
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet)
    if (Filter_on_RSD == TRUE){
      mSet<-FilterVariable(mSet, "rsd", "T", RSD_treshhold)
    }
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ref = names[1,])
    
    normalized_set <- mSet[["dataSet"]][["norm"]]
    
    # updating the variable_metadata and data_matrix
    
    data_matrix <- as.data.frame(t(normalized_set[order(row.names(normalized_set)), ]))
    variable_metadata <- variable_metadata[ (rownames(variable_metadata) %in% rownames(data_matrix)), ]
  }
  if (isTRUE(noise_reduction)){
    
    ################################################
    # filtering based on intensity of blank samples#
    ################################################
    # get names of blank samples
    
    blank_names <- rownames(sample_metadata[sample_metadata$sample_type == 'blank',])
    
    data <- (data_matrix) + log10(2)
    idx <- c()
    for (i in 1:(length(blank_names))){
      idx <- c(idx,grep(blank_names[[i]], colnames(data)))
    }
    blanks <- data[,idx] 
    blanks <- rowMeans(blanks)
    for (i in 1:length(rownames(data))){
      data_matrix[i,][data_matrix[i,] <= blanks[i]] <- NA   
    }
    
    #### Return the values to Original ###
    rm(data)
    
    
    # editing the variable metadata
    features_in_sample <- unique(sample_metadata$sample_type)
    for (i in 1:length(rownames(data_matrix))){
      for (ii in 1:length(features_in_sample)){
        data <- data_matrix[i,][rownames(sample_metadata[sample_metadata$sample_type == features_in_sample[ii],])]
        
        data <- sum(!is.na(data))
        variable_metadata[i, ][features_in_sample[ii]] <- data
        
      }
    }
    
    
    # filter on amount of hits #
    features_to_delete <- rownames(subset(variable_metadata,variable_metadata['sample'] < min_amoutn_of_sample_hit ))
    variable_metadata <- variable_metadata[ !(rownames(variable_metadata) %in% features_to_delete), ]
    data_matrix <- data_matrix[ !(rownames(data_matrix) %in% features_to_delete), ]
    
    
    
    ######### REplacing  NA with 1/5 of lowest found intensity
    
    for (i in 1:length(rownames(data_matrix))){
      data1 <- data_matrix[i,]
      data <- data1[!is.na(data1)]
      data <- min(data)-log10(5)
      data <- data1[is.na(data1)] <- data
      data_matrix[i,] <- data1
      
    }
  }
  
  
  ####################################
  # CREATE DATAFRAMES IN RIGHT FORMAT#
  ####################################
  
  #retrieving output name #
  
  tmp <- as.list(strsplit(gsub(".tsv*","",sample_metadata_path), '/')[[1]])
  tmp <- tmp[(length(tmp))]
  tmp <- paste( unlist(tmp), collapse='')
  tmp <- as.list(strsplit(tmp,"_")[[1]])
  
  tmp<- tmp[4:(length(tmp))]
  tmp <- paste( unlist(tmp), collapse='_')
  
  path <- paste0(outputfolder,tmp,'_pre_univariate_filterd_matrix.tsv')
  write.table(data_matrix,path , sep ='\t')
  line <- readLines(path)
  line[1] <- paste0('""\t',line[1])
  writeLines(line,path)
  
  
  path <- paste0(outputfolder,tmp,'_pre_univariate_filterd_samplemetadata.tsv')
  write.table(sample_metadata, path, sep ='\t')
  line <- readLines(path)
  line[1] <- paste0('""\t',line[1])
  writeLines(line,path)
  
  
  path <- paste0(outputfolder,tmp,'_pre_univariate_filterd_variable_metadata.tsv')
  write.table(variable_metadata, path, sep ='\t')
  line <- readLines(path)
  line[1] <- paste0('"Feature_ID"\t',line[1])
  writeLines(line,path)
  
########################
  # cleaning up#
######################

detach_package(dplyr)
detach_package( randomcoloR)
detach_package(reshape2)
detach_package(MetaboAnalystR)

rm(combined, mSet, normalized_set, data_matrix, original_feat, possible_istopes1,
   tmp_features_to_assign, feature_list, features_to_keep, features_to_delete, features_in_sample,
   intensities, metadata1, variable_metadata, line, blanks, sample_types, data1, sample_metadata1
   , sample_metadata, unique_groups, means2, means)

##############
# univariate #
##############



library(phenomis)
library(biosigner)
library(ropls)
library(dendextend)
library(data.table)
library(dplyr)
library(stringr)

###############################
# variables #
##########################

input_folder <- outputfolder
Corrected_data_matrix <- paste0(tmp,'_pre_univariate_filterd_matrix.tsv')
sample_meta_data <- paste0(tmp,'_pre_univariate_filterd_samplemetadata.tsv')
variable_meta_data <- paste0(tmp,'_pre_univariate_filterd_variable_metadata.tsv')


test <- c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman","limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[1]


variables_of_interest <- "gender,age,bmi"
correcting_data_set_according_variable <- "NULL"
factor_of_interest <- 'gender'
second_factor_of_interest <- 'age'

P_value_treshhold <- 0.05
max_features_output <- "NULL"
graph_title <- 'univariate testing is super cool'
pre_fix_of_report <- "this makes sure that it will not overwrite other reports"



pixelsize1 <- 20
pixelsize2 <- 12


######### filtering data ##########
# keeps data who is bigger then these numbers#
cutoff_hotelPval <-  0.001
cutoff_missPval <-  0.001
cutoffDecipval <- 0.001


######## heatmap ###########
cluster_groups_of_samples<- 5
cluster_groups_of_features<- 5

heat_map_statistics <- c("pearson", "kendall", "spearman")[1]
heat_map_algorithm <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
                        "1-cor", "1-abs(cor)")[7]






###############################
# correcting some dataformats #
###############################
if (correcting_data_set_according_variable == "NULL"){
  correcting_data_set_according_variable <- NULL
}

if (max_features_output == "NA"){
  max_features_output <- NA
}


try(variables_of_interest <- strsplit(variables_of_interest,","))

items_to_filter <- c('blank', 'QC')


##############################################
sampleType <- "sample_type"

output_folder  <- paste0(input_folder,'/',test,'_',factor_of_interest,'/')
dir.create(output_folder, showWarnings = T)
setwd(output_folder)

numerical_factor_of_interest <- second_factor_of_interest

# create output name#
result <- paste0(test,'_',factor_of_interest)

# remove file to resolve bug
try(file.remove(paste0(output_folder,result,'.txt')))
#retrieve best info
path_Corrected_data_matrix <- paste0(input_folder,'/',Corrected_data_matrix)
path_sample_meta_data <- paste0(input_folder,'/',sample_meta_data)
path_variable_meta_data <- paste0(input_folder,'/',variable_meta_data)




#############################################
# changing datasets for better graph display#
#############################################
sample_metadata <- (read.table(path_sample_meta_data, sep = '\t', header = TRUE, row.names = 1))
data_matrix <- (read.table(path_Corrected_data_matrix, sep = '\t', header = TRUE, row.names = 1))


#######################################################################
# edit datamatrix 1/5 of lowest feature in dataset will be replace NA#
#######################################################################

for (i in 1:length(rownames(data_matrix))){
  data1 <- data_matrix[i,]
  data <- data1[!is.na(data1)]
  data <- min(data)-log10(5)
  data <- data1[is.na(data1)] <- data
  data_matrix[i,] <- data1
}










factor <- sample_metadata[, factor_of_interest]
sample_names <- rownames(sample_metadata)
eset_names <- list()
for (i in 1:length(factor)){
  string <- paste0(sample_names[[i]],"-",factor[[i]])
  eset_names[i] <- string
  
}


colnames(data_matrix) <- eset_names
rownames(sample_metadata) <- eset_names

datamatrix_path <- paste0(output_folder,Corrected_data_matrix,'_tmp_',factor_of_interest,'.tsv')
write.table(data_matrix, datamatrix_path, sep ='\t')

sample_metadata_path <- paste0(output_folder,sample_meta_data,'_tmp_',factor_of_interest,'.tsv')
write.table(sample_metadata, sample_metadata_path, sep ='\t')


line <- readLines(datamatrix_path,)
line[1] <- paste0('""\t',line[1])
writeLines(line,datamatrix_path,)


line <- readLines(sample_metadata_path,)
line[1] <- paste0('""\t',line[1])
writeLines(line,sample_metadata_path, )

##################
# reading in data#
##################


set <- phenomis::reading(NA,
                         files.ls = list(dataMatrix = file.path(datamatrix_path),
                                         sampleMetadata = file.path(sample_metadata_path),
                                         variableMetadata = file.path(path_variable_meta_data)))


png(filename = paste0(output_folder,result,'overview.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)



sacurine.eset <- phenomis::inspecting(set,
                                      pool_as_pool1.l = FALSE,
                                      pool_cv.n = 0.3,
                                      span.n = 1,
                                      sample_intensity.c = c("median", "mean")[1],
                                      title.c = "Inital Overview",
                                      plot_dims.l = TRUE,
                                      figure.c = c("none", "interactive", "myfile.pdf")[2],
                                      
                                      report.c = c("none", "interactive", "myfile.txt")[3])

dev.off()

#######################
# Filtering of samples#
#######################


if (!is.null(items_to_filter)){
  for (i in 1:length(items_to_filter)){
    set <- set[, Biobase::pData(set)[, sampleType] != items_to_filter[i]]
    
  }
  png(filename = paste0(output_folder,result,'overview_after_sample_filtering.png'),
      width = 960, height = 960, units = "px", pointsize = pixelsize1,
      bg = "white",  res = NA)
  set <- sacurine.eset <- phenomis::inspecting(set,
                                               pool_as_pool1.l = FALSE,
                                               pool_cv.n = 0.3,
                                               span.n = 1,
                                               sample_intensity.c = c("median", "mean")[1],
                                               title.c = "Inital Overview",
                                               plot_dims.l = TRUE,
                                               figure.c = c("none", "interactive", "myfile.pdf")[2],
                                               
                                               report.c = c("none", "interactive", "myfile.txt")[3])
  dev.off()
}


####################################################
#correcting dataset according variable information##
####################################################
if(!(is.null(correcting_data_set_according_variable))){
  Biobase::exprs(set) <- sweep(Biobase::exprs(set),
                               2,
                               Biobase::pData(set)[, correcting_data_set_according_variable],
                               "/")
  set <- phenomis::inspecting(set)
}

#####################################################
######### Filtering according data statistics########
#####################################################

set <- set[, Biobase::pData(set)[, "hotel_pval"] >= cutoff_hotelPval &
             Biobase::pData(set)[, "miss_pval"] >= cutoff_missPval &
             Biobase::pData(set)[, "deci_pval"] >= cutoffDecipval]

phenomis::inspecting(set)


######################
# Univariate Testing #
######################

png(filename = paste0(output_folder,result,'signif_chart.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize2,
    bg = "white",  res = NA,
)


univariate_set <- phenomis::hypotesting(
  set,
  test.c = test,
  factor_of_interest,
  adjust_thresh.n = P_value_treshhold,
  signif_maxprint.i = max_features_output,
  title.c = 'Significant levels of univariate testing',
  prefix.c = pre_fix_of_report,
  adjust.c = "none",
  report.c = paste0(result,'.txt'),
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2])




dev.off()


################################################
##################### PCA ######################
################################################

png(filename = paste0(output_folder,result,'_Global_PCA.png'),
    width = 960, height = 960, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)
sacPca <- ropls::opls(univariate_set, info.txt = 'interactive',
                      y = NULL,
                      predI = NA,
                      orthoI = 0,
                      algoC = c("default", "nipals", "svd")[1],
                      crossvalI = 1,
                      log10L = FALSE,
                      permI = 7,
                      scaleC = c("none", "center", "pareto", "standard")[4],
                      subset = NULL,
                      plotSubC = 'Variance between Principal Components',
                      fig.pdfC = c("none", "interactive", "myfile.pdf")[2],
                      info.txtC = c("none", "interactive", "myfile.txt")[2],
                      printL = TRUE,
                      plotL = TRUE,
                      .sinkC = NULL,)

dev.off()

#update variables_of_interest
variables_of_interest <- unlist(variables_of_interest)
variables_of_interest <- variables_of_interest[variables_of_interest != factor_of_interest]


for (i in 1:length(variables_of_interest)){
  print(i)
  png(filename = paste0(output_folder,result,'_',variables_of_interest[i],'.png'),
      width = 960, height = 960, units = "px", pointsize = pixelsize1,
      bg = "white",  res = NA,
  )
  
  ropls::plot(sacPca,
              parAsColFcVn = Biobase::pData(univariate_set)[, variables_of_interest[i]],
              typeVc = "x-score",
              
              figure.c = paste0(result,'.pdf'),
              plotSubC = variables_of_interest[i])
  
  
  
  dev.off()
}








########################
##### heatmap ##########
########################



png(filename = paste0(output_folder,result,'_heatmap.png'),
    width = 1840, height = 1840, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)

sacurine.eset <- ropls::getEset(sacPca)
sacurine.eset <- phenomis::clustering(sacurine.eset, correl.c = heat_map_statistics,
                                      cex.vn = c(1,1),
                                      clusters.vi = c(cluster_groups_of_samples, cluster_groups_of_features),
                                      dissym.c = heat_map_algorithm ,
                                      agglo.c = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
                                                  "median", "centroid")[2],
                                      title.c = paste0('Clustering of Features and sampels ',heat_map_statistics))




dev.off()

###########################################
########### Supervised modeling############
########### (O)PLS(-DA) modeling###########
###########################################
png(filename = paste0(output_folder,result,'_Testing.png'),
    width = 1840, height = 1840, units = "px", pointsize = pixelsize1,
    bg = "white",  res = NA,
)

sacPlsda <- ropls::opls(sacurine.eset, predI = 5, factor_of_interest, seedI = 123)
sacurine.eset <- ropls::getEset(sacPlsda)




dev.off()


##############################################################
# rewriting the output for filtering#
##############################################################
######## removing tmp files ############


phenomis::writing(sacurine.eset, dir.c = getwd(),overwrite.l = TRUE)


if (!is.null(items_to_filter)){
  library(plyr)
  metadata <- read.table(path_sample_meta_data, sep = '\t', header = TRUE, row.names = 1)
  metadata1 <- read.table(paste0(output_folder,"/sampleMetadata.tsv"), sep = '\t', header = TRUE, row.names = 1)
  colnames1 <- setdiff(colnames(metadata1), colnames(metadata))
  
  metadata <- metadata[  metadata$sample_type %in% items_to_filter, ]
  
  
  for (i in 1:length(colnames1)){
    metadata[[colnames1[i]]] <- NA
  }
  
  
  metadata2 <- rbind(metadata1, metadata)
  write.table(metadata2, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  file.remove(paste0(output_folder,"/sampleMetadata.tsv"))
  
  sample_meta_data1 <- read.table(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), header = TRUE, row.names = 1)
  rownames(sample_meta_data1) <- gsub('(.*)-\\w+', '\\1', rownames(sample_meta_data1))
  write.table(sample_meta_data1, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  
  
  file.copy(path_Corrected_data_matrix,paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv') )
  file.copy(paste0(output_folder,"/dataMatrix.tsv"), paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  file.remove(paste0(output_folder,"/dataMatrix.tsv"))
  
  line <- readLines(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  
} else if(is.null(items_to_filter)) {
  file.rename(paste0(output_folder,"/sampleMetadata.tsv"), paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  
  sample_meta_data1 <- read.table(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), header = TRUE, row.names = 1)
  rownames(sample_meta_data1) <- gsub('(.*)-\\w+', '\\1', rownames(sample_meta_data1))
  write.table(sample_meta_data1, paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'), sep = '\t')
  
  file.rename(paste0(output_folder,"/dataMatrix.tsv"), paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  
  ####################
  # correcting output#
  ####################
  
  line <- readLines(paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line,paste0(output_folder,sample_meta_data,'_',test,"_",factor_of_interest,'.tsv'))
  
  line <- readLines(paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  line[1] <- paste0('""\t',line[1])
  writeLines(line,paste0(output_folder,Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv'))
  
}





file.rename(paste0(output_folder,"/variableMetadata.tsv"), paste0(output_folder,variable_meta_data,'_',test,"_",factor_of_interest,'.tsv'))


file.remove(datamatrix_path)
file.remove(sample_metadata_path)

##################
# cleaning up#
#################
detach_package(phenomis)
detach_package(biosigner)
detach_package(ropls)
detach_package(dendextend)
detach_package(data.table)
detach_package(dplyr)
detach_package(stringr)

rm(sacPca, sacPlsda, sacurine.eset, univariate_set, set, metadata2, line
   , sample_metadata1, metadata, metadata1, data1, eset_names, sample_metadata1, sample_names, totest)
############################
# post univariate filtering#
############################

#####################
#  filter the data#
####################
#install.packages("randomcoloR")
library(dplyr)
library( randomcoloR)
library(reshape2)


folder <- output_folder
variable_metadata <- paste0(variable_meta_data,'_',test,"_",factor_of_interest,'.tsv')
data_matrix <- paste0(Corrected_data_matrix,'_',test,"_",factor_of_interest,'.tsv')
sample_metadata <- paste0(sample_meta_data,'_',test,"_",factor_of_interest,'.tsv')


outputfolder_name <- "filterd"

filter_significant_hits = TRUE
filter_QC_SAMPLE_RATIO = TRUE
FILTER_blank_SAMPLE_RATIO = FALSE


minium_QC_sample_cutoff_ratio <- 0.1
maximum_QC_sample_cutoff_ratio <- 2.0

minium_QC_blank_cutoff_ratio <- 0.0
maximum_QC_blank_cutoff_ratio <- 1.0



#####################
# creating the paths#
#####################
outputfolder <- (paste0(folder,'/',outputfolder_name,'/'))
dir.create(outputfolder, showWarnings = T)




Variable_metadata_path <- paste0(folder,'/' ,variable_metadata)
data_matrix_path <- paste0(folder, '/' ,data_matrix)
sample_metadata_path <- paste0(folder, '/' ,sample_metadata)

####################
# reading in data###
####################


variable_metadata <- read.table(Variable_metadata_path, sep = '\t', header = TRUE, row.names = 1)
data_matrix <- read.table(data_matrix_path, sep = '\t', header = TRUE, row.names = 1)
sample_metadata <- (read.table(sample_metadata_path, sep = '\t', header = TRUE, row.names = 1))

##################
# debugging######
#################


colnames(data_matrix) <- rownames(sample_metadata)


###########################
#  ordering the dataframes#
###########################
data_matrix <- data_matrix[ , order(names(data_matrix))]
sample_metadata <- sample_metadata[order(row.names(sample_metadata)), ]

#######################################
# creating a boxplot of found features#
#######################################
features_in_sample <- unlist(unique(sample_metadata$sample_type))
features_in_sample <- variable_metadata[ , (names(variable_metadata) %in% features_in_sample)]

palette <-  palette(rainbow((length(colnames(features_in_sample)))))
png(filename = paste0(outputfolder,"Features count in each sample type before filtering.png"),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
barplot(t(as.matrix(features_in_sample)),beside=TRUE,
        main = "Features count in each sample type before filtering",
        col = c(palette)
)
legend("topleft",
       c(names(features_in_sample)),
       fill = c(palette)
)


dev.off()
################################
# filter based significant hits#
################################

if (isTRUE(filter_significant_hits)){
  signif_column_nr <- which( colnames(variable_metadata)==colnames(variable_metadata)[grepl('signif', colnames(variable_metadata))])
  variable_metadata <- variable_metadata[variable_metadata[,signif_column_nr] == 1, ] 
  
  filter_rowname <- rownames(variable_metadata)
  data_matrix <- data_matrix %>% filter(row.names(data_matrix) %in% filter_rowname)
}


######################################################
# filtering on prevelance QC/sample hits in features##
######################################################


if (isTRUE(filter_QC_SAMPLE_RATIO)){
  features_in_sample <- unique(sample_metadata$sample_type)
  features_in_sample <- variable_metadata[ , (names(variable_metadata) %in% features_in_sample)]
  
  
  features_in_sample1 <- features_in_sample %>% select('QC', 'sample')
  features_in_sample1$ratio <- (features_in_sample1[['QC']] / features_in_sample1[['sample']])
  features_to_delete <- rownames(features_in_sample1[features_in_sample1[,"ratio"] <= minium_QC_sample_cutoff_ratio, ])
  features_to_delete <- append(features_to_delete,rownames(features_in_sample1[features_in_sample1[,"ratio"] >= maximum_QC_sample_cutoff_ratio, ]))
  
  variable_metadata <- variable_metadata[ !(rownames(variable_metadata) %in% features_to_delete), ]
  data_matrix <- data_matrix[ !(rownames(data_matrix) %in% features_to_delete), ]
  
}

######################################################
# filtering on prevelance QC/blank hits in features##
######################################################
if(isTRUE(FILTER_blank_SAMPLE_RATIO)){
  features_in_sample <- features_in_sample %>% select('blank', 'QC')
  features_in_sample$ratio <- (features_in_sample[['blank']] / features_in_sample[['QC']])
  features_to_delete <- rownames(features_in_sample[features_in_sample[,"ratio"] <= minium_QC_blank_cutoff_ratio, ])
  features_to_delete <- append(features_to_delete,rownames(features_in_sample[features_in_sample[,"ratio"] >= maximum_QC_blank_cutoff_ratio, ]))
  
  
  variable_metadata <- variable_metadata[ !(rownames(variable_metadata) %in% features_to_delete), ]
  data_matrix <- data_matrix[ !(rownames(data_matrix) %in% features_to_delete), ]
}



#######################################
# creating a boxplot of found features#
#######################################
features_in_sample <- unique(sample_metadata$sample_type)
features_in_sample <- variable_metadata[ , (names(variable_metadata) %in% features_in_sample)]
palette <-  palette(rainbow((length(names(features_in_sample)))))
png(filename = paste0(outputfolder,"Features count in each sample type after filtering.png"),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
barplot(t(as.matrix(features_in_sample)),beside=TRUE,
        main = "Features count in each sample type after filtering",
        col = c(palette)
)
legend("topleft",
       c(names(features_in_sample)),
       fill = c(palette)
)


dev.off()



####################################
# CREATE DATAFRAMES IN RIGHT FORMAT#
####################################

#retrieving output name #

tmp <- as.list(strsplit(gsub(".tsv*","",data_matrix_path), '/')[[1]])
tmp <- tmp[(length(tmp))]
tmp <- paste( unlist(tmp), collapse='')
tmp <- as.list(strsplit(tmp,"_")[[1]])

tmp<- tmp[4:(length(tmp))]
tmp <- paste( unlist(tmp), collapse='_')

path <- paste0(outputfolder,tmp,'_filterd_data_matrix.tsv')
write.table(data_matrix,path , sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)

tmp <- as.list(strsplit(gsub(".tsv*","",sample_metadata_path), '/')[[1]])
tmp <- tmp[(length(tmp))]
tmp <- paste( unlist(tmp), collapse='')
tmp <- as.list(strsplit(tmp,"_")[[1]])

tmp<- tmp[4:(length(tmp))]
tmp <- paste( unlist(tmp), collapse='_')
path <- paste0(outputfolder,tmp,'_filterd_sample_metadata.tsv')
write.table(sample_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


tmp <- as.list(strsplit(gsub(".tsv*","",Variable_metadata_path), '/')[[1]])
tmp <- tmp[(length(tmp))]
tmp <- paste( unlist(tmp), collapse='')
tmp <- as.list(strsplit(tmp,"_")[[1]])

tmp<- tmp[4:(length(tmp))]
tmp <- paste( unlist(tmp), collapse='_')

path <- paste0(outputfolder,tmp,'_filterd_variable_metadata.tsv')
write.table(variable_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('"Feature_ID"\t',line[1])
writeLines(line,path)


  
  
  
