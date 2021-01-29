#####################
#  filter the noise#
####################
#install.packages("randomcoloR")
library(dplyr)
library( randomcoloR)
library(reshape2)
library(MetaboAnalystR)

folder <- '!@#$%^&1&^%$#@!'
outputfolder_name <- "!@#$%^&2&^%$#@!"



variable_metadata <- "!@#$%^&5&^%$#@!"
data_matrix <- "!@#$%^&3&^%$#@!"
sample_metadata <-"!@#$%^&4&^%$#@!"


Filter_on_isotopes <- !@#$%^&6&^%$#@!
Filter_on_RSD <- !@#$%^&7&^%$#@!
noise_reduction <- !@#$%^&8&^%$#@!

### isotope reduction #
ppm_mz = !@#$%^&9&^%$#@!
ppm_rt = !@#$%^&10&^%$#@!

### RSD filtering ###
RSD_treshhold <- !@#$%^&11&^%$#@!

### noise filtering ####
min_amoutn_of_sample_hit <- !@#$%^&12&^%$#@!



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

path <- paste0(outputfolder,tmp,'_pre_univriate_filterd_matrix.tsv')
write.table(data_matrix,path , sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


path <- paste0(outputfolder,tmp,'_pre_univriate_filterd_sample_metadata.tsv')
write.table(sample_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


path <- paste0(outputfolder,tmp,'_pre_univriate_filterd_variable_metadata.tsv')
write.table(variable_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('"Feature_ID"\t',line[1])
writeLines(line,path)





