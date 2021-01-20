#####################
#  filter the noise#
####################
#install.packages("randomcoloR")
library(dplyr)
library( randomcoloR)
library(reshape2)


folder <- '!@#$%^&1&^%$#@!'
outputfolder_name <- "!@#$%^&2&^%$#@!"



data_matrix<- "!@#$%^&3&^%$#@!"
sample_metadata <- "!@#$%^&4&^%$#@!"
variable_metadata<-"!@#$%^&5&^%$#@!"

min_amoutn_of_sample_hit <- !@#$%^&6&^%$#@!

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

path <- paste0(outputfolder,tmp,'_noice_reduced_matrix.tsv')
write.table(data_matrix,path , sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


path <- paste0(outputfolder,tmp,'_noice_reduced_sample_metadata.tsv')
write.table(sample_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


path <- paste0(outputfolder,tmp,'_noice_reduced_variable_metadata.tsv')
write.table(variable_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('"Feature_ID"\t',line[1])
writeLines(line,path)





