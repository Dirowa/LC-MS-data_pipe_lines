#####################
#  filter the data#
####################
#install.packages("randomcoloR")
library(dplyr)
library( randomcoloR)
library(reshape2)

Variable_metadata_path <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default/batch_correction/ttest_gender/variableMetadata.tsv"
data_matrix_path <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default/batch_correction/ttest_gender/Best_batchcorrected.tsv_ttest_gender.tsv"
sample_metadata_path <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default/batch_correction/ttest_gender/sample_meta_data_XCMS_default.tsv_batchcorrected.tsv_ttest_gender.tsv"

outputfolder <- "F:/avans/stage MM/Sherloktest_data_2/peakpick_output_XCMS_default/batch_correction/ttest_gender/filterd/"

sampletype <- "sampleType"

dir.create(outputfolder, showWarnings = T)

minium_QC_sample_cutoff_ratio <- 0.1
maximum_QC_sample_cutoff_ratio <- 2.0

minium_QC_blank_cutoff_ratio <- 0.0
maximum_QC_blank_cutoff_ratio <- 1.0

sample_in_sampleType <- 'sample'
blank_in_sampleType <- 'blank'
####################
# reading in data###
####################


variable_metadata <- read.table(Variable_metadata_path, sep = '\t', header = TRUE, row.names = 1)
data_matrix <- read.table(data_matrix_path, sep = '\t', header = TRUE, row.names = 1)
sample_metadata <- (read.table(sample_metadata_path, sep = '\t', header = TRUE, row.names = 1))


#this is tmp and needs to be removed
names(variable_metadata)[names(variable_metadata) == 'pool'] <- 'QC'




###########################
#  ordering the dataframes#
###########################
data_matrix <- data_matrix[ , order(names(data_matrix))]
sample_metadata <- sample_metadata[order(row.names(sample_metadata)), ]


################################
# filter based significant hits#
################################

signif_column_nr <- which( colnames(variable_metadata)==colnames(variable_metadata)[grepl('signif', colnames(variable_metadata))])
variable_metadata <- variable_metadata[variable_metadata[,signif_column_nr] == 1, ] 

filter_rowname <- rownames(variable_metadata)
data_matrix <- data_matrix %>% filter(row.names(data_matrix) %in% filter_rowname)

#######################################
# creating a boxplot of found features#
#######################################
features_in_sample <- unique(sample_metadata[[sampletype]])
features_in_sample <- variable_metadata[ , (names(variable_metadata) %in% features_in_sample)]
palette <-  palette(rainbow((length(names(features_in_sample)))))
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

################################################
# filtering based on intensity of blank samples#
################################################
# get names of blank samples
blank_names <- rownames(sample_metadata)[which(sample_metadata == 'blank', arr.ind=T)[, "col"]]
data <- (data_matrix) + log10(2)

idx <- match(blank_names, names(data))
idx <- sort(c(idx-1, idx))
blanks <- data[,idx] 
blanks <- rowMeans(blanks)
  

for (i in 1:length(rownames(data))){
  data[i,][data[i,] <= blanks[i]] <- NA   
}
  
data_matrix <- data

# editing the variable metadata
features_in_sample <- unique(sample_metadata[[sampletype]])
for (i in 1:length(rownames(data_matrix))){
  for (ii in 1:length(features_in_sample)){
    data <- data_matrix[i,][rownames(sample_metadata[sample_metadata[sampletype] == features_in_sample[ii],])]

    data <- sum(!is.na(data))
    variable_metadata[i, ][features_in_sample[ii]] <- data
    
  }
}





######################################################
# filtering on prevelance QC/sample hits in features##
######################################################
features_in_sample <- unique(sample_metadata[[sampletype]])
features_in_sample <- variable_metadata[ , (names(variable_metadata) %in% features_in_sample)]


features_in_sample1 <- features_in_sample %>% select('QC', sample_in_sampleType)
features_in_sample1$ratio <- (features_in_sample1[['QC']] / features_in_sample1[[sample_in_sampleType]])
features_to_delete <- rownames(features_in_sample1[features_in_sample1[,"ratio"] <= minium_QC_sample_cutoff_ratio, ])
features_to_delete <- append(features_to_delete,rownames(features_in_sample1[features_in_sample1[,"ratio"] >= maximum_QC_sample_cutoff_ratio, ]))

variable_metadata <- variable_metadata[-features_to_delete,]

variable_metadata <- variable_metadata[ !(rownames(variable_metadata) %in% features_to_delete), ]
data_matrix <- data_matrix[ !(rownames(data_matrix) %in% features_to_delete), ]



######################################################
# filtering on prevelance QC/blank hits in features##
######################################################
features_in_sample <- features_in_sample %>% select(blank_in_sampleType, 'QC')
features_in_sample$ratio <- (features_in_sample[[blank_in_sampleType]] / features_in_sample[['QC']])
features_to_delete <- rownames(features_in_sample[features_in_sample[,"ratio"] <= minium_QC_blank_cutoff_ratio, ])
features_to_delete <- append(features_to_delete,rownames(features_in_sample[features_in_sample[,"ratio"] >= maximum_QC_blank_cutoff_ratio, ]))

variable_metadata <- variable_metadata[-features_to_delete,]

variable_metadata <- variable_metadata[ !(rownames(variable_metadata) %in% features_to_delete), ]
data_matrix <- data_matrix[ !(rownames(data_matrix) %in% features_to_delete), ]




#######################################
# creating a boxplot of found features#
#######################################
features_in_sample <- unique(sample_metadata[[sampletype]])
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

tmp <- as.list(strsplit(gsub(".tsv*","",sample_metadata_path), '/')[[1]])
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


path <- paste0(outputfolder,tmp,'_filterd_sample_metadata.tsv')
write.table(sample_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('""\t',line[1])
writeLines(line,path)


path <- paste0(outputfolder,tmp,'_filterd_variable_metadata.tsv')
write.table(variable_metadata, path, sep ='\t')
line <- readLines(path)
line[1] <- paste0('"Feature_ID"\t',line[1])
writeLines(line,path)


