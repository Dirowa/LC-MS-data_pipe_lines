#####################
#  filter the data#
####################
#install.packages("randomcoloR")
library(dplyr)
library( randomcoloR)
library(reshape2)


folder <- "F:/avans/stage MM/Sherloktest_data_3/SherLOKdata_processed_XCMS_IPO_fitgaus_2L/batch_correction/Noice_reducted/ttest_sex"
variable_metadata <- "XCMS_IPO_fitgaus_2L_batchcorrected_noice_reduced_variable_metadata.tsv_ttest_sex.tsv"
data_matrix <- "XCMS_IPO_fitgaus_2L_batchcorrected_noice_reduced_matrix.tsv_ttest_sex.tsv"
sample_metadata <-"XCMS_IPO_fitgaus_2L_batchcorrected_noice_reduced_sample_metadata.tsv_ttest_sex.tsv"


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
   data_matrix <- data_matrix[row.names(data_matrix)%in%filter_rowname,]

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
  
  try(variable_metadata <- variable_metadata[-features_to_delete,])
  
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

