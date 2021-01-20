############
# Variables#
############

library( phenomis)
library(ropls)
library(dendextend)
library(data.table)
library(dplyr)
library(stringr)


input_folder <- "!@#$%^&1&^%$#@!"
Corrected_data_matrix <- "!@#$%^&2&^%$#@!"
sample_meta_data <- "!@#$%^&3&^%$#@!"
variable_meta_data <- "!@#$%^&4&^%$#@!"





variables_of_interest <- '!@#$%^&5&^%$#@!'
factor_of_interest <- '!@#$%^&6&^%$#@!'
second_factor_of_interest <- '!@#$%^&7&^%$#@!'
numerical_factor_of_interest <- '!@#$%^&8&^%$#@!'


knitr::opts_chunk$set(fig.width = 6,
                      fig.height = 6,
                      fig.path = 'figures/')

output_folder <- paste0(input_folder,'/','multivariate_',factor_of_interest,'_',second_factor_of_interest)
dir.create(output_folder, showWarnings = T)
output_folder<- paste0(output_folder,'/')
setwd(output_folder)


#####################
# editing variables #
####################

if (numerical_factor_of_interest == "NULL"){
  numerical_factor_of_interest <- NULL
}

try(variables_of_interest <- strsplit(variables_of_interest,",")[[1]])
##########
# import #
##########

path_Corrected_data_matrix <- paste0(input_folder,'/',Corrected_data_matrix)
path_sample_meta_data <- paste0(input_folder,'/',sample_meta_data)
path_variable_meta_data <- paste0(input_folder,'/',variable_meta_data)


######################
# retrieving dataset##
######################
dataMatrix <- t(read.table(path_Corrected_data_matrix, sep = '\t', header = T, row.names = 1))
sampleMetadata <- (read.table(path_sample_meta_data, sep = '\t', header = T, row.names = 1))
variableMetadata <- (read.table(path_variable_meta_data, sep = '\t', header = T, row.names = 1))




#####################################################
# correctiong dataset so that only samples are in it#
#####################################################

rownames_to_remove <- rownames(sampleMetadata %>% filter(sampleMetadata$sample_type == 'QC'| sampleMetadata$sample_type == 'blank'))
for (i in 1:length(rownames_to_remove)){
  sampleMetadata <- sampleMetadata[row.names(sampleMetadata) != rownames_to_remove[[i]], , drop = FALSE]
  dataMatrix <- dataMatrix[row.names(dataMatrix) != rownames_to_remove[[i]], , drop = FALSE]
  
}

#dropping uneeded columns#
sampleMetadata <- select(sampleMetadata, variables_of_interest)
#########################
# viewing the data #
########################







png(filename = paste0(output_folder,'intensities_across_samples','overview.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA
)
view(dataMatrix)

dev.off()


png(filename = paste0(output_folder,'samplemetadata', 'overview.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA
)
view(sampleMetadata)
dev.off()

png(filename = paste0(output_folder,'Variable_metadata','overview.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA
)
view(variableMetadata)
dev.off()
###################
# perfomring a PCA#
###################



PCA_data <- opls(dataMatrix, fig.pdfC = "none")
png(filename = paste0(output_folder,'overview.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA,
)
plot(PCA_data)
dev.off()




lagenda <- sampleMetadata[, factor_of_interest]



if (is.null(second_factor_of_interest)){
  
  png(filename = paste0(output_folder,'scoresPCA',factor_of_interest,'.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  
  plot(PCA_data,
       typeVc = "x-score",
       parAsColFcVn = lagenda)
  dev.off()
  
  
} else{
  png(filename = paste0(output_folder,factor_of_interest,'.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  
  plot(PCA_data,
       typeVc = "x-score",
       parAsColFcVn = lagenda)
  dev.off()
  
  png(filename = paste0(output_folder,factor_of_interest,'_',second_factor_of_interest,'.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  
  plot(PCA_data,
       typeVc = "x-score",
       parAsColFcVn = lagenda,
       parLabVc = as.character(sampleMetadata[, second_factor_of_interest]),
  )
  dev.off()
}

################################
# creating some different plots#
################################

sacurine.plsda <- opls(dataMatrix)

sacurine.plsda <- opls(dataMatrix, lagenda)


png(filename = paste0(output_folder,'scoresPCA',factor_of_interest,'_',second_factor_of_interest,'.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA
)
sacurine.oplsda <- opls(dataMatrix, lagenda,
                        predI = 1, orthoI = NA)
dev.off()



png(filename = paste0(output_folder,'scoresPCA',factor_of_interest,'_',second_factor_of_interest,'odd.png'),
    width = 960, height = 960, units = "px", pointsize = 12,
    bg = "white",  res = NA
)
sacurine.oplsda <- opls(dataMatrix, lagenda,
                        predI = 1, orthoI = NA,
                        subset = "odd")
dev.off()









try(for (i in 1){
## ----train--------------------------------------------------------------------
trainVi <- getSubsetVi(sacurine.oplsda)
table(lagenda[trainVi], fitted(sacurine.oplsda))

## ----test---------------------------------------------------------------------
table(lagenda[-trainVi],
      predict(sacurine.oplsda, dataMatrix[-trainVi, ]))

## ----overfit, echo = FALSE----------------------------------------------------
set.seed(123)
obsI <- 20
featVi <- c(2, 20, 200)
featMaxI <- max(featVi)
xRandMN <- matrix(runif(obsI * featMaxI), nrow = obsI)
yRandVn <- sample(c(rep(0, obsI / 2), rep(1, obsI / 2)))

layout(matrix(1:4, nrow = 2, byrow = TRUE))

})


try(for (featI in featVi) {
  
  png(filename = paste0(output_folder,'OPLS',featI,'.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  
  randPlsi <- opls(xRandMN[, 1:featI], yRandVn,
                   predI = 2,
                   permI = ifelse(featI == featMaxI, 100, 0),
                   fig.pdfC = "none",
                   info.txtC = "none")
  plot(randPlsi, typeVc = "x-score",
       parCexN = 1.3, parTitleL = FALSE,
       parCexMetricN = 0.5)
  mtext(featI/obsI, font = 2, line = 2)
  if (featI == featMaxI)
    plot(randPlsi,
         typeVc = "permutation",
         parCexN = 1.3)
  
  
  
  mtext(" obs./feat. ratio:",
        adj = 0, at = 0, font = 2,
        line = -2, outer = TRUE)
  
  
  
  
  
  dev.off()
}
)




## ----vip----------------------------------------------------------------------

if (!is.null(numerical_factor_of_interest)){
  ageVn <- sampleMetadata[, numerical_factor_of_interest]
  
  pvaVn <- apply(dataMatrix, 2,
                 function(feaVn) cor.test(ageVn, feaVn)[["p.value"]])
  
  vipVn <- getVipVn(opls(dataMatrix, ageVn,
                         predI = 1, orthoI = NA,
                         fig.pdfC = "none"))
  
  quantVn <- qnorm(1 - pvaVn / 2)
  rmsQuantN <- sqrt(mean(quantVn^2))
  
  opar <- par(font = 2, font.axis = 2, font.lab = 2,
              las = 1,
              mar = c(5.1, 4.6, 4.1, 2.1),
              lwd = 2, pch = 16)
  
  
  png(filename = paste0(output_folder,'significant_levels.png'),
      width = 960, height = 960, units = "px", pointsize = 12,
      bg = "white",  res = NA,
  )
  plot(pvaVn, vipVn,
       col = "red",
       pch = 16,
       xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")
  
  
  box(lwd = 2)
  
  curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)
  
  abline(h = 1, col = "blue")
  abline(v = 0.05, col = "blue")
  
  par(opar)
  dev.off()
  
  toW4M(sacSet, paste0(getwd(), "/out_"))
}
