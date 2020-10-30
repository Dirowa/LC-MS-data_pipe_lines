####################
# mait package     #
####################
# all files should be in one directory
# in this directory are child directorys of the class name
# in these child directory's are corrospendend mzML files

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAIT")


###################
# library's########
###################
library(MAIT)
library(faahKO)

#####################
#importing data path#
#####################

setwd( "F:/avans/stage MM/test_data_sherlok/dataxcms/")

####################################
# pre=processing  #################
##################################

MAIT <- sampleProcessing(dataDir = "F:/avans/stage MM/test_data_sherlok/dataxcms/", project = "MAIT_Demo_7", 
                         snThres = 5,
                         Sigma = 5/2.3548,
                         mzSlices = 0.3,
                         retcorrMethod = "loess",
                         groupMethod = "density",
                         bwGroup = 3,
                         mzWidGroup = 0.25,
                         filterMethod = "centWave",
                         prefilter = c(3,3000),
                         rtStep = 0.03,
                         nSlaves = 0,
                         minfrac = 0.5,
                         minsamp = 1,
                         peakwidth = c(5, 20),
                         ppm = 10,
                         family = c("gaussian", "symmetric"),
                         span = 0.2,
                         fwhm = 30
)

MAIT <- peakAnnotation(MAIT.object = MAIT,corrWithSamp = 0.7, 
                       corrBetSamp = 0.75,perfwhm = 0.6)

table <- getScoresTable


write.csv(table, "mait_output.csv")
