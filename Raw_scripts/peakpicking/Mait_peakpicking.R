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

#where subprocces are allowed to safe data
setwd( "F:/avans/stage MM/test_data_sherlok/dataxcms/")

# folder containg the data folders , 1 folder containging samples, 1 containing pooled
folder_containg_data <- 'F:/avans/stage MM/test_data_sherlok/dataxcms/data/'
output_folder_peaktable <- 'F:/avans/stage MM/peakpicking/MAIT_default.txt'
####################################
# pre=processing  #################
##################################
#data fles must not be in the same folder. every group ( sample, pooled, must be i different folders!)
MAIT <- sampleProcessing(dataDir = folder_containg_data, project = "M3")
                          
                         #snThres = 5,
                         #Sigma = 5/2.3548,
                         #mzSlices = 0.3,
                         #retcorrMethod = "loess",
                         #groupMethod = "density",
                         #bwGroup = 3,
                         #mzWidGroup = 0.25,
                         #filterMethod = "centWave",
                         #prefilter = c(3,3000),
                         #rtStep = 0.03,
                         #nSlaves = 0,
                         #minfrac = 0.5,
                         #minsamp = 1,
                         #peakwidth = c(5, 20),
                         #ppm = 10,
                         #family = c("gaussian", "symmetric"),
                         #span = 0.2,
                         #fwhm = 30
#)


#annotation
MAIT <- peakAnnotation(MAIT.object = MAIT)
#                       ,corrWithSamp = 0.7, 
                       #corrBetSamp = 0.75,perfwhm = 0.6)

# creating features
MAIT<- spectralSigFeatures(MAIT.object = MAIT,pvalue=0.05,
                           p.adj="none",scale=FALSE)


#creating the important tables
signTable <- sigPeaksTable(MAIT.object = MAIT, printCSVfile = FALSE)
write.csv(signTable, output_folder_peaktable)





#################### REST OF THE SCRIPT STILL NEEDS TO WORK OUT, not usefull for testing peakpicking#########


###################################################
### code chunk number 8: MAIT_Vignette.Rnw:282-283
###################################################
signTable <- sigPeaksTable(MAIT.object = MAIT, printCSVfile = FALSE)


###################################################
### code chunk number 9: MAIT_Vignette.Rnw:293-294
###################################################
MAIT


###################################################
### code chunk number 10: BoxHeatplots (eval = FALSE)
###################################################
## plotBoxplot(MAIT)
## plotHeatmap(MAIT)


###################################################
### code chunk number 11: PLSPCA
###################################################
MAIT<-plotPCA(MAIT,plot3d=FALSE)
MAIT<-plotPLS(MAIT,plot3d=FALSE)
PLSmodel <- model(MAIT, type = "PLS")
PCAmodel <- model(MAIT, type = "PCA")


###################################################
### code chunk number 12: showPLSmodel
###################################################
PLSmodel


###################################################
### code chunk number 13: showPLSscores
###################################################
pcaScores(MAIT)


###################################################
### code chunk number 14: resultsPath
###################################################
resultsPath(MAIT)


###################################################
### code chunk number 15: Biotransformations
###################################################
Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005) 


###################################################
### code chunk number 16: myBiotransf
###################################################
data(MAITtables)
myBiotransformation<-c("custom_biotrans",105.0)
myBiotable<-biotransformationsTable
myBiotable[,1]<-as.character(myBiotable[,1])
myBiotable<-rbind(myBiotable,myBiotransformation)
myBiotable[,1]<-as.factor(myBiotable[,1]) 
tail(myBiotable) 


###################################################
### code chunk number 17: identifyMetabolites
###################################################
MAIT <- identifyMetabolites(MAIT.object = MAIT, peakTolerance = 0.005)


###################################################
### code chunk number 18: metaboliteTable
###################################################
metTable<-metaboliteTable(MAIT)
head(metTable)


###################################################
### code chunk number 19: validation
###################################################
MAIT <- Validation(Iterations = 20, trainSamples= 3, 
                   MAIT.object = MAIT)


###################################################
### code chunk number 20: summaryMAIT2
###################################################
summary(MAIT)


###################################################
### code chunk number 21: classifRatioClasses
###################################################
classifRatioClasses(MAIT)


###################################################
### code chunk number 22: defData
###################################################
peaks <- scores(MAIT)
masses <- getPeaklist(MAIT)$mz
rt <- getPeaklist(MAIT)$rt/60


###################################################
### code chunk number 23: MAITbuilder
###################################################
importMAIT <- MAITbuilder(data = peaks, masses = masses, 
                          rt = rt,significantFeatures = TRUE, 
                          spectraEstimation = TRUE,rtRange=0.2,
                          corThresh=0.7)


###################################################
### code chunk number 24: BiotransformationsBuilder
###################################################
importMAIT <- Biotransformations(MAIT.object = importMAIT, 
                                 adductAnnotation = TRUE, 
                                 peakPrecision = 0.005, adductTable = NULL)


###################################################
### code chunk number 25: identifyMetabolitesBuilder
###################################################
importMAIT <- identifyMetabolites(MAIT.object = importMAIT, 
                                  peakTolerance=0.005,polarity="positive")







