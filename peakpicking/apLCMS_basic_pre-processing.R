###########################
# apLCMS pre-processing ##
#############################
# required installing R tools: https://cran.rstudio.com/bin/windows/Rtools/
  # after installen
    shell("R CMD build apLCMS")
    
    # restart R
    Sys.which("make")
    ## "C:\\rtools40\\usr\\bin\\make.exe"
    
    install.packages("jsonlite", type = "source")
    
##########################
# installing software    #
##########################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("devtools",'MASS', 'rgl', 'mzR' , 'splines', 'doParallel', 'foreach', 'iterators', 'snow','ROCR', 'ROCS','e1071','randomForest','gbm'))
BiocManager::install('Rtools')

#######################
# library's           #
#######################
library(c('apLCMS','MASS', 'rgl', 'mzR' , 'splines', 'doParallel', 'foreach', 'iterators', 'snow','ROCR', 'ROCS','e1071','randomForest','gbm'))
library('apLCMS')


######################
# variables          #
######################
#path to apLCMS zipfile
path_apLCMSzipfolder <- "F:/avans/stage MM/apLCMS/"
path_to_package <- "F:/avans/stage MM/apLCMS/apLCMS"

# path to data
path <- "F:/avans/stage MM/test_data_sherlok/"

#path to safe data
safe <-"F:/avans/stage MM/apLCMS/output.txt"


####################
# installing apLCMS#
####################
library('devtools')


setwd(path_apLCMSzipfolder)
shell("R CMD build apLCMS")
install.packages(path_apLCMSzip, repos = NULL)
install.packages(list.files(pattern="apLCMS*.gz"), repos = NULL)




#################################
# importing data & peakpicking  #
#################################

setwd(path)
cdf.files<-dir(pattern=".mzML")
cdf.files

registerDoParallel(10)

############## cant use this due bug #######################################
#Error in { : task 2 failed - cannot allocate vector of size 8 Kb
# tried to solve with limit.

# the loop
#require(foreach)
#features <- foreach(i=1:length(cdf.files),.packages = 'apLCMS'  ) %dopar%{
#  this.spectrum<-proc.cdf(filename=cdf.files[i], min.pres = 0.5, min.run = 10, tol = 1e-5)
#  this.peak.table<-prof.to.features(this.spectrum, bandwidth = 0.5, min.bw=10, max.bw=30)
#  feature <- this.peak.table
#  return(feature)
#  }
#############################################################################

features<-new("list")
for(i in 1:length(cdf.files))
{
lcms_data<-proc.cdf(filename=cdf.files[i],
                        min.pres=0.5,
                        min.run=12,
                        tol=1e-5,
                        baseline.correct=1e-5,
                        do.plot=FALSE,
                        intensity.weighted=FALSE
                        )

feature<-prof.to.features(lcms_data,
                                  bandwidth=0.5,
                                  min.bw=10,
                                  max.bw=30,
                                  sd.cut=c(1,60),
                                  sigma.ratio.lim=c(0.2, 5),
                                  shape.model="bi-Gaussian",
                                  estim.method="moment",
                                  do.plot=TRUE, power=2,
                                  component.eliminate=0.01,
                                  BIC.factor=2
                                  )

features[[i]]<-feature
 }

summary(features)


############################
# retention time correction#
############################
features2<-adjust.time(features,
                              mz.tol = NA,
                              chr.tol = NA,
                              colors=NA,
                              find.tol.max.d=1e-4,
                              max.align.mz.diff=0.01
                              )


############################
# allignment              #
###########################

aligned<-feature.align(features2,
                       min.exp = 2,
                       mz.tol = NA,
                       chr.tol = NA,
                       find.tol.max.d=1e-4,
                       max.align.mz.diff=0.01
                       )

head(aligned$aligned.ftrs)
############################
#filling data           #
##########################

post.recover<-aligned
for(i in 1:length(cdf.files))
  {
     this.recovered<-recover.weaker(filename=cdf.files[i], loc=i,
                                     aligned.ftrs=aligned$aligned.ftrs, pk.times=aligned$pk.times, align.mz.tol=aligned$mz.tol,
                                     align.chr.tol=aligned$chr.tol, this.f1=features[[i]], this.f2=features2[[i]], orig.tol=1e-5)
     post.recover$aligned.ftrs[,i+4]<-this.recovered$this.ftrs
     post.recover$pk.times[,i+4]<-this.recovered$this.times
     post.recover$features[[i]]<-this.recovered$this.f1
     post.recover$f2[[i]]<-this.recovered$this.f2
     }
summary(post.recover)

# MZ summary with intensity for each sample
head(post.recover$aligned.ftrs)
head(post.recover$f2)

mz_table <- post.recover$aligned.ftrs



write.csv(mz_table, safe)
