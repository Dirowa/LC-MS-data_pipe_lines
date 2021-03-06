############################
# installingof metaboshiny
############################


# RUN ON TERMINAL#
## run on terminal sudo apt install -y libudunits2-0 libudunits2-dev
 sudo apt install libgdal-dev
 sudo apt-get install libglu1-mesa-dev
 sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev
sudo apt-get install -y librsvg2-dev
sudo apt-get install -y default-jre
sudo apt-get install -y default-jdk
sudo R CMD javareconf

sudo add-apt-repository ppa:marutter/c2d4u3.5
sudo apt-get update
sudo apt-get install r-cran-rjava
sudo apt-get install libbz2-dev libpcre3-dev
sudo apt-get install liblzma-dev zlib1g-dev
sudo apt install libomp-dev -y
R CMD javareconf -e
sudo apt install openjdk-8-jre-headless
apt-get install libsqlite-dev
apt-get install libmariadbclient-dev
apt-get install libmariadbd-dev -f
 apt-get install libmariadb-client-lgpl-dev
sudo add-apt-repository ppa:webupd8team/java
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+
 sudo apt install gcc
sudo apt-get install libpcre2-dev - y
apt-get install -y libgit2-dev
apt install -y lsb-release
sudo apt-get install g++
sudo apt-get install libudunits2-dev
apt install gfortran gfortran-doc
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get install libgeos-dev
sudo apt-get install libproj-dev
sudo apt-get install gdal-bin 
sudo apt-get install proj-bin 
sudo apt-get install libignition-math4-dev -y
 sudo apt-get install libmysqlclient-dev -f
sudo apt-get install default-libmysqlclient-dev
sudo apt-get install default-libmysqlclient-dev
sudo apt-get install libgdal-dev -f




install.packages("devtools")
install.packages("BiocManager")
install.packages("rJava")


install.packages('orca')
BiocManager::install('orca')
install.packages("tidyverse", repo = 'https://mac.R-project.org')

BiocManager::install("org.Hs.eg.db")
install.packages("pacman")
BiocManager::install(' mutoss')
BiocManager::install(' metap')
BiocManager::install(' units')
BiocManager::install(' sf')
BiocManager::install('globaltest')
BiocManager::install('Rgraphviz')
BiocManager::install('SSPA')
BiocManager::install('KEGGgraph')
BiocManager::install('siggenes')
BiocManager::install('GlobalAncova')
BiocManager::install("mzR")
BiocManager::install("ctc")
BiocManager::install('fgsea')
BiocManager::install('glasso')
BiocManager::install('huge')
BiocManager::install('RBGL')
BiocManager::install('crmn')

BiocManager::install(gdata)
install.packages("gdata")
install.packages("glasso")
install.packages("huge")
install.packages("ppcor")
install.packages("plotly")
install.packages('rgl')

devtools::install_github("gaospecial/ggVennDiagram")

library(devtools)
BiocManager::install('rsvg')
install.packages('librsvg2-dev')
BiocManager::install("ChemmineR")
BiocManager::install("ChemmineR")
devtools::install_github("tgirke/ChemmineR")

BiocManager::install('MetaDBparse')
install.packages('MetaDBparse')
devtools::install_github("joannawolthuis/MetaDBparse")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") 
devtools::install_github("rwehrens/ChemometricsWithR")
devtools::install_github("rwehrens/BatchCorrMetabolomics")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")


devtools::install_github("xia-lab/MetaboAnalystR")



 #‘GlobalAncova’, ‘Rgraphviz’, ‘SSPA’, ‘KEGGgraph’, ‘siggenes’



install.packages("pacman")
pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea"))

devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
library(MetaboAnalystR)

devtools::install_github("xia-lab/MetaboAnalystR")


 BiocManager::install('multtest')
 BiocManager::install('MSnbase')
 BiocManager::install('siggenes')
 BiocManager::install('KEGGgraph')


 BiocManager::install('SSPA')
 BiocManager::install('preprocessCore')
 BiocManager::install('Rgraphviz')
 BiocManager::install('GlobalAncova')


 BiocManager::install('globaltest')
 BiocManager::install('pcaMethods')


#######################
# library's###########
######################

library(pacman)
# METADBPARSE REQUIREMENTS
pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table, 
               MetaDBparse, RSQLite, DBI, gsubfn, utils, RCurl, XML, base, 
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2, 
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist)

# METABOSHINY REQUIREMENTS
pacman::p_load(ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices, 
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools, 
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles, 
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel, 
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown, 
               enviPat, ROCR, tsne, e1071,
               pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker, 
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr, 
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,
               stats, pROC, car, doParallel, missForest)
#devtools::install_github("UMCUGenetics/MetaboShiny", ref="dev")
library(MetaboShiny)

devtools::install_github("joannawolthuis/MetaboShiny", force = TRUE)

###### dev branches#####
##devtools::install_github("joannawolthuis/MetaboShiny", ref = 'dev')
#devtools::install_github("joannawolthuis/MetaboShiny", ref = 'v1')
#devtools::install_github("joannawolthuis/MetaboShiny", ref = 'v2')


MetaboShiny::start_metshi(inBrowser=T)
