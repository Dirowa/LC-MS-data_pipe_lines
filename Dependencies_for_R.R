#################################
# R Dependencies for the Wrapper#
#################################

install.package('devtools')
nstall.packages("BiocManager")
BiocManager::install(c("xcms","RColorBrewer","pander","pheatmap","SummarizedExperiment","magrittr","MSnbase"))
install.packages('doParallel','stats')
BiocManager::install("IPO")
BiocManager::install("CAMERA")
BiocManager::install("ProteoMM")
remotes::install_github("ricoderks/Rcpm")

#metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T, force = TRUE)
install('dplyr')
BiocManager::install("ropls")
BiocManager::install('MultiDataSet')
devtools::install_github("SciDoPhenIA/phenomis")
BiocManager::install("biosigner")

