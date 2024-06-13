options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
go_list <- list(  EVT_PE_score =colnames(dtrain[,-1])
      )
gsva_mat <- gsva(expr=as.matrix(exprdt), 
                 gset.idx.list=go_list, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核
#########testdata#################

