library(Seurat)
library(dplyr)
#library(magrittr)
#library(hdf5r)
library(ggplot2)

################################################
QC<-function(data_10X){
  ot<-CreateSeuratObject(count= d, project = "tnf", min.cells = 3, min.features = 200)
  ot[["percent.mt"]] <- PercentageFeatureSet(ot, pattern = "^mt-")
  ot[["time"]]  <- time
  ot <- subset(ot, subset = nFeature_RNA > Rmin & nFeature_RNA < Rmax & percent.mt < mt.per)
  ot <- NormalizeData(ot, normalization.method = "LogNormalize", scale.factor = 10000)
  ot <- FindVariableFeatures(ot, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(ot) # get gene names
  ot<- ScaleData(ot, features = all.genes)
  ot <- RunPCA(ot,npcs = 60, features = VariableFeatures(object = ot))
  p1<-ElbowPlot(ot,ndims = npca)
  ot <- JackStraw(ot,dims = npca, num.replicate = 100) 
  ot <- ScoreJackStraw(ot, dims = 1:npca)
  JackStrawPlot(ot, dims = 1:npca)
  ot <- FindNeighbors(ot, dims = 1:ndim)
  ot <- FindClusters(ot, resolution = 0.4) 
  ot <- RunUMAP(ot, dims = 1:ndim)
  return(ot)
} 

# -- load 10X dataset
setwd("C:\\Users\\admin\\OneDrive - st.shou.edu.cn\\桌面\\sc")
##
Rmin <- 200
Rmax <- 2500
mt.per <- 5
npca <- 30
ndim <-50  
###
d <- Read10X(data.dir = ".\\P12\\P12_1")

time<-"P12"

P12_o<-QC(d)

DimPlot(P12_o, reduction = "umap",label = TRUE)
p1
