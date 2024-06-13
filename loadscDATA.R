library(Seurat)
library(patchwork)
library(ggplot2)
setwd('D:\\1.paperdata/placentascrna')
projectname <- 'PRJNA672658'
sample_name <- c('GSM5261695','GSM5261696','GSM5261699','GSM5261700')
# 创建一个函数，用于读取数据并创建 Seurat 对象
######1file#####
create_data <- function(sample_id) {
  pbmcdata <- read.delim(paste(projectname,'/',sample_id,'_counts.tsv.gz', sep = ''),
                          row.names=1)
  pbmc <- CreateSeuratObject(counts = pbmcdata , 
                             project = projectname, 
                             min.cells = 3, min.features = 300)
  # 可以添加其他处理步骤
  return(pbmc)
}
######3file#####
create_data <- function(sample_id) {
  file_path <- paste(projectname,'/',sample_id,'/outs/filtered_feature_bc_matrix', sep = '')
  pbmc.data <- Read10X(data.dir = file_path)
  pbmc <- CreateSeuratObject(counts = pbmc.data, 
                             project = projectname
                             , min.cells = 3, min.features = 300)
  # 可以添加其他处理步骤
  return(pbmc)
}
#######looom####
library(loomR)
library(SeuratDisk)
library(tidyverse)
library(data.table)
lfile <- connect(filename = "E:\\GSE156793_S3_gene_count.loom", mode = "r+")
matrix=lfile[["matrix"]][3999040:4028915,] 
matrix=t(matrix)
colnames(matrix)= GSE156793_S1_metadata$sample
rownames(matrix)= Metadata_genes$gene_short_name
symbol <- rownames(matrix)
table(duplicated(symbol))
matrix <- as.data.table(matrix)
matrix$symbol <- rownames(matrix)
matrix <- matrix[, lapply(.SD, sum), by = symbol]
matrix <- column_to_rownames(matrix,'symbol')
GSE156793<- CreateSeuratObject(counts = matrix, 
                           project = 'GSE156793',meta.data = GSE156793_S1_metadata[,-1]
                           , min.cells = 3, min.features = 300)
#####merge####
data_list <- lapply(sample_name, create_data)
datalist <- c(data_list1, data_list2, data_list3,data_list4,data_list5)
samplenames <- unlist(lapply(datalist,function(x) as.character(unique(x$orig.ident))))
names(datalist) <- samplenames
placenta_data = merge(datalist[[1]],y = datalist[-1], add.cell.ids =samplenames)
#######qc######
samplenames <- unlist(lapply(datalist,function(x) as.character(unique(x$orig.ident))))
placenta_data = merge(datalist[[1]],y = datalist[-1], add.cell.ids =samplenames)
HB.genes_total <- c("HBA1","HBA2","HBB") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(placenta_data))
HB.genes <- rownames(placenta_data@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
placenta_data[["percent.HB"]]<- PercentageFeatureSet(placenta_data,features =HB.genes)
placenta_data[["percent.mt"]] <- PercentageFeatureSet(placenta_data, pattern = "^MT-")
VlnPlot(placenta_data,group.by = 'orig.ident',pt.size = 0,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols = pal(25),ncol = 1,combine = T)+
  theme(axis.text.x.bottom =element_text(size = 3))
placenta_data <- subset(placenta_data, subset = nFeature_RNA > 500 & nFeature_RNA < 10000&
                          nCount_RNA < 100000  & percent.mt < 25 & percent.HB<1)
#######basic anlysis####
###placenta_data <- SplitObject(placenta_data,split.by ='orig.ident')
placenta_data <- NormalizeData(placenta_data)
placenta_data <- FindVariableFeatures(placenta_data)
placenta_data <- ScaleData(placenta_data)
placenta_data <- RunPCA(placenta_data)
placenta_data <- FindNeighbors(placenta_data, dims = 1:30, reduction = "pca")
placenta_data <- FindClusters(placenta_data, resolution = 2, cluster.name = "unintegrated_clusters")
placenta_data <- RunUMAP(placenta_data, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(placenta_data, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype"))
#######integrate####
gc()
placenta_data <- IntegrateLayers(
  object = placenta_data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
save.image(file = 'Placenta.RData')
# re-join layers after integration
placenta_data[["RNA"]] <- JoinLayers(placenta_data[["RNA"]])
placenta_data <- FindNeighbors(placenta_data, reduction = "harmony", dims = 1:30)
placenta_data <- RunUMAP(placenta_data, dims = 1:20, reduction = "harmony")
placenta_data <- FindClusters(placenta_data, resolution = 0.4)
DimPlot(placenta_data, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
saveRDS(placenta_data,'placenta_data.Rds')
