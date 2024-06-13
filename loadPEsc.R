library(Seurat)
library(patchwork)
library(ggplot2)
library(cols4all)
set.seed(4180)
######颜色#######
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
setwd('D:\\1.paperdata/placentascrna')
projectname <- 'GSE173193'
sample_name <- c('GSM5261695','GSM5261696','GSM5261699','GSM5261700')
# 创建一个函数，用于读取数据并创建 Seurat 对象
create_data <- function(sample_id) {
  file_path <- paste('1.rawdata/',projectname,'/',sample_id, sep = '')
  pbmc.data <- Read10X(data.dir = file_path)
  pbmc <- CreateSeuratObject(counts = pbmc.data, 
                             project = projectname
                             , min.cells = 3, min.features = 300)
  pbmc$orig.ident <- sample_id
  # 可以添加其他处理步骤
  return(pbmc)
}
data_list <- lapply(sample_name, create_data)
data_list[[4]]$group <- 'PE'
PE_placenta = merge(data_list[[1]],y = data_list[-1], add.cell.ids = sample_name )
HB.genes_total <- c("HBA1","HBA2","HBB") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(PE_placenta))
HB.genes <- rownames(PE_placenta@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
PE_placenta[["percent.HB"]]<- PercentageFeatureSet(PE_placenta,features =HB.genes)
PE_placenta[["percent.mt"]] <- PercentageFeatureSet(PE_placenta, pattern = "^MT-")
VlnPlot(PE_placenta,group.by = 'orig.ident',pt.size = 0,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols = pal(25),ncol = 1,combine = T)+
  theme(axis.text.x.bottom =element_text(size = 3))
PE_placenta <- subset(PE_placenta, subset = nFeature_RNA > 500 & nFeature_RNA < 10000&
                          nCount_RNA < 100000  & percent.mt < 25 & percent.HB<1)
#######basic anlysis####
###PE_placenta <- SplitObject(PE_placenta,split.by ='orig.ident')
PE_placenta <- NormalizeData(PE_placenta)
PE_placenta <- FindVariableFeatures(PE_placenta)
PE_placenta <- ScaleData(PE_placenta)
PE_placenta <- RunPCA(PE_placenta)
PE_placenta <- FindNeighbors(PE_placenta, dims = 1:30, reduction = "pca")
PE_placenta <- FindClusters(PE_placenta, resolution = 2, cluster.name = "unintegrated_clusters")
PE_placenta <- RunUMAP(PE_placenta, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
#######integrate####
gc()
PE_placenta<- IntegrateLayers(
  object = PE_placenta, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
save.image(file = 'PE_placenta.RData')
# re-join layers after integration
PE_placenta[["RNA"]] <- JoinLayers(PE_placenta[["RNA"]])
PE_placenta <- FindNeighbors(PE_placenta, reduction = "harmony", dims = 1:30)
PE_placenta <- RunUMAP(PE_placenta, dims = 1:20, reduction = "harmony")
PE_placenta <- FindClusters(PE_placenta, resolution = 0.4)
DimPlot(PE_placenta, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
saveRDS(PE_placenta,'PE_placenta.Rds')
