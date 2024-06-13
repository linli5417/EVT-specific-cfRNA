library(dplyr)
library(Seurat)
library(cols4all)
library(ggplot2)
setwd('D:\\1.paperdata/placentascrna')
set.seed(4180)
######颜色#######
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
###################
EVT_cells <- placenta_data %>%
  subset( celltype%in% "EVT")
EVT_cells <- EVT_cells[,sample(1:ncol(EVT_cells), 5000) ]
Other_cells <- placenta_data %>%
  subset( celltype!="EVT")
Other_cells <- Other_cells[, sample(1:ncol(Other_cells), 5000)]
final_sample <- merge(EVT_cells, Other_cells)
final_sample <- JoinLayers(final_sample)
gc()
metadata <- final_sample@meta.data
metadata$cellinfo <- ifelse(metadata$celltype %in% 'EVT','EVT','Other cells')
heatdt <- data.frame(final_sample@assays$RNA@layers$data)
colnames(heatdt) <- colnames(final_sample)
row.names(heatdt) <- rownames(final_sample@assays$RNA@features) 
EVT_marker <- subset(placenta.markers,cluster %in% 'EVT')
heatdt <- heatdt[EVT_marker$gene,]
library(ComplexHeatmap)
ha = HeatmapAnnotation(celltype= metadata$cellinfo)
show_rows <- c("HLA-G", "MMP2",'SOX4','KRT7')
index <- which(rownames(heatdt) %in% show_rows)
labs <- rownames(heatdt)[index]
lab2 = rowAnnotation(foo = anno_mark(at = index,
                                     side='left',
                                     labels = labs,
                                     labels_gp = gpar(fontsize = 8),
                                     lines_gp = gpar(linejoin='round')))
Heatmap(heatdt, top_annotation = ha,cluster_columns = F,show_row_dend = F,
        show_row_names = F,cluster_rows = T,row_names_gp = gpar(fontsize = 8),left_annotation = lab2 ,
        show_column_names = F,use_raster=TRUE,col =c('gray90','firebrick'))
