library(Seurat)
library(cols4all)
library(ggplot2)
setwd('D:\\1.paperdata/placentascrna')
set.seed(4180)
######颜色#######
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
#####UMAP########
DimPlot(placenta_data, reduction = "umap",label = T,
        group.by = c("celltype"),
        cols = pal(14))+
  tidydr::theme_dr(xlength = 0.4, 
                   ylength = 0.4,
                   arrow = arrow(length = unit(3, "mm"),type = "closed"))+
       theme(panel.grid = element_blank(), #移除背景网格线
             plot.title = element_text(hjust = 0.5,size = 16),
         axis.text.x = element_blank(), #x轴标签大小调整
         axis.text.y = element_blank(), #y轴标签大小调整
         axis.ticks = element_blank(),  #移除刻度
         axis.title.x = element_text(size = 14), #x轴标题大小调整
         axis.title.y = element_text(size = 14), #移除y轴标题
         legend.title = element_text(size = 12), #图例标题大小调整
         legend.text = element_text(size = 10),#图例标签大小调整
         legend.position="none")+
         guides(color = guide_legend(ncol =  1,override.aes = list(size = 4)))+ ##图例改成一行
        labs(title = "Cell Type",color= 'Sample Sources')

##
DimPlot(placenta_data, reduction = "umap", 
        group.by = c("seurat_clusters"),
        cols = pal(30),label = T)
#####DOTplot#####
markers.to.plot <- unique(marker$Cell.marker)
markers.to.plot <- split(marker$Cell.marker, marker$Cell.name) 
###绘图
DotPlot(placenta_data, features =markers.to.plot , 
        cols = pal(2), dot.scale = 3) +
        RotatedAxis()
##分页dot
DotPlot(placenta_data,dot.scale = 4,group.by = 'celltype',
        features = c(markers.to.plot),
        cols = c('white',pal(1))
) + RotatedAxis() + # 来自Seurat
  theme(legend.position = 'top',
        legend.title = element_text(size = 12),
    strip.text = element_blank(),#修改分页标签
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank()
  )+
  guides(colour = guide_colourbar(title.vjust = 0.1, title.hjust = 0,label.position = "top"),
         size = guide_legend(title.vjust = 0.1, title.hjust = 0,label.position = "top"))+
  labs(size = "Percent Expressed", color = "Average Expression")
########修改名字############
placenta_data@active.ident <- placenta_data$seurat_clusters
new.cluster.ids <- c("EVT", "DSC", "CTB", "EVT", "Mac", "FB",
                     "NK", "T cell",'HB','EVT','STB',
                     'EndoC','Mix','CTB','VSMC','EpiC','B cell','EVT',
                     'DSC','NK','T cell','HB','DSC','FB','EndoC')
names(new.cluster.ids) <- levels(placenta_data)
placenta_data <- RenameIdents(placenta_data,new.cluster.ids)
placenta_data$celltype <- placenta_data @active.ident
placenta_data$celltype <- factor(placenta_data$celltype,
                                 levels = c('EVT','CTB','STB',
                                            'DSC','VSMC','FB','EndoC','EpiC',
                                            'T cell','NK','B cell',
                                            'Mac','HB','Mix'))
#########feature Plot##################
FeaturePlot(placenta_data,features = c('HLA-G'),reduction = "umap",
            cols =c('white',pal(1)) )+
  tidydr::theme_dr(xlength = 0.4, 
                   ylength = 0.4,
                   arrow = arrow(length = unit(3, "mm"),type = "closed"))+
  theme(panel.grid = element_blank(), #移除背景网格线
        axis.text.x = element_blank(), #x轴标签大小调整
        axis.text.y = element_blank(), #y轴标签大小调整
        axis.ticks = element_blank(),  #移除刻度
        plot.title = element_text(hjust = 0.5,size = 16),
        axis.title.x = element_text(size = 14), #x轴标题大小调整
        axis.title.y = element_text(size = 14), #移除y轴标题
        legend.title = element_text(size = 12), #图例标题大小调整
        legend.text = element_text(size = 10),#图例标签大小调整
        legend.position="right")+labs(color= 'Expression')

########找marker#####
placenta.markers <- FindAllMarkers(placenta_data, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.25)
write.csv(placenta.markers,'placenta_maerkers.csv')
