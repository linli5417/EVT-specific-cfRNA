library(Seurat)
library(cols4all)
library(ggplot2)
library(ggunchained) 
library(reshape2) 
setwd('D:\\1.paperdata/placentascrna/5.PE')
set.seed(4180)
######颜色#######
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
#####UMAP########
DimPlot(PE_placenta, reduction = "umap",label = F,
        group.by =c("celltype"),
        cols = pal(13))+
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
        legend.position="right")+
  guides(color = guide_legend(ncol =  1,override.aes = list(size = 4)))+ ##图例改成一行
  labs(title = "",color= 'Celltype')
#########dot###################
markers.to.plot <- unique(marker$Cell.marker)
markers.to.plot <- split(marker$Cell.marker, marker$Cell.name) 
###绘图
DotPlot(PE_placenta,dot.scale = 4,group.by = 'celltype',
        features = c(markers.to.plot,'S100A8','S100A9','CEACAM8'),
        cols = c('white',pal(1))
) + RotatedAxis() + # 来自Seurat
  theme(legend.position = 'top',
        legend.title = element_text(size = 12),
        strip.text =  element_text(),#修改分页标签
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm"),
        axis.title = element_blank()
  )+
  guides(colour = guide_colourbar(title.vjust = 0.1, title.hjust = 0,label.position = "top"),
         size = guide_legend(title.vjust = 0.1, title.hjust = 0,label.position = "top"))+
  labs(size = "Percent Expressed", color = "Average Expression")
########修改名字############
PE_placenta@active.ident <- PE_placenta$seurat_clusters
new.cluster.ids <- c("CTB", "HB", "CTB", "Granu", "CTB", "EVT",
                     "Granu", "CTB",'T cell','NK','Mac',
                     'Granu','STB','CTB','STB','Mix','Mac','B cell',
                     'DSC','Mix','CTB','EndoC')
names(new.cluster.ids) <- levels(PE_placenta)
PE_placenta <- RenameIdents(PE_placenta,new.cluster.ids)
PE_placenta$celltype <- PE_placenta @active.ident
PE_placenta$celltype <- factor(PE_placenta$celltype,
                                 levels = c('EVT','CTB','STB',
                                            'DSC','EndoC',
                                            'T cell','NK','B cell',
                                            'Mac','HB','Granu','Mix'))
###########feature###############
FeaturePlot(PE_placenta,features = c('HBA'),reduction = "umap",
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
############VLN##################
VlnPlot(PE_placenta,features = c('TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D'),
        cols =c(pal(3)),group.by = 'celltype',
        pt.size = 0,split.by = 'group',split.plot = T)
#############findmarker##############
PE_marker<- FindMarkers(subset(PE_placenta,celltype%in%'EVT'),ident.1='PE',only.pos = F, min.pct = 0.25, 
                        logfc.threshold = 1,group.by = "group")
