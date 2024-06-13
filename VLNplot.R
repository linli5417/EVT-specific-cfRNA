library(cols4all)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggunchained) 
library(reshape2)  
set.seed(4180)
pal <- colorRampPalette(c4a('20',15))
########载入数据##########
##meta
dt <- PE_placenta@meta.data
dt <- dt %>% pivot_longer(c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
                    names_to = "gene",
                    values_to = "value")
dt$gene <- factor(dt$gene, levels = unique(dt$gene))
##gene
inter_gene <- intersect(subset(df_sig,module %in% c('EVT-M1','EVT-M3','EVT-M7','EVT-M9')
                         )$gene_name,cfgenes$x)
metadata <- subset(PE_placenta,celltype %in% 'EndoC')@meta.data
dt <- data.frame(subset(PE_placenta,celltype %in% 'EndoC')@assays$RNA@layers$data)
row.names(dt) <- row.names(subset(PE_placenta))
colnames(dt) <- row.names(metadata)
vlndt <- dt[c('TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D'), ]
vlndt <- data.frame(t(vlndt))
vlndt$group <- metadata$group
vlndt$celltype <- metadata$celltype
vlndt <- vlndt %>% pivot_longer(c('TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D'),
                          names_to = "gene",
                          values_to = "value")
vlndt$gene <- factor(vlndt$gene ,levels =c('TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D') )
###############绘图#######################
ggplot(data = vlndt,
       aes(x = value, y = celltype, fill = celltype)) +
       geom_violin(scale = 'width',
              draw_quantiles= c( 0.5),
              color= 'black',
              size= 0.45, #描边粗细
              alpha= 0.8) +
  facet_grid(cols = vars(gene), scales = 'free_x')+
  scale_fill_manual(values = pal(12)) + #填充色修改
  theme_bw()+
  theme(panel.grid = element_blank(), #移除背景网格线
    axis.text.x = element_text(size = 10), #x轴标签大小调整
    axis.text.y = element_text(size = 10), #y轴标签大小调整
    axis.title.x = element_text(size = 12), #x轴标题大小调整
    axis.title.y = element_blank(), #移除y轴标题
    strip.background = element_blank(), #移除分面外围边框
    strip.text.x = element_text(size = 10, angle = 0), #分面文本标签倾斜60°
    legend.title = element_text(size = 16), #图例标题大小调整
    legend.text = element_text(size = 15),#图例标签大小调整
    legend.position="none"
    ) + labs(x = 'Expression') #x轴标题本文内容修改

##########分列############

ggplot(data = vlndt,
       aes(x = gene, y = value, fill = group)) +
  geom_split_violin(scale = 'width',
              color= 'black',
              size= 0.45, #描边粗细
              alpha= 0.8) +
  scale_fill_manual(values = pal(3)) + #填充色修改
  theme_bw()+
  theme( #移除背景网格线
        axis.text.x = element_text(size = 10), #x轴标签大小调整
        axis.text.y = element_text(size = 10), #y轴标签大小调整
        axis.title.x = element_text(size = 16), #x轴标题大小调整
        strip.background = element_blank(), #移除分面外围边框
        strip.text.x = element_text(size = 16, angle = 0), #分面文本标签倾斜60°
        legend.title = element_blank(), #图例标题大小调整
        legend.text = element_text(size = 15),
        legend.background = element_blank(),#图例标签大小调整
        legend.position=c(.15,0.85)
  ) + labs(title='TNFSF10 receptors in endothelial cells',x='',y='Expression')+
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(width = 0.9),show.legend = F) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 2)), vjust = -2, position = position_dodge(width = 0.9))+
  stat_compare_means(method = ".test", label = "p.format")
##########box################
ggplot(vlndt,aes(gene,value))+
  geom_boxplot(aes(fill=group))+
  scale_fill_manual(values=pal(3))+
  theme(legend.position = "top")+
  ggplot2::labs(x=NULL,y="Expression")+
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position = "top",
                 axis.text.x = ggplot2::element_text(angle = 45,hjust = 1))+
  theme(axis.text.x = element_text(angle = 45,size =12,hjust = 1))+
  scale_color_manual(values =pal(2))+
  labs(title='',x='',y='expression')
