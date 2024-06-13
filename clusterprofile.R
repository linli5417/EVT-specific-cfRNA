library(clusterProfiler)
library(ggplot2)
setwd('D:\\1.paperdata/placentascrna/4.PE')
df_sig <- subset(hubdf,module %in% c('EVT-M1','EVT-M3','EVT-M9','EVT-M7'))
group <- data.frame(gene=df_sig$gene_name,
                    group=df_sig$module)
Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)
dotplot(data_GO_sim, showCategory=5,font.size = 10)
df_GO <- data_GO@compareClusterResult
library(forcats)
df_GO$Description <- as.factor(df_GO$Description)
df_GO$Description <- fct_inorder(df_GO$Description)
ggplot(df_GO, aes(Cluster, Description)) +
  geom_point(aes(fill=p.adjust, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="purple",high="yellow")+
  labs(x=NULL,y=NULL)+
  coord_flip()
