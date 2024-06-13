# single-cell analysis package
setwd('D:\\1.paperdata/placentascrna/4.PE')
.libPaths(c("~/SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(GlobalOptions)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(4180)
##############预处理############
sampled_cell<-PE_placenta
sampled_cell[["RNA"]] <- as(object =PE_placenta[["RNA"]],Class = 'Assay')
##sampled_cell[["RNA"]] <- as(object = sampled_cell[["RNA"]], Class = "Assay")
#num_cells <- ncol(sampled_cell)
#sample_indices <- sample(1:num_cells, 50000, replace = FALSE)
#sampled_cell<- sampled_cell[, sample_indices]
sampled_cell <- SetupForWGCNA(
  sampled_cell,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "placenta" # the name of the hdWGCNA experiment
)
length(sampled_cell@misc$placenta$wgcna_genes)
sampled_cell  <- MetacellsByGroups(
  seurat_obj = sampled_cell ,
  group.by = c("celltype", "orig.ident"), # specify the columns in sampled_cell@meta.data to group by
  k = 20, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)
######颜色#######
library(cols4all)
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
########
sampled_cell <- NormalizeMetacells(sampled_cell)
sampled_cell <- ScaleMetacells(sampled_cell, features=VariableFeatures(sampled_cell))
sampled_cell <- RunPCAMetacells(sampled_cell, features=VariableFeatures(sampled_cell))
sampled_cell <- RunHarmonyMetacells(sampled_cell, group.by.vars='orig.ident')
sampled_cell <- RunUMAPMetacells(sampled_cell, reduction='harmony', dims=1:15)
p1 <- DimPlotMetacells(sampled_cell, group.by='celltype') + umap_theme() + ggtitle("Cell Type")
p2 <- DimPlotMetacells(sampled_cell, group.by='orig.ident') + umap_theme() + ggtitle("Sample")
p1 | p2
#########sft##############
sampled_cell <- SetDatExpr(
  sampled_cell,
  group_name = c('EVT'),
  group.by='celltype'
)
# Test different soft powers:
sampled_cell <- TestSoftPowers(
  sampled_cell,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(sampled_cell)
wrap_plots(plot_list, ncol=2)
#########多线程开启hdwgcna#####
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
options(mc.cores = 4) # 设置线程数
sampled_cell <- ConstructNetwork(
  sampled_cell, soft_power=6,
  setDatExpr=FALSE,nThreads=4,parallelize = TRUE,
  tom_name = 'Placenta-M' # name of the topoligical overlap matrix written to disk
)
stopCluster(cl)
PlotDendrogram(sampled_cell, main='Placenta hdWGCNA Dendrogram')
# need to run ScaleData first or else harmony throws an error:
sampled_cell <- ScaleData(sampled_cell, features=VariableFeatures(sampled_cell))
# compute all MEs in the full single-cell dataset
sampled_cell <- ModuleEigengenes(
  sampled_cell
)
# harmonized module eigengenes:
hMEs <- GetMEs(sampled_cell)
# module eigengenes:
MEs <- GetMEs(sampled_cell, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
sampled_cell <- ModuleConnectivity(
  sampled_cell,
  group.by = 'celltype', group_name = 'EVT'
)
# rename the modules
sampled_cell <- ResetModuleNames(
  sampled_cell,
  new_name = "EVT-M"
)
modules <- GetModules(sampled_cell)
hub_df <- GetHubGenes(sampled_cell, n_hubs = 50)
# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  sampled_cell,raster = TRUE,raster_dpi = 100,
  features='hMEs', # plot the hMEs
  order=FALSE # order so the points with highest hMEs are on top
)
wrap_plots(plot_list, ncol=3)
# with UCell method
library(UCell)
sampled_cell <- ModuleExprScore(
  sampled_cell,
  n_genes = 50,
  method='UCell'
)
VlnPlot(
  sampled_cell,
  features = 'Placenta-M4',
  group.by = 'celltype',
  pt.size = 0 # don't show actual data points
)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
############vision################
group1 <- sampled_cell@meta.data %>% subset(celltype %in% 'EVT'&group%in% 'Normal') %>% rownames
group2 <- sampled_cell@meta.data %>% subset(celltype %in% 'EVT'& group%in% 'PE') %>% rownames
DMEs <- FindDMEs(
  sampled_cell,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='placenta',
  
)
library(ggrepel)
PlotDMEsVolcano(
  sampled_cell,
  DMEs,
  wgcna_name = 'placenta'
)+labs(title='PE vs Normal')+xlim(-5,5)
PlotDMEsLollipop(
  sampled_cell, 
  DMEs, 
  wgcna_name='placenta', 
  pvalue = "p_val_adj"
)+labs(title='PE vS Normal')+
theme(axis.title = element_text(size = 10))+theme(legend.position = 'right')
# network analysis & visualization package:#########
sampled_cell <- RunModuleUMAP(
  sampled_cell,
  genes_use = hub_df$gene_name,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
library(igraph)
ModuleUMAPPlot(
  sampled_cell,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
HubGeneNetworkPlot(
  sampled_cell,
  n_hubs = 3, n_other=10,
  edge_prop = 0.75,
  mods = 'all'
)+theme(axis.text = geom_text_repel(size=3))
