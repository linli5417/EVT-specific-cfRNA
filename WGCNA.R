library(dplyr)
library(GEOquery)
library(stringr)
library(WGCNA)
gset <- getGEO('GSE154377', destdir=".",
               AnnotGPL = T,     ## 注释文件
               getGPL = F)       ## 平台文件
setwd("D:/1.paperdata/cfrna/wgcna")
####testdata#####
file_list <- list.files(path = "GSE154377/",pattern = "Counts.txt.gz", full.names = TRUE)
dattest <- read.table(file_list[1], header = FALSE, skip = 4)
colnames(dattest) <- c('Geneid', substr(basename(file_list[1]), 1, 10))
for (file in file_list[-1]) {
  data <- read.table(file, header = FALSE, skip = 4)
  colnames(data) <- c('Geneid', substr(basename(file), 1, 10))
  # 使用 dplyr 的 bind_cols 进行横向合并
  dattest <- left_join(dattest, data, by = 'Geneid')
}
rownames(dattest) <- dattest[,1]
dattest <- dattest[,-1]
group_data <- gset[[1]]@phenoData@data
group_data <-  group_data %>% 
  mutate(group = str_extract(title, "^[^\\s]+"))
group_data <- group_data[,c(1,2,4)]
group_data <- subset(group_data,group %in% c('Normal','Preeclampsia/gestational') )
dattest <- dattest[,group_data$geo_accession]
########wgcnadata######
genedata <- expdata[,1:2]
expdata <- expdata[,-(1:2)]
colnames(expdata) <-gsub( 'X','',colnames(expdata))
row.names(expdata) <-genedata$gene_num
group_data <- gset[[1]]@phenoData@data
row.names(group_data) <- group_data[,1]
group_data <- subset(group_data,title %in%colnames(expdata))
expdata <- expdata[,group_data$title]
sampledata<- colnames(expdata)
########计算fpkm######
library(biomaRt)
ensembl <- useMart("ensembl") 
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
test <- getBM(attributes=c('ensembl_gene_id', 'start_position',
                           'end_position','ensembl_transcript_id',
                           'transcript_length','hgnc_symbol'),mart = ensembl)
test <- test[order(test$ensembl_gene_id,test$transcript_length,decreasing = T),]
g <- test[!duplicated(test$ensembl_gene_id),]
g <- g[,c(1,5)]
ng=intersect(rownames(dattest),g$ensembl_gene_id) 
lengths=g[match(ng,g$ensembl_gene_id),2]
names(lengths) <- g[match(ng,g$ensembl_gene_id),1]
dattest <- dattest[names(lengths),]
total_count <- colSums(dattest )
datExpr0  <- t(do.call( rbind,
                      lapply(1:length(dattest),
                             function(i){
                               10^9*dattest[,i]/lengths/total_count[i]
                               #lengths向量自动遍历
                             }) ))
colnames(datExpr0 ) <- group_data$geo_accession
datExpr0 = t(log2(datExpr0+1))
##############WGCNA预处理##########################
datExpr0 <- t(exprdt)
colnames(group_data) <- gsub('sampling time group:ch1','time',colnames(group_data))
datExpr0<-datExpr0[subset(group_data,time %in% c('≤12w','13-20w'))$title,]
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",hang=-1, cex.lab =0.1,
     cex.axis = 1.5, cex.main = 2,labels=FALSE)
abline(h = 210, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 210, minSize = 10)
table(clust)
keepSamples = (clust==1)  #保留非离群(clust==1)的样本
datExpr = datExpr0[keepSamples, ]  #去除离群值后的数据
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
########################sft##########
powers = c(c(1:20), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="black");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="black")
###############net######
sft$powerEstimate
cor<-WGCNA::cor
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "cfTOM",
                       verbose = 3)    
##############画图######
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)   
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
###########表型数据##############
allTraits<- group_data[,c(48,49)]
colnames(allTraits)<- c('disease','ges_week') 
allTraits<- allTraits %>%
              mutate(
  control = as.numeric(case_when(disease == "control" ~ 1, TRUE ~ 0)),
  pre_eclampsia = as.numeric(case_when(disease == "pre-eclampsia" ~ 1, TRUE ~ 0)),
  severe_PE = as.numeric(case_when(disease == "severe pre-eclampsia" ~ 1, TRUE ~ 0)))
allTraits<- allTraits %>%
  mutate(
    weeks_before_12= as.numeric(case_when(ges_week == "≤12w" ~ 1, TRUE ~ 0)),
    weeks_13_20 = as.numeric(case_when(ges_week == "13-20w" ~ 1, TRUE ~ 0)),
    weeks_after_20 = as.numeric(case_when(ges_week == "≥23w" ~ 1, TRUE ~ 0)),
    Post_partum= as.numeric(case_when(ges_week == "Post-partum" ~ 1, TRUE ~ 0))
  )
allTraits <- allTraits[,-(1:2)]
allTraits <- subset(allTraits,Post_partum != 1 )
allTraits <- subset(allTraits,weeks_after_20 != 1 )
allTraits <- allTraits[,-(4:7)]
##################相关性#####################
library(corrplot)
library(ggplot2)
library(ggpubr)
###########模块-表型数据关联####################
# 重新计算带有颜色标签的模块
datTraits <- allTraits[rownames(datExpr),]
datTraits <- data.frame(datTraits,t(gsva_mat))
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
library(ggplot2)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
                       moduleColor = moduleColors)
