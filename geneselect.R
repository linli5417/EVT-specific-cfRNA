library(dplyr)
library(tidyverse)
setwd('D:\\1.paperdata/placentascrna/4.ML')
set.seed(4180)
########数据准备#########
Cell_marker <- read_excel("Cell_marker_Human.xlsx")
normal_placenta <- GSE75010[,subset(GSE75010_class,diagnosis.ch1 %in% 'non-PE')$geo_accession]
normal_placenta <- data.frame(scale(normal_placenta))
tissue_data <- data.frame(scale(t(gtex_mrna_pheno[,-(1:2)])))
colnames(tissue_data) <- gtex_mrna_pheno$sample_id
tissue_data$gene <- rownames(tissue_data)
normal_placenta$gene <- rownames(normal_placenta)
common_genes <- dplyr::intersect(tissue_data$gene,normal_placenta$gene)
###############markerselect####################
dat1 <- merge(normal_placenta[normal_placenta$gene %in% common_genes,], 
              tissue_data[tissue_data$gene %in% common_genes,], by = "gene")
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
dat1 <-na.omit(dat1[EVT_markers$genes,])
groupdata$Sample <- ifelse(groupdata$group%in% 'Placenta','Placenta','Others')
groupdata$id <- gsub('-','.',groupdata$id)
colnames(dat1) <- gsub('-','.',colnames(dat1))
#########计算Ttest#############
dfClass <- groupdata[,-2]
df = dat1 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-1,names_to = "id",values_to = "value") %>%  # 转换成长数据
  left_join(dfClass,by=c("id" = "id")) # 与分组数据合并
#平均值
mean_expression <- df %>%
  group_by(Sample, Gene) %>%
  summarise(MeanExpression = mean(value))
mean_expression<- mean_expression %>%
  pivot_wider(names_from =Sample, values_from = MeanExpression)
#ttest
p_value_data <- df %>%
  group_by(Gene) %>%
  summarise(p_value = t.test(value ~ Sample)$p.value)
significant_genes <- p_value_data %>%
  filter(p_value < 0.01) 
####平均值差值
mean_expression$Relative_Difference <- (mean_expression$Placenta - mean_expression$Others)
select_gene <- (vo_data %>%
            filter(Relative_Difference > 0.5,P.Value < 0.01) )$Gene
