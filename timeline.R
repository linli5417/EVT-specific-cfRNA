library(dplyr)
library(ggplot2)
library(cols4all)
set.seed(4180)
######颜色#######
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
setwd("D:/1.paperdata/cfrna/")
#########数据载入################
exprdt <- datExpr0[group_data$title,]
ng = intersect(colnames(exprdt),genedata$gene_num) 
exprdt <- data.frame(t(exprdt[,ng ]))
genedata <- subset(genedata,gene_num %in%ng )
exprdt$gene <- genedata$gene_name
exprdt <- exprdt%>%
  group_by(gene) %>%
  summarise_all(mean)##取平均
exprdt <-as.data.frame(exprdt)
row.names(exprdt) <- exprdt[,1]
exprdt <- exprdt[,-1]
colnames(exprdt) <- gsub('X','',colnames(exprdt))
saveRDS(group_data,'group_data')
saveRDS(exprdt,'exprdt.Rds')
#######数据选择######

group_data$`sampling time group:ch1` <- gsub(' ','',group_data$`sampling time group:ch1` )
group_data<-readRDS(group_data,'group_data')
exprdt<-readRDS(exprdt,'exprdt.Rds')
########绘图数据#######
dat.fig <- data.frame(t(exprdt[inter_gene,]),time=factor(group_data$`sampling time group:ch1`,
                                                                      levels =c('≤12w','13-20w',
                                                                                '≥23w','Post-partum')))
dat.fig <- data.frame(EVT_score=t(gsva_mat),time=factor(group_data$`sampling time group:ch1`,
                                                           levels =c('≤12w','13-20w',
                                                                     '≥23w','Post-partum')))[group_data$title,]
dat.fig$Disease <- group_data$`disease:ch1`

#####绘图########
library(tidyverse)
library(dplyr)
dat.fig <- na.omit(dat.fig)
dat.fig <- gather(dat.fig, key = "gene", value = "expression", -c(time, Disease))
dat.fig$gene <- factor(dat.fig$gene,levels =inter_gene )
dat.sum <- dat.fig%>%
  group_by(time,Disease, gene) %>%
  summarise(
    GSVA_score = mean(expression),
    std_err = sd(expression) / sqrt(n())
  )

ggplot(dat.sum,aes(x=time,y=GSVA_score,group=Disease))+
  facet_wrap(gene~.)+
  geom_point(aes(color=Disease),size=1)+
  geom_line(aes(color=Disease),
            show.legend = FALSE,size=1.5)+
  geom_errorbar(aes(ymin=GSVA_score-std_err,
                    ymax=GSVA_score+std_err,color=Disease),
                width=0.1)+theme_bw() + 
  scale_color_manual(values=pal(3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=8,angle =0,hjust = .5,face = "bold"),
        strip.text = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 12))+labs(x='')
