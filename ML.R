library(caret)
library(Matrix)
library(dplyr)
setwd('D:\\1.paperdata/cfrna/1.ML')
set.seed(4180)
########数据准备#########
dat <- data.frame(t(datExpr0))
genedata <- test[,c(1,6)]
genedata <- subset(genedata,ensembl_gene_id %in% rownames(dat))
genedata<- genedata %>% 
            distinct()
genedata <- genedata[!duplicated(genedata$ensembl_gene_id), ]
genedata <- genedata[!duplicated(genedata[, 2]), ]
dat <- dat[genedata$ensembl_gene_id,]
dat$gene <- genedata$hgnc_symbol
dat <- na.omit(dat)
dat <- dat%>%
  group_by(gene) %>%
  summarise_all(mean)###相同基因取平均
dat <- data.frame(dat )
row.names(dat)<- dat$gene
dat <- dat[,-1]
colnames(dat) <- gsub("X",'',colnames(dat))
dattest <- data.frame(t(na.omit(dat[colnames(dtrain)[-1],])))
#########
select_gene <- read.csv("D:/1.paperdata/placentascrna/4.ML/select_gene.csv")
EVT_markers <- cfdeg[,1]
EVT_m <- na.omit(match(EVT_markers,row.names(exprdt)))
dat <- exprdt[EVT_m,]
EVT_markers <- cfgenes[,2]
EVT_m <- na.omit(match(EVT_markers,row.names(dat)))
dat <- dat[EVT_m,]
dat <- t(dat)
group_data$target <- ifelse(group_data$`disease:ch1` %in% 'control','0','1')
dat <- data.frame(target = group_data$target,dat)
colnames(group_data) <- gsub('sampling time group:ch1','time',colnames(group_data))
dtrain <-dat[subset(group_data,time %in% c('≤12w','13-20w'))$title,]
rownames(dtrain) <- gsub('X','',rownames(dtrain))
parts = createDataPartition(dtrain$target, p = .7, list = F)
train1 = dtrain[parts, ]
test1 = dtrain[-parts, ]
saveRDS(dtrain,'dtrain.Rds')
saveRDS(test1,'train1.Rds')
#############全机器学习###########
library(tidymodels)
library(discrim)
library(bonsai)
tidymodels_prefer()
df <- dtrain %>% 
  mutate(target =factor(target)) %>% 
  na.omit()
rec <- recipe(target ~ ., data = df )
dt_split<- initial_split(df, strata =target)
####设置模型####
library(rules)
library(baguette)
# 神经网络
nnet_spec<- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet") %>% 
  set_mode("classification")
# 多元自适应样条
mars_spec<- 
  mars(prod_degree = tune()) %>% 
  set_engine("earth") %>% 
  set_mode("classification")
# SVM_rbf
svm_r_spec<- 
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")
#SVM_ploy
svm_p_spec<- 
  svm_poly(cost = tune(), degree = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")
# KNN
knn_spec<- 
  nearest_neighbor(neighbors = tune(), dist_power = tune(), weight_func = tune()) %>% 
  set_engine("kknn") %>% 
  set_mode("classification")
#决策树
cart_spec<- 
  decision_tree(cost_complexity = tune(), min_n = tune()) %>% 
  set_engine("rpart") %>% 
  set_mode("classification")
# 装袋决策树
bag_cart_spec<- 
  bag_tree() %>% 
  set_engine("rpart", times = 50L) %>% 
  set_mode("classification")
#随机森林
rf_spec<- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")
# xgboost
xgb_spec<- 
  boost_tree(tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), 
             min_n = tune(), sample_size = tune(), trees = tune()) %>% 
  set_engine("xgboost") %>% 
  set_mode("classification")
# lgbm
lgbm_spec<- 
  boost_tree(tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), 
             min_n = tune(), trees = tune()) %>% 
  set_engine("lightgbm") %>% 
  set_mode("classification")
# 逻辑回归
logistic_spec <-          
  logistic_reg(penalty = tune()) %>%          
  set_engine('glmnet')  
####设置工作集####
all_workflows<- 
  workflow_set(
    preproc = list(rec), 
    models = list(neural_network = nnet_spec,
                  logistic_reg= logistic_spec,
                  SVM_radial = svm_r_spec, 
                  SVM_poly = svm_p_spec, 
                  KNN = knn_spec, 
                  MARS = mars_spec, 
                  CART = cart_spec, 
                  CART_bagged = bag_cart_spec,
                  RF = rf_spec, 
                  xgb = xgb_spec,
                  lgbm = lgbm_spec)
 ,cross = T ) %>% 
  mutate(wflow_id = gsub("recipe_", 
                         "", 
                         wflow_id))
# 设置重采样
folds <- vfold_cv(df,
                  repeats = 5,v=10)
# 控制条件，保存预测值      
ctr <- control_grid(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )
# 模型拟合 
library(doParallel) 
cl <- makePSOCKcluster(8) # 加速，用12个线程
registerDoParallel(cl)
grid_results<-
  all_workflows %>%
  workflow_map(
    seed = 4180,
    resamples = folds,
    grid = 10,
    control = ctr,
    verbose = TRUE
  )
stopCluster(cl)
#########筛选模型###########
grid_results%>% 
  rank_results() %>% 
  filter(.metric == "rmse") %>% 
  select(model, .config, rmse = mean, rank)
##################可视化########
library(cols4all)
library(ggplot2)
cols <- c4a('20',15)
pal <- colorRampPalette(cols)
####
rank_results(grid_results,rank_metric = "roc_auc") %>% 
  filter(.metric=="roc_auc") %>% 
  select(wflow_id,mean)
###ROC
collect_predictions(grid_results) %>%           
  group_by(wflow_id) %>%           
  roc_curve(target,.pred_0) %>%          
  ggplot(aes(x=1-specificity,y=sensitivity,color=wflow_id))+          
  geom_line(lwd=0.8)+theme_light()+          
  theme(legend.position = c(.98, .65),          
        legend.justification = c("right", "top"))+
   scale_color_manual(values=pal(11))+
  labs(title = "TrainSet ROC curve")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  coord_fixed()
####直方图
grid_results %>% autoplot(select_best = T)+theme_bw()
##########箱型图#####
rank_results(grid_results, select_best = TRUE,rank_metric = 'roc_auc') %>% 
  ggplot( aes(x = reorder(wflow_id, rank), y = mean, color= wflow_id)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err),
                position = position_dodge(width = 0.9), width = 0.25)+
  facet_wrap(~ .metric, scales = "free_y", nrow = 2) +
  labs(title = "Model Performance Comparison",
       x = "Model",
       y = "Mean Value") +
  theme_minimal() +scale_color_manual(values=pal(11))+
  theme(legend.position = "right") 

autoplot(grid_results, 
         id = "SVM_radial", 
         metric = "roc_auc")
########输出模型###############
library(pROC)  
best_results <- 
  grid_results %>% 
  extract_workflow_set_result("SVM_radial") %>% 
  select_best(metric = "roc_auc")
testdf <- test1 %>% 
  mutate(target =factor(target)) %>% 
  na.omit()
SVM_fit <- 
  grid_results %>% 
  extract_workflow("SVM_radial") %>% 
  finalize_workflow(best_results)%>% 
  last_fit(dt_split)
collect_metrics(SVM_fit)
collect_predictions(SVM_fit) %>%            
  roc_curve(target,.pred_0) %>%          
  ggplot(aes(x=1-specificity,y=sensitivity))+          
  geom_line(lwd=0.8)+theme_light()+          
  theme(legend.position = c(.98, .65),          
        legend.justification = c("right", "top"))+
  scale_color_manual(values=pal(11))+
  labs(title = "TestSet ROC curve(AUC=0.701)")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+  coord_fixed()
