library(limma)
countData <- as.data.frame(t(dtrain[, -1]))  # 忽略第一列基因名称
condition <- model.matrix(~dtrain$target)
data <- list(counts=countData , design=condition)
fit <- lmFit(data$counts, design=data$design)
fit <- eBayes(fit)
results <- topTable(fit, number=Inf)
