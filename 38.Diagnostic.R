library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(glmnet)

rm(list = ls()) 

####TCGA####
TCGA = read.table("PyroExp.txt", header = T, check.names = F, sep = "\t", row.names = 1)   #02
TCGA = TCGA["SEZ6L2",]
#TCGA = TCGA[rownames(TCGA) %in% c("SEZ6L2"),]
TCGA = as.data.frame(t(TCGA))
TCGA$SEZ6L2 = log2(TCGA$SEZ6L2+1)

# 创建新列Type
TCGA$group <- ifelse(substr(rownames(TCGA), nchar(rownames(TCGA)) - 2, nchar(rownames(TCGA)) - 2) == "1", "Normal", "Tumor")
#TCGA$group <- ifelse(sapply(strsplit(sapply(strsplit(rownames(TCGA),"\\-"), "[", 4),""), "[", 1) == "1", "Normal", "Tumor")


####GPL96#####
GPL96 = read.table("./GPL96/bindGEO_GLP96_remove.txt", header = T, check.names = F, sep = "\t", row.names = 1)
GPL96 = GPL96["SEZ6L2",]
GPL96 = as.data.frame(t(GPL96))
GPL96.group = read.csv("./GPL96/group.csv",row.names = 1)
GPL96.merge <- merge(GPL96, GPL96.group, by = "row.names", all = TRUE)
rownames(GPL96.merge) = GPL96.merge[,1]
GPL96.merge = GPL96.merge[,-1]


####GPL570#####
GPL570 = read.table("./GPL570/bindGEO_GLP570_ComBat.txt",header = T, check.names = F, sep = "\t", row.names = 1)
GPL570 = GPL570["SEZ6L2",]
GPL570 = as.data.frame(t(GPL570))
GPL570.group = read.csv("./GPL570/group.csv",row.names = 1)
GPL570.merge <- merge(GPL570, GPL570.group, by = "row.names", all = TRUE)
rownames(GPL570.merge) = GPL570.merge[,1]
GPL570.merge = GPL570.merge[,-1]

####临床PCR#####
PCR = read.table("PCR.txt",header = T, check.names = F, sep = "\t")


#### 绘制ROC曲线 #### 
roc1 <- roc(TCGA$group,TCGA$SEZ6L2, ci = T) 
roc2 <- roc(GPL96.merge$group,GPL96.merge$SEZ6L2, ci = T) 
roc3 <- roc(GPL570.merge$group,GPL570.merge$SEZ6L2, ci = T) 
roc4 <- roc(PCR$group,PCR$SEZ6L2, ci = T) 

mycol = c("#e11a0c","#337cba","#48af45", '#FDBF6F')
plot(roc1, title=FALSE, lwd=3, legacy.axes = F,
     xlab="1-Specificity", ylab="Sensitivity", col= mycol[1],
     cex.main=1, cex.lab=1, cex.axis=1, font=1,
     #auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     #auc.polygon.col= alpha(mycol[1],0.1),  # 设置多边形的填充颜色
     grid= T,   # 显示网格背景线
     )
plot.roc(roc2, add=TRUE, 
         col = mycol[2], lwd =3)  
plot.roc(roc3, add=TRUE,  
         col = mycol[3], lwd =3)  
plot.roc(roc4, add=TRUE, 
         col = mycol[4], lwd =3)  

legend('bottomright',
       c(paste0('TCGA (AUC=',sprintf("%.03f",roc1$auc[1]),')'),
         paste0('GPL96 (AUC=',sprintf("%.03f",roc2$auc[1]),')'),
         paste0('GPL570 (AUC=',sprintf("%.03f",roc3$auc[1]),')'),
         paste0('Clinical (AUC=',sprintf("%.03f",roc4$auc[1]),')')),
       col= mycol,
       lwd = 3, cex = 1, bty = 'n', 
       x.intersp = 0.1, 
       y.intersp = 1)

####依次计算诊断指标####
#install.packages("reportROC")
library(reportROC)

#计算
reportROC(PCR$group,
          PCR$SEZ6L2,
          important = "se",
          plot = TRUE)

#保存为excel可以打开的格式
ROC.info <- reportROC(TCGA$group,
                      TCGA$SEZ6L2,
                      important = "se",
                      plot = F)

write.csv(ROC.info, "ROC.info.TCGA(log2).csv")


#####其他#####
#先构建逻辑回归模型
GPL96.merge2 = GPL96.merge
GPL96.merge2$group = ifelse(GPL96.merge2$group =='PTC',1,0)

fit1 <- glm(group ~ SEZ6L2,
            data = GPL96.merge2,
            family = binomial())  

summary(fit1)

#使用predict()函数，可以观察某个预测变量在各个水平时对结果概率的影响。
GPL96.merge2$prob <- predict(fit1, 
                             newdata=GPL96.merge2, 
                             type="response")

