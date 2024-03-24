library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(tibble)
library(bitops)
library(rattle)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(randomForest)
library(venn)
library(dplyr)
library(tidyverse)
library(glmnet) 
library(xgboost)
library(DALEX)
library(gbm)
library(VennDiagram)
library(neuralnet)
library(NeuralNetTools)


#Data reading and processing
importantvlue=0
uniGenes <- 'train.uniSigExp.txt'
unidata <- read.table(uniGenes, header = T, row.names = 1, sep = '\t',check.names = F)
genedata <- t(unidata[,-grep('DSSEvent',colnames(unidata))])

#Set comparison group
con <- unidata[unidata$DSSEvent == '0', ]
treat <- unidata[unidata$DSSEvent == '1', ]

conNum = nrow(con)
treatNum = nrow(treat)

type = c(rep("con",conNum),rep("treat",treatNum))

merge2 <- rbind(con,treat)
merge2 <- merge2[,-grep('DSSEvent',colnames(merge2))]
merge2 <- t(merge2)

#Output
out=rbind(id=paste0(colnames(merge2),"_",type), merge2)
write.table(out, file="data.txt", sep="\t", quote=F, col.names=F)

data = read.table('data.txt',header = T, row.names = 1, sep = '\t',check.names =F)
data = t(data)
group=gsub("(.*)\\_(.*)", "\\2", rownames(data))


#####GBM#####
set.seed(1984)
metric <- "RMSE"
myControl <- trainControl(method="cv", number=5)
fitControl <- trainControl( method = "repeatedcv", number = 4, repeats = 4)
fit <- train(x=data,
             y=as.factor(group),  
             method = "gbm", trControl = fitControl,verbose = FALSE)
importances <- varImp(fit)
GBMimportance <- as.matrix(importances$importance)
GBMGenes=GBMimportance[order(GBMimportance[,"Overall"], decreasing = TRUE),]
GBMGenes=names(GBMGenes[GBMGenes>0])


#####随机森林#####
set.seed(1986)
rf=randomForest(as.factor(group)~., data=data, ntree=500)
optionTrees=which.min(rf$err.rate[,1])
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0]) 

#####决策树#####
data2=data
dim(data2)
colnames(data2)
aggr(data2)
set.seed(1988)
data2<-as.data.frame(data2)
#建模
mod1<-rpart(as.factor(group)~.,data = data2, method = "poisson")
#显示重要性
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>% arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])

#####lasso#####
set.seed(1990)
x=data
y=group
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGenes=lassoGene[-1]
lassoGenes


#####XGboost#####
set.seed(1990)
TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)
importance <- varImp(model)
XGBimportant <- as.matrix(importance$importance) 
XGBGenes=XGBimportant[order(XGBimportant[,"Overall"], decreasing = TRUE),]
XGBGenes=names(XGBGenes[XGBGenes>0])



####Venn画图####
library(venn)
venn(
  x = list(
    'XGBoost' = XGBGenes,
    'Random forest' = rfGenes,
    'Lasso' = lassoGenes,
    'Gradient boosting machine' = GBMGenes,
    'Decision tree' = JCSGenes),
  opacity = 0.3,  # 调整颜色透明度
  #zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，也可以自定义颜色
  zcolor = c("dodgerblue", "forestgreen", "darkorange1","yellow","darkorchid"),   #填充颜色
  col = 'gray10',
  #col =  c("dodgerblue", "forestgreen", "darkorange1","yellow","darkorchid"),   #边框颜色
  lwd = 0.8,   #边框粗细
  box = F,        # 是否添加边框
  ilcs = 0.9,     # 数字大小
  sncs = 0.8,        # 组名字体大小
  ilabels = F,  #使用表达式情况下是否展示counts数
  ellipse = F, #设置为椭圆。四、五个集合时可设置
)

#export:3x3

#####提取基因#####

#genes <- Reduce(intersect, list(A, B, C, D ...))   #可以试一试这个嵌套？https://www.jingege.wang/2022/12/22/r%e8%af%ad%e8%a8%80%e5%a4%9a%e4%b8%aa%e5%90%91%e9%87%8f%e5%8f%96%e4%ba%a4%e9%9b%86/

Gene1=as.data.frame(GBMGenes)
Gene2=as.data.frame(rfGenes)
Gene3=as.data.frame(JCSGenes)
Gene4=as.data.frame(lassoGenes)
Gene5=as.data.frame(XGBGenes)

sameSample=intersect(Gene1$GBMGenes, Gene3$JCSGenes)
sameSample=as.data.frame(sameSample)

sameSample=intersect(Gene5$XGBGenes,sameSample$sameSample)
sameSample=as.data.frame(sameSample)

sameSample=intersect(Gene4$lassoGenes,sameSample$sameSample)
sameSample=as.data.frame(sameSample)

sameSample=intersect(Gene2$rfGenes,sameSample$sameSample)

Genes=as.data.frame(sameSample)
colnames(Genes) <- 'Gene'
Genes=cbind(id=row.names(Genes),Genes)
write.table(Genes, file="5AI_Genes.txt",sep="\t",quote=F,row.names = F)
GenesExp=genedata[sameSample,,drop=F]
write.table(GenesExp, file="5AI_GenesExp.xls",sep="\t",quote=F)


####查看结果####
Genes



#####⛔️构建人工神经网络模型#####
#install.packages('neuralnet')
library(neuralnet)
library(NeuralNetTools)
rt=GenesExp
rt=as.matrix(rt)
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
diffRT=read.table("/Users/f/Library/CloudStorage/OneDrive-个人/Bioinformatics/A.Pyroptosis.RFS/03-2.newTCGA/014.model/AI4Genes/Count_vst_diff_PRGs226_FC1.5.xls", header=T, sep="\t", check.names=F, row.names=1)
diffRT=diffRT[row.names(data),]
dataUp=data[diffRT[,"log2FoldChange"]>0,]
dataDown=data[diffRT[,"log2FoldChange"]<0,]
dataUp2=t(apply(dataUp,1,function(x)ifelse(x>median(x),1,0)))
dataDown2=t(apply(dataDown,1,function(x)ifelse(x>median(x),0,1)))
outTab=rbind(dataUp2, dataDown2)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="GeneScore.txt", sep="\t", quote=F, col.names=F)

inputFile="GeneScore.txt"     
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data$con=ifelse(group=="con", 1, 0)
data$treat=ifelse(group=="treat", 1, 0)
fit=neuralnet(con+treat~., data, hidden=2)
fit$result.matrix
fit$weight
pdf(file="model.pdf",width=10, height=10)
plotnet(fit)
dev.off()

net.predict=compute(fit, data)$net.result
net.prediction=c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)
predict.table
conAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
treatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))
colnames(net.predict)=c("con", "treat")
outTab=rbind(id=colnames(net.predict), net.predict)
write.table(outTab, file="neural.predict.txt", sep="\t", quote=F, col.names=F)




