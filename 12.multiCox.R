library("glmnet")
library("survival")

##数据处理####
inputFile="train.expTime.txt"     
geneFile='5AI_Genes.txt'

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
rt$DSS=rt$DSS/365   #⭐️
gene=read.table(geneFile, header=T, sep="\t", check.names=F, row.names=1)    
sameGene=intersect(colnames(rt),gene$Gene)
data=rt[,c('DSS','DSEvent',sameGene)]

#### 多因素独立预后分析####
multiCox=coxph(Surv(DSS, DSEvent) ~ ., data = data)   #⭐️
#basehaz <- survfit(multiCox)   #查看基准风险函数
#plot(basehaz)

#multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=as.data.frame(cbind(id=row.names(multiTab),multiTab))
multiOutTab=multiTab[as.numeric(multiTab[,"pvalue"])<1,]  #⭐️⭐️⭐️
write.table(multiOutTab,"train.multiCox.txt",sep="\t",row.names=F,quote=F)

#输出多因素显著基因的表达量
multiSigExp=data[,c("DSS","DSEvent",row.names(multiOutTab))]
multiSigExp=cbind(id=row.names(multiSigExp),multiSigExp)
write.table(multiSigExp,file="train.multiSigExp.txt",sep="\t",row.names=F,quote=F)

#输出相关基因系数
coef=coef(multiCox, s=multiCox$coefficients)
index=which(coef != 0)
actCoef=coef[index]
coef=as.data.frame(coef)
multiGene=rownames(coef)
geneCoef=as.data.frame(cbind(Gene=multiGene, Coef=actCoef))
geneCoef=geneCoef[rownames(multiOutTab),]
write.table(geneCoef, file="multi_geneCoef.txt", sep="\t", quote=F, row.names=F)


####🐒用predict函数计算RiskScore（只能计算全部多因素基因的Score，即P<1）####
#输出train组风险值
trainScore=predict(multiCox, type="risk", newdata=data)
trainScore=log2(trainScore+1)   #⭐️
coxGene=rownames(multiOutTab)
#coxGene=gsub("`","",coxGene)
outCol=c("DSS","DSEvent",coxGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"High","Low"))
outTab=cbind(data[,outCol],RiskScore=as.vector(trainScore),Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="trainRisk.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
testFile='test.expTime.txt'
rt=read.table(testFile, header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
testScore=log2(testScore+1)   #⭐️
outCol=c("DSS","DSEvent",coxGene)
risk=as.vector(ifelse(testScore>median(trainScore),"High","Low"))
outTab=cbind(rt[,outCol],RiskScore=as.vector(testScore),Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRisk.txt",sep="\t",quote=F,row.names=F)

#输出total组风险值
rt=read.table("total.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
testScore=log2(testScore+1)   #⭐️
outCol=c("DSS","DSEvent",coxGene)
risk=as.vector(ifelse(testScore>median(trainScore),"High","Low"))
outTab=cbind(rt[,outCol],RiskScore=as.vector(testScore),Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="totalRisk.txt",sep="\t",quote=F,row.names=F)



####🦁手动计算加权分数####
gene = geneCoef
rownames(gene) <- gene$Gene
gene <- as.data.frame(t(gene))
gene <- gene[2,]

exp <- multiSigExp
exp <- multiSigExp[,3:ncol(multiSigExp)]

vec <- as.numeric(gene[1,])
mat <- matrix(rep(vec, nrow(exp)), nrow = nrow(exp), byrow = TRUE)
cal <- data.frame(mat, row.names = rownames(exp))
#cal <- data.frame(matrix(rep(as.numeric(gene[1,]), nrow(exp)), nrow = nrow(exp), byrow = TRUE), row.names = rownames(exp))

exp$RiskScore <- rowSums(exp[, colnames(exp) %in% colnames(gene)] * cal)
#exp$RiskScore <- rowSums(exp * cal)
mydata <- cbind(multiSigExp[,1:2],exp)


####🦁手动进行RiskGroup分组####
trainScore = mydata$RiskScore
#risk=as.vector(ifelse(trainScore>median(trainScore),"High","Low"))
risk <- as.vector(ifelse(trainScore > quantile(trainScore, 0.5), "High", "Low"))   #定义trainScore的前40%为High，后60%为Low
outTab=cbind(mydata,Risk=risk)
write.table(outTab,file="trainRisk_H.txt",sep="\t",quote=F,row.names=F)

#test
testFile='test.expTime.txt'
rt=read.table(testFile, header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365
coxGene=rownames(multiOutTab)
test=rt[,coxGene]
cal_test <- data.frame(matrix(rep(as.numeric(gene[1,]), nrow(test)), nrow = nrow(test), byrow = TRUE), row.names = rownames(test))
test$RiskScore <- rowSums(test[,colnames(test) %in% colnames(gene)] * cal_test)
mydata <- cbind(rt[,1:2],test)
testScore = mydata$RiskScore
#risk=as.vector(ifelse(testScore>median(trainScore),"High","Low"))
risk <- as.vector(ifelse(testScore > quantile(trainScore, 0.5), "High", "Low"))   #定义trainScore的前40%为High，后60%为Low
outTab=cbind(mydata,Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRisk_H.txt",sep="\t",quote=F,row.names=F)

#total
rt=read.table("total.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365
coxGene=rownames(multiOutTab)
total=rt[,coxGene]
cal_total <- data.frame(matrix(rep(as.numeric(gene[1,]), nrow(total)), nrow = nrow(total), byrow = TRUE), row.names = rownames(total))
total$RiskScore <- rowSums(total[,colnames(total) %in% colnames(gene)] * cal_total)
mydata <- cbind(rt[,1:2],total)
totalScore = mydata$RiskScore
#risk=as.vector(ifelse(totalScore>median(trainScore),"High","Low"))
risk <- as.vector(ifelse(totalScore > quantile(trainScore, 0.5), "High", "Low"))   #定义trainScore的前40%为High，后60%为Low
outTab=cbind(mydata,Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="totalRisk_H.txt",sep="\t",quote=F,row.names=F)


#### 求train组的cutoff####
rt=read.table('trainRisk.txt', header=T, sep="\t", check.names=F)
#计算ROC曲线
ROC_rt=timeROC(T=rt$DSS, delta=rt$DSEvent,    #⭐️
               marker=rt$RiskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
#计算cutoff
library(tidyverse)
library(survivalROC)
auc_text=c()
you_roc <- survivalROC(Stime=rt$DSS,
                       status = rt$DSEvent,
                       marker = rt$RiskScore,
                       predict.time = 5,   #⭐️
                       method = "KM")

cutoff_5years <- you_roc$cut.values[which.max(you_roc$TP-you_roc$FP)]
cutoff_5years
y1 <- you_roc$TP[you_roc$cut.values==cutoff_5years]   #0.7591272
x1 <- you_roc$FP[you_roc$cut.values==cutoff_5years]   #0.2710061

#绘图
pdf(file='ROC_train.pdf',width=5,height=5)
#新建画布
plot(you_roc$FP,you_roc$TP, xlab="", ylab="", col='white')
#背景参考线
abline(h = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
       v = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
       col = gray(0.99,0.95))
par(new=T)

#ROC曲线
plot(ROC_rt,time=1,col='green2',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)

#图例
legend('bottomright',   #x=0.4,y=0.3,
       c(paste0('1-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
         paste0('3-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[3]), ')'),
         paste0('5-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[5]), ')')
       ),
       col=c('green','blue','red'), lwd=2, bty = 'n')  

#加箭头
arrows(x0=0.22, y0=0.84,x1=x1-0.02,y1=y1+0.02,
       length = 0.08, angle = -20, code = 2,col = "red2", lwd = 1.5, lty = 1)
#加文本
text(0.2,0.92, labels = paste("Cutoff value: ",round(cutoff_5years,3)), col='red2', cex=0.9, font=2)

dev.off()


#### 根据cutoff分组 ####
rt=read.table('trainRisk.txt', header=T, sep="\t", check.names=F)
rt=select(rt,-Risk)
rt$Risk = as.vector(ifelse(rt$RiskScore>=cutoff_5years,"High","Low"))
write.table(rt,file="trainRisk_CO.txt",sep="\t",quote=F,row.names=F)

rt=read.table('testRisk.txt', header=T, sep="\t", check.names=F)
rt=select(rt,-Risk)
rt$Risk = as.vector(ifelse(rt$RiskScore>=cutoff_5years,"High","Low"))
write.table(rt,file="testRisk_CO.txt",sep="\t",quote=F,row.names=F)

rt=read.table('totalRisk.txt', header=T, sep="\t", check.names=F)
rt=select(rt,-Risk)
rt$Risk = as.vector(ifelse(rt$RiskScore>=cutoff_5years,"High","Low"))
write.table(rt,file="totalRisk_CO.txt",sep="\t",quote=F,row.names=F)


#### 🖤辅助方法：可用来快速查看cutoff数值 #####
library(pROC)
rt=read.table('testRisk.txt', header=T, sep="\t", check.names=F, row.names=1)  

# 生成roc曲线对象
rocobj <- roc(rt$DSEvent, rt$RiskScore)

# 绘制roc曲线
plot(rocobj,
     legacy.axes = F,
     main="ROC curve cutoff",
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best", # 在roc曲线上显示最佳阈值点(敏感度、特异度)
)



#### 绘制双向柱状图 ####
data <- read.table('multi_geneCoef.txt', header=T, sep="\t", check.names=F)

# 初步绘图：
ggplot(data)+
  geom_col(aes(reorder(Gene, Coef), Coef))+
  theme_classic()+
  ylim(-20,20)+
  coord_flip()

# 添加颜色，调整主题：
# 设置颜色变量：
color <- rep('red', 9)
color[which(data$Coef < 0)] <- "blue"
data$color <- color
mycol = c("#337cba","#e11a0c")

# 美化绘图：
ggplot(data)+
  geom_col(aes(reorder(Gene, Coef), Coef, fill = color),width = 0.6)+
  scale_fill_manual(values = mycol)+
  geom_segment(aes(y = 0, yend = 0,x = 0, xend = 9), lty=1, lwd=0.3)+
  coord_flip()+
  # 调整主题：
  theme_test(base_size = 10, base_line_size = 0.5,base_rect_size = 0.6)+
  theme(
    #坐标轴字体：
    axis.title = element_text(colour = 'black', size = 10),
    # 去除图例：
    legend.position = "none",
    # 标题居中：
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black", size=10),
    axis.text.x = element_text(color = 'black', size = 10)
  ) +
  ylab("Multivariate Cox regression coefficient")+
  # 添加label：
  geom_text(data=data[which(data$Coef > 0), ],aes(x = Gene, y = 0, label = sprintf("%.03f",Coef)), 
            hjust = 1.1, size = 4)+
  geom_text(data=data[which(data$Coef < 0), ],aes(x = Gene, y = 0, label = sprintf("%.03f",Coef)), 
            hjust = -0.1, size = 4)+
  
  scale_x_discrete(expand=expansion(add=c(0,0)))

dev.new()
print(p)

ggsave("barplot.pdf", height = 3.5, width = 3.8)

