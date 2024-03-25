#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(estimate)
inputFile="filterTPM.txt"       #表达数据文件01

#读取文件,并对输入文件进行整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=log2(data+1)    #⭐️TPM
 
#删除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#输出整理好的数据
out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol_TPM.txt",sep="\t",quote=F,col.names=F)

#运行estimate包
filterCommonGenes(input.f="uniq.symbol_TPM.txt", 
                  output.f="commonGenes_TPM.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes_TPM.gct",
              output.ds="estimateScore_TPM.gct")

#输出每个样品的打分
scores=read.table("estimateScore_TPM.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)

write.table(out, file="TMEscores_TPM.txt", sep="\t", quote=F, col.names=F)

#### 绘图 ####
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)
library(rstatix)   #显著性
library(gghalves)   #分面小提琴图

##数据处理和计算####
riskFile="totalRisk.txt"       #风险文件#12
scoreFile="TMEscores_TPM.txt"      #肿瘤微环境打分文件#24

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$Risk=factor(risk$Risk, levels=c("Low","High"))

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=as.matrix(score)
#row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))  #修改样品名称，保持一致
score=avereps(score)
score=score[,1:4]   #StromalScore\ImmueScore\ESTIMATEScore\TumorPurity   #⭐️计算3种score时，此行：score=score[,1:3]
#rt=rt[,c(1,5)]    #⭐️计算Risk和TumorPurity。

#样品取交集
sameSample=intersect(row.names(risk), row.names(score))
risk2=risk[sameSample,"Risk",drop=F]   #1列sameSample，⭐️1列“Risk”或者“RiskCO”
score=score[sameSample,,drop=F]  #1列sameSample，其余列保留“score”列数
rt=cbind(risk2, score)   #合并为sameSample+Risk+4列score
#rt=rt[,c(1:4)]    #⭐️计算Risk和Stroma
rt=rt[,c(1:5)]    #⭐️计算Risk和TumorPurity。

##离群值处理####
#定义函数，将离群值替换为NA
#离群值定义为为小于Q1(25%位置)－1.5IQR或大于Q3(75%位置)+1.5IQR的值。
re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   #⚠️⭐️  (0.25,0.75)
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}
#删除含有outliers(NA)的行
df_is <- rt%>%
  group_by(Risk)%>%
  mutate(ImmuneScore = re_outlier(ImmuneScore))
df_is<-df_is[complete.cases(df_is),]

##提取删除离群值后的样本名单
sam <- intersect(rt$TumorPurity, df_is$TumorPurity)
newdf<-rt[which(rt$TumorPurity%in%sam),]
newrow <- rownames(newdf)
#write.table(newrow, 'newrow.txt',sep="\t", quote=F, col.names=F, row.names = F)


##绘图前处理####
##将合并后的数据转换为ggplot2的输入文件
data_all=melt(df_is, id.vars=c("Risk"))   
colnames(data_all)=c("RiskGroup", "scoreType", "Score")

##自定义画图的ESTIMATE分数类型————(前3个score)
data <- data_all[data_all$scoreType == 'StromalScore'|
                   data_all$scoreType == 'ImmuneScore'| 
                   data_all$scoreType == 'ESTIMATEScore',]
data$scoreType = ifelse(data$scoreType=='StromalScore','Stromal score',
                        ifelse(data$scoreType=='ImmuneScore','Immune score','ESTIMATE score'))  ##修改单元格名称(前3个score)

##自定义画图的ESTIMATE分数类型—————(tumor)
data_tumor <- data_all[data_all$scoreType =='TumorPurity',]
data_tumor$scoreType = ifelse(data_tumor$scoreType=='TumorPurity','Tumor purity')##修改单元格名称(tumor)

##自定义画图的ESTIMATE分数类型————(4合1)
data_4in1 <- data_all[data_all$scoreType == 'StromalScore'|
                        data_all$scoreType == 'ImmuneScore'| 
                        data_all$scoreType == 'ESTIMATEScore'|
                        data_all$scoreType == 'TumorPurity',]
data_4in1$scoreType = ifelse(data_4in1$scoreType=='StromalScore','Stromal score',
                             ifelse(data_4in1$scoreType=='ImmuneScore','Immune score',
                                    ifelse(data_4in1$scoreType=='ESTIMATEScore','ESTIMATE score','Tumor purity')))   ##修改单元格名称————(4合1)
colnames(data_4in1)[colnames(data_4in1)=="RiskGroup"] <- "Risk group"   #修改列名

##计算显著性####
#data
data <- data %>%  group_by(scoreType)
data$Score <- as.numeric(data$Score)
stat.test_data <- wilcox_test(data = data,   #pairwise_
                              Score ~ RiskGroup, 
                              paired = F,   
                              p.adjust.method = "fdr",
                              #alternative = "less"
) %>%
  add_xy_position(x = "scoreType")
stat.test_data$p.scient <- format(stat.test_data$p, scientific = TRUE)
stat.test_data$p.round3 <- round(stat.test_data$p, 3)
stat.test_data$p.adj.signif <- ifelse(stat.test_data$p<0.001,"***",ifelse(stat.test_data$p<0.01,"**",ifelse(stat.test_data$p<0.05,"*","")))
stat.test_data$p.value <- ifelse(stat.test_data$p<0.001, format(stat.test_data$p, scientific = TRUE),round(stat.test_data$p, 3))
#data_tumor
data_tumor <- data_tumor %>%  group_by(scoreType)
data_tumor$Score <- as.numeric(data_tumor$Score)
stat.test_data_tumor <- wilcox_test(data = data_tumor,
                                    Score ~ RiskGroup, 
                                    paired = F,   
                                    p.adjust.method = "fdr",
                                    #alternative = "less"
) %>% add_xy_position(x = "scoreType")
stat.test_data_tumor$p.scient <- format(stat.test_data_tumor$p, scientific = TRUE)
stat.test_data_tumor$p.round3 <- round(stat.test_data_tumor$p,3)
stat.test_data_tumor$p.adj.signif <- ifelse(stat.test_data_tumor$p<0.001,"***",ifelse(stat.test_data_tumor$p<0.01,"**",ifelse(stat.test_data_tumor$p<0.05,"*","")))
stat.test_data_tumor$p.value <- ifelse(stat.test_data_tumor$p<0.001, format(stat.test_data_tumor$p, scientific = TRUE),round(stat.test_data_tumor$p, 3))


#### ggplot2画分半小提琴图（3score）####
mycol = c("#e11a0c","#337cba")    #c("#fb3e35",'#564efd')
p1 <- ggplot()+
  geom_half_violin(data = data %>% filter(RiskGroup == "High"),
                   aes(x = scoreType, y = Score), side = "l",size= 1,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)
  ) +
  geom_half_violin(data = data %>% filter(RiskGroup == "Low"),
                   aes(x = scoreType,y = Score), side = "r", size= 1,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)
  ) +
  #中间添加：分半箱线图
  geom_half_boxplot(data = data %>% filter(RiskGroup == "High"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[1], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[1],side = "l", alpha = 0.5, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data %>% filter(RiskGroup == "Low"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[2], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1)
  ) +
  #添加平均数白点
  #geom_line(data = data, aes(x = scoreType, y = Score, fill = RiskGroup),
  #           stat = 'summary', fun=mean, col= 'white', pch = 15, size =1,
  #           position = position_dodge(width = 0.2), show.legend = FALSE
  #           ) +
  #添加图例
  geom_line(data = data, aes(x = scoreType, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, size =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="PyroGroup",labels = c("High-risk","Low-risk"))+
  #显著性差异
  stat_pvalue_manual(
    stat.test_data, 
    label = 'P = {p.value}\n{p.adj.signif}',   #"p.adj.signif", 
    bracket.size = 0.3, # 粗细
    bracket.shorten = 0.15,   #宽度
    tip.length = 0.01,  # 两边竖线的长度
    size = 3.5
  ) + 
  xlab("") +
  ylab("Tumor microenvironment score") +
  scale_y_continuous(limits = c(-2500,4300)) +    #TPM
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5)+
  theme(
    #x轴标签
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black', margin = margin(0.2,0,0,0, 'cm')),    #x轴标签
    #轴刻度
    axis.ticks =element_line(linewidth = 0.3),
    #y轴标签
    axis.text.y = element_text(color = 'gray30', hjust = 1, # 左对齐
                               size = 9, lineheight = 1),
    #标题居中:
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    #图例：
    legend.position = "top",   #图例位置
    legend.key.size = unit(10,'pt'), #图例大小
    legend.justification = "centre")


#### ggplot2画分半小提琴图（tumor）####
p2 <- ggplot()+
  geom_half_violin(data = data_tumor %>% filter(RiskGroup == "High"),
                   aes(x = scoreType, y = Score), side = "l",size= 1,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)) +
  geom_half_violin(data = data_tumor %>% filter(RiskGroup == "Low"),
                   aes(x = scoreType,y = Score), side = "r", size= 1,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)) +
  #中间添加：分半箱线图
  geom_half_boxplot(data = data_tumor %>% filter(RiskGroup == "High"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[1], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[1],side = "l", alpha = 0.5, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data_tumor %>% filter(RiskGroup == "Low"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[2], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1)) +
  #添加平均数白点
  #geom_point(data = data_tumor, aes(x = scoreType, y = Score, fill = RiskGroup),
  #           stat = 'summary', fun=mean, col='white', pch = 15, size =1,
  #           position = position_dodge(width = 0.2),
  #           show.legend = FALSE) +
  #添加图例
  geom_line(data = data_tumor, aes(x = scoreType, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, size =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="PyroGroup",labels = c("High-risk","Low-risk"))+
  #显著性差异
  stat_pvalue_manual(
    stat.test_data_tumor, 
    label = 'P = {p.value}\n{p.adj.signif}', #"p.adj.signif", 
    bracket.size = 0.3, # 粗细
    bracket.shorten = 0.15,   #宽度
    tip.length = 0.01,  # 两边竖线的长度
    y.position = 1,
    size = 3.5
  ) + 
  xlab("") +
  ylab("Percentage") +
  scale_y_continuous(position = "right", limits = c(0.45,1.1)) +   #TPM
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5)+
  theme(
    #x轴标签
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black', margin = margin(0.2,0,0,0, 'cm')),    #x轴标签
    #轴刻度
    axis.ticks =element_line(linewidth = 0.3),
    #y轴标签
    axis.text.y = element_text(color = "gray30", size = 9, lineheight = 1, hjust = 1), # 左对齐
    axis.title.y.right = element_text(vjust = 2),
    #标题居中:
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    #图例：
    legend.position = 'top',   #图例位置
    legend.key.size = unit(10,'pt'), #图例大小
    legend.justification = "centre") 
#labs(fill = "Risk group") 

##拼接图片
library(patchwork)
p1 + p2 + 
  plot_annotation(#title = "Wilcox.test, ** P<0.01, * P<0.05"   #左上角标题
    #caption = "Wilcox.test, ** P<0.01, * P<0.05",  #右下角标题
    #tag_levels = "A",   #每张图的序号
  ) +
  plot_layout(ncol=2,   #图像设置为2列，默认按列填充    # nrow=  按行填充
              widths = c(3, 1),   # #两列之间相对宽度比为3：1   # heights=c(2,1)  相对高度
              guides='collect'   #合并图例为1个
  ) & theme(legend.position='top',
            legend.key.size = unit(0.3,"cm"),
            legend.text = element_text(size = 10)) 

ggsave("4score.pdf", height = 3.5, width = 6)   #new =4x6   new2=5x6




#### Spearman相关性 #####
#数据处理
riskScore = risk[sameSample,"RiskScore",drop=F]
riskScore = cbind(riskScore,score)

re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   #⚠️⭐️  (0.25,0.75)
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}

#删除含有outliers(NA)的行
riskScore2 <- riskScore%>%
  mutate(RiskScore = re_outlier(RiskScore))
riskScore2 <-riskScore2[complete.cases(riskScore2),]

#计算相关系数
SpearmanR <- cor(riskScore2$RiskScore, riskScore2$ESTIMATEScore,method="spearman",use="complete.obs")
SpearmanP <- cor.test(riskScore2$RiskScore, riskScore2$ESTIMATEScore, method="spearman", use="complete.obs")

#绘图
mycol2=c("#e11a0c","#48af45","#337cba","orange")
ggplot(riskScore2, aes(x=RiskScore, y=ESTIMATEScore)) +
  geom_point(color = mycol2[1], alpha=0.7, pch=20, cex=2) +
  geom_smooth(method=lm , formula = y ~ x, 
              #置信区间：
              color=mycol2[1], fill=mycol2[1], alpha = 0.3, se=TRUE) +
  theme_test(base_line_size = 0.3)+
  ylab("ESTIMATE score") + 
  xlab("PyroScore") +
  theme(
    # 去除网格线：
    panel.grid = element_blank(),
    # 修改坐标轴标签
    axis.title = element_text(size = 10),
    axis.title.y = element_text(vjust = 0),
    axis.text = element_text(color = "gray30",size = 9)
  ) + 
  #ESTIMATE:  total: x = 4.4, total2: x = 3
  annotate("text", x = 4.4, y = max(riskScore2$ESTIMATEScore), fontface = 1, hjust = 1, label = paste0("rho = ",round(SpearmanR,3), "\n","P = ", round(SpearmanP$p.value,3)))
ggsave("ESTIMATEScore.pdf", height = 3, width = 3.2)

