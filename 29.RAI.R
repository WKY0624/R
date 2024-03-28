#### TCGA队列 ####
# 输入文件
Riskfile = 'totalRisk.txt'  #12
RAIfile = "RAI.txt"  #DataSets，根据ATA定义将TCGA病例分为3个组别：non-RAI、RAIS、RAIR（具体方法详见论文章节3.2.8）

# 读取Score文件
risk=read.table(Riskfile, header = T, sep = '\t', check.names = F)
rownames(risk) = risk[,1]

# 读取二分类文件（RAI队列）
RAI = read.table(RAIfile, header = T, sep = '\t', check.names = F)
RAI = RAI[RAI$RAI %in% c("Refractory","Sensitive"),,drop=F]
rownames(RAI) = RAI[,1]

# 合并Socre与分类
sameSample=intersect(risk$id,RAI$ID)
RAI=RAI[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]

rt = merge(RAI, risk, by.x = 'ID', by.y = "id")   #将两个数据框按照ID列进行合并
#write.table(rt, file = 'merge.txt', sep="\t", quote=F, col.names= T, row.names = F)

#### GEO队列 ####
##下载平台文件
gset = getGEO('GSE151179', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]

## 读取平台中分组信息
pdata <- pData(gset[[1]])
group_list <- ifelse(str_detect(pdata$characteristics_ch1.1,"tissue type: Primary tumor"),"PTC","Other")
group_list = factor(group_list,levels = c("PTC","Other"))
group_list
pdata$PTCgroup = group_list

RAIlist <- ifelse(str_detect(pdata$characteristics_ch1.4,"patient rai responce: Refractory"),"RR","RS")
RAIlist = factor(RAIlist, levels = c("RR","RS"))
pdata$RAIgroup = RAIlist

pdata3 = pdata[,c('geo_accession', 'PTCgroup', 'RAIgroup')]
PTC_RAI = pdata3[pdata3$PTCgroup == "PTC",]
write.table(PTC_RAI, "PTC_RAI.txt", sep = "\t", row.names = F,quote = F)

##合并
sameSample = intersect(colnames(exprSet),rownames(PTC_RAI))
exp = exprSet[,sameSample,drop=F]
out <- cbind(id=row.names(exp),exp)
write.table(out,file="PTC_RAI_exp.txt",sep="\t",row.names=F,quote=F)

##从矩阵中提取出模型基因
#读取模型基因
GeneCoef = read.table("multi_geneCoef.txt", sep = "\t", header = T, check.names = F)   #12
modelGene = GeneCoef$Gene

exp2 <- exp
index <- which(rownames(exp2) == "HIST1H2BG")   #查找"HIST1H2BG"的位置
rownames(exp2)[index] <- "H2BC8"  #将"HIST1H2BG"替换为"H2BC8"

expGene = exp2[modelGene,,drop=F]
expGene = as.data.frame(t(expGene))
expGene$ID = rownames(expGene)
merge = merge(expGene, PTC_RAI, by.x = "ID", by.y = "geo_accession")
rownames(merge) = merge$ID
merge = merge[,-grep("ID|PTCgroup",colnames(merge)), drop=F]

colnames(merge)[colnames(merge)=="RAIgroup"] <- "DSEvent"   #修改列名
merge$DSEvent = ifelse(merge$DSEvent == "RR","1","0")

##数据矫正
exp3 <- exp2
exp3 = normalizeBetweenArrays(exp3)

##根据表达矩阵，计算RiskScore
library("glmnet")
library("survival")

inputFile = "train.expTime.txt"      #09
geneFile = '5AI_Genes.txt'    #11
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
rt$DSS=rt$DSS/365
gene=read.table(geneFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
sameGene=intersect(colnames(rt),gene$Gene)
data=rt[,c('DSS','DSEvent',sameGene)]
#多因素独立预后分析
multiCox=coxph(Surv(DSS, DSEvent) ~ ., data = data) 
#用predict函数计算RiskScore
trainScore=predict(multiCox, type="risk", newdata=merge)
trainScore=log2(trainScore+1)   
risk=as.vector(ifelse(trainScore>0.96,"High","Low"))
outTab=cbind(merge[,c("DSEvent",modelGene)],RiskScore=as.vector(trainScore),Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="RiskScore_predict.txt",sep="\t",quote=F,row.names=F)
GEOfile = 'RiskScore_predict.txt'

#读取GEO-Score文件
GEO=read.table(GEOfile, header = T, sep = '\t', check.names = F)
rownames(GEO) = GEO[,1]
GEO = GEO[-1]
GEO = GEO[,grep("DSEvent|RiskScore",colnames(GEO))]


#### 云雨图：PyroScore差异 ####
library(ggplot2) #绘图
library(ggsignif) #添加统计检验
library(ggdist) #云雨图
## 读入数据
RAIFile = 'RAI.txt'  #DataSets
scoreFile = 'totalRisk.txt'  #12
RAI <- read.table(RAIFile, header = T, sep = '\t', check.names = F)
score <- read.table(scoreFile, header = T,sep = '\t', check.names = F)
merge <- merge(RAI, score, by.x = 'ID', by.y = 'id')
merge <- merge[,c('RAI','RiskScore')]

##离群值处理
#定义函数，将离群值替换为NA
#离群值定义为为小于Q1(25%位置)－1.5IQR或大于Q3(75%位置)+1.5IQR的值。
re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}
#删除含有outliers(NA)的行
merge2 <- merge %>%
  group_by(RAI)%>%
  mutate(RiskScore2 = re_outlier(RiskScore))
merge2 <-merge2[complete.cases(merge2),]

##统计计算：计算三组间的P值
pvalue = kruskal.test(merge2$RiskScore2 ~ merge2$RAI, data = merge2)
data_long <- merge2
data_long$RAI <- factor(data_long$RAI, levels = unique(RAI$RAI))

##此处生成一个列表，用于ggsignif两两组间的检验
Vec1 <- c("No","Sensitive","Refractory")
comb_list <- list()
for(i in 1:(length(Vec1)-1)) {   #1:2
  for(j in (i+1):length(Vec1)) {   #2:3
    comb <- combn(c(Vec1[i], Vec1[j]), 2)   #1，2,2
    if(!any(comb[1,] == comb[2,])) {  
      comb_list[length(comb_list)+1] <- list(comb)
    }
  }
}


## 计算每组数据的中位数
medians <- data_long %>%
  group_by(RAI) %>%
  summarise(median_value = median(RiskScore))

## 绘图
# 绘图说明：因为geom_jitter散点图无法同时设置宽度和偏移，因此散点图必须对齐X轴刻度，而我希望让中间的图位于中间，所以此处隐藏了X轴刻度，并对X轴标签做了偏移，使其居中箱线图
mycol <- c("#48af45","#337cba","#e11a0c")
ggplot(data_long, aes(x = RAI, y = RiskScore, fill = RAI)) +
  # 绘制散点图：
  geom_jitter(mapping = aes(color = RAI), width = .05, alpha = 0.5,size=0.9) +
  # 绘制箱线图，并通过position设置偏移：
  geom_boxplot(position = position_nudge(x = 0.14), color = mycol, size=0.3, width=0.1, outlier.size = 0, outlier.alpha =0, notch = T) +
  # 绘制中位数线和标记：
  stat_summary(data = medians, fun = median, geom = "segment", aes(x = as.numeric(RAI) + 0.1, xend = as.numeric(RAI) + 0.18, y = median_value, yend = median_value),
               color = "white", size = 0.5) +
  # 绘制云雨图，并通过position设置偏移：
  stat_halfeye(mapping = aes(fill= RAI), width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + 
  scale_fill_manual(values = mycol) +   #映射云雨图和箱线图的颜色
  scale_color_manual(values = mycol) +  #映射散点的颜色
  expand_limits(x = c(1, 3.5))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("RAI response") +  #设置X轴标题
  ylab("PyroScore") +  #设置Y轴标题
  # 按照特定顺序绘制x轴标签：
  scale_x_discrete(labels = c("Refractory"="Refractory","Sensitive"="Sensitive","No"="Not-recevied")) + 
  scale_y_continuous(limits = c(-0.1,6), breaks = c(0,1,2,3,4,5)) +  
  theme_classic(base_line_size = 0.3) +
  theme(#axis.ticks.x = element_line(size = 0,color = "white"),  #自定义主题
    axis.ticks = element_line(lineend = 0.3, color = "#333333"),
    legend.position = "none", #隐藏图例
    axis.title.x = element_text(size = 8, vjust = -1),  #调整X轴标题字体大小
    axis.title.y = element_text(size = 8), #调整Y轴标题字体大小
    axis.text.x = element_text(size = 8, hjust = 0.2, color = "#333333"), #设置x轴刻度字体偏移
    axis.text.y = element_text(size = 8, color = "#333333"), #设置Y轴刻度字体大小
  ) +
  geom_signif(comparisons = list(c("No","Sensitive"), c("Sensitive","Refractory"),c( "No","Refractory")),
              stat = "signif",
              test = "wilcox.test",
              step_increase = 0.1,   #每个p值之间的高度距离
              map_signif_level = T,
              margin_top = 0,   #p值与图形之间的高度距离
              vjust = 0, hjust= 0.5,
              size = 0.3,   #线条粗细
              textsize = 3.5
  ) +
  # 添加星号注释：
  labs(tag = paste0("Kruskal-Wallis, P =", format(pvalue$p.value, scientific = T, digits = 3))) +
  theme(plot.tag.position = c(0.5,0.97),
        plot.tag = element_text(size = 9, color = '#333333',face='bold',
                                vjust = 0,   #垂直距离
                                hjust = 0.4,  #水平距离
                                lineheight = 0   #行间距
        ))
ggsave('RAIresponse.pdf', height = 3, width = 3.5)

#### ROC曲线：PyroScore预测RAI反应 ####
library(pROC)
library(ggplot2)

roc <- roc(rt$RAI, rt$RiskScore,levels = c("Sensitive", "Refractory"),direction = '<')
auc(roc)
roc2 <- roc(GEO$DSEvent, GEO$RiskScore, levels = c(0,1),direction = '>')
auc(roc2)

mycol = c("#fb3e35","#337cba")
plot(roc, 
     add = T,  # 增加曲线
     print.auc = TRUE, 
     print.auc.x = 0.4, print.auc.y=0.2, #图像上输出AUC值,坐标为(x，y)     
     auc.polygon = F,   #是否显示曲线下面积 
     col = mycol[1],  # 曲线颜色  
     lwd = 2,
     auc.polygon.col= NA, # 设置ROC曲线下填充色
     max.auc.polygon = F,  # 填充整个图像     
     max.auc.polygon.col = NA,
     #grid = c(0.2, 0.2), 
     #grid.col = c("black", "black"),  # 设置间距为0.1，0.2，线条颜色     
     #print.thres = T,  #是否显示截断值 
     #print.thres.cex = 1, # 图像上输出最佳截断值，字体缩放倍数     
     smooth = F, # 绘制不平滑曲线   
     axes = T,
     #xlim = c(0, 1),   # 设置X轴范围
     #ylim = c(0, 1),   # 设置Y轴范围
     #expand = c(0,0),
     #main = "Comparison of two ROC curves", # 添加标题 
     legacy.axes= F)   # 使横轴从0到1，表示为1-特异度

plot(roc2, 
     print.auc=TRUE, 
     print.auc.x=0.4, print.auc.y=0.1,   
     auc.polygon=F,  
     col = mycol[2], 
     lwd=2,
     auc.polygon.col= NA, 
     max.auc.polygon=F,  
     max.auc.polygon.col=NA,
     smooth = F,  
     axes = T,
     legacy.axes= F)   

legend(x=0.5,y=0.2,
       paste0('AUC = ',sprintf("%.03f",roc$auc)),
       text.col ="#fb3e35", 
       text.font = 2,  #加粗
       bty = "n") 

#export:4x4

#### 百分比柱状图 ####
library(ggplot2)
library(tidyverse)
library(reshape2)

##数据读取，分组设置
Riskfile = 'totalRisk.txt'   #12
RAIfile = 'RAI.txt'   #DataSets
risk = read.table(Riskfile, header = T, sep = '\t', check.names = F)
RAI = read.table(RAIfile, header = T, sep = '\t', check.names = F)
rownames(risk) = risk[,1]
rownames(RAI) = RAI[,1]
sameSample = intersect(row.names(RAI), row.names(risk)) 
RAI = RAI[sameSample,,drop=F] 
risk = risk[sameSample,,drop=F]
rt = merge(RAI, risk[,c(1,ncol(risk)-1,ncol(risk))], by.x = 'ID', by.y = "id")   #将两个数据框按照ID列进行合并  #⭐️risk的最后2列
high  = rt[rt$Risk=='High',]
low = rt[rt$Risk=='Low',] 

##计算分布差异
prop.table(table(high$RAI))  #占比
prop.table(table(low$RAI))   
table(high$RAI)  #个数
table(low$RAI)   

##正态性检验
#rt$RAI <- as.logical(rt$RAI)
#ks.test(scale(rt$RAI),'pnorm')   #非正态(D = 0.39889, p-value < 2.2e-16)，应用卡方检验

##转化为2x2表格
mytable <- table(rt$RAI, rt$Risk)
fisher = fisher.test(mytable, alternative = "two.sided",conf.int = T, simulate.p.value = TRUE)
fisher
#chisq.test(mytable)

##构建数据
data <- data.frame(HighRisk = c(sum(high$RAI=='No'), sum(high$RAI=='Sensitive'), sum(high$RAI=="Refractory")), 
                   LowRisk = c(sum(low$RAI=='No'), sum(low$RAI=='Sensitive'), sum(low$RAI=="Refractory")),   
                   group = c('Not-received','Sensitive',"Refractory")) %>%  
  pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group <- factor(data$group, levels = c("Not-received","Sensitive", "Refractory"))

# 计算占比并增加一列
data <- data %>%
  group_by(X) %>%
  mutate(prop = count / sum(count))

# 计算文本位置并增加一列
data <- data %>%
  group_by(X) %>%
  mutate(text_position = ifelse(group == "Sensitive", 0.5 * prop + (1 - prop), 0.5 * prop)) %>%
  ungroup()


##绘图：
mycol = c("#f3d2ce", "#e98b7b", "#ec5136")
ggplot(data) +
  geom_bar(aes(X, count, fill = group), color = "snow", 
           position = position_fill(), 
           #position = position_dodge2(width = 0.7, preserve = "single"),
           stat = "identity", 
           size = 0.5, width = 0.7) +
  scale_fill_manual(values = alpha(mycol, 1)) +
  #添加百分比显示
  geom_text(aes(X, text_position, label = paste0(round(prop * 100), "%")), size = 4, color = "white", fontface = 2) +
  # 添加星号注释：
  annotate("text", x = 1.5, y = 1.08, label = paste0("Fisher.test, P = ", format(fisher$p.value,scientific = T, digits = 3),"***"), size = 4, fontface = 2)+
  scale_x_discrete(labels = c("High risk", "Low risk")) +
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
  # 设置主题：
  #theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.4)+
  theme_classic(base_line_size = 0.3) +
  theme(panel.grid = element_blank(),
        # 修改背景色：
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = '#333333'), 
        axis.text.y = element_text(color = '#333333', hjust = 1,  size = 10, lineheight = 1),
        legend.position = "bottom", #legend.justification = "left", 
        #legend.direction = "vertical",
        legend.key.size = unit(8,'pt'),
        legend.margin = margin(t=-15,r=0,b=0,l=0,unit = "pt"),
        legend.text = element_text(size = 10)) +
  # 图例顺序：
  #guides(fill=guide_legend(reverse=T)) + 
  # 图例标题：
  labs(fill="")   #ICI therapy

ggsave('RAI.distribution.pdf', height = 4, width = 3.5)
