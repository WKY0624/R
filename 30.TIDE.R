library(limma)
library(sva)

#### TIDE数据标准化 ####
tcgaExpFile ="total.normalize.txt"  #08 #TCGA表达数据文件:477例tumor的TPM，不经log2处理

#按照网站说明，对于没有Control样本的测序数据，把表达量减去每个基因在所有样本中的平均值即可，即按照行计算平均值，再拿这一行所有表达量-该平均值。
rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]   
exp=rt[,2:ncol(rt)]   
dimnames=list(rownames(exp),colnames(exp))  
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga2 <- t(apply(tcga, 1, function(x)x-(mean(x)))) #TIDE标准化

#输出矫正后的数据
tcgaTab=rbind(ID=colnames(tcga2), tcga2)
write.table(tcgaTab, file="TIDE_total.normalize", sep="\t", quote=F, col.names=F)

#### TIDE计算 ####
#上传TIDE网站：http://tide.dfci.harvard.edu/
#获取结果文件：TIDE.csv
#Patient	No benefits	Responder	TIDE	IFNG	MSI Expr Sig	Merck18	CD274	CD8	CTL.flag	Dysfunction	Exclusion	MDSC	CAF	TAM M2
#TCGA-ET-A3BV	FALSE	FALSE	1.48	192.8	1	16.39	1.31	-4.57	FALSE	0.07	1.48	0.02	0.2	0


#### ROC曲线：PyroScore分数预测ICI反应的 ####
TIDEfile="TIDE.csv"  #29
Riskfile = 'totalRisk.txt'  #12

tide=read.table(TIDEfile, header=T, sep=",", check.names=F) 
risk=read.table(Riskfile, header = T, sep = '\t', check.names = F)
rownames(tide) = tide[,1]
rownames(risk) = risk[,1]
sameSample = intersect(row.names(tide), row.names(risk))
tide = tide[sameSample,,drop=F]
risk = risk[sameSample,,drop=F]
rt = merge(tide, risk, by.x = 'Patient', by.y = "id")   #将两个数据框按照ID列进行合并
#write.table(rt, file = 'merge.txt', sep="\t", quote=F, col.names= T, row.names = F)

#绘制ROC曲线
library(pROC)
library(ggplot2)

roc <- roc(rt$Responder, rt$RiskScore)
auc(roc)

plot(roc, 
     #add=T,  # 增加曲线
     print.auc=TRUE, 
     print.auc.x=0.4, print.auc.y=0.2, #图像上输出AUC值,坐标为(x，y)     
     auc.polygon=T,   #是否显示曲线下面积 
     col="#fb3e35",  # 曲线颜色  "#FF2E63"
     auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色     
     max.auc.polygon=FALSE,  # 填充整个图像     
     #grid=c(0.2, 0.2), 
     #grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色     
     #print.thres=T,  #是否显示截断值 
     #print.thres.cex=1, # 图像上输出最佳截断值，字体缩放倍数     
     smooth = F, # 绘制不平滑曲线     
     #main="Comparison of two ROC curves", # 添加标题 
     legacy.axes= F)   # 使横轴从0到1，表示为1-特异度

legend(x=0.5,y=0.2,
       paste0('AUC = ',sprintf("%.03f",roc$auc)),
       text.col ="#fb3e35", 
       bty = "n") 

#export:4x4


#### 百分比柱状图：fisher检验 ####
##数据处理，分组设置
TIDEfile="TIDE.csv"  #29
Riskfile = 'totalRisk.txt'  #12
tide = read.csv(TIDEfile, header=T, sep=",", check.names=F)
risk = read.table(Riskfile, header = T, sep = '\t', check.names = F)
rownames(tide) = tide[,1]
rownames(risk) = risk[,1]
sameSample = intersect(intersect(row.names(tide), row.names(risk)))   #RAI
tide=tide[sameSample,,drop=F] 
risk=risk[sameSample,,drop=F]
rt = merge(tide, risk[,c(1,13:14)], by.x = 'Patient', by.y = "id")   #将两个数据框按照ID列进行合并  #⭐️risk的最后2列
high  = rt[rt$Risk=='High',]
low = rt[rt$Risk=='Low',] 

##计算分布差异
prop.table(table(high$Responder))  #占比
prop.table(table(low$Responder))   
table(high$Responder)  #个数
table(low$Responder)   

## 正态性检验
rt$Responder <- as.logical(rt$Responder)
ks.test(scale(rt$Responder),'pnorm')   #非正态(D = 0.39889, p-value < 2.2e-16)，应用卡方检验

## 转化为2x2表格，统计分析
mytable <- table(rt$Responder, rt$Risk)
fisher = fisher.test(mytable, alternative = "two.sided",conf.int = T, simulate.p.value = TRUE)
fisher
#chisq.test(mytable)

## 构建绘图数据
data <- data.frame(HighRisk = c(sum(high$Responder=='FALSE'), sum(high$Responder=='TRUE')), 
                   LowRisk = c(sum(low$Responder=='FALSE'), sum(low$Responder=='TRUE')),    
                   group = c('Non-Responder','Responder')) %>%  
  pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group = sort(data$group, decreasing = F)

## 计算占比并增加一列
data <- data %>%
  group_by(X) %>%
  mutate(prop = count / sum(count))

## 计算文本位置并增加一列
data <- data %>%
  group_by(X) %>%
  mutate(text_position = ifelse(group == "Non-Responder", 0.5 * prop + (1 - prop), 0.5 * prop)) %>%
  ungroup()


## 绘图
mycol = c("#337cba","#de3124")
ggplot(data) +
  geom_bar(aes(X, count, fill = group), color = "#f3f4f4",
           position = position_fill(), 
           stat = "identity", 
           size = 0.5, width = 0.6) +
  scale_fill_manual(values = alpha(mycol, 0.9), labels = rev(c('Respond', 'Non-Respond'))) +
  geom_text(aes(X, text_position, label = paste0(round(prop * 100), "%")), 
            size = 4, color = "white", fontface = 2) +
  # 添加星号注释：
  annotate("text", x = 1.5, y = 1.08, label = paste0("Fisher.test, P = ", round(fisher$p.value,3)), size = 4, fontface = 2)+
  scale_x_discrete(labels = c("High risk", "Low risk")) +
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
  # 设置主题：
  #theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.4)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        #panel.background = element_rect(fill = "#f3f4f4")
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = '#333333'), 
        axis.text.y = element_text(color = '#333333', hjust = 1,  size = 10, lineheight = 1),
        legend.position = "top",
        legend.key.size = unit(10,'pt'),
        legend.text = element_text(size = 10)) +
  # 图例顺序：
  guides(fill=guide_legend(reverse=T)) + 
  # 图例标题：
  labs(fill="")  

ggsave('distribution.pdf', height = 4, width = 3)


#### 箱线图：TIDE其他结果 ####
library(limma)
library(reshape2)
library(ggpubr)
library(ggridges)

##计算组间差异
library(rstatix)  
library(ggpubr)

data_multi <- rt[,c("Patient",'Risk','TIDE','IFNG','MSI Expr Sig','Merck18','CD274','CD8','Dysfunction','Exclusion','MDSC','CAF','TAM M2')]
rownames(data_multi) = data_multi$Patient
data_multi = data_multi[,-1]
#因为MDSC、CAF、TAMM2三项的数值与其他参数不在一个数量级，统一*10以匹配其他参数标度。不影响差异分析结果。
data_multi$MDSC = data_multi$MDSC*10
data_multi$CAF = data_multi$CAF*10
data_multi$`TAM M2` = data_multi$`TAM M2`*10

## 正态性检验：经检验数据为正态分布，可以使用t.test
#shapiro.test(rt$TIDE)      #样本量小于50时选择W检验：W接近1，P值大于0.05，数据为正态分布
ks.test(scale(rt$TIDE),'pnorm')   #样本量大于50时选择D检验：D值越小，越接近0，表示样本数据越接近正态分布(D越小越好)；P<显著性水平α(0.05)，则拒绝H0（p越大越好）
ks.test(scale(low$TIDE),'pnorm')  

## 计算p值
#定义执行一次秩和检验的操作流程
wilcox.test(CD274 ~ Risk, data = data_multi, alternative = 'two.sided')
#能够执行一次那么就可以执行多次，批量执行秩和检验
stat_multi <- data.frame(Tides=colnames(data_multi)[2:ncol(data_multi)])
for (i in 2:ncol(data_multi)){
  print(i)
  stat_multi[i-1,2] <- wilcox.test(data_multi[,i] ~ Risk, data = data_multi, 
                                   alternative = 'two.sided', paired= F,
                                   exact = FALSE)[["p.value"]]}

## 计算fdr
stat_multi$fdr <- p.adjust(stat_multi$V2, method = "fdr")
stat_multi$bon <- p.adjust(stat_multi$V2, method = 'bonferroni')

## 输出结果
colnames(stat_multi) <- c("Tides","P.value","FDR","Bonferroni")

#stat_multi$P.scient <- format(stat_multi$P.value, scientific = TRUE, digits=4)
stat_multi$P.round3 <- round(stat_multi$P.value, 3)
#stat_multi$P.sig <- ifelse(stat_multi$P.value<0.001,"***",ifelse(stat_multi$P.value<0.01,"**",ifelse(stat_multi$P.value<0.05,"*","")))
#stat_multi$PP<- ifelse(stat_multi$P.value<0.001, sprintf(stat_multi$P.scient), sprintf("%.03f", stat_multi$P.value))


## 绘图数据
data2 = melt(data_multi, id.vars = "Risk")
colnames(data2) = c("Risk","Pathway","Score")

## 绘制箱线图
mycol = c("#e11a0c","#337cba")
ggboxplot(data2, x='Pathway', y= "Score", color = "Risk",
          notch = F, size = 0.4, width = 0.6, outlier.shape = 3, outlier.size = 1,
          xlab="",ylab="ssGSEA score", add = "none", 
          palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +    #转置xy
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "top", 
        legend.key.size = unit(10,'pt'),
        # 修改网格线：
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #x轴标签
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        # y轴标签：
        axis.text.y = element_text(color = "black", hjust = 1, size = 10, lineheight = 2),
        # 标题居中：
        plot.title = element_text(hjust = 0, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = stat_multi$P.round3,hjust = 0, size = unit(2.5,'mm'), color = 'black',
           x = unique(data2$Pathway), y=max(data2$Score)-0.3)   #文本注释还是在框里

ggsave("TIDEother.pdf", height = 4, width = 4.5)


