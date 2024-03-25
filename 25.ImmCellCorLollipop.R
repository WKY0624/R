#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("scales")
#install.packages("ggplot2")
#install.packages("ggtext")
#install.packages("tidyverse")
#install.packages("ggpubr")

#引用包
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)

riskFile="totalRisk.txt"      #风险文件#12
immFile="./DataSets/infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件


#读取风险的文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫细胞浸润的文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))  
immune=avereps(immune) 

#对风险文件和免疫细胞浸润文件取交集，得到交集样品
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "RiskScore"]
immune=immune[sameSample,]  

#对风险得分和免疫细胞进行相关性分析
x=as.numeric(risk)
x[x>quantile(x,0.99)]=quantile(x,0.99)
outTab=data.frame()  
for(i in colnames(immune)){    #对每个免疫细胞进行循环，计算每个免疫细胞的含量
  y=as.numeric(immune[,i])
  if(sd(y)<0.001){next}   
  corT=cor.test(x, y, method="spearman")   #进行相关性检验
  cor=corT$estimate    #得到相关系数和p值
  pvalue=corT$p.value
  if(pvalue<0.05){   
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))    #保留：免疫细胞名称，相关系数，p值
  }}
    #绘制相关性散点图
    #outFile=paste0("cor.", i, ".pdf")   #保存位置
    #outFile=gsub("/", "_", outFile)
    #df1=as.data.frame(cbind(x,y))
    #p1=ggplot(df1, aes(x, y)) +   #x轴：风险得分，y轴：免疫细胞名称
    #  xlab("Risk score") + ylab(i)+    #x轴和y轴的名称
    #  geom_point() +   #对样品打点
    #  geom_smooth(method="lm",formula = y ~ x) +   #根据打点模拟出线段
    #  theme_bw()+  
    #  stat_cor(method = 'spearman', aes(x =x, y =y)  #分析风险得分和免疫细胞之间的相关性
    #           )
    #相关性图形
    #pdf(file=outFile, width=5, height=4.7)
    #print(p1)
    #dev.off()
    #}}

#输出相关性结果，3列：免疫细胞名称、相关系数cor、p值
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)   #3列：免疫细胞名称，相关系数，相关性检验的p值


#绘图前处理####
data <- outTab
#获取算法信息
data$Algorithm = sapply(strsplit(data[,1],"_"), '[', 2)
data$Algorithm = factor(data$Algorithm,level=as.character(unique(data$Algorithm[rev(order(as.character(data$Algorithm)))])))
data$immune = sapply(strsplit(data[,1],"_"),'[', 1)
#data$immune = factor(data$immune, levels = rev(as.character(data$immune)))   
ysort = unique(data$immune)
data = data[order(data$immune, decreasing = F),]   #排序
colslabels = rep(hue_pal()(length(levels(data$Algorithm))),table(data$Algorithm))  #定义颜色
data$cor = as.numeric(data$cor)
data$pvalue = as.numeric(data$pvalue)

#绘图
mycol = c("tomato","yellow2","#337cba","#48af45","hotpink2","orange1","mediumpurple2")
ggplot(data, aes(x = immune, y = cor)) +
  geom_segment(aes(x = immune, xend = immune, y = 0, yend = cor),
               size = 0.3, linetype = 'solid' , color = 'gray30') +   #棒棒
  geom_point(data = data, aes(size = -log10(pvalue), color = Algorithm), shape = 19, alpha =0.9,
             #color = 'gray40'
             ) +    #糖  
  #scale_color_manual(values = mycol) +
  scale_size(range = c(2,5)) +
  geom_hline(yintercept = 0, linetype = 2, color = 'gray20', size = 0.3) +
  scale_x_discrete(limits=factor(sort(ysort, decreasing = T))) +
  labs(x = "Immune cell",y = 'Spearman correlation coefficient', title = '') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right", # 去掉图例
        legend.key.size = unit(8,'pt'),
        legend.title = element_text(size = 10, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 4),
        # 修改网格线：
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        #panel.grid.minor.y = element_line(linetype = "dashed"),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        # 去掉y轴刻度：
        #axis.ticks.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        # y轴标签：
        axis.text.y = element_text(color = 'black',
                                   hjust = 1, # 左对齐
                                   size = 8,
                                   lineheight = 2),
        axis.title = element_text(size = 10, color = 'black'),
  ) 

ggsave('ImmCellCorLollipop.pdf', height = 5, width = 5)

