#install.packages("ggplot2")
#install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)
library(dplyr)

clusterFile = "cluster.txt"      #04     
scoreFile="totalRisk.txt"   #12
clinicalFile = '06.clinical477.txt'

#data矩阵要求：row=sample、col=每一个层级
clusterRT = read.table(clusterFile, header = T, sep = '\t', check.names = F)
scoreRT = read.table(scoreFile, header = T, sep = '\t', check.names = F)
clinicalRT = read.table(clinicalFile, header = T, sep = '\t', check.names = F)

score = scoreRT[,c('id','Risk'),drop=F]
clinical = clinicalRT[,c('Sample','ATA','TNMStage'),drop = F]

merge2 = merge(score, clinical, by.x = 'id', by.y = "Sample")   #将两个数据框按照ID列进行合并
merge3 = merge(clusterRT,merge2, by.x = 'ID', by.y = "id")

data = merge3[,c(2:5),drop=F]
data = na.omit(data) #去除所有含有NA值的行

#修改名称
data$TNMStage = ifelse(data$TNMStage =='Stage I','Stage I',
                       ifelse(data$TNMStage =='Stage II','Stage II','Stage III&IV'))
data$Cluster = ifelse(data$Cluster == 'C1','Cluster1',
                      ifelse(data$Cluster == 'C2','Cluster2','Cluster3'))

colnames(data)[colnames(data)=="Cluster"] <- "PyroCluster"   #修改列名
colnames(data)[colnames(data)=="Risk"] <- "PyroGroup"   #修改列名
colnames(data)[colnames(data)=="TNMStage"] <- "TNM"   #修改列名

#整理绘图矩阵
corLodes = to_lodes_form(data, axes = 1:ncol(data), id = "Cohort")

#绘图前设置
corLodes$stratum = factor(corLodes$stratum, 
                          levels = c('Cluster1','Cluster2','Cluster3','High','Intermediate','Low',
                                     'Stage III&IV','Stage II','Stage I'))
#配色-深色
#mycol <- c('firebrick2','mediumseagreen','steelblue3','firebrick2','mediumseagreen','steelblue3','firebrick1','mediumseagreen','steelblue3')  #改进配色

#配色-浅色
mycol <- c('#ed766d','#91cf8f','#85b0d6','#ed766d','#91cf8f','#85b0d6','#ed766d','#91cf8f','#85b0d6')

#绘图
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 1)) + 
  geom_flow(width = 1/10, aes.flow = "forward") +    #空隙   #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
  geom_stratum(alpha = 1, width = 3.5/10,color = 'white',linetype=1, lwd=0.3) +  #柱子
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum", size = 3.8, color="black") + #文字大小、颜色
  xlab("") + ylab("") + 
  theme_bw() + #去除背景色
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10)) + #去除坐标轴
  theme(panel.grid =element_blank()) +   #去除网格线
  theme(panel.border = element_blank()) +   #去除外层边框
  theme(legend.position = "none") +  #隐藏图例
  ggtitle("") + 
  guides(fill = FALSE)  

ggsave('Sankey.pdf', height = 3.5, width = 4.5)

