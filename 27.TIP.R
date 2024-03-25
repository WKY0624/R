#### 直接官网获取TIP数据####
BiocManager::install("TCGAbiolinks") 
library(TCGAbiolinks)

project <- getGDCprojects()$project_id 
project <- project[grep("TCGA-", project)]
project <- "TGCA-THCA"

url = "http://biocc.hrbmu.edu.cn/TIP/RCodeAndData/pancancerData/"

dir.create ("TIP")

dir.create ("TIP/Immune activity scores/") 
dir.create("TIP/Immune cell infiltration/") 
dir.create("TIP/SignatureGenes.Expression/")

#proj="TCGA-LUSC" 
for(proj in project){
  message(proj)
  cancer = unlist(strsplit(proj, "-"))[2] 
  #Immune activity scores
  download.file(url = paste0(url,cancer,"/ssGSEA.normalized.score.txt"),
                 destfile = paste0("TIP/Immune activity scores/", proj, ".txt"))
  #Immune cell infiltration
  download.file(url = paste0(url,cancer,"/CIBER_",cancer,"_lm14_allsample_Result.txt"),
                destfile = paste0("TIP/Immune cell infiltration/", proj, ".txt"))
  #SignatureGenes. Expression
  download.file(url = paste0(url,cancer,"/SignatureGenes.Expression.txt"),
                destfile = paste0("TIP/SignatureGenes.Expression/",proj,".txt"))
}

#### 热图 ####
library(ggpubr)
rt=read.table("./TIPdownload/Immune.Activity.Scores.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件
#规范TCGA命名
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
colnames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(rt))

#读取分型文件
ClusterFile="cluster.txt"  #04
Cluster = read.table(ClusterFile,sep="\t",check.names=F,row.names=1,header=T)
Cluster = Cluster[,ncol(Cluster),drop=F]
RiskFile='totalRisk.txt'  #12
Risk = read.table(RiskFile,sep="\t",check.names=F,row.names=1,header=T)
Risk = Risk[,c(ncol(Risk)-1,ncol(Risk)),drop=F]

#合并
samSample = intersect(colnames(rt),row.names(Cluster))
data = rt[,samSample,drop=F]
Cluster = Cluster[samSample,,drop=F]
#data = t(data)
merge=cbind(t(data),Cluster,Risk) 

#将Step4合并
merge$Step4 <- rowMeans(merge[,4:20])
merge <- merge[,-c(4:20)]
merge <- merge[, c("Step1", "Step2", "Step3", "Step4", "Step5", "Step6", "Step7", "Cluster","RiskScore")]
merge <- merge[order(merge$RiskScore, decreasing = F),]

#merge.out <- cbind(id=row.names(merge),merge)
#write.table(merge.out,"merge.txt",sep="\t",row.names=F,quote=F)


## 构建绘图数据
C1 <- merge[merge$Cluster=="C1",,drop=F]
C1 <- C1[order(C1$RiskScore, decreasing = T),]
C2 <- merge[merge$Cluster=="C2",,drop=F]
C2 <- C2[order(C2$RiskScore, decreasing = T),]
C3 <- merge[merge$Cluster=="C3",,drop=F]
C3 <- C3[order(C3$RiskScore, decreasing = T),]

merge.order <- rbind(C1,C2,C3)

data_plot <- merge.order[,-c(8,9)]
data_plot <- as.data.frame(t(data_plot))

## 计算
outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(merge.order[,1:(ncol(merge.order)-2)])){
  rt1 = merge.order[,c(i,"Cluster")]
  colnames(rt1)= c("expression","Cluster")
  ksTest = kruskal.test(expression ~ Cluster, data = rt1)
  pValue = ksTest$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 3))
  if(pValue < 1){        #⭐️
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pSig$stepSig = paste0(rownames(pSig),pSig$pSig)
pNum = as.data.frame(pNum)

## 绘图前处理
# 修改行右侧文本信息，绘图矩阵的行名
#rownames(data_plot) <- pSig$stepSig
rownames(data_plot) <- pNum$pNum

# 构建列注释信息
annotation_col = data.frame(merge.order$Cluster)
colnames(annotation_col) = "PyroCluster"
rownames(annotation_col) <- rownames(merge.order)

# 构建行注释信息
annotation_row = data.frame(factor(c(1:7)))
colnames(annotation_row) = "Cycle"
rownames(annotation_row) = rownames(data_plot)
head(annotation_row)

# 自定注释信息的颜色列表
mycol = c("#e11a0c","#48af45","#337cba")

# 取渐变色：定义起始颜色和结束颜色
library(RColorBrewer)
# 定义要使用的颜色方案
color_scheme <- "Purples"
# 使用brewer.pal函数生成7个颜色
col7 <- brewer.pal(7, color_scheme)

ann_colors = list(
  PyroCluster = c(C1 = mycol[1], C2 = mycol[2], C3 = mycol[3]),
  Cycle = c("1" = col7[1], "2" = col7[2], "3" = col7[3], "4" = col7[4], "5" = col7[5], "6" = col7[6], "7" = col7[7]))

gap1 = nrow(C1)
gap2 = nrow(C2)


# 对每一行的数据进行归一化处理
#将每一行的数据减去该行数据的最小值，然后除以该行数据的极差（最大值减最小值），从而将每一行的数据都归一化到0到1之间。
#data_norm <- t(apply(data_plot, 1, function(x) (x - min(x)) / (max(x) - min(x))))

pheatmap(data_plot,    #归一化数据：data_norm
         # 去掉聚类树：
         cluster_cols = F,
         cluster_rows = F,
         # 加color bar：列注释信息；
         annotation_col = annotation_col,
         # 行注释信息：
         annotation_row = annotation_row,
         # color bar 颜色设定：
         annotation_colors = ann_colors,
         # 设置单元格颜色渐变；(100)表示分100段渐变；
         color = colorRampPalette(c(rep("#337cba",1.5), "white", rep("#e11a0c",1.5)))(100),
         # 行、列标签的字体大小
         fontsize_col = 8,
         fontsize_row = 10,
         # 是否显示行、列名
         show_colnames = F,
         gaps_col = c(gap1, gap1+gap2),
         # 设置每个单元格的宽度和高度
         cellwidth = 0.5, 
         cellheight = 20)

pdf("TIPheatmap.pdf", width = 6.5, height = 3)
dev.off()


#### Step4 雷达图 ####
#提取Step4
step4 = merge[,grep("Step4|Cluster|Risk", colnames(merge))]
#计算Cluster
outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(step4[,1:(ncol(step4)-2)])){
  rt1 = step4[,c(i,"Cluster")]
  colnames(rt1)= c("expression","Cluster")
  Test = kruskal.test(expression ~ Cluster, data = rt1)   #⭐️
  pValue = Test$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 1))
  if(pValue < 1){        #⭐️
    outTab=rbind(outTab,cbind(rt1,step=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pSig$stepSig = paste0(rownames(pSig),pSig$pSig)
pNum = as.data.frame(pNum)

#计算Risk
outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(step4[,1:(ncol(step4)-2)])){
  rt1 = step4[,c(i,"Risk")]
  colnames(rt1)= c("expression","Risk")
  Test = wilcox.test(expression ~ Risk, data = rt1)
  pValue = Test$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 1))
  if(pValue < 1){        #⭐️
    outTab=rbind(outTab,cbind(rt1,step=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pSig$stepSig = paste0(rownames(pSig),pSig$pSig)
pNum = as.data.frame(pNum)

## 绘图
#devtools::install_github("ricardo-bion/ggradar",dependencies = TRUE)
library(ggradar)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)

# 读取数据：
data <- outTab
table(data$step)

# 分组计算均值和中位数：
new_data <- data %>%
  group_by(step, Cluster) %>%
  mutate(Means = mean(expression),
         Median = median(expression))

# 筛选绘图数据,并实现长宽数据转换：
plot_data <- unique(new_data[,c(2,3,5)]) %>%     #⭐️
  pivot_wider(names_from = step, values_from = Median) %>%
  column_to_rownames("Cluster")

data2 <- plot_data
data2$group <- rownames(data2)

data3 <- data2 %>% 
  as_tibble(rownames = "Cluster") %>% 
  select(1:18)    

colnames(data3) = c("Cluster",pSig$stepSig)
#colnames(data3) <- str_wrap(colnames(data3), width = 15)   #列名的宽度超过15个字符时将被分成多行
colnames(data3) <- gsub(".recruiting", "", colnames(data3))  #删除列名中的“.recruiting”
colnames(data3) <- gsub("Step4.", "", colnames(data3))  #删除列名中的“.recruiting”

ggradar(data3,
        font.radar = 'sans',
        base.size = 2,
        
        ##------网格线的刻度范围
        #values.radar = c("0", "0.1", "0.2"),
        grid.min = min(plot_data),
        grid.mid = (max(plot_data) + min(plot_data))/2,
        grid.max = max(plot_data),
        
        ##------网格线的格式设置
        label.gridline.min = F,  #是否显示3层网格线的标签
        label.gridline.mid = F,
        label.gridline.max = F,
        grid.label.size	= 3,    #网格线label的大小
        grid.line.width	= 0.7,  #网格线的粗细 
        gridline.min.linetype = 2, #'longdash',  
        gridline.mid.linetype	= 2, #'longdash',
        gridline.max.linetype = 2, #'longdash',
        gridline.min.colour	= alpha('grey',0.7),
        gridline.mid.colour	= alpha('grey',0.7),
        gridline.max.colour	= alpha('grey',0.7),
        
        ##------辐射线
        axis.line.colour = alpha('grey',0.7), #辐射线的颜色
        axis.label.offset	= 1.1,   #外层label与圆圈的距离
        axis.label.size	= 3,    #外层label的字体大小
        
        ##------背景
        background.circle.colour = 'gray',
        background.circle.transparency	= 0.07,
        
        ##------点线
        group.colours = c("#e11a0c","#48af45","#337cba"),   #点的颜色   #48af45
        #group.colours = c("#e11a0c","#337cba"),   #点的颜色   #48af45
        group.point.size = 2,      #点的大小
        group.line.width = 0.8,      #线的粗细
        fill = F,    #是否填充 
        fill.alpha = 0.2,    #填充透明度
        
        ##------图例
        legend.title = "PyroCluster",
        legend.text.size = 1,
        legend.position = "top") +
  
  #-----主题设置
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.direction = 'horizontal',
        legend.position = 'top',
  )


ggsave('TIP.Step4.Radar.pdf', width = 5, height = 4)



