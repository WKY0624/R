rm(list = ls())
#install.packages("oncoPredict")
#BiocManager::install("GenomicFeatures") 
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("rtracklayer")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

#### 药敏分析 ####
dir='./DataSets/32.Drug/Training Data/'
dir(dir)

#表达数据
myexp = read.table("total.normalize.txt", sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = F)  #08
myexp = as.matrix(myexp)
#myexp = log2(myexp+1)
dim(myexp)

#分组数据
group = read.table("totalRisk.txt", sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)  #12
group = group[,ncol(group),drop=F]
head(group)

#选择药敏数据库
exp = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
exp[1:4,1:4]
dim(exp)    #17419个基因，805个样本
#drug = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
drug = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
#drug = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
drug <- exp(drug) #下载到的数据是被log转换过的，用这句代码逆转回去
#ggboxplot(melt(drug[ , 1:4]), "Var2", "value") 
drug[1:4,1:4]
dim(drug)   #805个样本，198个药物
identical(rownames(drug),colnames(exp))

#运行主代码
calcPhenotype(trainingExprData = exp,
                trainingPtype = drug,
                testExprData = myexp,
                batchCorrect = 'eb',  #   "eb" for array,standardize  for rnaseq  
              #batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData',
              ##Options are 'homogenizeData' and 'rawData'
              #homogenizeData is likely better if there is ComBat batch correction, raw data was used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data.
              )

calcPhenotype(trainingExprData = exp,
              trainingPtype = drug,
              testExprData = myexp,
              batchCorrect = "eb",
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              selection = 1,
              printOutput = TRUE,
              pcr = FALSE,
              removeLowVaringGenesFrom = "homogenizeData",
              report_pc = FALSE,
              cc = FALSE,
              percent = 80,
              rsq = FALSE)

#https://mirrors.sjtug.sjtu.edu.cn/cran/web/packages/oncoPredict/vignettes/calcPhenotype.html

#查看结果
library(data.table)
testPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', row.names = 1,check.names = F)
testPtype[1:4, 1:4]
dim(testPtype)
identical(colnames(testPtype),colnames(drug))  #198种药物IC50的预测结果就在这个表格里啦。


#初步画图看趋势
#可以画个图比较一下预测结果与真实数据，可以肉眼计算相关性系数基本是1，也就知道了计算的结果确实是IC50值，而且计算的还挺准。
#（当然准啦，因为数据是从矩阵里面截取的）

library(stringr)
p = str_remove(rownames(testPtype),"test")   #移除行名中的"test"字符
a = t(rbind(drug[p,],testPtype))
a = a[,c(1,5,2,6,3,7,4,8)]
par(mfrow = c(2,2))
plot(a[,1],a[,2])
plot(a[,3],a[,4])
plot(a[,5],a[,6])
plot(a[,7],a[,8])

#前面的箱线图，我们展现的是某个药物的八百多细胞系的IC50，这样可以看得出来有一些药物在很多癌症细胞系的表现就是废材，#比如 "Cisplatin_1005"和 "Cytarabine_1006"，
#当然了，因为我仅仅是展现了4个药物，所以说它们是废材仅仅是相当于 "Camptothecin_1003" 和"Vinblastine_1004"来说。
#还有另外一个展现方式，就是看针对具体的细胞系来说，那些药物有奇效那些药物是打酱油。
ggboxplot(melt(drug[ 1:4 ,]), "Var1", "value") 

#因为每个细胞系的箱线图里面都是约200个药物，所以这样的可视化看不出来具体 的药物表现，并没有太大的意义。我们应该是直接看top药物即可：
round(apply(drug[ 1:4 ,], 1, function(x){
  return(c(
    head(sort(x)),
    tail(sort(x))
  ))
}),2)

#可以看到， 每个细胞系都是有自己的特异性药物和废物药物，IC50接近于0的就是神药，那些大几千的就是辣鸡药物啦。
#但是，我们可能是更想看到的是药物名字啦！
#前面的6个药物是各自细胞系的神药，后面的6个是废物药物啦。
apply(drug[ 1:4 ,], 1, function(x){ 
  names(x)=gsub('_[0-9]*','',colnames(drug))
  return(c(
    names(head(sort(x))),
    names(tail(sort(x)))
  ))
})


#### 绘图 ####
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)

Drugfile = "./calcPhenotype_Output_CTRP2/DrugPredictions.csv"
Riskfile = 'totalRisk.txt'  #12
RAIfile = "RAI.txt"

## 数据读取
group = read.table(Riskfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
group = group[,ncol(group),drop=F]

RAI = read.table(RAIfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F )
sameSample = intersect(row.names(group), row.names(RAI))
RAIgroup = RAI[sameSample,,drop=F]
group = group[sameSample,,drop=F]
RAIgroup = RAIgroup[RAIgroup$RAI %in% c("Sensitive","Refractory"),,drop=F]   #"No","Sensitive","Refractory" ##【N+S+R=477，S+R=308，R=40】
RAIgroup$SampleID = rownames(RAIgroup)

drug <- read.csv(Drugfile, row.names = 1, check.names = F)   # 药敏结果
sameSample2 = intersect(row.names(RAIgroup), row.names(drug))
drug = drug[sameSample2,,drop=F]
group = group[sameSample2,,drop=F]

## 数据处理
# 风险文件和药物敏感性结果合并
sameSample=intersect(row.names(group), row.names(drug))
group = group[sameSample,"Risk",drop=F]
group$SampleID = rownames(group)
drug = drug[sameSample,,drop=F]
# 将drug和group矩阵进行合并，并整理为长格式
combined.data <- cbind(data.frame(SampleID = rownames(drug)), as.data.frame(drug))
combined.data <- pivot_longer(combined.data, -SampleID, names_to = "Drug", values_to = "IC50")
combined.data <- left_join(combined.data, group, by = "SampleID")

#### 1. 单药的差异分析 ####
##1.1 筛选出"单药"的数据 ####
drug.data <- combined.data %>% filter(grepl("AGK-2", Drug, ignore.case = TRUE))
#drug.data <- combined.data %>% filter(Drug == "selumetinib:GDC-0941 (4:1 mol/mol)")
#分子靶向药 #sorafenib   #lenvatinib   #pazopanib  #sunitinib  #axitinib  #Apatinib
#化疗药     #doxorubicin_133/1386  #paclitaxel  #carboplatin   #docetaxel
#预测的药物 #GDSC1：  #doramapimod  #elesclomol  #olaparib_1495
            #CTRP2：  #zebularine  #temozolomide  #azacitidine

## 处理离群值【可选】：定义函数，将离群值替换为NA
# 离群值定义为为小于Q1(25%位置)－1.5IQR或大于Q3(75%位置)+1.5IQR的值。
re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)  
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}
#删除含有outliers(NA)的行
drug.data  <- drug.data %>%
  group_by(Risk)%>%
  mutate(IC50 = re_outlier(IC50))
drug.data <- drug.data[complete.cases(drug.data),]

##1.2 统计分析 ####
drug.diff <- drug.data %>%
  group_by(Risk) %>%
  summarize(Mean = mean(IC50), SD = sd(IC50),
            Median = median(IC50), Q1 = quantile(IC50, 0.25), Q3 = quantile(IC50, 0.75), Count = n())
print(drug.diff)

## 使用wilcox.test计算差异
stat.data <- wilcox.test(drug.data$IC50[drug.data$Risk == "Low"],
                         drug.data$IC50[drug.data$Risk == "High"],
                         paired = F, alternative = 'two.sided')

## 1.3绘图 #### 
library(ggplot2)
mycol <- c("#e11a0c", "#337cba")

ggplot(drug.data, aes(x = Risk, y = IC50, fill = Risk))+
  scale_fill_manual(values = mycol) + 
  geom_violin(aes(color = Risk), width = 0.5, size = 1, alpha=0.1) +
  #geom_jitter(mapping = aes(color = Risk), shape = 16, width = 0.2, alpha = 0.6, size = 1.5)+ 
  scale_color_manual(values = mycol)+  #映射散点的颜色
  geom_boxplot(notch = F, outlier.size = -1, width = 0.5,lwd=0.5, color="black", alpha = 0.6)+ 
  #geom_point(shape = 21, size = 2, alpha = 1,color='black',position = position_jitterdodge())+ 
  #ylim(14,16) +
  ylab("IC50") +
  xlab("PyroGroup")  +
  theme_classic() +
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    axis.ticks.y = element_line(size=0.5, color="#333333"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(0.3,"cm"),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank()) +
  labs(tag = paste0(unique(drug.data$Drug), "\nP ", ifelse(stat.data$p.value<0.001, paste0("= ", format(stat.data$p.value, scientific = TRUE, digits = 3)), paste0("= ",round(stat.data$p.value,3)))))+
  theme(plot.tag.position = c(0.55,0.9),
        plot.tag = element_text(size = 10, color = 'black',
                                vjust = 0,   #垂直距离
                                hjust = 0.5,  #水平距离
                                lineheight = 1.48   #行间距
        )) 

##保存
unique(drug.data$Drug)
ggsave('CTRP2_AGK-2_outlier.pdf', height = 3.2, width = 3.5)

#### 2.棒棒糖相关性图：High vs Low组间差异最大的TOP药物 ####
## 2.1查看TOP药物 ####
# 创建一个空的数据框来存储每种药物的差异检验结果
results <- data.frame(Drug = character(0), p_value = numeric(0),
                      Median.High = numeric(0), Median.Low = numeric(0),
                      Mean.High = numeric(0), Mean.Low = numeric(0))

# 循环遍历每种药物
for (drug.name in unique(combined.data$Drug)) {
  # 提取当前药物的IC50数据
  drug.data <- combined.data %>% filter(Drug == drug.name)
  
  #离群值定义为为小于Q1(25%位置)－1.5IQR或大于Q3(75%位置)+1.5IQR的值。
  re_outlier<-function(x,na.rm = TRUE,...){
    qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   #⚠️⭐️  (0.25,0.75)
    h<-1.5*IQR(x,na.rm = na.rm)
    y<-x
    y[x<(qnt[1]-h)]<-NA
    y[x>(qnt[2]+h)]<-NA
    y}
  #删除含有outliers(NA)的行
  drug.data  <- drug.data %>%
    group_by(Risk)%>%
    mutate(IC50 = re_outlier(IC50))
  drug.data <- drug.data[complete.cases(drug.data),]
  
  #统计描述
  drug.diff <- drug.data %>%
    group_by(Risk) %>%
    summarize(Mean = mean(IC50), SD = sd(IC50),
              Median = median(IC50), Q1 = quantile(IC50, 0.25), Q3 = quantile(IC50, 0.75), Count = n())
  
  # 提取Low和High组的IC50数据
  drug.low <- drug.data$IC50[drug.data$Risk == "Low"]
  drug.high <- drug.data$IC50[drug.data$Risk == "High"]
  
  # 使用wilcox.test计算差异
  wilcox.result <- wilcox.test(drug.high,drug.low)
  
  
  # 将差异检验结果添加到数据框中
  results <- rbind(results, data.frame(Drug = drug.name, 
                                       Median.High = drug.diff$Median[drug.diff$Risk == "High"],
                                       Median.Low = drug.diff$Median[drug.diff$Risk == "Low"],
                                       Mean.High = drug.diff$Mean[drug.diff$Risk == "High"],
                                       Mean.Low = drug.diff$Mean[drug.diff$Risk == "Low"],
                                       p_value = wilcox.result$p.value))
}

# 根据p值从小到大排序
results <- results[order(results$p_value), ]
results <- results[order(results$Median.High), ]

# 提取差异最大的前10/30种药物
top_10_drugs <- head(results, 10)
top_30_drugs <- head(results, 30)



## 2.2计算PyroScore和IC50的相关性 ####
#PyroScore
score = read.table(Riskfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
score = score[,ncol(score)-1,drop=F]
score = score[sameSample,,drop=F]
score$ID = rownames(score)
drug2 <- drug
drug2$ID = rownames(drug2)
drug.score = merge(score, drug2, by.x = "ID", by.y = "ID")
rownames(drug.score) = drug.score$ID
drug.score = subset(drug.score,select = -c(ID)) 
score1 = score[,ncol(score)-1,drop=F]

#Spearman相关性分析
Spearman=data.frame()
for(agents in colnames(drug)){
  for(RiskScore in colnames(score1)){
    x=as.numeric(drug[,agents])
    y=as.numeric(score1[,RiskScore])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    Spearman = rbind(Spearman,cbind(Agents = agents, cor, text, pvalue))
  }
}

## 2.3提取差异最大的前10种药物 ####
data = merge(Spearman, results, by.x = "Agents", by.y = "Drug")
data <- data[order(data$p_value), ]
data30 <- head(data, 32)   #⭐️

data30$Agents = factor(data30$Agents, levels = data30$Agents)
data30$cor = as.numeric(data30$cor)
data30$p_value = as.numeric(data30$p_value)
data30$size <- with(data30, 2 + (5 - 2) * ((Mean.Low - min(Mean.Low)) / (max(Mean.Low) - min(Mean.Low))))
data30 = data30[order(data30$p_value,decreasing = T), ]
data30$log10P = as.numeric(-log10(data30$p_value))

## 2.4绘图数据 ####
#宽数据转为长数据
data30_long <- data30 %>%
  pivot_longer(c(Mean.High,Mean.Low),
               names_to = "group",
               values_to = "IC50")

#为了绘制横线，需要找出每组中log2FC的最小值和最大值：
data30_long_mm <- data30_long %>%
  group_by(Agents) %>%
  mutate(x_min = min(IC50),
         x_max = max(IC50))


## 2.5绘图 ####
ggplot(data30, aes(x=reorder(Agents, log10P), y=cor)) +
  geom_segment(aes(x=reorder(Agents, log10P), xend=Agents, y=0, yend=cor,  
                   #linetype=factor(data$tumortype,levels = c("Primary","Metastatic"))
                  ), 
               alpha = 1, linewidth = 0.3
               ) +
  geom_point(aes(size = data30$Mean.Low, color = data30$log10P),pch = 20) +
  geom_hline(yintercept = 0, linetype = 2, color = 'gray20', linewidth = 0.3) +
  scale_size(range = c(2,5)) +
  scale_colour_distiller(palette ='Reds', direction = 1, name = "-log10(pvalue)") +
  #scale_colour_gradient2(low = alpha(mycol[3],1), mid = alpha(mycol[1],0.3), high = alpha(mycol[1],1),   #"dodgerblue2"  "#fb3e35"
  #                       #breaks = c(min(data30$p_value), max(data30$p_value)), 
  #                       #labels = c(min(data30$p_value), max(data30$p_value)),   #图例标签
  #                       limits = c(min(data30$log10P), max(data30$log10P)),
  #                       breaks = c(5,6,7,8,9,10),
  #                       midpoint = 5,
  #                       name = "-log10(pvalue)") +  #单元格颜色
  ylim(-0.41,0.35) +
  coord_flip() +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.text = element_text(size = 8, color = '#333333'),
        axis.text.x = element_text(vjust = 0, hjust = 0.5),   #angle = 60, 
  ) +
  xlab("") +
  ylab("Spearman correlation") +
  #scale_y_continuous(limits=c(0,7),breaks=c(0,2,4,6)) +
  theme(legend.position = 'right',
        legend.key.size = unit(10,'pt'),
        legend.title = element_text(size = 8, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 5),)  +
  guides(size = guide_legend(title = "IC50", order = 0, keyheight=1)) 

ggsave('GSDC1_top30.pdf', width = 5, height = 5)

