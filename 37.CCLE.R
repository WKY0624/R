library(ggplot2)
library(limma)
library(dplyr)
library(tidyverse)
library(readr)

#### 读取数据 ####
## 在CCLE下载所有细胞株转录组数据：https://sites.broadinstitute.org/ccle/datasets
# TPM
ccleExp <- read.table('./DataSets/OmicsExpressionProteinCodingGenesTPMLogp1(240223).csv',header = T, row.names = 1, sep = ',', check.names = F)

# 甲状腺癌的细胞株
thyroid <- read.table('37.cell lines in Thyroid.csv',header = T, row.names = 1, sep = ',', check.names = F, quote = "")

#### 处理ccle数据（TPM）####
same <- intersect(rownames(ccleExp),rownames(thyroid))
ccleTHCA <- ccleExp[same,,drop = F]
thyroid <- thyroid[same,,drop = F]
colnames(ccleTHCA) = gsub("\\(.*?\\)","",colnames(ccleTHCA))    #规范细胞株名称：删除第一个左括号和第一个有括号，及其之间的内容
colnames(ccleTHCA) = as.character(colnames(ccleTHCA))
ccleTHCA = as.data.frame(ccleTHCA)

#### 保存数据 ####
merge2 = cbind(thyroid,ccleTHCA)
merge2 <- cbind(id=row.names(merge2),merge2)
write.table(merge2, 'CCLE_THCA_TPM(240223).txt',col.names = T, row.names = F, quote = F, sep = '\t')

######### 读取处理CCLE+THCA之后的数据直接画 图###########
#Counts：expression有时会很大，而且部分低表达基因可能不全）
#TPM:
ccleTHCA <- read.table("CCLE_THCA_TPM(240223).txt",header = T,sep = '\t', check.names = F)

## 选择想画图的基因
filterGene <- colnames(ccleTHCA) %in% c("Cell Line","Tumor Type","Primary Disease",
                                        'CXCL8','H2BC8','IFI27','GJA1','TRAF6','PYCARD','SEZ6L2','SIGLEC15','PRDM1')  
merge <- ccleTHCA[, filterGene]
merge$Abbr <- ifelse(merge$`Primary Disease`=='Well-Differentiated Thyroid Cancer','Well-DTC',
                     ifelse(merge$`Primary Disease`=='Anaplastic Thyroid Cancer','ATC',
                            ifelse(merge$`Primary Disease`=='Medullary Thyroid Cancer','MTC','Poorly-DTC')))
merge <- merge[order(merge$Abbr, decreasing =F),]


## 构建画图数据
data <- data.frame(x=merge$`Cell Line`,y=merge$CXCL8, tumortype=merge$`Tumor Type`, Abbr=merge$Abbr)
data$tumortype = as.factor(data$tumortype)
data$x=as.factor(data$x)

#data <- data[order(data$tumortype,decreasing = F),]  #根据disease排序
data <- data[order(data$y, decreasing = F),]   #根据基因表达量排序
data$x <- factor(data$x,levels = data$x)


## 绘图
ggplot(data, aes(x=x, y=y, color=data$Abbr)) +
  geom_segment(aes(x=x, xend=x, y=0, yend=y, linetype=factor(data$tumortype,levels = c("Primary","Metastatic"))), 
               color=ifelse(data$tumortype %in% c("Primary"), "black", "black"),
               #size=ifelse(data$tumortype %in% c("Primary"), 0.3,0.3), alpha=1,
               #linetype=ifelse(data$tumortype %in% c("Primary"),1,5)
               ) +
  geom_point(#color=ifelse(data$Abbr %in% c("Well-DTC"), "#48af45",ifelse(data$Abbr %in% c('ATC'),"#e11a0c",ifelse(data$Abbr %in% c('MTC'), 'orange','#337cba'))), 
             #size=ifelse(data$tumortype %in% c("Primary"), 2.5, 2.5),
             #size=data$y,
             size=2.5,
             alpha=0.9, pch=19) + 
  scale_color_manual(values = c("#e11a0c",'orange','#337cba',"#48af45"), name="Primary disease") +
  scale_linetype_manual(values = c(1,2), name ='Tumor type')+
  coord_flip() +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.text = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(vjust = 3),
        axis.title.x =element_text(vjust = -2),
        plot.margin = margin(10,10,10,10)) +
  xlab("Cell line") +
  ylab("CXCL8 expression") +
  #scale_y_continuous(limits=c(0,7),breaks=c(0,2,4,6)) +
  theme(legend.position = 'right',
        legend.key.size = unit(10,'pt'),
        legend.title = element_text(size = 8, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 5),) +
  guides(linetype=guide_legend(order = 0, keyheight=0.8),
         color=guide_legend(order = 1, keyheight=0.8, override.aes = list(size=1.5))) 


ggsave('CXCL8.pdf', height = 3, width = 3.7)




#### 【略】或者可以使用Count数据，但是需要清洗数据矩阵####
## F注：Counts的expression有时会很大，而且部分低表达基因可能不全）
## 处理Count原数据的细胞系名称：由PR-改为ACH-（有对应的Cell lines名称）
# 读取ID对应数据
ref <- read.table("OmicsProfiles.csv", header = T, sep = ",", check.names = F)
ref <- ref[,c('ProfileID','ModelID')]
rownames(ref) <- ref[,1]

#读取Count原数据
count <- read.table("OmicsExpressionGenesExpectedCountProfile.csv", header = T, sep = ",", check.names = F)
rownames(count) <- count[,1]
colnames(count)[1] <- "ProfileID"   #修改第1列的列名

#取交集名称，再分别以此提取2个数据集
intersect <- intersect(rownames(count),rownames(ref))
count2 <- count[intersect,,drop=F]
ref2 <- ref[intersect,,drop=F]

#合并提取后的数据集，包括有Pr-和ACH-
data <- merge(ref2, count2, by.x = 'ProfileID', by.y = "ProfileID")   #将两个数据框按照ID列进行合并
#删除多余的Pr-列
data2<- subset(data,select=-c(ProfileID)) 
write.table(data2,'CCLE_counts.txt',sep="\t",quote=F, col.names=T, row.names = F)

ccleExp <- read.table('CCLE_count.txt',header = T, sep = '\t', quote = "",check.names = F)
ccleExp_mt <- as.matrix(ccleExp)
rownames(ccleExp_mt)=ccleExp_mt[,1]    #将rt的第一列（id）设置为rt的行名   
exp=ccleExp_mt[,2:ncol(ccleExp_mt)]    #rt的第2列至最后1列为表达数据
dimnames=list(rownames(exp),colnames(exp))    #分别提取exp的行名和列名，组成list
ccleExp_mt2=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
same <- intersect(rownames(ccleExp_mt2),rownames(thyroid))
ccleTHCA <- ccleExp_mt2[same,,drop = F]
thyroid <- thyroid[same,,drop = F]
colnames(ccleTHCA) = gsub("\\(.*?\\)","",colnames(ccleTHCA))    #删除第一个左括号和第一个有括号，及其之间的内容
colnames(ccleTHCA) = as.character(colnames(ccleTHCA))
ccleTHCA = as.data.frame(ccleTHCA)

