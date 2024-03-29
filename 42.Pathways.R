library(GSVA)
library(limma)
library(GSEABase)

#### ssGSEA计算通路Score ####
#定义ssGSEA的函数
immuneScore=function(expFile=null, gmtFile=null, project=null){
  ##读取表达输入文件，并对输入文件处理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  mat=avereps(mat)
  mat=mat[rowMeans(mat)>0,]
  
  ##读取数据集文件
  geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
  
  ##ssgsea分析
  ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
    #method=用于估计每个样本的基因集富集分数的方法。默认设置为:gsva，其他选项为:ssgsea、zscore\或plage。
    #后两种方法将第一个表达谱标准化为样本的z分数，在zscore的情况下，它将它们组合在一起，作为它们的总和除以基因集大小的平方根，而在plage的情况下，它们用于计算基因集合中基因的奇异值分解（SVD）
    #kcdf=在非参数估计样本间表达式级别的累积分布函数时使用。
    #默认情况下，kcdf=“Gaussian”，当输入表达值连续时，如对数尺度的微阵列荧光单位、RNA seq log CPMs、log RPKMs或log TPM时适用。
    #当输入表达式值是整数计数时，例如从RNA-seq实验得出的值，则该参数应设置为kcdf=“Poisson”。
  
  ##定义ssGSEA score矫正函数
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  ##对ssGSEA score进行矫正
  #ssgseaOut = ssgseaScore
  ssgseaOut=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
  write.table(ssgseaOut, file=paste0(project, ".score.txt"), sep="\t", quote=F, col.names=F)
}

#调用函数,进行ssGSEA分析
immuneScore(
  ## ----- 表达矩阵
  expFile = "total.normalize.txt",   #08 TCGA组织表达矩阵
  #expFile = "CCLE_THCA_TPM(240223).txt",    #CCLE细胞株表达矩阵
  
  ## ----- 标志物          
  gmtFile = "22.h.all.v2023.1.Hs.symbols.gmt",  #50个Cancer Hallmark
  #gmtFile = "23.allTCpathway",    #甲状腺相关标志物
  #gmtFile = "42.18PCD.gmt", #18种细胞程序性死亡   
  #gmtFile = "42.Proliferation.gmt",  #细胞增殖生长相关标志物
  #gmtFile = "42.Metastasis.gmt",  #侵袭转移相关标志物


  ## ----- 保存文件
  project = "TCGA.Hall"  #18PCD #thyroid  #Proliferation #Metastasis
  )

#### TCGA山脊图 ####
library(limma)
library(reshape2)
library(ggpubr)

## 读取ssGSEA结果文件
ssgsea = read.table("TCGA.18PCD.score.txt", header=T, sep="\t", check.names=F,row.names = 1)
ssgsea = t(ssgsea)

## 分组为oe、ue、me
group = read.table("sample_group.txt", header=T, sep="\t", check.names=F, row.names=1) 
#行名是你的两个矩阵的行名（行索引）
common_rows <- intersect(rownames(ssgsea), rownames(group))
#根据相同的行名筛选出矩阵的子集
subset_matrix1 <- ssgsea[common_rows, ]
SEZ6L2 <- group[common_rows, ]
# 将两个子集矩阵合并
merged <- cbind(subset_matrix1, SEZ6L2)
merged = as.data.frame(merged)
# 将merged矩阵的第1到17列转换为数字格式
merged[, 1:17] <- apply(merged[, 1:17], 2, as.numeric)
# 宽数据 转化 长数据
data = melt(merged, id.vars = "SEZ6L2")
colnames(data) = c("SEZ6L2","Pathway","Score")
# 转置Y轴顺序
data$Pathway = sort(data$Pathway, decreasing = F)   #xy转置后，Z→A。为了正序，需要把数据按照倒序排列；之后在x的label再把名称倒序排列


## Mann-Whitney U 检验循环
# 创建一个空的数据框来存储结果
Wilcox_data <- data.frame(Pathway = colnames(merged2)[1:(ncol(merged2) - 1)])
# 针对每个Pathway执行 Mann-Whitney U 检验
for (i in 1:17) {
  print(i)
  p_value <- wilcox.test(merged2[, i] ~ SEZ6L2, data = merged2)$p.value
  Wilcox_data[i, 2] <- p_value
}
# 计算fdr
Wilcox_data$fdr <- p.adjust(Wilcox_data$V2, method = "fdr")
Wilcox_data$BH <- p.adjust(Wilcox_data$V2, method = "BH")
# 输出结果
colnames(Wilcox_data) <- c("Pathway","p.value","adj.p.value")    #,"logFC")
Wilcox_data$p.value3 = ifelse(Wilcox_data$p.value<0.001,"<0.001",sprintf("%.03f", Wilcox_data$p.value))
Wilcox_data$p.signif = ifelse(Wilcox_data$adj.p.value<0.001,"***",
                              ifelse(Wilcox_data$adj.p.value<0.01,"**",
                                     ifelse(Wilcox_data$adj.p.value<0.05,"*","NS")))
Wilcox_data$p.sci = ifelse(Wilcox_data$p.value>0.001,sprintf("%.03f", Wilcox_data$p.value),
                           format(Wilcox_data$p.value, scientific = TRUE, digits = 1))

## Kruskal循环
kruskal_data <- data.frame(Pathway = colnames(merged)[1:ncol(merged)-1])
for (i in 1:17){
  print(i)
  kruskal_data[i,2] <- kruskal.test(merged[,i] ~ SEZ6L2, data = merged)[["p.value"]]}
# 计算fdr
kruskal_data$fdr <- p.adjust(kruskal_data$V2, method = "fdr")
kruskal_data$BH <- p.adjust(kruskal_data$V2, method = "BH")
# 输出结果
colnames(kruskal_data) <- c("Pathway","p.value","adj.p.value")    #,"logFC")
kruskal_data$p.value3 = ifelse(kruskal_data$p.value<0.001,"<0.001",sprintf("%.03f", kruskal_data$p.value))
kruskal_data$p.signif = ifelse(kruskal_data$adj.p.value<0.001,"***",
                               ifelse(kruskal_data$adj.p.value<0.01,"**",
                                      ifelse(kruskal_data$adj.p.value<0.05,"*","NS")))
kruskal_data$p.sci = ifelse(kruskal_data$p.value>0.001,sprintf("%.03f", kruskal_data$p.value),
                            format(kruskal_data$p.value, scientific = TRUE, digits = 1))


## ANOVA循环
# 将数据从宽格式转换为长格式
library(tidyr)
data_long <- rt2 %>%
  pivot_longer(cols = colnames(rt2)[1:26], names_to = "Pathway", values_to = "Score")
# 创建一个空的结果数据框来存储ANOVA结果
anova_results <- data.frame(Pathway = character(0), P_value = numeric(0))
# 遍历每条路径
for (p in unique(data_long$Pathway)) {
  pathway_data <- data_long %>% filter(Pathway == p)
  anova_result <- aov(Score ~ Cluster, data = pathway_data)
  p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
  anova_results <- rbind(anova_results, data.frame(Pathway = p, P_value = p_value))
}
# 显示每条路径的ANOVA结果
print(anova_results)


## 绘图
# 设置X轴标签颜色
mycol = c("#e11a0c","#48af45","#337cba")

data2$SEZ6L2 <- factor(data2$SEZ6L2, levels = c("ue","oe"))  #调整山脊图的图层顺序
data2$Score = as.numeric(data2$Score)
data2$Pathway = as.factor(data2$Pathway)

# 绘制山脊图
library(ggridges)
ggplot(data2, aes(y = Pathway, x = Score, fill = SEZ6L2)) +
  geom_density_ridges(scale = 1.2, alpha = 0.6, size = 0.2, #from = -10, to = 10,
                      jittered_points = F , #rel_min_height = 0.01,
  ) +
  scale_y_discrete(limits = rev(levels(data2$Pathway))) +
  xlim(0,1.1) +
  scale_fill_manual(values = c("oe"=alpha(mycol[1],0.7),"ue"=mycol[3])) +
  #scale_color_manual(values = c("oe"= alpha(mycol[1],0.7),"me"=mycol[2],"ue"=mycol[3])) +
  labs(title = "", x = "Score", y = "") +
  theme_ridges() +
  theme(legend.position = "top") +
  annotate(geom = 'text', label = Wilcox_data$p.sci, hjust = 0, size = unit(3,'mm'), color = 'black', 
           y = unique(data2$Pathway), x = 1.01) 

ggsave("18PCD_wilcox.pdf",height = 7, width = 8)


#### CCLE相关性趋势图 ####
rt = read.table("CCLE_THCA_TPM(240223).txt", header = T, check.names = F, sep = "\t", row.names = 1)
SEZ6L2 <- rt["SEZ6L2", ]
Score = read.table("CCLE.18PCD.score.txt",header = T, check.names = F, sep = "\t", row.names = 1)
data = rbind(SEZ6L2, Score)
data = as.data.frame(t(data))

#计算相关系数
SpearmanR <- cor(data$SEZ6L2, data$Entosis,method="spearman",use="complete.obs")
SpearmanP <- cor.test(data$SEZ6L2, data$Entosis,method="spearman", use="complete.obs")
p.sci = ifelse(SpearmanP$p.value > 0.001, sprintf("%.03f", SpearmanP$p.value),
               format(SpearmanP$p.value,scientific = TRUE, digits = 1))


## 绘图
mycol2=c("#e11a0c","#48af45","#337cba","orange")

ggplot(data, aes(x=data$Entosis,y=SEZ6L2)) +
  geom_point(color = "#333333", alpha=0.7, pch=20, cex=2) +
  geom_text(aes(label = rownames(data)), size = 1.5, vjust = 2, hjust = 0.5, alpha = 1) +
  geom_smooth(method=lm , formula = y ~ x, 
              #置信区间：
              color=mycol2[2], 
              fill=mycol2[2], alpha = 0.2, se=T) +
  
  theme_test(base_line_size = 0.3)+
  ylab("SEZ6L2") + 
  xlab("Entosis") +
  theme(
    # 去除网格线：
    panel.grid = element_blank(),
    # 修改坐标轴标签
    axis.title = element_text(size = 10),
    axis.title.y = element_text(vjust = 0),
    axis.text = element_text(color = "gray30",size = 9)
  ) + 
  
  annotate("text", x = median(data$Entosis), y = max(data$SEZ6L2), 
           fontface = 1, hjust = 1, size = 3,
           label = paste0("rho = ",round(SpearmanR,3), ", ","P = ", p.sci))

ggsave("Entosis.pdf", height = 3, width = 3.2)
