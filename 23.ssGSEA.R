#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSEABase")


#引用包
library(GSVA)
library(limma)
library(GSEABase)

####ssGSEA分析####
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
immuneScore(expFile="total.normalize.txt",  #08
            gmtFile="./DataSets/allTCpathway.gmt",  #甲状腺相关通路标志物
            project="ssgsea")



#### 分组间的差异分析 ####
#读取ssGSEA结果文件
ssgsea = read.table("ssgsea.score.txt", header=T, sep="\t", check.names=F, row.names=1)
ssgsea = t(ssgsea)

## 1.Cluster####
cluster=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)  #Cluster数据04，或者Group数据12
rt2 = cbind(ssgsea, cluster)
data2 = melt(rt2, id.vars = "Cluster")
colnames(data2) = c("Cluster","Pathway","Score")

#转置Y轴顺序
data2$Pathway = sort(data2$Pathway, decreasing = F)   #xy转置后，Z→A。为了正序，需要把数据按照倒序排列；之后在x的label再把名称倒序排列

##统计分析
#单个通路
#aa <- aov(`Signaling by WNT`~ Cluster, data=rt2)
#summary(aa)[[1]][["Pr(>F)"]]
#kruskal.test(`Signaling by WNT`~Cluster,data=rt2)

##Kruskal循环####
kruskal_data <- data.frame(Pathway = colnames(rt2)[1:ncol(rt2)-1])
for (i in 1:24){   #⭐️
  print(i)
  kruskal_data[i,2] <- kruskal.test(rt2[,i] ~ Cluster, data = rt2)[["p.value"]]}
##计算fdr
kruskal_data$fdr <- p.adjust(kruskal_data$V2, method = "fdr")
kruskal_data$BH <- p.adjust(kruskal_data$V2, method = "BH")
##输出结果
colnames(kruskal_data) <- c("Pathway","p.value","adj.p.value")    #,"logFC")
kruskal_data$p.value3 = ifelse(kruskal_data$p.value<0.001,"<0.001",sprintf("%.03f", kruskal_data$p.value))
kruskal_data$p.signif = ifelse(kruskal_data$adj.p.value<0.001,"***",
                               ifelse(kruskal_data$adj.p.value<0.01,"**",
                                      ifelse(kruskal_data$adj.p.value<0.05,"*","NS")))
kruskal_data$p.sci = ifelse(kruskal_data$p.value>0.001,sprintf("%.03f", kruskal_data$p.value),
                            format(kruskal_data$p.value, scientific = TRUE, digits = 1))



##⏹ANOVA循环#####
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

###绘图####
#设置X轴标签颜色
mycol = c("#e11a0c","#48af45","#337cba")
data2$Cluster <- factor(data2$Cluster, levels = c("C3","C1","C2"))  #调整山脊图的图层顺序

# 绘制山脊图
library(ggridges)
ggplot(data2, aes(y = Pathway, x = Score, fill = Cluster, color = Cluster)) +
  geom_density_ridges(scale = 1.2, alpha = 0.7, size = 0.1, #from = -10, to = 10,
                      jittered_points = F , rel_min_height = 0.01,
                      #point_shape = "|", 
                      #point_size = 1,
                      position = position_points_jitter(height = 0.2, width = 0.1), #color = "snow"
                      #quantile_lines = TRUE, 
                      #quantiles = 2,
                      #stat = "binline",bins = 15, draw_baseline = F, 
                      #alpha = 0.7, scale = 1, size = 1, fill =NA
  ) +
  scale_fill_manual(values = c("C1"=alpha(mycol[1],0.7),"C2"=mycol[2],"C3"=mycol[3])) +
  scale_color_manual(values = c("C1"= alpha(mycol[1],0.7),"C2"=mycol[2],"C3"=mycol[3])) +
  scale_y_discrete(limits = rev(levels(data2$Pathway))) +
  xlim(0,1.1) + 
  labs(title = "", x = "Score", y = "") +
  theme_ridges() +
  theme(legend.position = "top") +
  annotate(geom = 'text', label = kruskal_data$p.sci, hjust = 0, size = unit(3,'mm'), color = 'black', #fontface = faceCell,
           y = unique(data2$Pathway), x=1.01)   #文本注释还是在框里


ggsave("Cluster.ssGSEA.pdf",height = 9, width = 8)





#### 2.Group ####
risk=read.table("totalRisk.txt", header=T, sep="\t", check.names=F, row.names=1)  #12

#合并数据
sameSample=intersect(row.names(ssgsea),row.names(risk))
ssgsea = ssgsea[sameSample,,drop=F]
risk = risk[sameSample,,drop=F]
rt=cbind(ssgsea,risk[,c("RiskScore","Risk")])   #⭐️
rt=rt[,-(ncol(rt)-1)]
data=melt(rt,id.vars=c("Risk"))    #data有3列：'Risk'、variable、value   #⭐️
colnames(data)=c("PyroGroup","Pathway","Score")  

#转置Y轴顺序
data$Pathway = sort(data$Pathway, decreasing = T)   #xy转置后，Z→A。为了正序，需要把数据按照倒序排列；之后在x的label再把名称倒序排列

##【单独计算组间差异】
library(rstatix)  
library(ggpubr)
##计算p值
#定义执行一次秩和检验的操作流程
#wilcox_test(aDCs ~ Risk, data = rt1, alternative = 'two.sided')
#能够执行一次那么就可以执行多次，批量执行秩和检验
wilcox_data <- data.frame(Pathways=colnames(rt)[1:ncol(rt)-1])
for (i in 1:24){
  print(i)
  wilcox_data[i,2] <- wilcox.test(rt[,i] ~ Risk, data = rt, 
                                  alternative = 'two.sided', 
                                  exact = FALSE)[["p.value"]]
}

##计算fdr
wilcox_data$fdr <- p.adjust(wilcox_data$V2, method = "fdr")

##输出结果
colnames(wilcox_data) <- c("Pathway","p.value","adj.p.value")    #,"logFC")
wilcox_data$p.value3 = ifelse(wilcox_data$p.value<0.001,"<0.001",sprintf("%.03f", wilcox_data$p.value))
wilcox_data$p.signif = ifelse(wilcox_data$adj.p.value<0.001," ***",
                              ifelse(wilcox_data$adj.p.value<0.01," **",
                                     ifelse(wilcox_data$adj.p.value<0.05," *","")))
wilcox_data$p.sci = ifelse(wilcox_data$p.value>0.001,sprintf("%.03f", wilcox_data$p.value),
                           format(wilcox_data$p.value, scientific = TRUE, digits = 3))
wilcox_data$p.plus = ifelse(wilcox_data$p.value>0.001,paste0(sprintf("%.03f", wilcox_data$p.value),wilcox_data$p.signif),
                            paste0(format(wilcox_data$p.value, scientific = TRUE, digits = 3),wilcox_data$p.signif))

#设置X轴标签颜色
mycol = c("#e11a0c","#337cba")
data$PyroGroup=factor(data$PyroGroup, levels=c("High","Low"))#调整山脊图的图层顺序

# 绘制山脊图
ggplot(data, aes(y = Pathway, x = Score, fill = PyroGroup, color = PyroGroup)) +
  geom_density_ridges(scale = 1.5, size = 0.2, from = 0, to = 1,
                      jittered_points = F, rel_min_height = 0.05,
                      point_shape = "|", point_size = 2,
                      position = position_points_jitter(height = 0)
  ) +
  scale_fill_manual(values = c("High"=alpha(mycol[2],0.5),"Low"=alpha(mycol[1],0.5))) +
  scale_color_manual(values = c("High"= alpha(mycol[2],0.5),"Low"=alpha(mycol[1],0.3))) +
  scale_y_discrete(limits = rev(levels(data$Pathway))) +
  #xlim(-10,10) + 
  labs(title = "", x = "ssGSEA score",y = "") +
  theme_ridges() +
  theme(legend.position = "top") +
  annotate(geom = 'text', label = rev(wilcox_data$p.plus), hjust = 0, size = unit(3,'mm'), color = 'black', #fontface = faceCell,
           y = unique(data$Pathway), x=1)   #文本注释还是在框里

ggsave("Group.ssGSEA.pdf", height = 5, width = 8)
