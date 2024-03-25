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
    #method=用于估计每个样本的基因集富集分数的方法。默认设置为gsva，其他选项为ssgsea、zscore或plage。
    #后两种方法将第一个表达谱标准化为样本的z分数，在zscore的情况下，它将它们组合在一起，作为它们的总和除以基因集大小的平方根，而在plage的情况下，它们用于计算基因集合中基因的奇异值分解（SVD）
    #kcdf=在非参数估计样本间表达式级别的累积分布函数时使用。
    #默认情况下，kcdf=“Gaussian”，当输入表达值连续时，如对数尺度的微阵列荧光单位、RNA seq log CPMs、log RPKMs或log TPM时适用。
    #当输入表达式值是整数计数时，例如从RNA-seq实验得出的值，则该参数应设置为kcdf=“Poisson”。
  
  ##定义ssGSEA score矫正函数
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  ##对ssGSEA score进行矫正
  ssgseaOut=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
  write.table(ssgseaOut, file=paste0(project, ".score.txt"), sep="\t", quote=F, col.names=F)
}

#调用函数,进行ssGSEA分析
immuneScore(expFile="total.normalize.txt", #08
            gmtFile="./DataSets/immune.gmt", #免疫标记基因集
            project="total")


#### 绘图 ####
#引用包
library(limma)
library(reshape2)
library(ggpubr)

#读取ssGSEA结果文件
data=read.table("TPM100_total.score.txt", header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#读取风险文件
risk=read.table("totalRisk.txt", header=T, sep="\t", check.names=F, row.names=1)  #12

#合并数据
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
rt=cbind(data,risk[,c("RiskScore","Risk")])   #⭐️
rt=rt[,-(ncol(rt)-1)]

#绘图
mycol=c("#e11a0c","#337cba")

##（转置版）16种免疫细胞：箱线图####
immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c("Risk",immCell)]    #⭐️
data=melt(rt,id.vars=c("Risk"))    #data有3列：'Risk'、variable、value   #⭐️
colnames(data)=c("PyroGroup","Type","Score")   #重命名3列名称
data$PyroGroup=factor(data$PyroGroup, levels=c("Low","High"))

#转置顺序
data$Type = sort(data$Type, decreasing = T)   #xy转置后，Z→A。为了正序，需要把数据按照倒序排列；之后在x的label再把名称倒序排列

##【单独计算组间差异】##
library(rstatix)  
library(ggpubr)
##计算p值
#定义执行一次秩和检验的操作流程
wilcox_test(aDCs ~ Risk, data = rt1, alternative = 'two.sided')
#能够执行一次那么就可以执行多次，批量执行秩和检验
wilcox_data <- data.frame(Genesymbol=colnames(rt)[1:ncol(rt)-1])
for (i in 1:ncol(rt)){
  print(i)
  wilcox_data[i,2] <- wilcox.test(rt[,i] ~ Risk, data = rt, 
                                  alternative = 'two.sided', 
                                  exact = FALSE)[["p.value"]]
}

##计算fdr
wilcox_data$fdr <- p.adjust(wilcox_data$V2, method = "fdr")

##输出结果
colnames(wilcox_data) <- c("Genesymbol","p.value","adj.p.value")    #,"logFC")
wilcox_data$p.value3 = ifelse(wilcox_data$p.value<0.001,"<0.001",sprintf("%.03f", wilcox_data$p.value))
wilcox_data$p.signif = ifelse(wilcox_data$adj.p.value<0.001,"***",
                              ifelse(wilcox_data$adj.p.value<0.01,"**",
                                     ifelse(wilcox_data$adj.p.value<0.05,"*","NS")))

#提取Type分类
#Cell = unique(data$Type)
#immCell2=paste0(immCell,"(",wilcox_data$p.signif,")")  ##paste0连接cli和Pvalue和Star 

#设置X轴（Immune Type）标签颜色
colorsCell = ifelse(wilcox_data$p.value<0.05,"#fb3e35","black")
faceCell = ifelse(wilcox_data$p.value<0.05,2,1)

#绘图
ggboxplot(data, x='Type', y= "Score", color = "PyroGroup",
          notch = F, size = 0.4, width = 0.6, outlier.shape = 3,
          xlab="Immune cell",ylab="ssGSEA score", add = "none", 
          palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +     #转置xy
  #scale_x_discrete(labels = sort(immCell,decreasing = T),   #label倒序排列   #左=immCell；右=immCell2
  #                 position = ) +  #左=空白；右='top'
  scale_y_continuous(limits = c(0,1.05)) +
  rotate_x_text(50) +
  #theme_bw(base_size = 10) +
  #theme_classic() + 
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "top", # 去掉图例
        legend.key.size = unit(10,'pt'),
        # 修改网格线：
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #panel.grid.minor.y = element_line(linetype = "dashed"),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        # 去掉y轴刻度：
        #axis.ticks.y = element_blank(),
        #　拼图时去掉ｘ轴的标题
        #axis.title.y = element_blank(),
        #x轴标签
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        # y轴标签：
        axis.text.y = element_text(face = rev(faceCell),
                                   #color = rev(colorsCell),
                                   hjust = 1, # 左对齐
                                   size = 10,
                                   lineheight = 2),
        # 标题居中：
        plot.title = element_text(hjust = 0, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = wilcox_data$p.signif,hjust = 0, size = unit(3,'mm'), color = 'black',fontface = faceCell,
           x = unique(data$Type), y=1.03)   #文本注释还是在框里

ggsave('Cell.pdf',height = 4, width = 5)



##（转置版）13种免疫相关功能：箱线图#####
immFunc=c("APC_co_inhibition","APC_co_stimulation","CCR",
          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
rt2=rt[,c("Risk",immFunc)]
data2=melt(rt2,id.vars=c("Risk"))
colnames(data2)=c("PyroGroup","Type","Score")
data2$PyroGroup=factor(data2$PyroGroup, levels=c("Low","High"))

#转置顺序
data2$Type = sort(data2$Type, decreasing = T)   #xy转置后，Z→A。为了正序，需要把数据按照倒序排列；之后在x的label再把名称倒序排列

##【单独计算组间差异】##
library(rstatix)  
library(ggpubr)
##计算p值
#定义执行一次秩和检验的操作流程
wilcox_test(APC_co_inhibition ~ Risk, data = rt2)
#能够执行一次那么就可以执行多次，批量执行秩和检验
wilcox_data2 <- data.frame(Genesymbol=colnames(rt2)[2:ncol(rt2)])
for (i in 2:ncol(rt2)){
  print(i)
  wilcox_data2[i-1,2] <- wilcox.test(rt2[,i] ~ Risk, data = rt2, exact = FALSE)[["p.value"]]
}

##计算fdr
wilcox_data2$fdr <- p.adjust(wilcox_data2$V2,method = "fdr")

##计算logFC
data_control2 <- rt2[c(which(rt2$Risk == "Low")),]
data_trt2 <- rt2[c(which(rt2$Risk == "High")),]
wilcox_data2$foldchange <- colMeans(data_trt2[2:ncol(data_trt2)])-colMeans(data_control2[2:ncol(data_control2)])

##输出结果
colnames(wilcox_data2) <- c("Genesymbol","p.value","adj.p.value","logFC")
wilcox_data2$p.value3 = ifelse(wilcox_data2$p.value<0.001,"<0.001",sprintf("%.03f", wilcox_data2$p.value))
wilcox_data2$p.signif = ifelse(wilcox_data2$adj.p.value<0.001,"***",
                               ifelse(wilcox_data2$adj.p.value<0.01,"**",
                                      ifelse(wilcox_data2$adj.p.value<0.05,"*","NS")))

#设置X轴标签：提取Type分类+加上p值
#Cell = unique(data$Type)
immFunc2=paste0(immFunc,"(",wilcox_data2$p.value,")")  ##paste0连接cli和Pvalue和Star 

#设置X轴（Immune Type）标签颜色
colorsFunc = ifelse(wilcox_data2$p.value<0.05,"#fb3e35","black")
faceFunc = ifelse(wilcox_data2$p.value<0.05,2,1)



#绘图
p2 <- ggboxplot(data2, x='Type', y= "Score", color = "PyroGroup",
                notch = F, size = 0.4, width = 0.6, outlier.shape = 3,
                xlab="Immune function",ylab="ssGSEA score", add = "none",
                palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +     #转置xy
  scale_x_discrete(labels = sort(immFunc,decreasing = T),   #label倒序排列    #左=immFunc；右=immFunc2
                   position = ) +    #左=空白；右='top'
  scale_y_continuous(limits = c(0.25,1.05)) +
  rotate_x_text(50) +
  #theme_classic() +
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "none", # 去掉图例
        legend.key.size = unit(10,'pt'),
        # 修改网格线：
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        #panel.grid.minor.y = element_line(linetype = "dashed"),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        # 去掉y轴刻度：
        #axis.ticks.y = element_blank(),
        #x轴标签
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        # y轴标签：
        axis.text.y = element_text(
          face = rev(faceFunc),
          #color = rev(colorsFunc),
          hjust = 1, # 左对齐
          size = 10,
          lineheight = 2),
        # 标题居中：
        plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = wilcox_data2$p.signif,hjust = 0, size = unit(3,'mm'), color = 'black',fontface = faceFunc,
           x = unique(data2$Type), y=1.03)   #文本注释还是在框里
#是否显示统计值
stat_compare_means(aes(group=`Risk group`),
                   symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                    symbols = c("***", "**", "*", "")),
                   method="wilcox.test",   #默认为“wilcox.test”（非参数检验），可指定method = “t.test”，表示T检验（参数检验）
                   #method.args = list(alternative = "greater"),    # “two.sided”, “less”, “greater” 
                   label.y = 1.05,
                   label.x = ,
                   label.y.npc = 'centre',
                   label = "p.format")   #"p.signif"    "p.format" 


ggsave('Function.pdf',height = 4, width = 5)

#ggtitle("ssGSEA : Immune cell", ) +
#theme(plot.title = element_text(size =12))

##图片组合####
library(patchwork)
p1 + p2 +
  plot_layout(guides = "collect", #合并图例为1个
              nrow = 2,   #图像设置为2行 # nrow=按行填充   ncol=按列填充
              heights = c(1.1, 1)) &     #两列之间相对宽度比为3：1
  theme(legend.position='top')  #图例位置


plot_annotation(#title = "Risk group",   #左上角标题
  caption = "Wilcox.test, ** P<0.001, * P<0.01",  #右下角标题
  #tag_levels = "A",   #每张图的序号
)

#输出图片文件
ggsave(file="Cell+Func.pdf", height=8, width=6)
