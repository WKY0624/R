#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggplot2")
#BiocManager::install('ComplexHeatmap')

#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(dplyr)
library(grid)

expFile="total.normalize.txt"    ##08
riskFile="totalRisk.txt"        #12
geneFile = "./DataSets/ImmunomodulatorsList.txt"   #免疫调节分子列表文件


#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件，获取免疫检查点相关基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data1=t(data[sameGene,])

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data1), row.names(risk))
data1=data1[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

#相关性分析
outTab=data.frame()
for(checkpiont in colnames(data1)){   #对免疫检查点基因进行循坏
  for(gene in colnames(risk)){        #对模型基因进行循环
    x=as.numeric(data1[,checkpiont])    #获取各自的表达量
    y=as.numeric(risk[,gene])         
    corT=cor.test(x,y,method="spearman")   #相关性检验
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, checkpiont=checkpiont, cor, text, pvalue))
  }
}


#### 绘图 ####
#绘图前处理
outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$checkpiont=factor(outTab$checkpiont, levels=colnames(data1))
outTab$cor=as.numeric(outTab$cor)
outTab$pvalue=as.numeric(outTab$pvalue)
ysort = unique(outTab$checkpiont)
xface = c('plain','plain','plain','plain','plain','plain','plain','plain','plain','bold')
mycol = c('#cb6351','#2d8dc5','#499d71','#e9a04b','#9390c4')
mycol2 = c("#e11a0c","#337cba")

##热图
library(aplot)
library(tidyr)
p <- ggplot(outTab, aes(Gene, checkpiont)) + #设置xy轴
  geom_tile(aes(fill = cor), size = 0)+   #标题颜色
  scale_fill_gradient2(low = mycol2[2] , mid = "white", high = mycol2[1],   #"dodgerblue2"  "#fb3e35"
                       breaks=c(-0.5,0,0.5), labels = c(-0.5,0,0.5),   #图例标签
                       limits=c(-0.85,0.85),
                       name = "rho") +  #单元格颜色
  geom_text(aes(label=text),col ="#333333",size = 3) +   #字体颜色和大小
  theme_void() +    #去掉背景
  scale_y_discrete(limits=factor(sort(ysort, decreasing = T)),position = 'left')+
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) +
  theme(axis.title = element_blank(), 
        axis.ticks.x = element_line(linewidth = 0.3),    #⭐️前4张修改X轴显示与否
        axis.ticks.y = element_line(linewidth = 0.3), 
        #axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = xface, color = "#333333", margin = margin(0.2,0,0,0, 'cm')),   #x轴字体
        #axis.text.y = element_text(size = 10, color = "black")
        ) +       #y轴字体
  #labs(fill =paste0("***P<0.001","\n", "** P<0.01","\n", " * P<0.05","\n", "\n","Correlation")) +   #设置图例
  theme(legend.position = "bottom", legend.key.height = unit(2,"mm"), legend.key.width = unit(10,"mm"),
        legend.direction = "horizontal") +
  scale_x_discrete(position = "bottom") +        #X轴名称显示位置
  guides(fill='none')    #删除图例 +

table(geneRT$Immunomodulators)

left <- geneRT$Gene %>%  as.data.frame() %>% 
  mutate(group = rep(unique(geneRT$Immunomodulators), times = c(48,24,35,24))) %>%     #mod1: 48,30,35,24  #24,46,21
  mutate(p="")  %>%   #增加1列p，空值
  ggplot(aes(x=p, y=geneRT$Gene, fill = group)) + 
  geom_tile() +
  scale_fill_manual(values = mycol[1:4]) +
  #geom_tile(fill = mycol[1:4]) + 
  #geom_tile(fill = '#cb6352')+    #红
  #geom_tile(fill ='#e9a04b')+    #黄
  #geom_tile(fill ='#2d8dc5')+    #蓝
  #geom_tile(fill ='#499d71')+    #绿
  #geom_tile(fill ='#9390c4')+    #红
  #color = c('#cb6352','#e9a04b','#2d8dc5','#499d71','#9390c4')
  scale_y_discrete(position="right") +  
  theme_void() +
  theme(axis.text.y = element_text(colour = '#333333', size= 10, hjust = 1, margin = margin(0,0.1,0,0, 'cm'))) +
  theme(legend.position = "none") +
  xlab(NULL) + ylab(NULL) +
  labs(fill='group')

p %>% insert_left(left, width = .08) 

ggsave('ImmunomodulatorsCor.pdf', width = 3.5, height = 13)  #mod1（137）:3.5x20，mod2（78）：3.5x13\11\12


##单独画图例####
#连续性图例
dt = data.frame(Gene = c('A','B'), Geneset = c('Receptor','MHC','Immunoinhibitor','Immunostimulator','Chemokine'), cor = 1:10)  
for_legend <- ggplot(dt, aes(Gene, Geneset)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.5)+  #内部边框颜色
  scale_fill_gradient2(low = mycol2[2], mid = "white", high = mycol[1],  #图形颜色
                       breaks=c(-0.5,0,0.5),labels=NULL,   #图例标签   c(-0.5,0,0.5)
                       limits=c(-1,1),
                       name = "rho") +
  guides(fill = guide_colorbar(#title = paste0("*** P<0.001","\n", " ** P<0.01","\n", "  * P<0.05","\n","\n","Correlation Coefficient"),
                               title.position = 'left',
                               title.theme = element_text(size = 8,face = "bold",colour = "black"),
                               label = T,   #图例的标签
                               label.theme = element_text(size = 8,face = "plain",colour = "black"),
                               raster = T,
                               frame.colour = NULL,    #图例边框颜色
                               barwidth = unit(20,"mm"),
                               barheight = unit(3,"mm"),
                               nbin = 50,   #指定绘制颜色的纸槽数，较大的值可以使色标更平滑
                               ticks = T,   #指定色条上的刻度是否可见,FALSE则为一个完全连续的图例
                               draw.ulim = T,   #指定上限刻度线是否可见
                               draw.llim = T,   #指定上限刻度线是否可见
                               )) +
  theme(legend.direction = "horizontal")

legend1 <- cowplot::get_legend(for_legend)
grid.newpage()
grid.draw(legend1)

#export: 3x3  0-legend-continus.pdf

# 注释图例
library(ComplexHeatmap)
color = c('#cb6352','#e9a04b','#2d8dc5','#499d71','#9390c4')
legend2 <- Legend(labels = c('Receptor','MHC','Immunoinhibitor','Immunostimulator','Chemokine'), 
                  title = "Immunomodulators", title_gap = unit(3,'mm'),
                  row_gap = unit(2,"mm"),
                  legend_gp = gpar(fill = color))
draw(legend2)

#export: 3x3 