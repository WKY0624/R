#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")
#install.packages("reshape2")
#install.packages("ggplot2")
#BiocManager::install("HDF5Array")

#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(reshape2)
library(ggplot2)

expFile="filterTPM.txt"      #01
riskFile="totalRisk.txt"     #12
gmtFile="./DataSets/h.all.v2023.1.Hs.symbols.gmt"    #基因集文件

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
#gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
#write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#删除正常样品
data=t(gsvaResult)
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
#row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
#data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]    #-2删除RiskScore
#xsort=colnames(risk)

#相关性分析
outTab=data.frame()
for(Geneset in colnames(data)){
  for(gene in colnames(risk)){
    x=as.numeric(data[,Geneset])
    y=as.numeric(risk[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Geneset=Geneset, cor, text, pvalue))
  }
}

##绘制相关性热图####
#前处理
outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
ysort = unique(outTab$Geneset)
xcolor = c("#564efd","#564efd","#fb3e35","#fb3e35","#fb3e35",
           "#fb3e35","#fb3e35","#564efd",'black')
xface = c('plain','plain','plain','plain','plain',
          'plain','plain','plain','plain','bold')


#绘图
ggplot(outTab, aes(Gene, Geneset)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.5)+  #内部边框颜色
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "#fb3e35",  #图形颜色
                       breaks=c(-0.5,0,0.5),labels=c(-0.5,0,0.5),   #图例标签
                       limits=c(-0.9,0.9)) +    
  geom_text(aes(label=text),col ="black",size = 3) +   #字体颜色和大小
  #theme_minimal() +    #去掉背景
  scale_y_discrete(limits=factor(sort(ysort, decreasing = T)),position = 'left') +
  scale_x_discrete(limits=c('H2BC8','CXCL8','PYCARD',"SEZ6L2","PRDM1","IFI27","GJA1","SIGLEC15","TRAF6","RiskScore")) +
  #scale_x_discrete(limits=factor(sort(colnames(risk),decreasing =F)),position = "bottom") +      #X轴名称显示顺序
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.title.x=element_blank(), 
        #axis.ticks.x=element_blank(), 
        axis.ticks = element_line(linewidth = 0.3),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   colour = 'black',
                                   #colour = xcolor,
                                   face = xface,
                                   size = 8),    #x轴字体
        axis.text.y = element_text(size = 8,
                                   colour = 'black')) +  # face = "bold")) +       #y轴字体
  ggtitle("Hallmarks GSVA correlation", ) +
  theme(plot.title = element_text(size =10, hjust = 0.5)) + #标题居中
  guides(fill = guide_colorbar(title = paste0("*** P<0.001","\n", " ** P<0.01","\n", "  * P<0.05","\n","\n","Correlation\nCoefficient"),
                               title.position = 'top',
                               title.theme = element_text(size = 8,face = "plain",colour = "black"),
                               label = T,   #图例的标签
                               label.theme = element_text(size = 8,face = "plain",colour = "black"),
                               raster = T,
                               frame.colour = NULL,    #图例边框颜色
                               barwidth = unit(3,"mm"),
                               barheight = unit(18,"mm"),
                               nbin = 50,   #指定绘制颜色的纸槽数，较大的值可以使色标更平滑
                               ticks = T,   #指定色条上的刻度是否可见,FALSE则为一个完全连续的图例
                               draw.ulim = T,   #指定上限刻度线是否可见
                               draw.llim = T,   #指定上限刻度线是否可见
                               )
         )
  #labs(fill = paste0("*** P<0.001","\n", " ** P<0.01","\n", "  * P<0.05","\n","\n","Correlation\nCoefficient")) +   #设置图例
  
ggsave("GSVAcor.pdf", height = 8.5, width = 5)



