#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")   #人的比对基因
library("enrichplot")
library("ggplot2")
library("stringr")

pvalueFilter=0.05           #p值过滤条件
qvalueFilter=0.05           #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #19

#### 基因名字转换为基因id ####
genes=as.vector(rt[,1])
entrezIDs=mget(genes, 
               org.Hs.egSYMBOL2EG, 
               ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#### GO富集分析 #####
kk=enrichGO(gene = gene,
            OrgDb = org.Hs.eg.db,    #参照基因：人类
            pvalueCutoff =0.05,   #p值过滤条件，一般选择0.01或0.05
            qvalueCutoff = 0.05,    ##矫正后的p值过滤条件
            ont="all",   #all=BP\MF\CC
            readable =T)
GO=as.data.frame(kk)

#保存富集结果
write.table(GO,file="GO.total.txt",sep="\t",quote=F,row.names = F)


#### 绘图 ####
#定义显示Term数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}


##气泡图
pdf(file="bubble.pdf",width = 10,height =7)
bub=dotplot(kk,
            x = "GeneRatio",
            color = colorSel,   #"p.adjust",
            size = "Count", 
            showCategory = showNum,
            orderBy = "GeneRatio",   #Y轴顺序
            label_format = 100,   #GO条目的字数长度
            split="ONTOLOGY",
            font.size = 10) + facet_grid(ONTOLOGY~., scale='free')  

print(bub)
dev.off()


##柱状图
pdf(file="barplot.pdf",width = 10,height =7)
bar=barplot(kk, 
            drop = TRUE, 
            showCategory =showNum,
            split="ONTOLOGY",
            color = colorSel,
            label_format = 100, #GO条目的字数长度
            font.size = 10) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()


