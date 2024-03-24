library(Rtsne)
library(scatterplot3d)
library(ggplot2)

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
  #读取输入文件,提取数据
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  data=rt[c(3:(ncol(rt)-2))]
  Risk=rt[,"Risk"]  #⭐️
  
  #### PCA ####
  data.pca=prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  #绘图前处理
  group = levels(factor(Risk))
  mycol = c("#e11a0c","#337cba")
  col=mycol[match(Risk,group)]
  #绘制PCA图（三维）
  PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2], PC3 = pcaPredict[,3],
                   PC4 = pcaPredict[,4], PC5 = pcaPredict[,5], PC6 = pcaPredict[,6],
                   PC7 = pcaPredict[,7], PC8 = pcaPredict[,8], PC9 = pcaPredict[,9],Risk=Risk)	
  pdf(file=pcaFile, height = 4.5, width = 4.5)
  p <- scatterplot3d(PCA[,c(3,5,7)], grid = T, box=T, angle = 45, pch = 19, color = alpha(col,0.5),   #train:c(7,3,5)
                col.axis="gray30",   #盒子线条的颜色
                xlab = "PC1", ylab = "PC2", zlab = "PC3"
                #col.grid="lightblue"   #背景网格线的颜色
                )
  print(p)
  dev.off()  
  
  #绘制PCA图(二维)
  #pdf(file=pcaFile, height=4.5, width=5.5)     
  #ggplot(data = PCA, aes(PC1, PC7)) + geom_point(aes(color = Risk)) +
  #  scale_colour_manual(name="Risk",  values =c("#e11a0c","#337cba"))+    #"#fb3e35","#564efd"
  #  theme_bw()+
  #  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #print(p)
  #dev.off()
  
  #### t-SNE ####
  tsneOut=Rtsne(data, dims=3, perplexity=10, verbose=F, max_iter=500, check_duplicates=F)
  tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2], tSNE3 = tsneOut$Y[,3], Risk=Risk)	
  #绘制tSNE图（二维）
  pdf(file=tsneFile, height=4.5, width=4.5)      
  p <- scatterplot3d(tsne[,c(2,3,1)], grid = T, box=T, angle = 45, pch = 19, color = alpha(col,0.5),  #train:c(3,1,2)
                col.axis="gray30",   #盒子线条的颜色
                xlab = "tSNE1", ylab = "tSNE2", zlab = "tSNE3"
                #col.grid="lightblue"   #背景网格线的颜色
  )
  print(p)
  dev.off()
  
  
  #绘制tSNE图（二维）
  #pdf(file=tsneFile, height=4.5, width=5.5)     
  #ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Risk)) +
  #  scale_colour_manual(name="Risk",  values =c("#e11a0c","#337cba"))+   #"#fb3e35","#564efd"
  #  theme_bw()+
  #  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #print(p)
  #dev.off()
}


bioPCA(inputFile="trainRisk.txt", pcaFile="PCA.train.pdf", tsneFile="t-SNE.train.pdf")
bioPCA(inputFile="testRisk.txt", pcaFile="PCA.test.pdf", tsneFile="t-SNE.test.pdf")
bioPCA(inputFile="totalRisk.txt", pcaFile="PCA.total.pdf", tsneFile="t-SNE.total.pdf")
