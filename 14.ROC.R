
library(survival)
library(survminer)
library(timeROC)

mycol <- c("#48af45","#337cba","#e11a0c")   

#定义绘制ROC曲线函数
bioROC=function(inputFile=null, rocFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #ROC曲线
  ROC_rt=timeROC(T=rt$DSS, delta=rt$DSEvent,    #⭐️
                 marker=rt$RiskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col=mycol[1],title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col=mycol[2],add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col=mycol[3],add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',   #x=0.4,y=0.3,
         c(paste0('1-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
           paste0('3-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[2]), ')'),
           paste0('5-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[3]), ')')),
         col=mycol, lwd=2, bty = 'n')   
  dev.off()
}

bioROC(inputFile="trainRisk.txt", rocFile="trainROC.pdf")
bioROC(inputFile="testRisk.txt", rocFile="testROC.pdf")
bioROC(inputFile="totalRisk.txt", rocFile="totalROC.pdf")



