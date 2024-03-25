
library(limma)
logFCfilter=0         #logFC过滤条件(logFC=1: Fold change=2     logFC=0.585: Fold change=1.5)
fdrFilter=0.05        

RiskDiff=function(expFile=null, riskFile=null, outFile=null){
  #读取文件,并对输入文件整理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  
  #读取risk文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  sameSample=intersect(colnames(data),row.names(risk))
  data=data[,sameSample]
  risk=risk[sameSample,]
  
  #提取low risk和high risk样品
  riskLow=risk[risk$Risk=="Low",]
  riskHigh=risk[risk$Risk=="High",]
  dataLow=data[,row.names(riskLow)]
  dataHigh=data[,row.names(riskHigh)]
  data=cbind(dataLow,dataHigh)
  data=data[rowMeans(data)>0,]   #原0.5
  conNum=ncol(dataLow)
  treatNum=ncol(dataHigh)
  Type=c(rep(1,conNum),rep(2,treatNum))
  
  
  #差异分析
  outTab=data.frame()
  for(i in row.names(data)){
    geneName=unlist(strsplit(i,"\\|",))[1]
    geneName=gsub("\\/", "_", geneName)
    rt=rbind(expression=data[i,],Type=Type)
    rt=as.matrix(t(rt))
    wilcoxTest=wilcox.test(expression ~ Type, data=rt)
    conGeneMeans=mean(data[i,1:conNum])
    treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    logFC=treatGeneMeans-conGeneMeans
    pvalue=wilcoxTest$p.value
    conMed=median(data[i,1:conNum])
    treatMed=median(data[i,(conNum+1):ncol(data)])
    diffMed=treatMed-conMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
      outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    }
  }
  pValue=outTab[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
  outTab=cbind(outTab,fdr=fdr)
  
  #输出差异表格
  outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
  write.table(outDiff, file=outFile, sep="\t", row.names=F, quote=F)
}

#调用函数，进行差异分析
RiskDiff(expFile="total.normalize.txt",  #08
         riskFile="totalRisk.txt",  #12
         outFile="riskDiff.txt")
