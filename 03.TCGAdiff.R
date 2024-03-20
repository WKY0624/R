

####Method 1: Extract from DESeq2 results####
#Read DESeq2 result
DESeq=read.table("DESeq2.diff.txt", header=T, sep="\t", check.names=F)
DESeq=as.matrix(DESeq)
rownames(DESeq)=DESeq[,1]  
exp=DESeq[,2:ncol(DESeq)]  
dimnames=list(rownames(exp), colnames(exp)) 
DESeq=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames) 

#Read expression matrix
Expfile='filterTPM.txt'
Exp <- read.table(Expfile, header=T, sep="\t", check.names=F)
Exp=as.matrix(Exp)
rownames(Exp)=Exp[,1] 
exp_Exp=Exp[,2:ncol(Exp)]  
dimnames_Exp=list(rownames(exp_Exp), colnames(exp_Exp))  
Exp=matrix(as.numeric(as.matrix(exp_Exp)), nrow=nrow(exp_Exp), dimnames=dimnames_Exp) 

#Extract PRGs results
gene=read.table("PyroExp.txt", header=F, sep="\t", check.names=F)  
sameGene=intersect(as.vector(gene[,1]), rownames(DESeq))  
geneDiff=DESeq[sameGene,]    

#Output all difference analysis results
out=rbind(ID=colnames(geneDiff),geneDiff)    
write.table(out,file="stat_PRGs.txt",sep="\t",quote=F,col.names=F)

#Set filtering thresholds for logFC and P values
geneDiff2 = as.data.frame(geneDiff)
logFCfilter = 0    #0=1x；0.379=1.3x；0.585=1.5x；1=2倍；1.58=3x
fdrFilter = 0.05      
geneFilter=geneDiff2[(abs(as.numeric(as.vector(geneDiff2$log2FoldChange)))>logFCfilter & as.numeric(as.vector(geneDiff2$padj))<fdrFilter),]  

#Output the difference analysis results that meet the filtering conditions
outFilter=rbind(ID=colnames(geneFilter),geneFilter) 
write.table(outFilter, file="diff_PRGs.xls", sep="\t", row.names=T, col.names=F,quote=F)

#Output differential gene expression files
ExpDiff=Exp[as.vector(rownames(geneFilter)),]
expOut=rbind(ID=colnames(ExpDiff), ExpDiff)
write.table(expOut, file="diff_Exp_PRGs.txt", sep="\t", col.names=F, quote=F)



####Method 2: Calculate wilcox.test using TPM/Count data after filtering####
library(limma)
library(pheatmap)
expFile="Counts_vst.txt"

#Read input file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)]  
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  


#Clarify the number of normal and tumor cases   （# 4th digit of TCGA：0=tumor（01=primary；06=metastasis），1=Solid Tissue Normal，2=Control Analyte(CELLC)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)   #group = Divide the sample ID and obtain the 4th place, for example: "11A" in the original ID
group=sapply(strsplit(group,""), "[", 1)    #Divide the 4th digit (group) again (without a separator, so it is ""), take the first digit, for example: "1" in "11A"
group=gsub("2", "1", group)  
NormalNum=length(group[group==1])     
TumorNum=length(group[group==0])   
sampleType=c(rep(1,NormalNum), rep(2,TumorNum))


#Differential analysis
sigVec=c()
outTab=data.frame()
for(i in rownames(data)){
  if(sd(data[i,])<0.001){next}     
  wilcoxTest = wilcox.test(data[i,] ~ sampleType)   
  #In R, the wilcox. test() function can be used for Wilcoxon rank sum test and Mann Whitney U-test.
  #When the parameter is a single sample, or two samples are subtracted, or two parameters are paired=True, it is the Wilcoxon rank sum test.
  #When paired = FALSE (independent sample), it is the Mann Whitney U test.
  pvalue = wilcoxTest$p.value  
  #pvalue = p.adjust(pvalue, method = "fdr")
  if(pvalue<1){  
    Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    sigVec=c(sigVec, paste0(i, Sig))
    NormalGeneMeans=mean(data[i,1:NormalNum]) 
    TumorGeneMeans=mean(data[i,(NormalNum+1):ncol(data)])
    logFC=log2(TumorGeneMeans)-log2(NormalGeneMeans)  #The ratio of the difference in mean between two groups is fold change
    outTab=rbind(outTab,cbind(gene=i,
                              NormalMean=NormalGeneMeans,
                              TumorMean=TumorGeneMeans,
                              logFC=logFC,
                              pvalue=pvalue))  
  }
}

#Calculate FDR
pvalue=outTab[,"pvalue"]
fdr=p.adjust(as.numeric(as.vector(pvalue)), method="fdr")   
outTab=cbind(outTab, FDR=fdr)   
write.table(outTab, file = "stat_PRGs(Count_vst).xls", sep = "\t", row.names = F, quote = F)

#Set filtering thresholds for logFC and P values
logFCfilter = 0       #0=1x；0.58=1.5x；1=2x；1.58=3x
fdrFilter = 0.05       
outDiff = outTab[(abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$FDR))<fdrFilter),]  

#Output difference analysis results
write.table(outDiff, file="diff_PRGs(Count_vst).xls", sep = "\t", row.names = F, quote = F)

#Output differential gene expression files
exp = data[as.vector(outDiff[,1]),]
expOut = rbind(ID=colnames(exp), exp)
write.table(expOut, file = "diff_Exp_PRGs(Count_vst).txt", sep = "\t", col.names = F, quote = F)


