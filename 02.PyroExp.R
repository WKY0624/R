
library(limma)  
setwd("")     

####Read input files and process data ####
rt=read.table("Counts_vst.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)] 
dimnames=list(rownames(exp), colnames(exp))  
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames) 

####Extracting the expression level of cell pyroptosis genes⭐️####
gene=read.table("Datasets/PRGs_226.txt", header=F, sep="\t", check.names=F)  #
sameGene=intersect(as.vector(gene[,1]), rownames(data))     
geneExp=data[sameGene,]  

out=rbind(ID=colnames(geneExp),geneExp) 
write.table(out,file="PyroExp.txt",sep="\t",quote=F,col.names=F)
