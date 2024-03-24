library(limma)              
expFile = "diff_Exp_PRGs(Count_vst).txt"  #03
#expFile = "diff_Exp_PRGs.txt"   #03
cliFile = "./DataSets/06.time_TCFi.txt"

#Read the expression data file and organize the input file
rt=read.table(expFile, header=T, sep="\t", check.names=F)  
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)] 
dimnames=list(rownames(exp), colnames(exp))  
data=matrix(data=as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)   
group=sapply(strsplit(group,""), "[", 1)    
group=gsub("2", "1", group)   
data=data[,group==0]  
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data)) 
data=avereps(data)   
data=data[rowMeans(data)>0,]  
#data=log2(data+1)   
data=t(data)

cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)  

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="total.expTime.txt",sep="\t",row.names=F,quote=F)   

#Cut total into train and test

total = out
trainFile <- read.table("trainSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)
testFile <- read.table("testSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)

train <- total[row.names(trainFile),]
write.table(train,"train.expTime.txt", sep = "\t", quote = F, row.names = F, col.names = T)

test <- total[row.names(testFile),]
write.table(test,"test.expTime.txt", sep = "\t", quote = F, row.names = F, col.names = T)



