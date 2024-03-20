#BiocManager::install("ConsensusClusterPlus")
library(limma)
library(survival)
library(ConsensusClusterPlus)

expFile="diff_Exp_PRGs.txt"   
cliFile="./DataSets/04.time_TCFi.txt" 

#Read expression matrix data
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
#data=log2(data+1)   #⭐️TPM
#data=t(data)

#Read clinical prognostic data
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)  
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli,data)

##Selecting Cluster Genes: Manually####
#data=data[,grep("H2BC13|EZH2|MKI67",colnames(data))]
#data=t(data)

##Selecting Cluster Genes: Single Factor COX Analysis####
sigGenes=c()
for(i in colnames(rt)[3:ncol(rt)]){
  cox=coxph(Surv(DSS, DSEvent) ~ rt[,i], data = rt)   #⭐️Need to modify the header to be consistent with time
  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){ sigGenes=c(sigGenes,i) }
}
data2=t(data[,sigGenes])
sigGenes


##clustering####
results=ConsensusClusterPlus(
  data2,
  maxK = 9,    #Maximum number of categories evaluated (k)
  reps = 100,      #The number of iterations for each k
  pItem = 0.8,    #The proportion of resampled samples is set to 80%
  pFeature = 1,   #The proportion of resampling features, set to 100% or 80%
  title = 'PDF_cluster',    #Location for storing images and files
  clusterAlg = "km",  #Clustering algorithm, "km", hierarchical clustering "pam", "hc"
  distance = "euclidean",  #Select "euclidean" for clustering distance km, and Pearson correlation distance for hc
  seed = 666,  
  tmyPal = c("white","#F7FBFF","#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C2","#4292C6", "#2171B5","#2171B9"),
  plot = "pdf",
  writeTable = F,  #If it is true, output the consistency matrix, ICL, and log to a CSV file
  verbose = F  #If it is true, progress information can be output on the screen
  )
##If the prompt "10 iterations still not aggregated" indicates that the number of iterations is not enough, the resulting partition is unstable, and the algorithm has not converged to the optimal solution


##Select the optimal k based on the algorithm####
#See details：https://www.bioinfo-scrounger.com/archives/734/
Kvec = 2:9
x1 = 0.1; x2 = 0.9          # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")          # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}                                 
optK = Kvec[which.min(PAC)]  # optimal k
optK



##Output Typing Results####
clusterNum = optK     #Based on the best K, output the typing result
Cluster = results[[clusterNum]][["consensusClass"]]
Cluster = as.data.frame(Cluster)
Cluster[,1]=paste0("C", Cluster[,1])
ClusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(ClusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)

