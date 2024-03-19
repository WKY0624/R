
#BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)
library(latex2exp)   


####Count data processing + DESeq2 differential expression analysis####
rt <- read.table("newCounts_477+59.txt",header=T,sep="\t",comment.char="",check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp),colnames(exp))
data = matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data = avereps(data)    #duplicate genes take the average
data = data[rowMeans(data)>100,]   #Filter and retain Counts with an average greater than 100（14826）
data = data[apply(data, 1, sum) > 0 , ]    ##Preserve rows with sum>0, exclude genes that are 0 in all samples
data = as.data.frame(data)
data <- floor(data)   #Downward retrieval
#data <- ceiling(data)   #Upward retrieval

####Read sample grouping information####
group1 = sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1 = sapply(strsplit(group1,""), "[", 1)
group1 = gsub("2", "1", group1)
conNum = length(group1[group1==1])       
treatNum = length(group1[group1==0])     

Sample = colnames(data)
Type = c(rep(1,conNum), rep(2,treatNum))
exp = cbind(Sample, Type)
exp = as.data.frame(exp)
colnames(exp) = c("id", "Type")
exp$Type = ifelse(exp$Type==1, "Normal", "Tumor")
exp = as.matrix(exp)
rownames(exp) = exp[,1]
exp = exp[,2:ncol(exp)]

countData = as.data.frame(exp)
colnames(countData) = c("condition")
colData$condition <- as.factor(countData$condition)
write.table(countData,"Counts.txt",sep="\t",quote=F,col.names = NA)

####Running DESeq2####
# Building objects in DESeq2
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
# Designate which group as the control group
dds$condition <- relevel(dds$condition, ref = "Normal")
# Calculate the normalization coefficient for each sample
dds <- estimateSizeFactors(dds)
# Estimating the dispersion of genes
dds <- estimateDispersions(dds)
# Differential analysis
dds <- nbinomWaldTest(dds)
dds <- DESeq(dds)
res <- results(dds)
write.table(res,"DESeq2.diff.txt",sep="\t",quote=F,col.names = NA)

#### Normalized (reducing dispersion) ####
vsd <- vst(dds, blind = FALSE)
countData_old<- assay(dds)
countData_new<- assay(vsd)
n.sample=ncol(countData)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
write.table(countData_new,"Counts_vst.txt",sep="\t",quote=F,col.names = NA)

# Pre normalized plot
pdf(file="rawBox.pdf")
boxplot(countData_old, col = cols,main="expression value",xaxt = "n")
dev.off()
# Normalized plot
pdf(file="normalBox.pdf")
boxplot(countData_new, col = cols,main="expression value",xaxt = "n")
dev.off()
pdf(file="histold.pdf")
hist(countData_old)
dev.off()
pdf(file="histnew.pdf")
hist(countData_new)
dev.off()



####Filter low expression genes based on Count and save TPM matrix####
count <- read.table('Counts.txt', header=T,sep="\t",comment.char="",check.names=F)
tpm <- read.table('newTPM_477.txt', header=T,sep="\t",comment.char="",check.names=F)
count=as.matrix(count)
rownames(count)=count[,1]
exp=count[,2:ncol(count)]
dimnames=list(rownames(exp),colnames(exp))
count=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
count=as.data.frame(count)

tpm=as.matrix(tpm)
rownames(tpm)=tpm[,1]
exp2=tpm[,2:ncol(tpm)]
dimnames2=list(rownames(exp2),colnames(exp2))
tpm=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
tpm=as.data.frame(tpm)
sameGene <- intersect(rownames(count), rownames(tpm))
filterTPM <- tpm[sameGene,]
out=rbind(ID=colnames(filterTPM),filterTPM)  #行整合=新增1行ID-样品名 + filterTPM   

write.table(out,file="filterTPM.txt",sep="\t",quote=F,col.names=F)


####🐰Visualization of single gene expression####
#If using countData to draw a graph, some parts cannot be annotated. Draw a graph using standardized Counts_vst
read.countData <- read.table('Counts_vst.txt',header=T,sep="\t",comment.char="",check.names=F)
rt2 =as.matrix(read.countData)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
countData2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)

group1=sapply(strsplit(colnames(countData2),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])   
treatNum=length(group1[group1==0])     


#Extraction of expression levels of single genes
gene="IFI27"   
data=t(countData2[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

mycol <-  c("lightsteelblue1",'mistyrose')

ggboxplot(exp, 
          x="Type", y="gene", 
          color= "Type",  
          fill = 'snow',
          linetype = 1,
          lwd = 0.6,
          xlab="",
          ylab=paste0(gene, " expression"),
          #ylim= c(3,4),
          width = 0.3,   
          error.plot = "errorbar",
          bxp.errorbar = T,
          bxp.errorbar.width = 0.15,
          outlier.shape = 3,
          legend.title="Group", 
          palette = "aaas",  
          notch = T) +  
  stat_compare_means(comparisons=my_comparisons, label = 'p.format') +
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black'),
        legend.position = "none")

ggsave(file=paste0(gene,".diff.pdf"), width=2.5, height=4)

