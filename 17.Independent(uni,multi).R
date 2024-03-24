library(survival)   

####数据处理####
#输入数据
riskFile = "totalRisk.txt"
cliFile = "06.clinical477.txt"  #需要赋值

#数据读取
risk = read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    
cli = read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      
cli = cli[,grep("Age|Gender|Tstage|Nstage|Mstage",colnames(cli))]  

#数据合并
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(RFS=risk[,1],    #⭐️
         Recurrence=risk[,2], cli, RiskScore=risk[,(ncol(risk)-1)])

####uniCox####
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(RFS, Recurrence) ~ rt[,i], data = rt)   #⭐️
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
uniTab
write.table(uniTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)



####multiCox####
uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
rt1=rt[,c("RFS", "Recurrence", as.vector(uniTab[,"id"]))]   #⭐
multiCox=coxph(Surv(RFS, Recurrence) ~ ., data = rt1)   #⭐️
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
multiTab
write.table(multiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)



####森林图####
#install.packages("forplo")
library(forplo)

###单变量Cox绘图###
#数据读取
uni = read.table('uniCox.txt', header = T, row.names = 1, check.names = F, sep = '\t')  #含p值
uni2 = uni[,1:3]   #仅HR、LL、HL 3列 
uni3 = uni2[c(1,2,7,8,9,10),]   # 选定特定的行标签
uni4 = uni[c(1,2,7,8,9,10),]    # 含p值
#绘图
pdf('uniCox.pdf', height = 3, width = 7.5)
forplo(uni3,     #仅HR、LL、HL 3列 
       em = "HR",
       row.labels = c('Age','Gender', 'T Stage','N Stage','M Stage','PyroScore'),
       pval = sprintf("%.03f",uni4$pvalue),    #加1列p值
       xlim = c(0.5,12),
       ci.sep = '-',
       ci.lwd = 2,
       ci.edge = T,   #置信区间首位的边缘线
       char = 15,    #设置标记的符号，类似基础R绘图函数中的pch
       size = 1.2,   #设置标记的大小 
       insig.col = 'gray60',
       fill.colors = c("black","gray60","black","gray60","black","black"),
       right.bar = T,   #添加右垂直线
       rightbar.ticks = T,   #添加右垂直线的刻度
       left.bar = T,
       leftbar.ticks = T,
       left.align = F,  #居右对齐
       margin.top = 10,
       margin.bottom = 5,
       margin.right = 1,
       title = ' ')
dev.off()

###多变量Cox绘图###
multi = read.table('multiCox_total_6.txt', header = T, row.names = 1, check.names = F, sep = '\t')
multi2 = multi[,1:3]
#绘图
pdf('multiCox.pdf', height = 3, width = 7.5)
forplo(multi2,
       em = "HR",
       row.labels = c('Age', 'Gender','T Stage','N Stage','M Stage','PyroScore'),
       pval = sprintf("%.03f",multi$pvalue),
       xlim = c(0.4,7),
       ci.sep = '-',
       ci.lwd = 2,
       ci.edge = T,   
       char = 15,  
       size = 1.2,  
       insig.col = 'gray60',
       right.bar = T,   
       rightbar.ticks = T,  
       left.bar = T,
       leftbar.ticks = T,
       left.align = F,  
       margin.top = 10,
       margin.bottom = 5,
       margin.right = 1,
       title = ' ')
dev.off()