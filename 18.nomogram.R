library(foreign)
library(rms)

#输入文件
riskFile="totalRisk.txt"   
cliFile="06.clinical477.txt"    #提取所需行列整理
#Sample	Age	Gender	Tstage	Nstage	Mstage
#TCGA-4C-A93U	74	Female	T3-T4	N1	M0

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)
cli$Tstage=as.factor(cli$Tstage)   #将分类变量处理成factor，二分类不影响结果，多分类必须处理成factor


#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("RFS", "Recurrence", "Risk")], cli)
rt=rt[,grep("RFS|Recurrence|Risk|Age|Gender|Tstage|Nstage|Mstage",colnames(rt))] 
rt<-na.omit(rt)
colnames(rt)[colnames(rt)=="Risk"] <- "RiskGroup"   #修改列名
colnames(rt)[colnames(rt)=="Tstage"] <- "TStage"   #修改列名
colnames(rt)[colnames(rt)=="Nstage"] <- "NStage"   #修改列名
colnames(rt)[colnames(rt)=="Mstage"] <- "MStage"   #修改列名

#定义终点事件
#rt$Recurrence <- ifelse(rt$Recurrence=="Recurrence",1,0)

#设置参考组
##rt$Recurrence <- relevel(rt$Recurrence ,ref='Recurrence')


#数据打包
dd<-datadist(rt)
options(datadist='dd') 

######构建Cox回归方程##########
##拟合模型:lrm=logistic回归；glm=广义线性模型；psm=倾向得分匹配；【cph/coxph=Cox比例风险回归模型】。
#psm是用的rms包中的参数回归模型，假设生存函数是lognormal分布。
#【cph是cox等比例风险模型，不注重生存函数的分布，只关注暴露因素的效应，流行病学中主要采用f2拟合方程】。

coxm <-cph(Surv(RFS,Recurrence==1)~RiskGroup+TStage+NStage+MStage+Gender+Age, data=rt, x=T, y=T, surv=T)
surv<- Survival(coxm) # 建立生存函数
surv1<- function(x)surv(1,lp=x) # 定义time.inc,1年RFS
surv2<- function(x)surv(3,lp=x) # 定义time.inc,3年RFS
surv3<- function(x)surv(5,lp=x) # 定义time.inc,5年RFS


###### 图1.绘制列线图 #########
plot(nomogram(coxm, fun=list(surv1,surv2,surv3),
              lp=F,  #是否显示线性预测坐标（不知道怎么解读，所以就不显示了吧）
              #conf.int = c(0.1,0.3),  #是否显示置信区间
              funlabel=c("1-Year RFS",'3-Year RFS','5-Year RFS'),
              maxscale=100,   #maxscale参数指定最高分数，一般设置为100或者10分
              #fun.at = c('0.05','seq(0.1,0.9,by=0.05)','0.95'), 
              #fun.at = seq(0.05,0.95,0.05)
              fun.at=c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  #fun.at 设置生存率的刻度
              ),
     xfrac=.25,  #xfrac 设置数值轴与最左边标签的距离
     #fun.side = c(3,1,3,1,3,1,3,1,3)   #fun的刻度显示在上=3或下=1
     label.every = 1,  #变量之间的间隔
     col.grid = gray(c(0.8,0.95)),   #参考线的颜色
     #conf.space = c(0.1,0.3),
     #col.conf = c('red','green')
     )   

#Export：9x5
#Export(cutoff)：8x5
dev.off()


#####输出列线图的风险打分文件########
nomoRisk=predict(coxm, data=rt, type="lp")
rt2=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt2), rt2)  #有个报错：因子层次有错，产生了NA（是因为输出文件增加了1行，可忽略此报错）
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


######计算C-index########
f<-coxph(Surv(RFS,Recurrence==1)~Age+Gender+TStage+NStage+MStage+RiskGroup,data=rt)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
write.table(c_index,file = "cindex.txt", sep = '\t', col.names = F,quote = F)
#此处再次解释一下上述R语言代码中计算的C-Index统计量的含义，
#此值可以参考ROC曲线下面积去理解，取值范围0~1，
#越接近1说明我们构建的Cox回归模型预测价值越大，
#一般来说C-Index=0.7即表明模型具有良好的预测价值。

####计算C-index的补数（可选）
library(Hmisc)
S<-Surv(rt$RFS,rt$Recurrence==1)
rcorrcens(S~predict(coxm),outx=TRUE)

##### 图2.校准曲线（三合一）#######
#计算
#1年校准曲线
f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=1)
cal1 <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt2)/3), B=1000)
#3年校准曲线
f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=3)
cal3 <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt2)/3), B=1000)
#5年校准曲线
f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=5)
cal5 <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt2)/3), B=1000)

#绘图
pdf(file="calibration.pdf", width=5, height=5)
plot(cal1, xlim=c(0.6,1), ylim=c(0.6,1),
     xlab="Nomogram-predicted RFS", 
     ylab="Observed RFS", 
     lwd=1,
     lty=5,
     col="green",
     sub=F)
plot(cal3, xlim=c(0.6,1), ylim=c(0.6,1), xlab="", ylab="", 
     lwd=1, lty=5, col="blue", sub=F, add=T)
plot(cal5, xlim=c(0.6,1), ylim=c(0.6,1), xlab="", ylab="",  
     lwd=1, lty=5, col="red", sub=F, add=T)

#加粗线条
lines(cal5[,c("mean.predicted","KM")],type="b",lwd=3,col='red', pch=4)
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=3,col='blue', pch=4)
lines(cal1[,c("mean.predicted","KM")],type="b",lwd=3,col='green', pch=4)
#美化参考对角线
#abline(0,1,lty=1,lwd=2,col='gray')

#图例
legend('bottomright', c('1-Year', '3-Year', '5-Year'),
       col=c("green","blue","red"), lwd=2, bty = 'n')

dev.off()


#基于训练集采用bootstrap冲抽样的方法验证列线图的预测准确性。
#结果解读：横坐标为根据列线图预测的每个患者的10年生存概率，
#纵坐标为每个患者实际的10年生存概率，
#如果图中的红色线条恰好与蓝色虚线完全重合时最为理想。



#### 图3.C指数条形图 ####
modelFile="5models.txt"       #模型数据文件
#整理各模型分类数据
#Sample	PRNomogram	TNM	ATA	MACIS	EORTC
#TCGA-4C-A93U	1.891853008	3	3	4	5

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)
cli$Tstage=as.factor(cli$Tstage)  

#读取模型文件
model=read.table(modelFile,header=T, sep="\t", check.names=F, row.names=1)
model$PRNomogram=as.numeric(model$PRNomogram)
model$TNM=as.factor(model$TNM)   
model$ATA=as.factor(model$ATA)
model$MACIS=as.factor(model$MACIS)
model$EORTC=as.factor(model$EORTC) 


#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
model = model[samSample,,drop=F]
rt=cbind(risk1[,c("RFS", "Recurrence", "RiskScore")], cli, model)
rt<-na.omit(rt)

#数据打包
dd<-datadist(rt) # 将数据整合
options(datadist='dd') # 将数据整合


######计算C-index
f<-coxph(Surv(RFS,Recurrence==1)~ PRNomogram, data=rt)    #PRNomogram\ATA等
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index   #将结果保存到5cindex.txt

######绘制条形图+误差值
cindex <- read.table('5cindex.txt',header=T, sep="\t", check.names=F)
#model	Cindex	SE
#PR Nomogram	0.775	0.038
#TNM	0.631	0.041
#ATA	0.648	0.038
#MACIS	0.657	0.038
#EORTC	0.63	0.046


library(ggprism)

pdf('Cindex.pdf', height = 3, width = 4)   #方=3x3；长=3x4
ggplot(cindex, aes(x=factor(model, levels = c(unique(cindex$model))), y=Cindex)) +
  geom_bar(stat = 'identity', 
           fill = 'white',   # 改变填充颜色
           color = c('red','green2',"blue",'orange','purple'),   #边框颜色,
           position= "dodge",
           width = 0.7  # 柱子宽度
  ) + 
  coord_cartesian(ylim = c(0.5, 0.8)) +
  geom_errorbar(aes(ymin=Cindex-SE, ymax=Cindex+SE), width=.2,
                color = c('red','green2',"blue",'orange','purple')) + 
  labs(y = 'C-index', x="") + # 调整横纵坐标轴标题字号
  theme_prism(base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_size = 12,  # 图形的字体大小
              base_line_size = 0.5, # 坐标轴的粗细
              axis_text_angle = 45) + # 可选值有 0，45，90，270
  theme(plot.margin = margin(t = 10,  # 顶部边缘距离
                             r = 10,  # 右边边缘距离
                             b = 0,  # 底部边缘距离
                             l = 10)) # 左边边缘距离

dev.off()


#### 图4.DCA曲线 ####
#options(unzip ='internal')  #如果下行报错了，再运行此行
#devtools::install_github('yikeshu0611/ggDCA')
library(rms)
library(ggDCA)
library(survival)  
library(ggprism)

rm(list = ls()) 

riskFile = "totalRisk.txt"     #风险文件
cliFile="clinical_total.txt"   #临床数据文件
#Sample	Age	Gender	Tstage	Nstage	Mstage
#TCGA-4C-A93U	74	Female	T3-T4	N1	M0
modelFile="5models.txt"       #模型数据文件

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)
cli$Tstage=as.factor(cli$Tstage)   #将分类变量处理成factor，二分类不影响结果，多分类必须处理成factor

#读取模型文件
model=read.table(modelFile,header=T, sep="\t", check.names=F, row.names=1)
model$PRNomogram=as.numeric(model$PRNomogram)
model$TNM=as.factor(model$TNM)  
model$ATA=as.factor(model$ATA)
model$MACIS=as.factor(model$MACIS)
model$EORTC=as.factor(model$EORTC) 

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
model=model[samSample,,drop=F]
rt=cbind(risk1[,c("RFS", "Recurrence", "RiskScore")], cli, model)
#rt=rt[,grep("RFS|Recurrence|Risk|Age|Gender|Tstage|Nstage|Mstage",colnames(rt))] 
rt<-na.omit(rt)

#数据打包
dd<-datadist(rt) # 将数据整合
options(datadist='dd') # 将数据整合

##### model：构建多因素Cox回归模型 ######
PRNomogram <- coxph(Surv(RFS,Recurrence==1) ~ Age+Gender+Tstage+Nstage+Mstage+RiskScore, data=rt)
ATA <- coxph(Surv(RFS,Recurrence==1) ~ ATA, data=rt)
TNM <- coxph(Surv(RFS,Recurrence==1) ~ TNM, data=rt)
MACIS <- coxph(Surv(RFS,Recurrence==1) ~ MACIS, data=rt)
EORTC <- coxph(Surv(RFS,Recurrence==1) ~ EORTC, data=rt)

###### model：计算DCA######
dca_1<-dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=1)
dca_3<-dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=3)
dca_5<-dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=5)
dca_10 <- dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=10)
dca_13510<- dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=c(1,3,5,10))
dca_3510<- dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=c(3,5,10))
dca_35<- dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=c(3,5))
dca_135<- dca(PRNomogram,ATA,TNM,MACIS,EORTC, times=c(1,3,5))


###### model：绘制DCA曲线######  
pdf('DCA.pdf', height = 5, width = 5)  
ggplot(dca_135,aes(x=thresholds, y=NB, group=model),
       lwd = 0.6) + 
  #scale_x_continuous(limits = c(0,1), guide = 'prism_minor') +    #3年：limits = c(0,0.42)；5年：limits = c(0,0.65)；1年：limits = c(0,0.157),
  #scale_y_continuous(limits = c(-0.01,0.03), guide = 'prism_minor') +   #3年：limits = c(-0.03,0.095)；5年：limits = c(-0.04,0.15)；1年：limits = c(-0.01,0.03),
  scale_color_manual(values = c('red','blue','green','orange','purple','gray','black')) + 
  #scale_size_manual(values=c(1,2,3,3,3,3,3)) +
  #scale_linetype_manual(values = c('twodash','twodash','twodash','twodash','twodash','solid','solid')) +  #自定义线条的类型
  scale_linetype_manual(values = c('solid','solid','solid','solid','solid','twodash','twodash')) +
  theme(legend.title = element_blank()) +   #删除图例标题) +    #更改图例位置：'top''bottom''left''right'或放置于图里c()坐标 
  annotate('text', x = 0.22, y = 0.03, label = "1-Year", colour="black", size = 5) +   #3年：x = 0.21 , y = 0.095；5年：x = 0.32 , y = 0.15；1年：x = 0.075, y = 0.03,
  theme_classic() +
  theme(axis.title.x=element_text(vjust=-5, size=15),  # X axis title
        axis.title.y=element_text(size=15, vjust=5),  # Y axis title
        axis.text.x=element_text(size=15,
                                 #angle = 45,
                                 #color = "red",
                                 vjust= -1),  # X axis text
        axis.text.y=element_text(size=15, hjust = 0)) +  # Y axis text
  theme(plot.margin = margin(t = 15, r = 10, b = 30, l = 30))  # 顶部边缘距离
                             
dev.off()  

