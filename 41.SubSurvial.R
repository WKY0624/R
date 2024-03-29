library(survival)
library(survminer)

riskFile ='totalRisk.txt'   #12

#risk=read.table(riskFile, header = T, sep='\t', check.names = F)

####TCGA####
TCGA = read.table("PyroExp.txt", header = T, check.names = F, sep = "\t", row.names = 1)
TCGA = TCGA["SEZ6L2",]
group=sapply(strsplit(colnames(TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
TCGA = TCGA[,group==0]   #删除正常样本
colnames(TCGA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(TCGA))   #规范TCGA样本名称
TCGA = as.data.frame(t(TCGA))
TCGA$SEZ6L2 = log2(TCGA$SEZ6L2+1)
TCGA$group = ifelse(TCGA$SEZ6L2 > median(TCGA$SEZ6L2),"High","Low")


#### 临床数据 ####
cliFile = '06.clinical477.txt'
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names = 1)
cli$RFS = cli$RFS/365

#### 合并数据 ####
samSample = intersect(row.names(TCGA),row.names(cli))
rt = cbind(SEZ = TCGA[,"group"],cli)
#rt<-na.omit(rt)   #删除NA行

#### 临床参数赋值 ####
#新增1列，分为Age>=55和<55
rt$Age2 = ifelse(rt$Age < "55", "<55",">=55")

#新增1列，将肿瘤最大径分为微小癌1和非微小癌2
rt$TD1 <- rt$TumorSize
rt$TD1 = ifelse(rt$TD1 =='0.1~1.0',1,2)

rt$TD4 <- rt$TumorSize
rt$TD4 = ifelse(rt$TD4 =='4.1~max',2,1)

#新增1列，分为有HT vs 无HT
rt$HT <- rt$`Coexist Disease`
rt$HT = ifelse(rt$HT =='Lymphocytic thyroiditis',1,0)

#新增一列，分为HT、Nodular、None
rt$`Combined Disease2` <- rt$`Combined Disease`
rt$`Combined Disease2`=ifelse(rt$`Combined Disease2`=='None',0,
                              ifelse(rt$`Combined Disease2`=='Nodular hyperplasia',1,2))    #2=HT；1=Nodular；0=None
#新增一列，分为Stage12，Stage34
rt$TNM2 <- rt$TNMStage
rt$TNM2 <- ifelse(rt$TNM2 == 'Stage I',12,
                   ifelse(rt$TNM2 =='Stage II',12,34)) 

#新增一列，分为Stage1、Stage2、Stage3
rt$TNM3 <- rt$TNMStage
rt$TNM3 <- ifelse(rt$TNM3 == 'Stage I',1,
                  ifelse(rt$TNM3 =='Stage II',2,3)) 

#新增一列，分为ETE\None
rt$ETE <- rt$ExtrathyroidalExtension
rt$ETE <- ifelse(rt$ETE == 'None',"None","ETE") 

#### 临床参数的亚组提取 ####
#Age
ageless55 <- rt[rt$Age2=='<55',c('RFS','Recurrence','SEZ')]   #把age<55的行提取出来，生存曲线(RFS, Recurrence) ~ RiskGroup
agemore55 <- rt[rt$Age2=='>=55',c('RFS','Recurrence','SEZ')]

#Gender
female <- rt[rt$Gender=='Female',c('RFS','Recurrence','SEZ')]
male <- rt[rt$Gender=='Male',c('RFS','Recurrence','SEZ')]

#Tstage
T4 <- rt[rt$`T Stage`=='T4',c('RFS','Recurrence','Risk')]

#Nstage
N0 <- rt[rt$NStage=='N0',c('RFS','Recurrence','SEZ')]
N1 <- rt[rt$NStage=='N1',c('RFS','Recurrence','SEZ')]

#Mstage
M1 <- rt[rt$`M Stage`=='M1',c('RFS','Recurrence','RiskCO')]

#TumorDiamater
TD4_less <- rt[rt$TD4=='1', c('RFS','Recurrence','SEZ')]
TD4_more <- rt[rt$TD4=='2', c('RFS','Recurrence','SEZ')]

TD1_less <- rt[rt$TD1=='1', c('RFS','Recurrence','SEZ')]
TD1_more <- rt[rt$TD1=='2', c('RFS','Recurrence','SEZ')]

# ETE
ETE0 <- rt[rt$ETE=='None',c('RFS','Recurrence','SEZ')]
ETE1 <- rt[rt$ETE=='ETE',c('RFS','Recurrence','SEZ')]

#Histopathology
histo1 <- rt[rt$`Histological Type`=='Classical',c('RFS','Recurrence','SEZ')]
histo2 <- rt[rt$`Histological Type`=='Follicular variant',c('RFS','Recurrence','SEZ')]
histo3 <- rt[rt$`Histological Type`=='Aggressive variants',c('RFS','Recurrence','SEZ')]

#Combined Disease
CD0 <- rt[rt$`Combined Disease2`=='0',c('RFS','Recurrence','RiskCO')]
CD1 <- rt[rt$`Combined Disease2`=='1',c('RFS','Recurrence','RiskCO')]
CD2 <- rt[rt$`Combined Disease2`=='2',c('RFS','Recurrence','RiskCO')]

HT0 <- rt[rt$HT=='0',c('RFS','Recurrence','Risk')]

#Multifocality
Focal0 <-rt[rt$FocusType=='Unifocal',c('RFS','Recurrence','SEZ')]
Focal1 <-rt[rt$FocusType=='Multifocal',c('RFS','Recurrence','SEZ')]

#BRAF
BRAF0 <- rt[rt$BRAF=='Wild-Type',c('RFS','Recurrence','SEZ')]
BRAF1 <- rt[rt$BRAF=='Mutation',c('RFS','Recurrence','SEZ')]

#TERT
TERT0 <- rt[rt$TERT=='Wild-Type',c('RFS','Recurrence','SEZ')]
TERT1 <- rt[rt$TERT=='Mutation',c('RFS','Recurrence','SEZ')]

#RAS
RAS0 <- rt[rt$RAS=='Wild-Type',c('RFS','Recurrence','SEZ')]
RAS1 <- rt[rt$RAS=='Mutation',c('RFS','Recurrence','SEZ')]

#Gene Fusion
GF0 <- rt[rt$GeneFusion=='Absence',c('RFS','Recurrence','SEZ')]
GF1 <- rt[rt$GeneFusion=='Presence',c('RFS','Recurrence','SEZ')]

#TNM stage
TNM12 <- rt[rt$TNM2 =='12', c('RFS','Recurrence','SEZ')]
TNM34 <- rt[rt$TNM2 =='34', c('RFS','Recurrence','SEZ')]

#TNM stage = 3组
TNM1 <- rt[rt$TNM3 =='1', c('RFS','Recurrence','SEZ')]
TNM2 <- rt[rt$TNM3 =='2', c('RFS','Recurrence','SEZ')]
TNM3_34 <- rt[rt$TNM3 =='3', c('RFS','Recurrence','SEZ')]

#ATA
ATA1 <- rt[rt$ATA=='Low', c('RFS','Recurrence','SEZ')]
ATA2 <- rt[rt$ATA=='Intermediate', c('RFS','Recurrence','SEZ')]
ATA3 <- rt[rt$ATA=='High', c('RFS','Recurrence','SEZ')]


#### 生存曲线分析（亚组间）####
fit <- survfit(Surv(RFS, Recurrence) ~ SEZ, data = N1)  #⭐️  #拟合生存曲线，创建生存对象

# 统计值
diff=survdiff(Surv(RFS, Recurrence) ~ SEZ, data = N1)  #⭐️
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){pValue=paste0("P= ",format(pValue, scientific = TRUE))}else{pValue=paste0("P= ",sprintf("%.03f",pValue))}    #pValue="p<0.001"
HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
HR = sprintf("%.02f",HR)
up95 = sprintf("%.02f",up95)
low95 = sprintf("%.02f",low95)
printHR = c(HR, low95, up95)
printHR 


######❤️绘制生存曲线#####
mycol = c("#337cba","#e11a0c","#48af45")
mycol0 = c("#e11a0c","#337cba")

mycol2 = c("#e11a0c","gray")
mycol3 = c("#337cba","gray")
mycol4= c("#48af45","gray")


pdf('ETEnew_.pdf',width =3.5, height =3.5, onefile=F)
ggsurvplot(fit, 
           data = rt,   #⭐️
           conf.int = F,   #增加置信区间
           pval = pValue,
           pval.size= 5,
           pval.coord = c(0.2,0.7),
           ggtheme = theme_test(base_size = 12, base_line_size = 0.2,base_rect_size = 0.8) +
             theme(plot.title=element_text(hjust=0.5)),  #theme_pubr()
           title = paste0("ATA,","HR = ",HR,"[",low95,"-",up95,"]"),   #⭐️
           #legend.title = 'ETE',   # 设置图例标题
           #legend.labs=c("ETE","None"),   # 指定图例分组标签   #⭐️     #'High risk','Low risk'   
           palette = mycol,  #⭐️    #'blue2','red2'
           legend = c(0.7, 0.3),  # 指定图例位置:"top"(默认),"bottom","left","right","none"，也可以用数字向量c()指定
           font.legend = 12,
           font.main = c(12, "plain","black"),
           font.x = c(12,"plain", "black"),
           font.y = c(12,"plain", "black"),
           font.tickslab = c(12, "plain", 'black'),
           xlab="Time(years)",
           ylab="TCFi probability",
           ylim = c(0.5,1),
           break.time.by = 2
           ) 
  
  
dev.off()  


