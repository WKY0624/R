library(survival)
library(survminer)
library(timeROC)

inputFile="totalRisk.txt"  #12
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

#### ROC曲线 ####
ROC_rt=timeROC(T=rt$DSS, delta=rt$DSEvent,    #⭐️
               marker=rt$SEZ6L2, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

mycol = c("#48af45","#337cba","#e11a0c")

pdf(file="TCGA.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=mycol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=mycol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=mycol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',   #x=0.4,y=0.3,
       c(paste0('1-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
         paste0('3-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[2]), ')'),
         paste0('5-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[3]), ')')),
       col=mycol, lwd=2, bty = 'n')   
dev.off()

#### KM曲线 ####
## SEZ6L2表达数据
TCGA = read.table("PyroExp.txt", header = T, check.names = F, sep = "\t", row.names = 1)  #02
TCGA = TCGA["SEZ6L2",]
group=sapply(strsplit(colnames(TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
TCGA = TCGA[,group==0]   #删除正常样本
colnames(TCGA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(TCGA))   #规范TCGA样本名称
TCGA = as.data.frame(t(TCGA))
TCGA$SEZ6L2 = log2(TCGA$SEZ6L2+1)
#根据SEZ6L2中位数分为高表达组和低表达组
TCGA$group = ifelse(TCGA$SEZ6L2 > median(TCGA$SEZ6L2),"High","Low")

## 预后数据
cli = read.table("04.time_TCFi.txt", header=T, sep="\t", check.names=F, row.names = 1)
cli$DSS = cli$DSS/365

## 合并数据
samSample = intersect(row.names(TCGA),row.names(cli))
rt = cbind(SEZ = TCGA[,"group"],cli)

## 绘制生存曲线
fit <- survfit(Surv(DSS, Event) ~ rt$SEZ, data = rt, #拟合生存曲线，创建生存对象
               type = "kaplan-meier", error = "greenwood",
               conf.type = "plain", na.action = na.exclude) 

## 卡方检验
diff = survdiff(Surv(DSS, Event) ~ rt$SEZ, data = rt, na.action = na.exclude)  #⭐️
pValue = 1 - pchisq(diff$chisq, df = 1)
if(pValue < 0.001){pValue = paste0("P= ",format(pValue, scientific = TRUE))}else{pValue = paste0("P= ",sprintf("%.03f",pValue))}

HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
HR = sprintf("%.02f",HR)
up95 = sprintf("%.02f",up95)
low95 = sprintf("%.02f",low95)
print(c(HR, low95,up95))

## 绘制曲线
mycol = c("#e11a0c","#337cba")
pdf('KMplot.pdf',width =3.5, height =3.5, onefile=F)
ggsurvplot(fit, 
           data = rt,   
           conf.int = F,   #是否增加置信区间
           pval = pValue,
           pval.size= 5,
           pval.coord = c(0.2,0.7),
           ggtheme = theme_test(base_size = 12, base_line_size = 0.2,base_rect_size = 0.8) +
             theme(plot.title=element_text(hjust=0.5)),  
           title = "SEZ6L2 for TCFi", 
           legend.title = 'SEZ6L2',  
           legend.labs=c("High","Low"),  
           palette = mycol, 
           legend = c(0.7, 0.3),  # 指定图例位置:"top"(默认),"bottom","left","right","none"，也可以用数字向量c()指定
           font.legend = 12,
           font.main = c(12, "plain","black"),
           font.x = c(12,"plain", "black"),
           font.y = c(12,"plain", "black"),
           font.tickslab = c(12, "plain", 'black'),
           xlab="Time(years)",
           ylab="TCFi probability",
           ylim = c(0.5,1),
           break.time.by = 2,  # 设置x轴刻度间距
           ## 风险表
           risk.table= T,   #"absolute"、"percentage"、"abs_pct"",
           #risk.table.col = 'strata',   #按组更改风险表颜色
           risk.table.y.text.col = T,    #颜色风险表文本注释（按层）
           risk.table.y.text = T, #在风险表图例中的文本注释中显示条形而不是名称
           risk.table.height = 0.2,
           fontsize = 3.5,   #风险表字体大小
           ##积累事件表
           cumevents=T,
           cumevents.height = 0.2,  #占20%
           tables.theme = theme_test(base_size = 8, base_line_size = 0.2, base_rect_size = 0.4)+
             theme(axis.text.y = element_text(size = 10, face = 'bold'),
                   axis.text.x = element_text(size = 10, color = 'black'),
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks = element_line(size = 0.3))
           )

dev.off()  