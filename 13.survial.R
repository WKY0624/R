library(survival)
library(survminer)

#绘制生存曲线函数
bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(DSS, DSEvent) ~Risk,data = rt)  #⭐️
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue= paste0("P = ", format(pValue, scientific = TRUE, digits = 3))
  }else{
    pValue=paste0("P = ",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(DSS, DSEvent) ~ Risk, data = rt)  #⭐️  #拟合生存曲线，创建生存对象
  
  HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
  HR = sprintf("%.02f",HR)
  up95 = sprintf("%.02f",up95)
  low95 = sprintf("%.02f",low95)
  
  #绘制生存曲线
  mycol = c("#e11a0c","#337cba")   
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,   #增加置信区间
                     conf.int.alpha = 0.1,
                     conf.int.style = 'ribbon',
                     palette = mycol,   #palette=c("red", "blue"),   #自定义调色板可选调色板有 "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons"和"rickandmorty"
                     ### P值
                     pval = F,  
                     pval.size = 4,   
                     pval.coord=c(2,0.25),   
                     pval.method = TRUE,  
                     pval.method.size = 4, 
                     pval.method.coord=c(2,0.35),   
                     ### 主题
                     ggtheme = theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
                       theme(axis.title.y = element_text(size = 10),
                             axis.title.x = element_text(size = 10),
                             axis.text.y = element_text(size = 8),
                             axis.text.x = element_text(size = 9),
                             legend.title = element_text(size = 10),
                             legend.text = element_text(size = 10),
                             axis.ticks = element_line(size = 0.3)
                             ),
                     ### 图例
                     legend.title="",
                     legend.labs=c("High-risk", "Low-risk"),   # 指定图例分组标签
                     legend = c(0.8, 0.2),   #图示的位置
                     font.legend = 10,
                     ### 坐标轴
                     xlab="Follow-up time (year)",
                     ylab="Disease-specific survival probability",  #⭐️
                     break.time.by = 1,
                     ylim = c(0.5,1),
                     #surv.scale = "percent", #纵坐标以百分数形式呈现
                     )
  surPlot$plot <- surPlot$plot + 
    annotate("text",x = 0.2, y = 0.60, hjust = 0, fontface = 1, label = paste0("HR = ",HR,"[",low95,"-",up95,"]")) +
    annotate("text",x = 0.2, y = 0.55, hjust = 0, fontface = 2, label = pValue) 
    
  
  
  pdf(file=outFile,onefile = FALSE,width = 3.5, height = 3.5)
  print(surPlot)
  dev.off()
}


#AI4Genes：
bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")
bioSurvival(inputFile="totalRisk.txt", outFile="totalSurv.pdf")
