#install.packages("ggrisk")
library(ggrisk)
library(rms)    
library(pheatmap)

#####🐰方法3：分开画（ggplot2）####
inputFile = "totalRisk.txt"
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
rt=rt[order(rt$RiskScore),]      #根据病人风险得分对样品进行排序
rt$patient_id = seq(nrow(rt))  # 添加一个表示病例顺序的新列
riskClass=rt[,"Risk"]
lowLength=length(riskClass[riskClass=="Low"])
highLength=length(riskClass[riskClass=="High"])
lowMax=max(rt$RiskScore[riskClass=="Low"])
rt$RiskScore2 <- rt$RiskScore
rt$RiskScore2[rt$RiskScore2 > 5] = 5
rt$color <- c(rep("#337cba", lowLength), rep("#e11a0c", nrow(rt) - lowLength))
## 绘制风险曲线
p1 <- ggplot(rt, aes(x = rt$patient_id, y = rt$RiskScore2)) +
  geom_point(aes(color = I(color)), shape = 20, size = 1.5) +
  scale_color_manual(values = c("#337cba", "#e11a0c"), labels = c("Low-risk", "High-risk")) +
  labs(x = "Patients (increasing PyroScore)", y = "PyroScore") +
  theme_bw() +
  #theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  # 添加一条水平的虚线，表示lowMax
  geom_hline(yintercept = lowMax, linetype = "dashed") +
  # 添加一条垂直的虚线，表示lowLength
  geom_vline(xintercept = lowLength, linetype = "dashed") + 
  # 将图例放置在主图的右侧，且距离不太远
  theme(legend.position = "right", legend.box.spacing = unit(0.2, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))+
  # 修改图例标题和标签
  guides(color = guide_legend(title = "PyroGroup", override.aes = list(shape = 20, size=2.5)
                              )) +
  # 设置背景
  theme(axis.title = element_text(size=10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 10,face="bold")) +
  ggtitle("Entire cohort")    #⭐️
  
## 绘制生存状态图
p2 <- ggplot(rt, aes(x = rt$patient_id)) +
  geom_point(aes(y = rt$DSS, color = as.factor(rt$DSEvent)),shape = 20, size = 2.5,alpha = 1) +
  scale_color_manual(values = c("#337cba","#e11a0c"), labels = c("Censored","Event")) +
  labs(x = "Patient (increasing PyroScore)", y = "Follow-up time (year)") +
  geom_vline(xintercept = lowLength, linetype = "dashed") +
  theme_bw() +
  #theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  guides(color = guide_legend(title = "Status",override.aes = list(shape = 20, size=2.5))) +
  theme(legend.position = "right", legend.box.spacing = unit(0.2, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        axis.title = element_text(size=10)) 
    
## 拼接图片
library(patchwork)
p1 + p2 + 
  plot_layout(nrow = 2,   #图像设置为2列，默认按列填充    # nrow=  按行填充
              heights = c(1, 1),   #两列之间相对宽度比为3：1   # heights=c(2,1)  相对高度
              guides = "keep") &
  theme(legend.position = 'right',
        legend.box.spacing = unit(0.1,"cm"),
        legend.spacing = unit(0,'cm'),
        legend.justification = "centre",
        legend.key.size = unit(0.4,"cm")) 

ggsave("total.pdf",height = 4,width = 6)  #⭐️


#####方法2：分开画（ggrisk）####
mydata = read.table("totalRisk.txt",header=T, sep="\t", check.names=F, row.names=1)
mydata = mydata[,-grep("RiskScore|Risk",colnames(mydata))]

colnames(mydata)

fit <- cph(Surv(DSS,DSEvent)~SEZ6L2+PRDM1+CXCL8+GJA1+TRAF6+H2BC8+PYCARD+IFI27+SIGLEC15, mydata)

ggrisk(fit,   #模型名称
       cutoff.x = 145,    #cutoff标签位置
       cutoff.y = -0.8,   #cutoff标签位置
       cutoff.show = F,
       cutoff.value = "median",   #cutoff的选择：median、roc、cutoff(最小p值方法选择最优切点)、自定义
       code.0 = "No",
       code.1 = "Yes",
       title.A.ylab = "Risk score",
       title.B.ylab = "DSS (years)",
       title.A.legend = "Risk group",
       title.B.legend = "Disease-specific event",
       title.C.legend = "Expression",
       color.A = c(low = "#337cba", high = "#e11a0c"),
       color.B = c(code.0 ="#337cba" ,code.1 = "#e11a0c"),
       color.C = c(low = "blue", median = "white", high = "red"),   #热图的颜色
       size.Ctext = 11,
       relative_heights = c(1,1,0.05,1),   #ABC图的高度比例
       #expand.x = 3,
)  


#####方法1：3张一起画####
#定义风险曲线的函数
bioRiskPlot = function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
  rt=rt[order(rt$RiskScore),]      #根据病人风险得分对样品进行排序
  
  #绘制风险曲线
  riskClass=rt[,"Risk"]
  lowLength=length(riskClass[riskClass=="Low"])
  highLength=length(riskClass[riskClass=="High"])
  lowMax=max(rt$RiskScore[riskClass=="Low"])
  line=rt[,"RiskScore"]
  line[line>10]=10
  pdf(file=paste0(project, ".RiskScore.pdf"), width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing PyroScore)",
       ylab="PyroScore",
       col=c(rep("#337cba",lowLength),rep("#e11a0c",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("right", 
         c("High","Low"),
         bty="n",pch=20,col=c("#e11a0c","#337cba"),cex=0.8)  #⭐️
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$DSEvent)
  color[color==1]="#e11a0c"
  color[color==0]="#337cba"
  pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
  plot(rt$DSS, pch=20,    #⭐️
       xlab="Patients (increasing PyroScore)",
       ylab="Follow-up time (years)",
       col=color)
  legend("right", 
         c("Yes", "No"),   #⭐️
         bty="n",pch=20,col=c("#e11a0c","#337cba"),cex=0.8)  #⭐️
  abline(v=lowLength,lty=2)
  dev.off()
  
  #定义热图注释的颜色
  ann_colors=list()
  bioCol=c("#337cba","#e11a0c")
  names(bioCol)=c("Low", "High")
  ann_colors[["Risk"]]=bioCol
  
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(Risk=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation,
           annotation_colors = ann_colors, 
           cluster_cols = F,
           cluster_rows = T,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("#337cba",3), "white", rep("#e11a0c",3)))(100),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

#调用函数，绘制风险曲线
#bioRiskPlot(inputFile="trainRisk.txt", project="train")
#bioRiskPlot(inputFile="testRisk.txt", project="test")
bioRiskPlot(inputFile="totalRisk.txt", project="total")




