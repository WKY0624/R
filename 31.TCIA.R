library(ggpubr) 
library(ggplot2)
library(gghalves)
library(reshape)  #melt
library(rstatix)   # add_xy_position

tciaFile="./DataSets/TCIA-THCA_230310_tidy.tsv"     #TCIA打分文件：https://tcia.at/ （整理：仅保留4行，删除NA样品，新建TCIA.txt）
riskFile="totalRisk.txt"     #12

## 读取数据
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(ips), row.names(Risk))
ips=ips[sameSample, , drop=F]
Risk=Risk[sameSample, "Risk", drop=F]
data=cbind(Risk,ips)
high  = data[data$Risk=='High',]
low = data[data$Risk=='Low',] 

##计算组间差异
##正态性检验：经检验数据为非正态分布(P<0.05)，后续使用wilcox.test
#shapiro.test(rt$TIDE)      #样本量小于50时选择W检验：W接近1，P值大于0.05，数据为正态分布
ks.test(scale(data$ips_ctla4_neg_pd1_pos),'pnorm')   #样本量大于50时选择D检验：P小于显著性水平α(0.05)，则拒绝H0（p越大越好）；D值越小，越接近0，越接近正态分布(D越小越好)
ks.test(scale(low$ips_ctla4_pos_pd1_neg),'pnorm')  

##计算p值
#定义执行一次秩和检验的操作流程
wilcox.test(ips_ctla4_neg_pd1_pos ~ Risk, data = data, alternative = 'two.side')
#能够执行一次那么就可以执行多次，批量执行秩和检验
stat <- data.frame(TCIA=colnames(data)[2:ncol(data)])
for (i in 2:ncol(data)){
  print(i)
  stat[i-1,2] <- wilcox.test(data[,i] ~ Risk, data = data, 
                              alternative = 'two.side', paired= F,
                              exact = FALSE)[["p.value"]]}

##计算fdr
stat$fdr <- p.adjust(stat$V2, method = "fdr")
stat$bon <- p.adjust(stat$V2, method = 'bonferroni')

##输出结果
colnames(stat) <- c("TCIA","P.value","FDR","Bonferroni")

##构建绘图数据
data_long = melt(data, id.vars = ("Risk"))
colnames(data_long) = c("RiskGroup", "TCIA", "Score")

##获取统计数值和位置
data_long <- data_long %>% group_by(TCIA)
data_long$Score <- as.numeric(data_long$Score)

stat.data <- pairwise_wilcox_test(data = data_long, Score~RiskGroup, paired = F, alternative = 'two.side') %>% 
  add_xy_position(x = 'TCIA')
stat.data$p.scient <- format(stat.data$p.adj, scientific = TRUE)
stat.data$p.round3 <- round(stat.data$p.adj,3)
stat.data$Pvalue <- ifelse(stat.data$p<0.001, sprintf(stat.data$p.scient), sprintf("%.03f", stat.data$p))

##绘图
mycol =  c("#e11a0c","#337cba")
ggplot()+
  geom_half_violin(data = data_long %>% filter(RiskGroup == "High"),
                   aes(x = TCIA, y = Score), side = "l",size= 0,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,nudge = 0.01,    #"#fb3e35"
                   position = position_dodge(width = 0)
  ) +
  geom_half_violin(data = data_long %>% filter(RiskGroup == "Low"),
                   aes(x = TCIA,y = Score), side = "r", size= 0,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1.1, nudge = 0.01,   #'#564efd'
                   position = position_dodge(width = 0)
  ) +
  #中间添加：分半箱线图
  geom_half_boxplot(data = data_long %>% filter(RiskGroup == "High"),
                    aes(x = TCIA, y = Score), width = 0.2, 
                    colour=mycol[1], #lwd= 0.4, 
                    outlier.shape = NA, outlier.size = 0.8, outlier.stroke = T,
                    fill=mycol[1],side = "l", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data_long %>% filter(RiskGroup == "Low"),
                    aes(x = TCIA, y = Score), width = 0.2, 
                    colour=mycol[2], #lwd= 0.4, 
                    outlier.shape = NA, outlier.size = 0.8, outlier.stroke = T,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1)
  ) +
  #添加平均数白点
  #geom_point(data = data2, aes(x = scoreType, y = Score, fill = RiskGroup), stat = 'summary', fun=mean, col= 'white', pch = 3, size =1, position = position_dodge(width = 0.2), show.legend = FALSE) +
  #添加竖线
  #geom_vline(xintercept = c(1.5, 2.5), color = "#bcbdbf", alpha = 0.8, size = 0.3) +
  #添加图例
  geom_line(data = data_long, aes(x = TCIA, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, size =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="PyroGroup",labels = c("High","Low")) +
  #显著性差异
  stat_pvalue_manual(
    stat.data, 
    y.position = c(10.5,10.5,10.5,10.5),
    label = 'P = {Pvalue}\n{p.adj.signif}',   #"p.adj.signif", 
    bracket.size = 0.3, # 粗细
    bracket.shorten = 0.15,   #宽度
    tip.length = 0.01,  # 两边竖线的长度
  ) +
  xlab("") +
  ylab("TCIA Score") +
  ylim(4,11) + 
  scale_x_discrete(labels = c('IPS','PD1 blocker','CTLA4 blocker','PD1+CTLA4 blocker'))+
  theme_classic(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5)+
  theme(
    #x轴标签
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),    
    #轴刻度
    axis.ticks =element_line(linewidth = 0.3),
    #y轴标签
    axis.text.y = element_text(color = 'black', hjust = 1, # 左对齐
                               size = 10, lineheight = 1),
    #标题居中:
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    #图例：
    legend.position = "top",   #图例位置
    legend.key.size = unit(10,'pt'), #图例大小
    legend.justification = "centre")

ggsave('TCIA.pdf', height = 4.5, width = 6)
