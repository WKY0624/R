#### 读取合并数据 ####
## TCGA表达数据
TCGA = read.table("PyroExp.txt", header = T, check.names = F, sep = "\t", row.names = 1)  #02
TCGA = TCGA["SEZ6L2",] 
group=sapply(strsplit(colnames(TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
TCGA = TCGA[,group==0]   #删除正常样本
colnames(TCGA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(TCGA))   #规范TCGA样本名称
TCGA = as.data.frame(t(TCGA))
TCGA$SEZ6L2 = log2(TCGA$SEZ6L2+1)

## 临床数据
clin = read.table("06.clinical477.txt",header=T, sep="\t", check.names=F, row.names = 1)
merge = merge(TCGA, clin, by = "row.names", all = TRUE)
rownames(merge) = merge[,1]
merge = merge[,-1]
#write.table(merge, "mergeSEZ.txt", sep = "\t")

#### 提取参数 ####
# Demography
rt1 = merge[,c("SEZ6L2","Age","Gender","RadiationHistory","CombinedDisease")]

# Pathology
rt2 = merge[,c("SEZ6L2","Tcategory2","Ncategory","Mcategory","HistologicalType","ETE","Multifocality")]
rt2$ETE = ifelse(rt2$ETE == "None","None","ETE")

# Molecular
rt3 = merge[,c("SEZ6L2","BRAF","RAS","TERT","GeneFusion")]
rt3$BRAF = ifelse(rt3$BRAF == "Mutation","BRAF","non-BRAF")
rt3$RAS = ifelse(rt3$RAS == "Mutation","RAS","non-RAS")
rt3$TERT = ifelse(rt3$TERT == "Mutation","TERT","non-TERT")

# Prognosis
rt4 = merge[,c("SEZ6L2","RAIResponse","NewTumorEvent","TNM","ATA")]
rt4$NewTumorEvent = ifelse(rt4$NewTumorEvent == "None","None","TCFi")


#### 处理长数据 ####
library(tidyverse)
# 使用pivot_longer函数转换为长数据
#long_data <- pivot_longer(rt1, cols = c(Age, Gender,RadiationHistory,CombinedDisease), names_to = "class", values_to = "class_temp")
#long_data <- pivot_longer(rt2, cols = c(Tcategory2, Ncategory,Mcategory,HistologicalType,ETE,Multifocality), names_to = "class", values_to = "class_temp")
#long_data <- pivot_longer(rt3, cols = c(BRAF,RAS,TERT,GeneFusion), names_to = "class", values_to = "class_temp")
long_data <- pivot_longer(rt4, cols = c(RAIResponse,NewTumorEvent,TNM,ATA), names_to = "class", values_to = "class_temp")

# 合并class和PyroScore1列
result_data <- cbind(long_data$SEZ6L2, long_data$class_temp)  

# 重新命名列
colnames(result_data) <- c("SEZ6L2", "class")
# 打印结果
result_data = as.data.frame(result_data)
result_data$SEZ6L2 = as.numeric(result_data$SEZ6L2)


#### 绘图 #####
library(ggplot2)
library(ggpubr)   #统计值

## 设置期望的 X 轴顺序
#custom_order <- c("T1~T2","T3~T4", "N0", "N1", "M0", "M1","Classical","Follicular", "Aggressive","None","ETE","Unifocal","Multifocal")
#custom_order = c("non-BRAF","BRAF","non-RAS","RAS","non-TERT","TERT","Absence","Presence")
custom_order = c("Not received","Sensitive","Refractory","None","TCFi","Stage I","Stage II","Stage III~IV","Low","Intermediate","High")

result_data$class = factor(result_data$class, levels = custom_order)
result_data= na.omit(result_data)    #去除所有含有NA值的行


## 设置颜色
mycol <- c("tomato","tomato","tomato",
           "goldenrod1", "goldenrod1",
           "forestgreen", "forestgreen", "forestgreen",
           "steelblue", "steelblue","steelblue",
           "mediumpurple","mediumpurple","mediumpurple",
           "hotpink","hotpink",
           "orange","orange","orange",
           "lightpink","lightpink","lightpink","lightpink")

## 蜂群图+箱线图
library(ggbeeswarm)
ggplot(result_data, aes(x = class, y=SEZ6L2, color = class, fill = class)) + 
  geom_beeswarm(cex = 1, size = 2 , pch = 16, alpha = 0.3, stroke = 0, na.rm = TRUE,
                priority = "none",  #"descending", "density", "random" and "none"
  ) + 
  geom_boxplot(outlier.shape = NA, width = 0.3, notch = T, color="#333333", size = 0.7) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = alpha(mycol,0)) +
  
  # --------背景主题
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black',angle = 45, hjust = 1),
        legend.position = "none",
        #plot.margin = unit(c(2, 1,  1,  1), 'cm')
        ) +
  
  #---------设置坐标轴
  labs(y = "SEZ6L2", x = "") +  
  ylim(0,10) +
  
  # --------添加组间比较的统计值
  stat_compare_means(
    method = "wilcox.test", # 使用 Wilcoxon 秩和检验
    method.args = list(alternative = "two.side"),
    comparisons = list(
      #c("T1~T2","T3~T4"), c("N0","N1"), c("M0", "M1"), c("Classical","Follicular"),c("Classical","Aggressive"),c("Follicular","Aggressive"),c("None","ETE"), c("Unifocal", "Multifocal")
      #c("BRAF","non-BRAF"),c("RAS","non-RAS"),c("TERT","non-TERT"),c("Absence","Presence")
      c("Not received","Sensitive"),c("Not received","Refractory"),c("Sensitive","Refractory"),
      c("None","TCFi"),c("Stage I","Stage II"),c("Stage I","Stage III~IV"),c("Stage II","Stage III~IV"),
      c("Low","Intermediate"),c("Low","High"),c("Intermediate","High")
      #c("Not received","Sensitive","Refractory"),c("None","TCFi"),c("Stage I","Stage II","Stage III~IV"),c("Low","Intermediate","High")
      ),
    label = "p.format",
    #label = "p.signif", # 添加显著性标记
    hide.ns = F,
    step.increase = F
  )

ggsave("Prog.pdf",height = 3, width = 5)

