
library(survival)
library(survminer)
library(ggsignif) 
library(ggpubr)
ClusterFile="cluster.txt"  
cliFile="./DataSets/04.time_OS.txt"  #TCFi/PFS/OS

#Read input file
Cluster=read.table(ClusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("RFS", "Recurrence")   #⭐️
cli$RFS=cli$RFS/365  #⭐️

#Data merging
sameSample=intersect(row.names(Cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], Cluster[sameSample,,drop=F])


##Survival analysis ####
fitd <- survdiff(Surv(RFS, Recurrence) ~ Cluster, data = rt, na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(RFS, Recurrence) ~ Cluster, data = rt,
               type = "kaplan-meier", error = "greenwood",
               conf.type = "plain", na.action = na.exclude)

##Paired survival analysis####
ps <- pairwise_survdiff(Surv(RFS, Recurrence) ~ Cluster, data = rt, p.adjust.method = "none") # Correction is not used here. If correction is needed, none can be replaced with BH（“holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”）

##Draw survival curve####
mycol = c("#e11a0c","#48af45","#337cba")
length=length(levels(factor(rt$Cluster)))
mycol = mycol[1:length]

p <- ggsurvplot(fit, 
           data=rt,
           conf.int=F,
           conf.int.alpha = 0.1,
           conf.int.style = 'ribbon',
           linetype = 1,
           palette = mycol,
           #### P value
           pval = F,  
           pval.size = 4,  
           pval.coord=c(10,0.25),   
           pval.method = TRUE,    
           pval.method.size = 4, 
           pval.method.coord=c(10,0.35),  
           #### Theme
           ggtheme = theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
             theme(axis.title.y = element_text(size = 10),
                   axis.title.x = element_blank(),
                   axis.text.y = element_text(size = 9, color='black'),
                   axis.text.x = element_text(size = 0),
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   axis.ticks = element_line(size = 0.3)
                   ),
           #### Legend
           legend.title="PyroCluster",
           legend.labs=levels(factor(rt[,"Cluster"])),
           legend = "none",
           font.legend = 0,
           #### Coordinate axis
           xlab="Follow-up years",
           ylab = 'OS probability',  #⭐️
           break.time.by = 3,
           xlim = c(0,15.5),
           ylim = c(0.5,1),
           #### Risk table
           risk.table= T,   #"absolute"、"percentage"、"abs_pct"",
           #risk.table.col = 'strata',  #Change risk table color by group
           risk.table.y.text.col = T,   #Color Risk Table Text Annotations (by Layer)
           risk.table.y.text = T, 
           risk.table.height = 0.3,
           fontsize = 3.5,   
           #### Accumulate event table
           cumevents=F,
           cumevents.height = 0.3,
           tables.theme = theme_test(base_size = 8, base_line_size = 0.2, base_rect_size = 0.4)+
             theme(axis.text.y = element_text(size = 10, face = 'bold'),
                   axis.text.x = element_text(size = 10, color = 'black'),
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks = element_line(size = 0.3),
                   )) 

##Add overall p value####
p.lab <- paste0("P", ifelse(p.val < 0.001, paste0(" = ", format(p.val, scientific = TRUE, digits = 3)), paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + 
  #total P
  annotate("text",x = 0.2, y = 0.60, hjust = 0, fontface = 1, label = paste0("Log-rank test")) +
  annotate("text",x = 0.2, y = 0.55, hjust = 0, fontface = 2, label = p.lab)   
  #separate P
  annotate("text",x = 11, y = 0.84,  #PFS：y = 0.77, #DSS：y = 0.84,  #⭐️
           hjust = 0, size = 3, color = mycol[3], label = paste0("P (C2 vs. C3) = ", round(ps$p.value[2,2],3))) + 
  annotate("text",x = 11, y = 0.72,  #PFS：y = 0.68, #DSS:y = 0.72,   #⭐️
           hjust = 0, size = 3, color = mycol[2], label = paste0("P (C1 vs. C2) = ", round(ps$p.value[1,1], 3))) +
  annotate("text",x = 11, y = 0.68, #PFS：y = 0.64,#DSS:y = 0.68,    #⭐️
           hjust = 0, size = 3, color = mycol[3], label = paste0("P (C1 vs. C3) = ", round(ps$p.value[1,2], 3)))

p


## Output Image####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("cluster_OS.pdf", width = 3.5, height = 4)
print(p)
dev.off()


