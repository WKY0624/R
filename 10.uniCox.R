library(survival)  

#### UniCox regression analysis ####

coxPfilter=0.2 # ⭐ The significance filtering standard can be changed to 0.01 if there are too many results; If there are too few results, it can be changed to 0.1 or 0.2
inputFile="train.expTime.txt"   
inputFile2="test.expTime.txt" 

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365   #⭐️

#Cycle genes to identify prognostic related genes
outTab=data.frame()
sigGenes=c("DSS","DSEvent")   #⭐️
for(i in colnames(rt[,3:ncol(rt)])){
  #Cox analysis
  cox <- coxph(Surv(DSS, DSEvent) ~ rt[,i], data = rt)   #⭐️
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  #Preserve prognostic related genes
  if(coxP<coxPfilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}


#Output the results of single factor Cox regression analysis
write.table(outTab,file="train.uniCox.txt",sep="\t",row.names=F,quote=F)  

#Output the expression levels of significant genes in single factor Cox regression analysis
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="train.uniSigExp.txt",sep="\t",row.names=F,quote=F)  


#### Draw a forest plot ####
library(forestplot) 
library(readxl) 

options(scipen = 0)   
options(digits = 3)  

rs_forest1 <- read.table('total.uniCox.txt', header = T, check.names = F, sep = '\t')
rs_forest1<-data.frame(rs_forest1)
rs_forest1 <- rs_forest1[order(rs_forest1$pvalue, decreasing =F),]
rs_forest1$`95%CI` <- paste0('[', round(rs_forest1$HR.95L,3) ,'-', round(rs_forest1$HR.95H,3), ']')
rs_forest1$P.value <- ifelse(rs_forest1$pvalue<0.001,format(rs_forest1$pvalue, scientific = T, digits = 3),
                             sprintf("%.03f", rs_forest1$pvalue))

#Add a header to the data column
label<-cbind(c("Gene", rs_forest1$id), 
             c("HR",round(rs_forest1$HR,3)),
             c("95% CI", rs_forest1$`95%CI`),
             #c("LL",rs_forest1$HR.95L),
             #c("UL",rs_forest1$HR.95H),
             c("P value", rs_forest1$P.value))

forestplot(labeltext = label, 
           mean = c(NA,rs_forest1$HR),
           lower = c(NA,rs_forest1$HR.95L), 
           upper = c(NA,rs_forest1$HR.95H), 
           is.summary=c(F,F,F,F,F
                        #This parameter takes a logical vector to define whether each row in the data is a summary value. If so, it is set to True at the corresponding position. If not, it is set to FALSE; Rows set to True appear in bold
           ), 
           
           ### Set Table           
           lineheight = unit(4,'mm'),    
           colgap = unit(5,'mm'),        
           line.margin=unit(6, 'mm'),
           graphwidth = unit(0.3,"npc"), 
           hrzl_lines = list("2" = gpar(lty = 1, lwd = 1.5, col="black")),
           
           ### Reference line
           zero = NA, 
           lwd.zero = 1, 
           grid = structure(1, gp = gpar(col = "black", lty=2, lwd=1)), 
           
           ### Set Forest Plot
           graph.pos = 2,     
           lty.zero = '',
           lty.ci = 'solid',  
           lwd.ci = 1.3,    
           ci.vertices = T,   
           ci.vertices.height = 0.15,   
           fn.ci_norm="fpDrawDiamondCI",    #The shape of the box: fpDrawNormalCI、fpDrawDiamondCI、fpDrawCircleCI、fpDrawPointCI、fpDrawSummaryCI、fpDrawBarCI
           boxsize = 0.3,   
           
           ### Set X-axis
           clip=c(0,7),  #After the X-axis exceeds, the right end becomes an arrow
           lwd.xaxis = 1.5,     
           lwd.xticks = 0.5,
           xticks = c(0,1,2,3,4),   
           xlab="Hazard ratio", 
           xlog= F, 
           
           ### Set Color
           col=fpColors(box = 'firebrick1',  #箱体的颜色：  蓝 "#337cba",  红 "#e11a0c",
                        lines = "black",   #置信区间端点线条的颜色
                        summary = "yellow",   #文字的颜色
                        zero = "grey",   #无效线的颜色
                        text = "black",    #文本的颜色
                        axes = "black",    #x轴的颜色
                        #vrtcl_lines = "green",
                        hrz_lines = "black",   #设置线条的颜色
           ),
           
           ### Set Font
           txt_gp=fpTxtGp(label=gpar(fontfamily="",cex=1),    # 整个表格
                          ticks=gpar(fontfamily="",cex=0.8),    # X轴刻度标签
                          xlab=gpar(fontfamily="",cex=1),     # X轴
                          title=gpar(fontfamily="",cex=1)     # 标题
           ), 
)

#export: 4 x 7  


