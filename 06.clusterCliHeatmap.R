library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#install.packages("gtools")
library(gtools)


#### Data reading and processing ####
#Reading data files
Cluster=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)  
cli=read.table("./DataSets/clinical477.txt",header = T,sep = "\t",check.names=F, row.names=1 )  

set.seed(666)
#Merge data
samSample=intersect(row.names(Cluster), row.names(cli)) 
cli=cli[samSample,,drop=F]  
Cluster=Cluster[samSample,,drop=F]   
Type=cbind(Cluster, cli)  
#Type=Type[,grep("Cluster|Age|Gender|Coexist Disease|Radiation History|Tstage|Nstage|Mstage|Histological Type|Multifocality|Extrathyroid Extension|BRAF|RAS|TERT|Tumor Size|Lymph Node Count|TNM|ATA|MACIS|EORTC",colnames(Type))] 
Type=Type[,-grep("T category_2|AJCC 8th TNM|ATA|MACIS|EORTC",colnames(Type))] 
Type=Type[order(Type$Cluster),,drop=F]  


#### Statistical analysis ####
# Chi square test, calculating the p-value and significance asterisk for each clinical trait
#p_values = c()
p_num=c()
p_sci=c()
for (cli in colnames(Type[,2:ncol(Type)])) {
  data=Type[c("Cluster", cli)]
  colnames(data)=c("Cluster", "cli")
  data=data[(data[,"cli"]!="NA"),]
  tableStat=table(data)  
  stat=fisher.test(tableStat,simulate.p.value=TRUE)    #fisher.test() or chisq.test()
  p <- stat[["p.value"]]
  p1 <- ifelse(p<0.001,format(p, scientific = TRUE,digits = 3), sprintf("%.03f", p))
  p2 <- ifelse(p<0.001,"<0.001", sprintf("%.03f", p)) 
  p_num[cli] <- p2
  p_sci[cli] <- p1
}

star <- ifelse(p_num<0.001,"***",ifelse(p_num<0.01,"**",ifelse(p_num<0.05,"*","")))
p_star <- paste0(p_sci," ", star)


# top_annotation
col_group <- Type$Cluster

# Row segmentation
row_group <- rep(c(""," ","   ", "    "), c(4, 8, 4, 3))

#### Set color scheme ####
## 1.The color scheme of the main image
m <- apply(Type[,-1], 2, function(x){length(unique(x))})    #Calculate the number of categories for each column
qz <- sort(unique(as.character(as.matrix(Type[,-1]))))   #Extract elements from data and sort them, qz=names of various elements
length(qz)  #View how many unique elements are in the data

colors <- c("#f9e5e9",    #1:Age<55
            "#95353bFF",  #2:Age>55
            "ghostwhite", #3:Size-1
            "ghostwhite", #4:LNCount-1
            "#c7d79bFF",  #5:LNCount-2   #edf0e5FF→ghostwhite
            
            "#c7d79bFF",  #6:Size-2
            "#75be70FF",  #7:LNCount-3
            "#237c23FF",  #8:LNCount-4
            "#237c23FF",  #9:Size-3
            "aliceblue",  #10:GeneFusion-1Absence
            
            "#298f29FF",  #11:Histo-1 Aggressive
            "ghostwhite",   #12:Vital-1-Alive
            "yellow3",  #13:Vital-2 Any cause death
            "ghostwhite", #14:Histo-2 Classical    #edf0e5FF
            "purple3",  #15：Vital-3=Disease specific death + NTEevent-1
            
            "blue3",  #16：⭐️NTEvent-1 Distant
            "#f9e5e9",    #17:Gender-1-Female
            "#c7d79bFF",  #18：Histo-3(Follicular variant) 
            "#298f29FF",  #31：ETE-2=Moderate

            "#d1979aFF",  #26：Coexist-1(Lymphocytic Thyroiditis)
            "ghostwhite", #27：Mstage-1=M0
            
            "#237c23FF",  #28：Mstage-2=M1 
            "#95353bFF",  #29:Gender-2 Male
            "#c7d79bFF",  #30：ETE-1=Minimal
            "#237c23FF",  #32：Multifocality-1
            
            "#386692FF",  #33：Mutation-1
            "ghostwhite", #34：Nstage-1=N0
            "#237c23FF",  #35：Nstage-2=N1
            "green4",  #36：⭐️NTEvent-2 = New primary tumor
            "#f9e5e9",    #37:Radiation-1 No  
            
            "#95353bFF",  #38:Coexist-3(Nodular)   #dfbbbeFF
            "ghostwhite", #39:【None】ETE-3 or Coexist-4  or Event-4  ##ebe8d3FF
            'ghostwhite', #RAI-3 = Not received
            "#386692FF",  #40:GeneFusion-2【Presence】  
            'red3',              #41：⭐️NTEvent-3 = Recurrence
            
            "#4e3394FF",  #42：Refractory
            "#d2d4f5",  #43：Sensitive
            "ghostwhite","#c7d79bFF","#75be70FF","#237c23FF",   #47-50:TStage
            "#f9e5e9",    #51：Coexist-1(Dysfuntion)
            "ghostwhite", #52：Multifocality-2          
            "aliceblue",  #53：Mutation-2【Wild-Type】
            "#95353bFF"   #54：Radiation-2【Yes】
            )

colors_na <- c("gray95")

names(colors) <-  qz   #Assign qz=the names of various elements to each color name in the colors section
head(colors)

## 2.The color of the comment bar
color_an <- c("#e11a0c","#48af45","#337cba") 
names(color_an) <- unique(col_group)
color_an

#### Plotting ####
# Main image
data <- Type[,-1]
p <- Heatmap(t(data), 
             col = colors,    #⭐️
             na_col = colors_na,   #缺失值的颜色
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#e11a0c","#48af45","#337cba"), col = c("#e11a0c","#48af45","#337cba")),
                                                                 #labels = c("BRFA-like with gene mutation", "BRFA-like with gene fusion", "RAS-like"),
                                                                 labels = c("Cluster 1 (N=188)","Cluster 2 (N=96)","Cluster 3 (N=193)"),
                                                                 labels_gp = gpar(col = "white", fontsize = 10, fontface = 'bold'),
                                                                 height = unit(5,"mm"))),
             row_split = row_group,
             column_split = col_group,
             show_heatmap_legend = F,
             border = T,
             show_column_names = F,
             row_names_side = "left",
             # P value and star：
             right_annotation = rowAnnotation(
               "P value" = anno_text(p_star, show_name = TRUE),
               annotation_name_side = "top",
               annotation_name_rot = 0),
             heatmap_width = unit(27,"cm"),
             heatmap_height = unit(13,"cm"),
             column_title = NULL,
             column_title_gp = gpar(fontsize =11)
             )
draw(p)

# Custom legend
k = 1
lgd = list()
for(i in 1:ncol(data)){
  un = sort(unique(data[,i]))
  ti = colnames(data)[i]
  lgd[[i]] = Legend(seq(0, 1, length.out = length(un)), labels = un,
                    title = ti, legend_gp = gpar(fill = colors[un]))  #⭐️
  k = k + length(un)
}

# Legend for missing values
lgd_na = Legend(
  at = 1, 
  labels = "Unknown",
  legend_gp = gpar(fill = colors_na)
)


# Merge various legends together
pd = packLegend(list = c(lgd, lgd_na), 
                max_width = unit(11, "cm"),
                direction = "horizontal",
                column_gap = unit(2, "mm"), 
                row_gap = unit(4, "mm")
                )

pd2 = packLegend(list = c(lgd, lgd_na), 
                max_height = unit(5.5, "cm"),    
                max_width = unit(18, "cm"),
                direction = "horizontal",
                column_gap = unit(2, "mm"), 
                row_gap = unit(2, "mm"))

# Draw the main body and legend together
pdf('cluster_CliHeatmapp.pdf', height = 7, width = 16 )
draw(p, heatmap_legend_list = pd, heatmap_legend_side = "right")   #vertical
#draw(p, heatmap_legend_list = pd2, heatmap_legend_side = "bottom")   #horizontal
dev.off()

#Export:horizontal=9x11    vertical= 16x7

