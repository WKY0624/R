
#install.packages("recipes")
#install.packages("caret")
library(recipes)
library(caret)


#### Random cutting ####
data <- read.table("./DataSets/06.clinical477.txt",header=T, sep="\t", check.names=F)
set.seed(5659)  
sam <- createDataPartition(data$Recurrence, p = 0.7, list = F)   
trainSet <- data[sam,]
testSet <- data[-sam,]

#View the cutting ratio of some parameters between two groups
prop.table(table(trainSet$Recurrence))
prop.table(table(testSet$Recurrence))
prop.table(table(data$Recurrence))
table(trainSet$Recurrence)
table(testSet$Recurrence)
table(data$Recurrence)

#Output file
write.table(trainSet, file = "trainSet.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(testSet, file = "testSet.txt",sep="\t",quote=F,row.names=F,col.names=T)


#### table 1 ####
library(table1)
library(dplyr)

#Data import and preprocessing
train <- read.table('trainSet.txt', header = T, sep = '\t', check.names = F)
test <- read.table('testSet.txt',  header = T, sep = '\t', check.names = F)
total <- read.table('./DataSets/06.clinical477.txt', header = T, sep = '\t', check.names = F)

train$Group = c('Training cohort')  
test$Group = c('Testing cohort')    
total2 <- rbind(train,test)          

total2$RFS <- total2$RFS/365
total2$OS <- total2$OS/365

#Assigning variables
total2$Recurrence = ifelse(total2$Recurrence=='1','Yes','No')
total2$Age = ifelse(total2$Age<55, '<55','>=55')

#Designated Unit
units(total2$RFS) <- "years"　
units(total2$OS) <- "years"　　

#Specify variable name
label(total2$OS) <- "Follow-up time"   #The follow-up time is based on OS, as RFS with recurrent events will be smaller than its OS
label(total2$Recurrence) <- 'Progression events'
label(total2$TNMStage) <- 'AJCC 8th TNM'
label(total2$TStage) <- 'T category'
label(total2$NStage) <- 'N category'
label(total2$MStage) <- 'M category'

##Three line table (single variable column)
table1(~ Age + Gender + TStage + NStage + MStage + TNMStage + Recurrence + OS | Group, data = total2, overall = "Entire cohort")

##Add statistical values
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value    #Continuous variables: parameter t. test, non parameter wilcox. test (if independent, paid=F, if paired, paid=T)
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value   #Categorical variables: In theory, fisher.test can be used for more accuracy, requiring computer power but without statistical values# ⭐ The cardholder will provide the statistical value chisq.test
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~ Age + Gender + TStage + NStage + MStage + TNMStage + Recurrence + OS | Group, data=total2, 
       extra.col=list(`P-value`=pvalue),  
       overall=F)



############ Other tips

##Convert to 2x2 table
mytable <- table(total2$Group, total2$Gender)

fisher = fisher.test(mytable, alternative = "two.sided",conf.int = T, simulate.p.value = TRUE)
fisher$estimate
fisher

chisq = chisq.test(mytable,correct = T)   #Whether to perform continuity correction
chisq$statistic
chisq$p.value


num = t.test(total2$OS ~ total2$Group)
num

##File export
#Save as HTML, copy to Word to modify table format


##Modify Table Style
#There are some built-in table styles in this package, such as:

#Zebra: alternating rows with and without shadows (zebra crossing).
#Grid: Display all gridlines
#Shade: Apply gray shading to the title line
#Times: Use serif fonts
#Center: Center all columns, including the first column containing row labels.
#These styles can be selected through the topclass parameter in table1. Here are some examples.

#table1(~ age + sex + wt | treat, data=dat, topclass="Rtable1-zebra")
#table1(~ age + sex + wt | treat, data=dat, topclass="Rtable1-zebra Rtable1-grid")
