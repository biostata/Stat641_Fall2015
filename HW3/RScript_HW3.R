rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 3, Stat 641, Fall 2015
########################################################################
library(xtable)
########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW3/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

Data <- read.table(file = paste0(RScriptPath, 'data3.csv'), sep=',', header = T)
str(Data)
Data$trt <- as.factor(Data$trt)
Data$sex <- factor(Data$sex, labels=c('Male', 'Female')) 
str(Data)

########################################################################
## Part a
########################################################################
Model1 <- survfit(Surv(days,status) ~ trt, data = Data)
summary(Model1)

pdf("Plota.pdf", width=8,height=5)
plot(Model1,lty=1:2, lwd=2, xlab="Time (Days)", ylab="Survival Probability", mark.time=F)
legend("bottomright", c("Treatment: 0","Treatment: 1"),lwd=2, lty=1:2, bty="n")
dev.off()

########################################################################
## Part b
########################################################################
Model2 <- coxph(Surv(days,status)~trt, data=Data)
summary(Model2)

########################################################################
## Part d
########################################################################
Model3 <- coxph(Surv(days,status)~trt + age + sex, data=Data)
summary(Model3)

########################################################################
## Part e
########################################################################
pdf("Plote.pdf", width=8,height=5)
plot(survfit(Surv(days,status)~trt, data=Data), lty=1:2, fun="cloglog")
legend("bottomright", c("Treatment: 0","Treatment: 1"),lwd=2, lty=1:2, bty="n")
dev.off()
