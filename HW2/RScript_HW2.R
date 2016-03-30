rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 2, Stat 641, Fall 2015
########################################################################
library(xtable)
########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW2/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

########################################################################
## Problem 1
########################################################################
Data <- read.csv(paste0(RScriptPath, 'data2.csv'))

## (a)
Model1 <- lm(x ~ as.factor(z), data = Data)
summary(Model1)

xtable(Model1)

## (b)
Model2 <- lm(x ~ as.factor(z) + as.factor(w), data = Data)
summary(Model2)

########################################################################
## Problem 3
########################################################################
library(survival)

Data <- as.data.frame(cbind(
  time = c(8, 11, 16, 18, 23, 24, 26, 28, 30, 31, 
           9, 12, 13, 14, 14, 16, 19, 22, 23, 29), 
  status = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 
             1, 1, 1, 1, 1, 1, 0, 0, 0, 0), 
  grp = c(rep('A', 10), rep('B', 10))
))

Data$time <- as.numeric(as.vector(Data$time))
Data$status <- as.numeric(as.vector(Data$status))
str(Data)

Model2 <- survfit(Surv(time, status) ~ 1, data = Data)
summary(Model2)

Table <- cbind(Model2$time, Model2$n.risk, Model2$n.event, Model2$n.censor, Model2$surv)
colnames(Table) <- c('time', 'n.risk', 'n.event', 'n.censor', 'survival')


print(xtable(Table, digits = c(0, 0, 0, 0, 0, 4)), include.rownames=FALSE)

filename <- paste0(RScriptPath, 'Plot3a.pdf')
pdf(file = filename)
plot(Model2)
title("Kaplan-Meier estimate for both groups combined") 
dev.off()

Partb <- cbind(time = summary(Model2)[['time']], 
               St = summary(Model2)[['surv']],
               Var = (summary(Model2)[['std.err']])^2 ,
               Lower = summary(Model2)[['lower']],
               Upper = summary(Model2)[['upper']])

print(xtable(Partb, digits=c(0, 0, 4, 4, 4, 4)), include.rownames = FALSE)

options(digits=9)
Diff <- survdiff(Surv(time, status) ~ grp, data = Data)
Diff
        
DiffTable <- cbind(Diff[['n']],
      Diff[['obs']],
      Diff[['exp']],
      Diff[['chisq']]
      )
colnames(DiffTable) <- c('N', 'Observed', 'Expected', '')

xtable(DiffTable, digits = c(0, 0, 0, 4, 4))


Diff1 <- survdiff(Surv(time, status) ~ grp, data = Data, rho = 1)
Diff1

DiffTable1 <- cbind(Diff1[['n']],
      Diff1[['obs']],
      Diff1[['exp']],
      Diff1[['chisq']]
      )
colnames(DiffTable1) <- c('N', 'Observed', 'Expected', '')

xtable(DiffTable1, digits = c(0, 0, 0, 4, 4))

Model3 <- survfit(Surv(time, status) ~ grp, data = Data)
summary(Model3)
filename <- paste0(RScriptPath, 'Plot3c.pdf')
pdf(file = filename)
plot(Model3, lty = c(1, 2))
title("Kaplan-Meier estimate for treatment groups A and B") 
dev.off()
