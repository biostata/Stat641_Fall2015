rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 6, Stat 641, Fall 2015
########################################################################
library(xtable)

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW6/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

########################################################################
## Problem 1
########################################################################
Data <- read.table(file = paste0(RScriptPath, 'data6.csv'), sep=',', header = T)
str(Data)

## (a)
Means <- with(Data, tapply(X = y, INDEX = list(z, period), FUN = mean))
xtable(x = Means, digits = c(0, 4, 4))

## (b)
Model1 <- lm(y ~ z + period + id, data = Data)
summary(Model1)
xtable(summary(Model1))

## (c)
Model2 <- lm(y ~ z, data = subset(Data, period==1))
summary(Model2)
xtable(summary(Model2))
