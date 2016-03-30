rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 3, Stat 641, Fall 2015
########################################################################
library(xtable)
library(ggplot2)
########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW5/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

Data1 <- read.table(file = paste0(RScriptPath, 'data5.csv'), sep=',', header = T)

Data <- Data1

Data$z <- factor(Data$z, labels = c('Trt0', 'Trt1'))
Data$dead <- factor(Data$dead, labels = c('Alive', 'Dead'))
str(Data)

########################################################################
## Part a
########################################################################
Table_a <- table(Data$z, Data$dead)
chisq.test(Table_a, correct = F)
xtable(Table_a)

########################################################################
## Part b
########################################################################
Data.Alive <- subset(Data, dead == 'Alive')
Plotb <- qplot() + geom_boxplot(aes(x = z, y = y), data = Data.Alive) +
  ylab(label = 'y') + ggtitle(label = 'For people who are alive')
pdf(file = paste0(RScriptPath, 'Plotb.pdf'))
Plotb
dev.off()

Modelb <- lm(y ~ z, data = Data.Alive)
summary(Modelb)
xtable(summary(Modelb))

########################################################################
## Part c
########################################################################
summary(Data.Alive$y)
Data.c <- Data
Data.c$y[Data.c$dead == 'Dead'] <- 2.0

Plotc <- qplot() + geom_boxplot(aes(x = z, y = y), data = Data.c) +
  ylab(label = 'y') + ggtitle(label = 'Difference of two treatments')
pdf(file = paste0(RScriptPath, 'Plotc.pdf'))
Plotc
dev.off()

Modelc <- lm(y ~ z, data = Data.c)
summary(Modelc)
xtable(summary(Modelc))

median(Data.c$y[Data.c$z=='Trt0'])
median(Data.c$y[Data.c$z=='Trt1'])

y0 <- Data.c$y[Data.c$z=='Trt0']
y1 <- Data.c$y[Data.c$z=='Trt1']

wilcox.test(y0, y1, alternative='greater')

wilcox.test(y0,y1+2.25, alternative="greater")



