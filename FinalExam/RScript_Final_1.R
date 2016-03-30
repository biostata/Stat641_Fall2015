rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is Prob 5, Final Exam, Stat 641, Fall 2015
########################################################################
library(xtable)
library(survival)
library(ggplot2)
library(reshape2)
########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/FinalExam/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

Filename <- paste0(RScriptPath, 'dataFinal1.csv')
Data <- read.csv(Filename)

########################################################################
## Part a
########################################################################
Data_Long <- melt(data = Data, measure.vars = c('x1', 'x2'))
colnames(Data_Long) <- c('x0', 'z', 'x12', 'y')
Plota <- qplot() + geom_point(aes(x = x0, y = y, shape = factor(z), col = factor(z)), size = 4, data = Data_Long) +
  facet_grid(. ~ x12) +
  ggtitle(label = 'Reponse vs Baseline') + 
  ylab(label = 'x1, x2') +
  theme(strip.text.x = element_text(size = 14), 
        axis.title = element_text(size = 14),
        legend.position = 'top')

Plotb1 <- qplot() + geom_boxplot(aes(x = factor(z), y = x0, fill = factor(z)), data = Data) +
  ggtitle(label = 'Baseline vs treatment') +
  ylab(label = 'x0') + xlab(label = 'z') +
  theme(axis.title = element_text(size = 14),
        legend.position = 'top')
  
Plotb2 <- qplot() + geom_boxplot(aes(x = factor(z), y = y, fill = factor(z)), data = Data_Long) +
  facet_grid(. ~ x12) +
  ggtitle(label = 'Response vs treatment') +
  ylab(label = 'Response') + xlab(label = 'z') +
  theme(strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = 'top')


pdf(file = paste0(RScriptPath, 'Plot1a.pdf'), onefile = T)
Plota
Plotb1
Plotb2
dev.off()

########################################################################
## Part b
########################################################################
Modelb1 <- lm(x1 ~ z, data = Data)
summary(Modelb1)
xtable(summary(Modelb1))

Modelb2 <- lm((x1 - x0) ~ z, data = Data)
summary(Modelb2)
xtable(summary(Modelb2))

Modelb3 <- lm(x1 ~ z + x0, data = Data)
summary(Modelb3)
xtable(summary(Modelb3))

########################################################################
## Part c
########################################################################
Modelc1 <- lm(x2 ~ z, data = Data)
summary(Modelc1)
xtable(summary(Modelc1))

Modelc2 <- lm((x2 - x0) ~ z, data = Data)
summary(Modelc2)
xtable(summary(Modelc2))

Modelc3 <- lm(x2 ~ z + x0, data = Data)
summary(Modelc3)
xtable(summary(Modelc3))

Modelc4 <- lm(x2 ~ z + x0 + x1, data = Data)
summary(Modelc4)
xtable(summary(Modelc4))
