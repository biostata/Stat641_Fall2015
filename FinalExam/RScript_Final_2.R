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

########################################################################
## Data
########################################################################
Data_Obs <- as.data.frame(cbind(
                                time = c(8, 11, 21, 35, 38,
                                  37, 41, 76, 83, 119),
                                status = c(1, 1, 0, 1, 1,
                                  1, 0, 1, 0, 1),
                                group = c(rep(1, 5), rep(2, 5))
                          ))

Model1 <- survfit(Surv(time, status) ~ 1, data = Data_Obs)
summary(Model1)

Data <- Data_Obs
Data <- Data[order(Data$time),]
Data$n.risk <- c(10:1)
Data$n.risk.inv <- 1/Data$n.risk
Data$n.risk.inv.cumsum <- cumsum(Data$n.risk.inv * Data$status)

## Data$Score_C <- cumsum(1/Data$n.risk)*(Data$status) + cumsum(1/Data$n.risk - 1)*(1 - Data$status)
Data$Score_C <-  Data$n.risk.inv.cumsum - Data$status

xtable(Data, digits = c(rep(0, 5), rep(4, 3)))

T_Obs <- sum(Data$Score_C[Data$group == 2])
T_Obs

fn_return_grp2T <- function(Data){
  Data <- Data[order(Data$time),]
  Data$n.risk <- c(10:1)
  Data$n.risk.inv <- 1/Data$n.risk
  Data$n.risk.inv.cumsum <- cumsum(Data$n.risk.inv * Data$status)
  Data$Score_C <-  Data$n.risk.inv.cumsum - Data$status
  #print(round(sum(Data$Score_C), 3))

  T_grp2 <- sum(Data$Score_C[Data$group == 2])
  return(T_grp2)
}

T_Obs <- fn_return_grp2T(Data = Data)

fn_permute <- function(Data = Data_Obs, N = 10){
  T_Permute <- vector(mode = 'numeric', length = N)
  for(i in 1:N){
    Grp1 <- sample(x = 1:10, size = 5, replace = FALSE)
    Grp2 <- c(1:10) %w/o% Grp1
    Data_Permute <- Data[c(Grp1, Grp2),]
    Data_Permute$group <- c(rep(1, 5), rep(2, 5))
    T_Permute[i] <- fn_return_grp2T(Data = Data_Permute)
  }
  return(T_Permute)
}

T_Permute <- fn_permute(Data = Data_Obs, N = 10000)

summary(T_Permute)
pValue <- round(1 - mean(T_Obs > abs(T_Permute)), 4)
#pValue <- round(1 - mean(T_Obs > T_Permute), 4)

Plot2a <- qplot() + geom_histogram(aes(x = T_Permute), bindiwth = 0.5, fill = 'gray40') +
  ggtitle('Randomization distribution of T') +
  ylab(label = '') + xlab(paste('p-value:', pValue)) + 
  geom_vline(xintercept = T_Obs, size = 1)

Filename <- paste0(RScriptPath, 'Plot2.pdf')
pdf(file = Filename)
Plot2a
dev.off()

var(T_Permute)
T_Stat <- T_Obs^2/var(T_Permute)
pchisq(q = T_Stat, df = 1, lower.tail = FALSE)

########################################################################
## Part b
########################################################################
Diff1 <- survdiff(Surv(time, status) ~ group, data = Data_Obs, rho = 1)
Diff1

DiffTable1 <- cbind(Diff1[['n']],
      Diff1[['obs']],
      Diff1[['exp']],
      Diff1[['chisq']]
      )
colnames(DiffTable1) <- c('N', 'Observed', 'Expected', '')

xtable(DiffTable1, digits = c(0, 0, 0, 4, 4))

