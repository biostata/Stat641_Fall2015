rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is Prob 5, Final Exam, Stat 641, Fall 2015
########################################################################
library(xtable)
library(survival)

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/FinalExam/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

Filename <- paste0(RScriptPath, 'dataFinal2.csv')
Data <- read.csv(Filename)
source('ldbounds.R')
########################################################################
## Part a
########################################################################
mod1 <- survdiff(Surv(time1, dead1) ~ z, data = Data)
mod1

mod2 <- survdiff(Surv(time2, dead2) ~ z, data = Data)
mod2

mod3 <- survdiff(Surv(time3, dead3) ~ z, data = Data)
mod3

mod4 <- survdiff(Surv(time4, dead4) ~ z, data = Data)
mod4

mod5 <- survdiff(Surv(time5, dead5) ~ z, data = Data)
mod5

Table_a <- cbind(Period = c(1:5),
      Chi = c(round(mod1$chisq, 4), round(mod2$chisq, 4), round(mod3$chisq, 4), round(mod4$chisq, 4), round(mod5$chisq, 4)),
      Z = c(round(sqrt(mod1$chisq), 4), round(sqrt(mod2$chisq), 4), round(sqrt(mod3$chisq), 4), round(sqrt(mod4$chisq), 4), round(sqrt(mod5$chisq), 4))
      )
xtable(Table_a, digits = c(0, 0, 4, 4))

########################################################################
## Part b
########################################################################
N2 <- 1272
IRatio1 <- c(sum(mod1$obs), sum(mod2$obs), sum(mod3$obs), sum(mod4$obs), N2)/N2 
Bound1 <- bounds(IRatio1, iuse=3, alpha=0.025, phi = 2)
xtable(summary(Bound1)$bounds, digits = c(0, 4, 4, 4, 4, 4, 4))

########################################################################
## Part c
########################################################################
sum(mod5$obs)
IRatio2 <- c(sum(mod1$obs), sum(mod2$obs), sum(mod3$obs), sum(mod4$obs), sum(mod5$obs))/sum(mod5$obs)
round(IRatio2, 4)
Bound2 <- bounds(IRatio2, iuse=3, alpha=0.025, phi = 2)
xtable(summary(Bound2)$bounds, digits = c(0, 4, 4, 4, 4, 4, 4))

########################################################################
## Part e
########################################################################
DriftObj_e <- drift(
      t = IRatio,
      zb = c(Bound1$upper[1:4], Bound2$upper[5]),
      ## zb = c(Bound2$upper),
      za = rep(-10, 5),
      drft = 0
      )
DriftObj_e$power

lr1 <- drift(zb = c(Bound1$upper[1]), za = rep(-10,1), t = IRatio2[1], drft = 0)
lr2 <- drift(zb = c(Bound1$upper[1], Bound2$upper[5]), za = rep(-10,2), t = IRatio2[1:2], drft = 0)
lr3 <- drift(zb = c(Bound1$upper[1:2], Bound2$upper[5]), za = rep(-10,3), t = IRatio2[1:3], drft = 0)
lr4 <- drift(zb = c(Bound1$upper[1:3], Bound2$upper[5]), za = rep(-10,4), t = IRatio2[1:4], drft = 0)
lr5 <- drift(zb = c(Bound1$upper[1:4], Bound2$upper[5]), za = rep(-10,5), t = IRatio2[1:5], drft = 0)

lr1$exit
lr2$exit
lr3$exit
lr4$exit
lr5$exit

########################################################################
## Part f
########################################################################
bK <- Bound2$upper.bound[5]
T <- IRatio2
Z <- Table_a[,'Z']
Z[1] <- -Z[1]
B <- sqrt(T)*Z
Theta <- -log(0.8)*sqrt(1140/4)

1 - pnorm((bK - B[1]-(1-T[1])*Theta)/(sqrt(1-T[1])))

1 - pnorm((bK - B[2]-(1-T[2])*Theta)/(sqrt(1-T[2])))

1 - pnorm((bK - B[3]-(1-T[3])*Theta)/(sqrt(1-T[3])))

1 - pnorm((bK - B[4]-(1-T[4])*Theta)/(sqrt(1-T[4])))

