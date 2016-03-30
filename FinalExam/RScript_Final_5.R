rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is Prob 5, Final Exam, Stat 641, Fall 2015
########################################################################
library(xtable)

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW8/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

########################################################################
## Part (a)
########################################################################
qnorm(p = 0.975)
qnorm(p = 0.9)
N <- (qnorm(p = 0.975) + qnorm(p = 0.9))^2 / (0.5 * 0.5 * (log(0.8)^2))
N
Rho <- 1 + (exp(-0.36) - exp(-0.18))/(0.18)

N/Rho

((qnorm(p = 0.975) + qnorm(p = 0.9))^2)*2*0.09 / ( (2*0.09 - exp(-0.09*2) + exp(-0.09*4))*(0.5 * 0.5 * (log(0.8)^2)) )

########################################################################
## Part (b)
########################################################################
T <- 2
T <- seq(from = 1.45, to = 1.5, by = 0.01)
R <- 3583/2
Lambda <- 0.09
X <- (R/Lambda)*(1 - exp(-Lambda*T))
Events <- R*T - X
cbind(T, Events)

X2 <- 3279
E2 <- 304

X22 <- 3279 - 34; T22 <- (log(X2) - log(X22))/Lambda; T22
X33 <- 3279 - 34 - 169; T33 <- (log(X2) - log(X33))/Lambda; T33
X44 <- 3279 - 34 - 2*169; T44 <- (log(X2) - log(X44))/Lambda; T44
X55 <- 3279 - 34 - 3*169; T55 <- (log(X2) - log(X55))/Lambda; T55

########################################################################
## Part (c)
########################################################################
source('ldbounds.R')
IRatio <- c(0.2, 0.4, 0.6, 0.8, 1)
Bound <- bounds(IRatio, iuse=3, alpha=0.025, phi = 2)
summary(Bound)
xtable(summary(Bound)$bounds, digits = c(0, 4, 4, 4, 4, 4, 4))

########################################################################
## Part (d)
########################################################################
I5 <- 844*0.5*0.5
Theta <- abs(sqrt(I5)*log(0.8))
EZ <- sqrt(IRatio)*Theta
Bound_d <- bounds(IRatio, iuse=3, alpha=0.15, phi = 2)

Upper <- Bound_d$upper.bounds
ak <- EZ - Upper
round(ak, 4)

########################################################################
## Part (e)
########################################################################
#DriftObj <- drift(t = IRatio, zb = Bound$upper.bounds, za = rep(-10,5), drft = 3.24)
DriftObj <- drift(t = IRatio, zb = Bound$upper.bounds, za = ak, drft = 3.24)

DriftObj$power + sum(Bound$nom.alpha)
plot(DriftObj)

########################################################################
## Part (f)
########################################################################
pdf(file = 'Plot5f.pdf')
plot(DriftObj)
dev.off()

########################################################################
## Part (g)
########################################################################
Theta2 <- EZ[5] + Bound$upper.bounds[5] - ak[5]
DriftObj2 <- drift(t = IRatio, zb = Bound$upper.bounds, za = ak, drft = Theta2)
#plot(DriftObj2)
DriftObj2$power 

########################################################################
## Part (h)
########################################################################
N2 <- (Theta2/log(0.8))^2 / (0.5^2)

N2
N2/Rho
