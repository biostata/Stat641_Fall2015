rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 8, Stat 641, Fall 2015
########################################################################
library(xtable)

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW8/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

Data <- read.table(file = paste0(RScriptPath, 'data8.csv'), sep=',', header = T)
str(Data)

########################################################################
## Part (a)
########################################################################
Data1 <- subset(Data, t == 1)
n1 <- sum(Data1$z)
n0 <- nrow(Data1) - n1
Data1$rank <- rank(Data1$y)
W <- sum(Data1$rank[Data1$z == 1])
Mu <- n1*(n1 + n0 + 1)/2
Sigma2 <- n1*n0*(n1 + n0 + 1)/12
Sigma <- sqrt(Sigma2)
pnorm(q = W, mean = Mu, sd = Sigma, lower.tail = FALSE)
wilcox.test(x = Data1$y[Data1$z == 1], y = Data1$y[Data1$z == 0], alternative = 'greater', paired = FALSE)

Data2 <- subset(Data, t != 3)
n1 <- sum(Data2$z)
n0 <- nrow(Data2) - n1
Data2$rank <- rank(Data2$y)
W <- sum(Data2$rank[Data2$z == 1])
Mu <- n1*(n1 + n0 + 1)/2
Sigma2 <- n1*n0*(n1 + n0 + 1)/12
Sigma <- sqrt(Sigma2)
pnorm(q = W, mean = Mu, sd = Sigma, lower.tail = FALSE)
wilcox.test(x = Data2$y[Data2$z == 1], y = Data2$y[Data2$z == 0], alternative = 'greater', paired = FALSE)

Data3 <- Data
n1 <- sum(Data3$z)
n0 <- nrow(Data3) - n1
Data3$rank <- rank(Data3$y)
W <- sum(Data3$rank[Data3$z == 1])
Mu <- n1*(n1 + n0 + 1)/2
Sigma2 <- n1*n0*(n1 + n0 + 1)/12
Sigma <- sqrt(Sigma2)
pnorm(q = W, mean = Mu, sd = Sigma, lower.tail = FALSE)
wilcox.test(x = Data3$y[Data3$z == 1], y = Data3$y[Data3$z == 0], alternative = 'greater', paired = FALSE)

########################################################################
## Part (b)
########################################################################
Mu1 <- 20372.5; Sigma1 <- 677.0386
Mu2 <- 105120; Sigma2 <- 2426.256
Mu3 <- 239239; Sigma3 <- 4562.214

fn_TT <- function(Mu1, Mu2, Mu3, Sigma1, Sigma2, Sigma3){
T1 <- rnorm(n = 1, mean = Mu1, sd = Sigma1)
T2 <- rnorm(n = 1, mean = Mu2, sd = Sigma2)
T3 <- rnorm(n = 1, mean = Mu3, sd = Sigma3)
  return(c(T1 - Mu1, T2 - Mu2, T3 - Mu3))
}  

Matrix_b <- replicate(n = 10000, fn_TT(Mu1, Mu2, Mu3, Sigma1, Sigma2, Sigma3))
CovMatrix_b <- cov(t(Matrix_b))
CorMatrix_b<- cor(t(Matrix_b))
xtable(CovMatrix_b, digits = c(0,0, 0, 0))
xtable(CorMatrix_b, digits = c(0,4, 4, 4))

########################################################################
## Part (c)
########################################################################
Mu1 <- 20372.5/135; Sigma1 <- 677.0386/135
Mu2 <- 105120/320; Sigma2 <- 2426.256/320
Mu3 <- 239239/522; Sigma3 <- 4562.214/522

Matrix_c <- replicate(n = 10000, fn_TT(Mu1, Mu2, Mu3, Sigma1, Sigma2, Sigma3))
CovMatrix_c <- cov(t(Matrix_c))
CorMatrix_c <- cor(t(Matrix_c))
xtable(CovMatrix_c, digits = c(0,2, 2, 2))
xtable(CorMatrix_c, digits = c(0,4, 4, 4))

########################################################################
## Part (d)
########################################################################
Sigma1^2
Sigma2^2
Sigma3^2

I1 <- Sigma1^2/Sigma3^2
I2 <- Sigma2^2/Sigma3^2

########################################################################
## Part (e)
########################################################################
source('ldbounds.R')
Bound <- bounds(c(I1, I2, 1), iuse=3, alpha=0.025, phi = 1)
summary(Bound)
xtable(summary(Bound)$bounds, digits = c(0, 4, 4, 4, 4, 4, 4))

Z1 <- 2.6424

Crit1 <- Z1*Sigma1 + Mu1
pnorm(q = Crit1, mean = Mu1, sd = Sigma1, lower.tail = FALSE)

Z2 <- 2.4959
Crit2 <- Z2*Sigma2 + Mu2
pnorm(q = Crit2, mean = Mu2, sd = Sigma2, lower.tail = FALSE)

Z3 <- 2.5096
Crit3 <- Z3*Sigma3 + Mu3
pnorm(q = Crit3, mean = Mu3, sd = Sigma3, lower.tail = FALSE)

########################################################################
## Part (f)
########################################################################
source('ldbounds.R')
Bound_f <- bounds(c(0.3288, 0.7574, 1), iuse=3, alpha=0.025, phi = 1)
summary(Bound_f)
xtable(summary(Bound_f)$bounds, digits = c(0, 4, 4, 4, 4, 4, 4))
