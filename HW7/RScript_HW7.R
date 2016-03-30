rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 7, Stat 641, Fall 2015
########################################################################
library(xtable)

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
#source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW7/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

########################################################################
## Problem 1
########################################################################
A <- c(3, 1.7, 4.6, 3.6, 1.2, 3.1)
B <- c(3.9, 6.7, 7.5, 3.4, 8.8, 9.6)
D_Obs <- mean(B) - mean(A)

Pop <- c(A, B)
N <- 10000
D_Rands <- vector(mode = 'numeric', length = N)
for(i in 1:N){
  A_Rand <- sample(x = Pop, size = 6, replace = F)
  B_Rand <- Pop[!(Pop %in% A_Rand)]
  D_Rand <- mean(B_Rand) - mean(A_Rand)
  D_Rands[i] <- D_Rand  
}
hist(D_Rands, breaks = 20)
round(mean(D_Obs < D_Rands), 4)

wilcox.test(x = A, y = B, alternative = 'less', paired = FALSE)

library(ggplot2)
Plot1 <- qplot() + geom_histogram(aes(D_Rands), binwidth = 0.5) +
  ylab(label = '') + xlab(label = expression(paste('Randomized ', hat(beta)))) +
  ggtitle(label = 'Histogram of randomized test') +
  geom_vline(xintercept = D_Obs, col = 'blue')
pdf(file = paste0(RScriptPath, 'Plot1.pdf'))
Plot1
dev.off()


########################################################################
## Problem 3
########################################################################
Data <- read.table(file = paste0(RScriptPath, 'data7.csv'), sep=',', header = T)
str(Data)

Model1 <- lm( y ~ z, data = Data)
summary(Model1)
xtable(summary( Model1 ))

Model2 <- lm( y ~ z + w, data = Data)
summary(Model2)
xtable(summary( Model2 ))

## Complete Randomization
modCR <- lm( y ~ z, data = Data)
modCR$coef[2]   ## Real value

betaSimCR <- replicate(10000, {Data$z <- sample(x = 0:1, size = 60, replace = T) ; lm(y ~ z, data = Data)$coef[2]})
var(betaSimCR)
sd(betaSimCR)
round(mean(abs(betaSimCR)>=abs(modCR$coef[2])), 4)

## Random ALlocation 30/30
betaSimRA <- replicate(10000,{Data$z <- sample(x = rep(0:1,30), size = 60, replace = T) ; lm( y ~ z, data = Data)$coef[2]})
var(betaSimRA)
sd(betaSimRA)
round(mean(abs(betaSimRA)>=abs(modCR$coef[2])), 4)

## Permuted Block (size 4)
betaSimPB <- replicate(n = 10000, expr = {Data$z <- c(replicate(15, sample(rep(0:1,2),4))) ; lm( y ~ z, data = Data)$coef[2]})
var(betaSimPB)
sd(betaSimPB)
round(mean(abs(betaSimPB)>=abs(modCR$coef[2])), 4)
