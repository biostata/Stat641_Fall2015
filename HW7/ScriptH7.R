setwd("/home/cony/Documents/FALL 2014/STAT 641/HWK 7")
library(xtable)

data7<-read.csv("data7.csv")
names(data7)
head(data7)

data7$w<-as.factor(data7$w)
str(data7)

mod1<-lm(y~z, data=data7)
xtable(summary(mod1))

mod2<-lm(y~z+w, data=data7)
xtable(summary(mod2))


## Complete Randomization

modCR<-lm(y~z, data=data7)
modCR$coef[2]   ## Real value

betaSimCR<-replicate(10000,{data7$z<-sample(0:1,60, rep=T) ; lm(y~z, data=data7)$coef[2]})
var(betaSimCR)
sd(betaSimCR)
mean(abs(betaSimCR)>=abs(modCR$coef[2]))

hist(betaSimCR)

## Random ALlocation 30/30

sample(rep(0:1,30),60, rep=F)

betaSimRA<-replicate(10000,{data7$z<-sample(rep(0:1,30),60, repl=T) ; lm(y~z, data=data7)$coef[2]})
hist(betaSimRA)

var(betaSimRA)
sd(betaSimRA)
mean(abs(betaSimRA)>=abs(modCR$coef[2]))

## Permuted Block (size 4)


replicate(15,sample(rep(0:1,2),4))

c(replicate(15,sample(rep(0:1,2),4)))

betaSimPB<-replicate(10000,{data7$z<-c(replicate(15,sample(rep(0:1,2),4))) ; lm(y~z, data=data7)$coef[2]})
hist(betaSimPB)

var(betaSimPB)
sd(betaSimPB)
mean(abs(betaSimPB)>=abs(modCR$coef[2]))

pdf("betas.pdf", width=5, height=3)
par(mfrow=c(1,3))
hist(betaSimCR, freq=F, xlab=expression(beta) , main="")
hist(betaSimRA, freq=F, xlab=expression(beta) , main="")
hist(betaSimPB, freq=F, xlab=expression(beta) , main="")
dev.off()
