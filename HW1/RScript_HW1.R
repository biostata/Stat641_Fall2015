rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 1, Stat 641, Fall 2015
########################################################################

########################################################################
## Load header files and source functions
########################################################################
#source('~/Courses/Stat641_Fall2015/HeaderFile_Nandi.R')
source('~/RScripts/HeaderFile_Nandi.R')

RScriptPath <- '~/Courses/Stat641_Fall2015/HW1/'
source('~/Courses/Stat641_Fall2015/fn_Library_Stat641.R')

########################################################################
## Load Data
########################################################################
Data <- read.csv(paste0(RScriptPath, 'data1.csv'))

# a)
t.test(w ~ z, data = Data, var.equal = T)

# b)
round(aggregate(x~z, function(y) c(mean=mean(y),sd=sd(y)),data=Data),3)
round(aggregate(w~z, function(x) c(mean=mean(x),sd=sd(x)),data=Data),3)

# c)
t.test(x ~ z, data = Data, var.equal = T)

# d)
pdf("Plot_d.pdf", width=6, height=4, paper='special') 
qplot() + geom_point(aes(x = w, y = x, col = factor(z), shape = factor(z)), data = Data, size=3) +
  xlab(label='w') + ylab(label='x')
dev.off()

# e)
model1 <- lm(x ~ as.factor(z) + w, data = Data)
summary(model1)
xtable(summary(model1))

model0 <- lm(x ~  as.factor(z), data = Data)
summary(model0)
xtable(summary(model0))

## pdf("e.pdf",width=6,height=4,paper='special') 
## plot(x~w, pch=(z=="0")+16,col=(z=="0")+160)
## mA <- lm(x~w,data=data, subset=z=="0")
## mB <- lm(x~w,data=data, subset=z=="1")
## abline(mA, col=161,lwd=2)
## abline(mB, col=160,lwd=2)
## legend("bottomright",pch=c(16,17),col=c(160,161),c("z=1","z=0"), bty="n")
## dev.off()

# g)
t.test(y ~ z, data = Data, var.equal = T)

# h)
model2 <- lm(y ~ as.factor(z) + w , data = Data)
summary(model2)

xtable(summary(model12))
# i)

model3 <- lm(y ~ as.factor(z) + w + x, data = Data)
summary(model3)

model30 <- lm(y ~ as.factor(z) + w , data = Data)
summary(model30)

xtable(summary(model3))
