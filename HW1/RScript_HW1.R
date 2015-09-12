rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script is HW 1, Stat 641, Fall 2015
########################################################################

########################################################################
## Load header files and source functions
########################################################################
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
round(aggregate(x~z, function(x) c(mean=mean(x),sd=sd(x)),data=Data),3)
round(aggregate(w~z, function(x) c(mean=mean(x),sd=sd(x)),data=Data),3)

# c)
t.test(x ~ z, data = Data, var.equal = T)

# d)
ggplot() + geom_point(aes(x = w, y = x), data = Data)

#pdf("d.pdf",width=6,height=4,paper='special') 
plot(x~w, data = Data, pch=(z=="0")+16,col=(z=="0")+160)
legend("bottomright",pch=c(16,17),col=c(160,161),c("z=1","z=0"), bty="n")
#dev.off()
# e)
summary(lm(x ~ as.factor(z) , data=data))
model1<-lm(x ~ as.factor(z) + w, data=data)
summary(model1)

model10<-lm(x ~  w+as.factor(z), data=data)
summary(model10)

model100<-lm(x-w ~  as.factor(z), data=data)
summary(model100)

xtable(summary(model1))

pdf("e.pdf",width=6,height=4,paper='special') 
plot(x~w, pch=(z=="0")+16,col=(z=="0")+160)
mA <- lm(x~w,data=data, subset=z=="0")
mB <- lm(x~w,data=data, subset=z=="1")
abline(mA, col=161,lwd=2)
abline(mB, col=160,lwd=2)
legend("bottomright",pch=c(16,17),col=c(160,161),c("z=1","z=0"), bty="n")
dev.off()

# g)

t.test(y ~ z, data = data,var.equal = T)

# h)

model12<-lm(y ~ as.factor(z) + w , data=data)
summary(model12)
xtable(summary(model12))
# i)

model2<-lm(y ~ as.factor(z) + w + x, data=data)
summary(model2)
model20<-lm(y ~ as.factor(z) + w , data=data)
summary(model20)
xtable(summary(model2))
