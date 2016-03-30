rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

library(ggplot2)

RScriptPath <- '/ua/snandi/Courses/Stat641_Fall2015/HW1/'

Data <- read.csv(paste0(RScriptPath, 'data1.csv'))

qplot() + geom_point(aes(x = w, y = x, col = factor(z), shape = factor(z)), data = Data, size=3) +
  xlab(label='w') + ylab(label='x')
