library(nlme)
library(lme4)
library(splines)
library(ggplot2)


setwd('/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3')
Df.138=read.csv('Extracted_138.csv')
#

#
ggplot(Df.138,aes(y=X138max,x=I.WHO,group=I.WHO))+
  geom_boxplot()+
  scale_y_continuous(trans='log2')
#
ggplot(Df.138,aes(x=I.WHO,y=X138max))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y~x^2) +
  geom_point()