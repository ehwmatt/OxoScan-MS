library(nlme)
library(lme4)
library(splines)
library(ggplot2)
library(psych)
library(MASS)

setwd('/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3')
Df.138=read.csv('Extracted_138.csv')
Df.138$X138max<-scale(Df.138$X138max)
#
psych::describeBy(Df.138$X138max, Df.138$I.WHO,mat=TRUE)
#
ggplot(Df.138,aes(y=X138max,x=I.WHO,group=I.WHO))+
  geom_boxplot()
  #scale_y_continuous(trans='log2')
#
ggplot(Df.138,aes(x=I.WHO,y=X138max))+
  geom_smooth(method = "lm", color="black", formula = y~ns(x,df=2)) +
  geom_point()
#
#
lmnull<-lm(X138max~1,data=Df.138)
lm0<-lm(X138max~I.WHO,data=Df.138)
lmS1<-lm(X138max~ns(I.WHO,2),data=Df.138)
lmS<-lm(X138max~ns(I.WHO,3),data=Df.138)
#
anova(lmnull,lm0,lmS1,lmS)
###
#############################
###
# Ordinal Logistic regression 
###
#############################
###
Df.138$I.WHO<-factor(Df.138$I.WHO)
ord.lm<-polr(I.WHO ~ X138max, data = Df.138, Hess=TRUE)
ctable <- coef(summary(ord.lm))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p) ## combined table
#
ci <- confint(ord.lm)


