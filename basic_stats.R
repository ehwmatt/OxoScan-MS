library(nlme)
library(lme4)
library(splines)
library(ggplot2)
library(psych)
library(MASS)
library(quantmod)
library(class)

setwd('/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3')
Df.138=read.csv('Extracted_138.csv')
Df.138$X138max<-scale(Df.138$RTint1381)
#
IDs=unique(Df.138$pID)
##### Peak finding ######
dfpeaks=data.frame()
n.peaks=30
for (i in IDs){
  df<-Df.138[Df.138$pID==i,]
  df2<-df[seq(1,nrow(df),3),]
  p=which(diff(sign(diff(df2$RTint1381)))==-2)
  p=p[which(diff(p)>5)]
    #findPeaks(df2$RTint1381,3)
  df3<-df2[p,] %>%
    arrange(desc(RTint1381))
  df3<-df3[1:min(nrow(df3),n.peaks),] %>%
    arrange(MZ)
  dfpeaks<-rbind(dfpeaks,df3)
}
#
dfpeaks$F.Severity<-factor(dfpeaks$F.Severity,c("healthy","mild","severe","critical"))
# example plot
ggplot(df2, aes(x=MZ,y=RTint1381))+geom_line()+geom_point(data=df3,color='red')+
  labs(x='M/Z',y='RT integrated intensity for 138')
  
#############################
###
# KNN Clustering
###
#############################
###
train.N<-sample(1:nrow(dfpeaks),size=nrow(dfpeaks)*0.7,replace = FALSE) #random selection of 70% data

train<-dfpeaks[train.N,c('pID','MZ','X138max')]
test<-dfpeaks[-train.N,c('pID','MZ','X138max')] # MASS Spec data 
Labels.train<-dfpeaks[train.N,'F.Severity'] 
Labels.test<-dfpeaks[-train.N,'F.Severity'] # Disease severity, can be F.Severity or I.WHO

########################
K<-sqrt(nrow(train)) # optimal starting point for cluster #
K<-seq(3,ceiling(K),1)

KNN.fits<-rep(0,length(K))
KNN.mod<-vector('list',length(K))
i=1
for (k in K){
  tmp<-knn(train=train, test=test, cl=Labels.train, k=k)
  KNN.mod[[i]]<-tmp
  KNN.fits[i]<-100 * sum(Labels.test == tmp)/NROW(Labels.test) # accuracy rate
  i=i+1
}
plot(K,KNN.fits)
k.opt<-K[which(max(KNN.fits)==KNN.fits)[1]]
KNN.mod.opt<-KNN.mod[[which(max(KNN.fits)==KNN.fits)[1]]]
table(KNN.mod.opt,Labels.test)
Main.res<- confusionMatrix(table(KNN.mod.opt,Labels.test))
Main.res
#
### Results by disease severity
Res.byclass<-as.data.frame(Main.res$byClass)
Res.byclass$severity<-rownames(Res.byclass)
Res.byclass$severity<-factor(Res.byclass$severity)
Res.byclass$severity<-relevel(Res.byclass$severity,ref="Class: healthy")
Res.byclass$severity<-factor(Res.byclass$severity,c("Class: healthy","Class: mild","Class: severe","Class: critical"))

ggplot(Res.byclass,aes(x=severity,y=Precision))+geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 2)+
  labs(x='COVID19 disease severity',y='Prediction accuracy rate')+
  ylim(0,1)
