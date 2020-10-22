
library(stringr)
library(dplyr)

setwd('/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3')
folder="/camp/project/proj-data-challenge/project3/"
Files=list.files(paste(folder,'Data/',sep=''),pattern='.txt')
#
Meta.df=read.csv(paste(folder,'metatable.tsv',sep=''),sep='\t')
patient.df=Meta.df[!is.na(Meta.df$anonym_id),] ## SPlit patient data from commercial controls
QC.df=Meta.df[is.na(Meta.df$anonym_id),]
#
patient.df=patient.df[,c("File.Name","F.Severity","I.WHO","anonym_id")] ## only keep important patient info
#
#Get file names as on CAMP
Files=as.data.frame(Files)
key=unlist(lapply(Files[,1], function(f) str_match(f,"Glyco15min_(.*).wiff.")[,2]))
Files$key=key
patient.df=merge(patient.df,Files,by.x ="File.Name",by.y="key") ## merge file names with patient info
#
IDs=unique(patient.df$anonym_id)
#i=IDs[1]
Df.138<-data.frame('pID'=IDs,'X138max'=rep(0,length(IDs)))
for (i in IDs){
  df.sub=patient.df[patient.df$anonym_id==i,]
  MASS_SPEC=vector('list', nrow(df.sub)) ## collate data from all 3 replicates for each patient 
  names(MASS_SPEC)<-df.sub$Files
  IONS=vector('list', nrow(df.sub))
  names(IONS)<-df.sub$Files
  for (j in df.sub$Files){
    df<-read.csv(paste(folder,'Data/',j,sep=''),sep='')
    df$Window.mid<-rowMeans(df[,c("Window.Low","Window.High")]) # take mid between low and high window
    Ions<-"X138.055" #colnames(df)[5:length(colnames(df))-1]
    MASS_SPEC[[j]]<-df[,c('RT',"Window.mid",Ions)]
    tmp<-sum(df[,c(Ions)])
    IONS[[j]]<-tmp
  }
  IONS<-mean(unlist(IONS))
  Df.138[Df.138$pID==i,'X138max']=IONS
}
Df.138=merge(Df.138,patient.df[,c('anonym_id','I.WHO')],by.x='pID',by.y='anonym_id',all.x=TRUE)
Df.138=unique(Df.138)
###
write.csv(Df.138,'/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3/Extracted_138.csv')
#
