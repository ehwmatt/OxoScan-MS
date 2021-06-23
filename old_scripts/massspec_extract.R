
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
#
Df.138<-vector('list',length(IDs))
for (i in IDs){
  df.sub=patient.df[patient.df$anonym_id==i,]
  MASS_SPEC=vector('list', nrow(df.sub)) ## collate data from all 3 replicates for each patient 
  names(MASS_SPEC)<-df.sub$Files
  n<-1
  for (j in df.sub$Files){
    df<-read.csv(paste(folder,'Data/',j,sep=''),sep='')
    df$Window.mid<-rowMeans(df[,c("Window.Low","Window.High")]) # take mid between low and high window
    Ions<-"X138.055" #colnames(df)[5:length(colnames(df))-1]
    df.keep<-df[,c('RT',"Window.mid",Ions)]
    newcol<-paste('RTint138',n,sep='')
    df.new=data.frame('MZ'=unique(df.keep$Window.mid))
    for( mz in unique(df.keep$Window.mid)){
      dfs<-df.keep[df.keep$Window.mid==mz,]
      tmp<-integrate.xy(dfs$RT, dfs$X138.055)#
      df.new[df.new$MZ==mz,newcol]<-tmp
    }
    MASS_SPEC[[j]]<-df.new
    n<-n+1
  }
  IONS<-merge(MASS_SPEC[[1]],MASS_SPEC[[2]],by='MZ',all=TRUE)
  IONS<-merge(IONS,MASS_SPEC[[3]],by='MZ',all=TRUE)
  IONS$RTint1381[is.na(IONS$RTint1381)]<-IONS$RTint1382[is.na(IONS$RTint1381)]
  IONS$RTint1381[is.na(IONS$RTint1381)]<-IONS$RTint1383[is.na(IONS$RTint1381)]
  IONS$pID<-i
  Df.138[[i]]=IONS[,c('MZ',"RTint1381","pID")]
}
#
Df.138.all<-do.call(rbind, Df.138)
#
Df.138.all=merge(Df.138.all,patient.df[,c('anonym_id','I.WHO',"F.Severity")],by.x='pID',by.y='anonym_id',all.x=TRUE)
Df.138.all=unique(Df.138.all)
###
write.csv(Df.138.all,'/camp/stp/babs/working/schneid/Data_Challenge_2020/2020-Project3/Extracted_138.csv')
#
