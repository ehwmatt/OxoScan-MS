
library(stringr)
library(dplyr)

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
i=IDs[1]
df.sub=patient.df[patient.df$anonym_id==i,]
MASS_SPEC=vector('list', nrow(df.sub)) ## collate data from all 3 replicates for each patient 
names(MASS_SPEC)<-names(df.sub$Files)
for (j in df.sub$Files){
  df<-read.csv(paste(folder,'Data/',j,sep=''),sep='')
  df$Window.mid<-rowMeans(df[,c("Window.Low","Window.High")]) # take mid between low and high window
  Ions<-colnames(df)[5:length(colnames(df))-1]
  MASS_SPEC[[j]]<-df[,c('RT',"Window.mid",Ions)]
}