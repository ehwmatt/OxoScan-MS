
library(stringr)
library(dplyr)

folder="/camp/project/proj-data-challenge/project3/"
Files=list.files(paste(folder,'Data/',sep=''),pattern='.txt')
#
Meta.df=read.csv(paste(folder,'metatable.tsv',sep=''),sep='\t')
patient.df=Meta.df[!is.na(Meta.df$anonym_id),]
QC.df=Meta.df[is.na(Meta.df$anonym_id),]
#
patient.df=patient.df[,c("File.Name","F.Severity","I.WHO","anonym_id")]
#
Files=as.data.frame(Files)
key=unlist(lapply(Files[,1], function(f) str_match(f,"Glyco15min_(.*).wiff.")[,2]))
Files$key=key
patient.df=merge(patient.df,Files,by.x ="File.Name",by.y="key")
#
IDs=unique(patient.df$anonym_id)
i=IDs[1]
df.sub=patient.df[patient.df$anonym_id==i,]
MASS_SPEC=vector('list', nrow(df.sub))
names(MASS_SPEC)<-names(df.sub$Files)
for (j in df.sub$Files){
  df<-read.csv(paste(folder,'Data/',j,sep=''),sep='')
}