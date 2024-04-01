library(data.table)

files<-list.files(pattern="*parameters")
for (i in 1:length(files)){
  if(i==1){
    df<-fread(files[i])
    df$trait<-gsub("_bestmodel_parameters","",files[i])
  } else{
    df2<-fread(files[i])
    df2$trait<-gsub("_bestmodel_parameters","",files[i])
    df<-rbind(df,df2)
  }}

fwrite(df,"flagship_best_model_parameters.tsv",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
