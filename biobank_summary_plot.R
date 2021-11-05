library(googlesheets4)
library(googledrive)
library(readxl)
library(rcartocolor)
library(dplyr)
library(ggplot2)
library(stringi)
library(tidyverse)
library(reshape2)
library(ggrepel)

#I think this doesn't work because it's a .xlsx
#file<-"https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit?usp=sharing&ouid=115040676596308698055&rtpof=true&sd=true"

#download file to use 
file<-drive_find(pattern="Intervene_flagship_endpoint_collection.xlsx") #takes a bit of time..would be better to know path
drive_download(file$name,overwrite=TRUE)
#will have to authorize your google account, check option to see accounts, wil copy and apste an authorization code 

#read in names
names<-read_xlsx(file$name,sheet=1)
names2<-names[c("Endpoint","Endpoint family","PhenoID")]
my_colors = carto_pal(length(unique(names2$`Endpoint family`)), "Safe")
names2$Group<-names2$`Endpoint family`
names2$`Endpoint family`<-as.factor(names2$`Endpoint family`)
levels(names2$`Endpoint family`)<-my_colors #set colorblind frienldy colors

#list sheets
sheets<-excel_sheets(file$name)
indices<-grep("Metrics",sheets) #select sheets with metrics

##### TO DO: fix IBD, IBS, osteo and hip in the sheets

drop_rows_all_na <- function(x, pct=0.90) x[!rowSums(is.na(x)) >= ncol(x)*pct,]

### function to read and parse each metric sheet into consistent formatting
my_read_func<-function(sheet){
  df<-read_xlsx(file$name,sheet=sheet,na=c("NA","na"))
  df<-drop_rows_all_na(df) #handle finngen formatting with merged columns
  df<-df[c(1:39),c(1:10)] #cut extra columns that some sheets have, manually checked that base 10 are the same (moved extra to the end)
  names(df)<-c("endpoint","def","cases","controls","age_baseline","follow_dist","age_onset_dist","age_corr","sex_corr","female")
  
  print(sheet)
  
  df$cases<-as.numeric(df$cases) 
  df$controls<-as.numeric(df$controls)
  
  #format the IQR columns
  df$age_baseline_median<-trimws(sapply(stri_split_regex(as.character(df$age_baseline),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_baseline_iqr<- trimws(sapply(stri_split_regex(as.character(df$age_baseline),pattern="\\s|\n|\\(",n=2),"[",2))
  
  df$age_onset_dist_median<-trimws(sapply(stri_split_regex(as.character(df$age_onset_dist),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_onset_dist_iqr<-trimws(sapply(stri_split_regex(as.character(df$age_onset_dist),pattern="\\s|\n|\\(",n=2),"[",2))
 
  df$follow_dist_median<-trimws(sapply(stri_split_regex(as.character(df$follow_dist),pattern="\\s|\n|\\(",n=2),"[",1))
  df$follow_dist_iqr<-trimws(sapply(stri_split_regex(as.character(df$follow_dist),pattern="\\s|\n|\\(",n=2),"[",2))
  
  df$sex_corr<-trimws(sapply(stri_split_regex(as.character(df$sex_corr),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_corr<-trimws(sapply(stri_split_regex(as.character(df$age_corr),pattern="\\s|\n|\\(",n=2),"[",1)) #pearson?
  
  df$female_est<-as.numeric(trimws(sapply(stri_split_regex(as.character(df$female),pattern="\\s|\n|\\(",n=2),"[",1)))
  df$female_ci<-trimws(sapply(stri_split_regex(as.character(df$female),pattern="\\s|\n|\\d\\(",n=2),"[",2))
  
  df$biobank<-unlist(strsplit(sheets[sheet],"-"))[[2]] #label biobank
  
  df<-df %>% select("endpoint","cases","controls","age_baseline_median","age_baseline_iqr","age_onset_dist_median","age_onset_dist_iqr","follow_dist_median","follow_dist_iqr","age_corr","sex_corr","female_est","female_ci","biobank")
  return(df)
}
#loop over sheets with lapply
list_of_biobanks<-lapply(indices,my_read_func)

#do some cleaning of the dataframe
dat<-bind_rows(list_of_biobanks) %>% left_join(names2,by=c("endpoint"="Endpoint")) 
dat<-dat[!is.na(dat$endpoint),] #remove when endpoint is NA
dat<-dat[!is.na(dat$Group),]
dat<-dat[dat$endpoint!="BMI (obesity)",]#remove BMI since it is quant
dat$endpoint<-as.factor(dat$endpoint)
dat$Group<-as.factor(dat$Group)

######### number of cases 
options(scipen=10)
pdf(file="cases.pdf",height=10,width=12)
ggplot(dat,aes(y=endpoint,x=cases)) + geom_bar(aes(fill = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_fill_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + theme_bw() + 
  theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

#prevalence 
dat$prev<-dat$cases/(dat$cases+dat$controls)
pdf(file="prevalence.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$prev),],aes(y=endpoint,x=prev)) + geom_bar(aes(fill = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_fill_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + theme_bw() + labs(xlab="Prevalence") +
  theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

###leave one biobank vs rest for prevalence
prev_dat<-dat[c("endpoint","biobank","Group","Endpoint family","prev")]
prev_dat2<-prev_dat %>% left_join(prev_dat,by="endpoint") %>% subset(biobank.x!=biobank.y) %>% drop_na %>% mutate(label=paste(sep="\n",endpoint,biobank.y))
pdf(file="prevalence_LOBO.pdf",height=10,width=12)
ggplot(prev_dat2,aes(y=prev.y,x=prev.x,label=label)) + geom_point(aes(color=`Endpoint family.y`)) + facet_wrap(~biobank.x,nrow=2) + 
  labs(x="Leave one biobank out",y="Remaining biobanks",main="Prevalence") +  theme_bw() + geom_text_repel(size=2,max.overlaps=12) +
  theme(legend.position="bottom")  +  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") +   
  scale_color_manual(values = levels(prev_dat2$`Endpoint family.y`),labels=levels(prev_dat2$Group.y))
dev.off()

pdf(file="prevalence_matrix.pdf",height=10,width=12)
ggplot(prev_dat2,aes(y=prev.y,x=prev.x,label=endpoint)) + geom_point(aes(color=`Endpoint family.y`)) + 
  facet_grid(biobank.x~biobank.y,scales="free") +
  labs(main="Prevalence",x=" ",y=" ") +  theme_bw() + geom_text_repel(size=2,max.overlaps=15) +
  theme(legend.position="bottom")  +  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") +   
  scale_color_manual(values = levels(prev_dat2$`Endpoint family.y`),labels=levels(prev_dat2$Group.y))
dev.off()

pdf(file="prevalence_matrix_log.pdf",height=10,width=12)
ggplot(prev_dat2,aes(y=prev.y,x=prev.x,label=endpoint)) + geom_point(aes(color=`Endpoint family.y`)) + 
  facet_grid(biobank.x~biobank.y,scales="free")+ scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
  labs(main="Prevalence",x="-log10(prevalence of biobank y)",y="-log10(prevalence of biobank x") +  theme_bw() + geom_text_repel(size=2,max.overlaps=10) +
  theme(legend.position="bottom")  +  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") +   
  scale_color_manual(values = levels(prev_dat2$`Endpoint family.y`),labels=levels(prev_dat2$Group.y))
dev.off()

pdf(file="prevalence_matrix_log_fit.pdf",height=10,width=12)
ggplot(prev_dat2,aes(y=prev.y,x=prev.x,label=endpoint)) + geom_point(aes(color=`Endpoint family.y`)) + 
  facet_grid(biobank.x~biobank.y,scales="free")+ scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
  labs(main="Prevalence",x="-log10(prevalence of biobank y)",y="-log10(prevalence of biobank x") +  theme_bw() + geom_text_repel(size=2,max.overlaps=10) +
  theme(legend.position="bottom")  +  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") +   geom_smooth(method = "lm",alpha=0.5)  +
  scale_color_manual(values = levels(prev_dat2$`Endpoint family.y`),labels=levels(prev_dat2$Group.y))
dev.off()


## median follow up 
dat$follow_dist_low<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$follow_dist_iqr),pattern="-",n=2),"[",1)))
dat$follow_dist_high<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$follow_dist_iqr),pattern="-",n=2),"[",2)))
pdf(file="followup.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$follow_dist_median),],aes(y=endpoint,x=as.numeric(follow_dist_median))) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + 
  theme_bw() + labs(x="Median and interquartile range of follow-up from baseline (yrs)") + geom_linerange(aes(xmin=follow_dist_low,xmax=follow_dist_high,color=`Endpoint family`)) +
  theme(legend.position="bottom")  + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

##median age of onset
dat$age_onset_dist_low<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$age_onset_dist_iqr),pattern="-",n=2),"[",1)))
dat$age_onset_dist_high<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$age_onset_dist_iqr),pattern="-",n=2),"[",2)))
pdf(file="onset.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$age_onset_dist_median),],aes(y=endpoint,x=as.numeric(age_onset_dist_median))) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) +   theme_bw() + labs(x="Median age of onset and interquartile range") +
  geom_linerange(aes(xmin=age_onset_dist_low,xmax=age_onset_dist_high,color=`Endpoint family`)) +
theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

##### % female
dat$female_est_perc<-as.numeric(ifelse(dat$female_est<=1,dat$female_est*100,dat$female_est))
pdf(file="female.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$female_est_perc),],aes(y=endpoint,x=female_est_perc)) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + theme_bw() + labs(x="% Female") + geom_vline(xintercept=50,linetype="dashed",color="red") +
  theme(legend.position="bottom") +  scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

###correlations
dat_long<-melt(dat[,c("endpoint","age_corr","sex_corr","biobank","Endpoint family","Group")],id.vars=c("endpoint","biobank","Endpoint family","Group"))
dat_long$value<-as.numeric(dat_long$value)
dat_long$value<-ifelse(dat_long$biobank=="GNH" & dat_long$variable=="sex_corr",dat_long$value*1,dat_long$value*-1) #fix GNH because correlation is flipped for sex compared to the rest 
pdf(file="corr.pdf",height=10,width=12)
ggplot(dat_long[!is.na(dat_long$value),],aes(y=endpoint,x=value,shape=variable)) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) + geom_vline(linetype="dashed",color="red",xintercept=0) +
  scale_color_manual(values = levels(dat_long$`Endpoint family`), labels = levels(dat_long$Group)) + theme_bw()  + theme(legend.position="bottom") +  scale_y_discrete(limits=rev(levels(dat_long$endpoint)))
dev.off()

#### leave one biobank vs rest for age correlations
dat_long_age<-dat_long[dat_long$variable=="age_corr",]
corr_age<-dat_long_age %>% left_join(dat_long_age,by="endpoint") %>% subset(biobank.x!=biobank.y) %>% drop_na %>% mutate(label=paste(sep="\n",endpoint,biobank.y))
pdf(file="corr_age.pdf",height=12,width=12)
ggplot(corr_age,aes(x=value.x,y=value.y,shape=biobank.y,color=Group.x,label=label)) + geom_point() + theme_bw() + 
  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") + facet_wrap(~biobank.x) + 
  labs(title="Endpoint & Age correlations",x="Leave one out biobank correlation",y="Remaining biobank correlations") +
  geom_text_repel(size=3) + theme(legend.position="bottom") +  scale_shape_discrete(name="Biobank") +
  geom_vline(xintercept=0,color="black",alpha=0.5) +  geom_hline(yintercept=0,color="black",alpha=0.5) +
  scale_color_manual(values = levels(corr_age$`Endpoint family.y`), labels = levels(corr_age$Group.y),name="Endpoint Family") 
dev.off()

#### matrix for age correlations
pdf(file="corr_age_matrix.pdf",height=12,width=12)
ggplot(corr_age,aes(x=value.x,y=value.y,shape=biobank.y,color=Group.x,label=endpoint)) + geom_point() + theme_bw() + 
  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") + facet_grid(biobank.x~biobank.y) +
  labs(title="Endpoint & Age correlations") +
  geom_text_repel(size=2,max.overlaps=10) + theme(legend.position="bottom") + scale_shape_discrete(name="Biobank on x-axis") +
  geom_vline(xintercept=0,color="black",alpha=0.5) +  geom_hline(yintercept=0,color="black",alpha=0.5) +
  scale_color_manual(values = levels(corr_age$`Endpoint family.y`), labels = levels(corr_age$Group.y),name="Endpoint Family") 
dev.off()

#### leave one biobank vs rest for sex correlations
dat_long_sex<-dat_long[dat_long$variable=="sex_corr",]
corr_sex<-dat_long_sex %>% left_join(dat_long_sex,by="endpoint") %>% subset(biobank.x!=biobank.y) %>% drop_na  %>% mutate(label=paste(sep="\n",endpoint,biobank.y))
pdf(file="corr_sex.pdf",height=12,width=12)
ggplot(corr_sex,aes(x=value.x,y=value.y,shape=biobank.y,color=Group.x,label=label)) + geom_point() + theme_bw() + 
  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") + facet_wrap(~biobank.x) + 
  labs(title="Endpoint & Sex correlations",x="Leave one out biobank correlation",y="Remaining biobank correlations") +
  geom_text_repel(size=3) + theme(legend.position="bottom") + scale_shape_discrete(name="Biobank") +
  geom_vline(xintercept=0,color="black",alpha=0.5) +  geom_hline(yintercept=0,color="black",alpha=0.5) +
  scale_color_manual(values = levels(corr_age$`Endpoint family.y`), labels = levels(corr_age$Group.y),name="Endpoint Family") 
dev.off()

#### matrix for sex correlations
pdf(file="corr_sex_matrix.pdf",height=12,width=12)
ggplot(corr_sex,aes(x=value.x,y=value.y,shape=biobank.y,color=Group.x,label=endpoint)) + geom_point() + theme_bw() + 
  geom_abline(slope=1,intercept=0, color="black",linetype="dashed") + facet_grid(biobank.x~biobank.y) +
  labs(title="Endpoint & Sex correlations") +
  geom_text_repel(size=2,max.overlaps=10) + theme(legend.position="bottom") + scale_shape_discrete(name="Biobank on x-axis") +
  geom_vline(xintercept=0,color="black",alpha=0.5) +  geom_hline(yintercept=0,color="black",alpha=0.5) +
  scale_color_manual(values = levels(corr_age$`Endpoint family.y`), labels = levels(corr_age$Group.y),name="Endpoint Family") 
dev.off()


##### total cases across biobanks
tot<-read_xlsx(file$name,sheet=sheets[sheets=="Case-Totals"],na=c("NA","na"))
tot<-tot[tot$Endpoints!="BMI (obesity)",]
tot<-tot %>% left_join(names2,by=c("Endpoints"="Endpoint")) 
tot$Endpoints<-as.factor(tot$Endpoints)
tot$Group<-as.factor(tot$Group)
pdf(file="all_biobanks_cases.pdf",height=6,width=10)
ggplot(tot,aes(y=Endpoints,x=`N cases`)) + geom_bar(aes(fill = `Endpoint family`),stat="identity") +
  scale_fill_manual(values = levels(tot$`Endpoint family`), labels = levels(tot$Group)) + theme_bw() + 
  theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(tot$Endpoints)))
dev.off()



#would be good to plot the disease prevalence across biobank one vs the others. Like scatter plots.
#and also the correlation with age and sex
#to understand the consistency across endpoints

