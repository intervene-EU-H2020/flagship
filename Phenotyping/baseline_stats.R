
###### libraries ######
library(data.table)
library(tidyverse)
library(dplyr)
options(scipen=999)
library(lubridate)

####### variables ######
file_name<-"your/path/to/pheno/file.csv"
biobank<-"biobank_abbrev"
output_dir<-"path/to/output/directory/" #keeping the final "/" is required
#output_dir<-getwd()
drop<-c(NA) #if you are missing a phenotype, list the name here
birth_as_baseline<-FALSE #if you are using birth as the baseline, change to TRUE
ancestry<-"ancestry" #add ancestry EUR, SAS, EAS, AFR which the phenotype file corresponds to

#read in data 
df<-fread(file_name)

#if your SEX column is male/female/other you can comment this out, otherwise, replace 2 or 1 with the values you use for male/female
df$SEX<-recode(df$SEX, `2`="female", `1`="male", .default = NA_character_)

####################
#all 38 phenotypes
#p<-c("C3_CANCER","C3_COLORECTAL","C3_BREAST","T2D","C3_PROSTATE","I9_CHD","I9_SAH","C3_MELANOMA_SKIN","J10_ASTHMA","I9_HEARTFAIL_NS","I9_STR","G6_AD_WIDE","T1D","I9_AF","N14_CHRONKIDNEYDIS","COVID","F5_DEPRESSIO","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","K11_IBD_STRICT","I9_VTE","I9_THAORTANEUR","I9_ABAORTANEUR","COX_ARTHROSIS","KNEE_ARTHROSIS","M13_OSTEOPOROSIS","AUD_SWEDISH","E4_HYTHYNAS","E4_THYTOXGOITDIF","G6_SLEEPAPNO","IPF","ILD","GOUT","H7_GLAUCOMA","G6_EPLEPSY","GE_STRICT","FE_STRICT","K11_APPENDACUT")
###19 phenotypes of interest

p<- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", 
               "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
p<-p[!(p %in% drop)] #drop any phenotypes you don't have to avoid errors later

p<- c("C3_CANCER", "K11_APPENDACUT")

#phenotype file must have "SEX", DATE_OF_BIRTH", "START_OF_FOLLOWUP","END_OF_FOLLOWUP" in YYYY-MM-DD format
#each of the phenotypes above must have a matching column which is the same name with "_DATE" appended

for (idx in 1:length(p)){
  
  if(sum(is.na(pull(df,p[idx])))==nrow(df)){ #if entire case/control designation is NA then the phenotype is missing 
    next
  } else {
    date_col<-paste0(p[idx],"_DATE")
    print(p[idx])
    
    df$DATE_OF_BIRTH<-as.POSIXct(df$DATE_OF_BIRTH)
    df$START_OF_FOLLOWUP<-as.POSIXct(df$START_OF_FOLLOWUP)
    df$END_OF_FOLLOWUP<-as.POSIXct(df$END_OF_FOLLOWUP)
    #set proper baseline
    if (birth_as_baseline==TRUE){
      df$BASELINE<-df$DATE_OF_BIRTH
    } else{
      df$BASELINE<-df$START_OF_FOLLOWUP
    }
    
    #interquartile range for age at onset    
    tmp<-df %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(time_length(difftime(as.POSIXct(get(date_col)), DATE_OF_BIRTH), 'years')))
    if (idx==1){
      age_df<-data.frame(p[idx],t(as.data.frame(quantile(tmp$age, probs = c(.25, .5, .75)))))
    } else{
      age_df<-rbind(age_df,data.frame(p[idx],t(as.data.frame(quantile(tmp$age, probs = c(.25, .5, .75))))))
    }
    
    # of cases and controls
    n_case<-df %>% filter(get(p[idx])==1) %>% nrow()
    n_control<-df %>% filter(get(p[idx])==0) %>% nrow()
    prev<-n_case/(n_case+n_control)*100
    
    #age distribution at recruitment/baseline in CASES (yrs)
    tmp<-df %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric((BASELINE-DATE_OF_BIRTH)/365.5))
    age_recruitment_median_cases<-median(tmp$age,na.rm=TRUE)
    age_recruitment_IQR_cases<-IQR(tmp$age,na.rm=TRUE)
    
    #age distribution at recruitment/baseline in CONTROLS
    tmp<-df %>% filter(get(p[idx])==0) %>% mutate(age=as.numeric((BASELINE-DATE_OF_BIRTH)/365.5))
    age_recruitment_median_controls<-median(tmp$age,na.rm=TRUE)
    age_recruitment_IQR_controls<-IQR(tmp$age,na.rm=TRUE)
    
    #Age distribution at time of recruitment/baseline (yrs) (median, IQR)	
    tmp<-df %>% mutate(age=as.numeric((BASELINE-DATE_OF_BIRTH)/365.5))
    age_recruitment_median<-median(tmp$age,na.rm=TRUE)
    age_recruitment_IQR<-IQR(tmp$age,na.rm=TRUE)
    #we use this tmp later
    
    #Age of onset distribution (median, IQR)	ONLY CASES
    age_onset_median<-df %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(difftime(as.POSIXct(get(date_col)),DATE_OF_BIRTH,units="days"))/365.5) %>% summarize(median(age))
    age_onset_IQR<-df %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(difftime(as.POSIXct(get(date_col)),DATE_OF_BIRTH,units="days"))/365.5)  %>% summarize(IQR(age))
    
    #Distribution of time of follow-up in cohort since baseline (median, IQR)	in years 
    follow_up_median<-df %>% filter(!is.na(get(p[idx]))) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5)) %>% summarize(median(follow,na.rm=TRUE))
    follow_up_IQR<-df %>% filter(!is.na(get(p[idx]))) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5))  %>% summarize(IQR(follow,na.rm=TRUE))
    
    #time of follow up in cases
    follow_up_median_cases<-df %>% filter(get(p[idx])==1) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5))  %>% summarize(median(follow,na.rm=TRUE))
    follow_up_IQR_cases<-df %>% filter(get(p[idx])==1) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5)) %>% summarize(IQR(follow,na.rm=TRUE))
    
    #time of follow up in controls
    follow_up_median_controls<-df %>% filter(get(p[idx])==0) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5))  %>% summarize(median(follow,na.rm=TRUE))
    follow_up_IQR_controls<-df %>% filter(get(p[idx])==0) %>% mutate(follow=as.numeric(difftime(END_OF_FOLLOWUP,BASELINE,units="days")/365.5))  %>% summarize(IQR(follow,na.rm=TRUE))
    
    #correlations 
    age_corr<-cor.test(pull(df,p[idx]),tmp$age)$estimate #has to be age of recruitment because age of onset wouldn't have controls
    age_cor_ci<-cor.test(pull(df,p[idx]),tmp$age,method="pearson")$conf.int
    
    df$SEX_NUM<-recode(df$SEX, female=0, male=1)
    if (p[idx]=="C3_PROSTATE"|p[idx]=="C3_BREAST"){ #logic for sex-specific conditions
      sex_corr<-NA
      sex_corr_ci<-NA
    } else{
      sex_corr<-cor.test(pull(df,p[idx]),df$SEX_NUM,method="pearson")$estimate #positive is correlated with male
      sex_cor_ci<-cor.test(pull(df,p[idx]),df$SEX_NUM,method="pearson")$conf.int
    }
    
    #female percentage
    female<-df[df$SEX_NUM==0,]
    n_female_case<-nrow(df %>% filter(SEX_NUM==0 & get(p[idx])==1))
    n_female_control<-nrow(df %>% filter(SEX_NUM==0 & get(p[idx])==0))
    female_perc_cases<-n_female_case/n_case*100
    female_perc_controls<-n_female_control/n_control*100
    
    ##### proportion of total cases in age of diagnosis bins
    #all 
    tmp<-df %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(difftime(as.POSIXct(get(date_col)),DATE_OF_BIRTH,units="days"))/365.5)
    #lifetimerisks <- data.frame(age=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    tmp$bins<-cut_width(tmp$age,width=5, center=2.5,closed="left") #cut into lifetime risk bins
    total<-as.numeric(tmp %>% summarize(sum(get(p[idx]))))
    counts<- tmp %>% group_by(bins) %>% count(p[idx])
    counts$prop<-counts$n/total*100
    counts$sample<-"all"
    
    #female
    tmp<-df %>% filter(get(p[idx])==1 & SEX_NUM==0) %>% mutate(age=as.numeric(difftime(as.POSIXct(get(date_col)),DATE_OF_BIRTH,units="days"))/365.5)
    tmp$bins<-cut_width(tmp$age,width=5, center=2.5,closed="left") #cut into lifetime risk bins
    total<-as.numeric(tmp %>% summarize(sum(get(p[idx]))))
    counts_f<-tmp %>% group_by(bins) %>% count(p[idx])
    counts_f$prop<-counts_f$n/total*100
    counts_f$sample<-"female"
    
    #male
    tmp<-df %>% filter(get(p[idx])==1 & SEX_NUM==1) %>% mutate(age=as.numeric(difftime(as.POSIXct(get(date_col)),DATE_OF_BIRTH,units="days"))/365.5)
    tmp$bins<-cut_width(tmp$age,width=5, center=2.5,closed="left") #cut into lifetime risk bins
    total<-as.numeric(tmp %>% summarize(sum(get(p[idx]))))
    counts_m<-tmp %>% group_by(bins) %>% count(p[idx])
    counts_m$prop<-counts_m$n/total*100
    counts_m$sample<-"male"
    
    counts<-rbind(counts,counts_f,counts_m)
    
    
    #coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group")))
    
    #] is inclusive, ) is exclusive
    
  }
  if (idx==1){
    summary_stats_df<-data.frame(p[idx],n_case,n_control,prev,age_recruitment_median,age_recruitment_IQR,
                                 age_recruitment_median_cases,age_recruitment_IQR_cases,
                                 age_recruitment_median_controls,age_recruitment_IQR_controls,
                                 age_onset_median,age_onset_IQR,follow_up_median,follow_up_IQR,
                                 age_corr,sex_corr,n_female_case,n_female_control,female_perc_cases,female_perc_controls)
    
    counts_df<-counts
    
  } else{
    summary_stats_df<-rbind(summary_stats_df,data.frame(p[idx],n_case,n_control,prev,age_recruitment_median,age_recruitment_IQR,
                                                              age_recruitment_median_cases,age_recruitment_IQR_cases,
                                                              age_recruitment_median_controls,age_recruitment_IQR_controls,
                                                              age_onset_median,age_onset_IQR,follow_up_median,follow_up_IQR,
                                                              age_corr,sex_corr,n_female_case,n_female_control,female_perc_cases,female_perc_controls))
  
    counts_df<-rbind(counts_df,counts)
    
    }
}

#write age data frame
names(age_df)<-c("trait","X25","X50","X75")
write.csv(format(age_df,digits=3),paste0(output_dir,biobank,"_age_quartiles.csv"),row.names=FALSE,quote=FALSE)

#write counts data frame
names(counts_df)<-c("bins","pheno","count","prop","sample")
write.csv(format(counts_df,digits=3),paste0(output_dir,biobank,"_prevalences.csv"),row.names=FALSE,quote=FALSE)

#plot 
pdf(file=paste0(output_dir,biobank,"proportions.pdf"),height=4,width=12)
ggplot(counts_df,aes(x=bins,y=prop,color=pheno)) + geom_point() + facet_wrap(~sample) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#write summary stats data frame
names(summary_stats_df)<-c("trait","cases","controls","prevalence","age_recruitment_median","age_recruitment_IQR",
                                 "age_recruitment_median_cases","age_recruitment_IQR_cases",
                                 "age_recruitment_median_controls","age_recruitment_IQR_controls",
                                 "age_onset_median","age_onset_IQR","follow_up_median","follow_up_IQR","age_corr",
                                 "sex_corr","n_cases_female","n_controls_female",
                                 "female_prev_cases","female_prev_controls")
write.csv(format(summary_stats_df,digits=3),paste0(output_dir,biobank,"_baseline_summary_stats.csv"),row.names=FALSE,quote=FALSE)


