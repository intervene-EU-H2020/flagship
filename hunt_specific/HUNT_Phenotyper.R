
#/apps/statistics2/R-4.0.0/bin/R

setwd("~/intervene")

###### libraries ######
library(data.table)
library(tidyverse)
library(foreign)
library(stringr)
library(dplyr)
library(readstata13)
options(scipen=999)
library(haven)
#remotes::install_github("bayesiandemography/demprep")
library(demprep)
library(lubridate)

###### variables #######
endpoint_file<-"../flagship/Phenotyping/UKBB_definitions_demo_TEST.csv"
master_file<-"/mnt/work/master/DATASET_20170512/SAMPLE_QC/Masterkey_DATASET.20170512.txt.gz"
bridge_file<-"/mnt/work/bridge/allin-phecode-2018_41492/PID@108485-PID@105118.sav"
fam_file<-"/mnt/scratch/brooke/bcf/all.log.fam"

files<-list.files("/mnt/work/phenotypes/allin-phecode-2018_41492/kilde/hnt/",full.names=TRUE)
files<- files[!grepl("etter",files)] #what is issue with that sav file?

#/home/benb/work/phenotypes/allin-phecode-2018_41492/kilde/hnt

#Read in the information for each phenotype
endpoints <- fread(endpoint_file, data.table=FALSE)
head(endpoints)

#Transforms input into regex format for detecting ICD_10, ICD_9 and ATC codes that will be recognised by the following code using grepl.
endpoints$ICD_10 <- ifelse(!is.na(endpoints$COD_ICD_10), paste0("^10x(",endpoints$COD_ICD_10,")"), NA)
endpoints$ICD_10_EXC <- ifelse(!is.na(endpoints$COD_ICD_10_EXC), paste0("^10x(",endpoints$COD_ICD_10_EXC,")"), NA)
endpoints$ICD_9 <- ifelse(!is.na(endpoints$COD_ICD_9), paste0("^9x(",endpoints$COD_ICD_9,")"), NA)  
endpoints$ICD_9_EXC <- ifelse(!is.na(endpoints$COD_ICD_9_EXC), paste0("^9x(",endpoints$COD_ICD_9_EXC,")"), NA) 
endpoints$ATC <- ifelse(!is.na(endpoints$KELA_ATC), paste0("^ATCx(",endpoints$KELA_ATC,")"),NA) 
endpoints$INCLUDE <- ifelse(!is.na(endpoints$INCLUDE), paste0("(", endpoints$INCLUDE, ")"), NA)

#Creates a column 'CODE' which combines the regex format of the three previous separate regex. 
endpoints$CODE <- NA
endpoints$CODE_EXC <- NA

head(endpoints)

for(i in 1:nrow(endpoints)){
  codes <- c(endpoints$ICD_10[i], endpoints$ICD_9[i], endpoints$ATC[i])
  codes <- codes[!is.na(codes)]
  codes <- paste(codes, collapse="|")
  endpoints$CODE[i] <- codes
  
  exccodes <- c(endpoints$ICD_10_EXC[i], endpoints$ICD_9_EXC[i])
  exccodes <- exccodes[!is.na(exccodes)]
  exccodes <- paste(exccodes, collapse="|")
  endpoints$CODE_EXC[i] <- exccodes
}

head(endpoints)

#Data file contains the ID of the person, the date of the diagnosis (in numeric format) and the ICD CODE (CODE1). 
#dat <- fread("zcat ../tuomo/UKBB_longitudinal.txt", data.table=FALSE)
#head(dat)

#function to format file as above
make_file<-function(file){
  #loop over .sav file s
  data<-read_sav(file)
  print(file)
  names(data)[1]<-"ID"
  data$ID<-as.character(data$ID)
  print(names(data))
  #diagnosis date in POSIXct format (does this need to be numeric?)
  if (TRUE %in% grepl("diagnosedato",names(data))){
    data$dx<-as.POSIXct(as.Date(paste0(data$YYYYMM_diagnosedato,"01"),format='%Y%m%d'),origin="1970-01-01")
  }else if (TRUE %in% grepl("diagnosekode",names(data))){
    data$dx<-as.POSIXct(as.Date(paste0(data$YYYYMM_diagnosekode,"01"),format='%Y%m%d'),origin="1970-01-01")
  }else if (TRUE %in% grepl("diagnose",names(data))){
    data$dx<-as.POSIXct(as.Date(paste0(data$YYYYMM_diagnose,"01"),format='%Y%m%d'),origin="1970-01-01")
  }
  #trim white space on ICD code entries and add prefix
  data$ICD9<-paste0("9x",trimws(data$ICD9))
  data$ICD10<-paste0("10x",trimws(data$ICD10))
  #long form of data 
  datl<-pivot_longer(data,cols=c("ICD9","ICD10"))[c("ID","dx","value")]

  datl<-datl[!datl$value=="10x"|datl$value=="9x",] #drop empty columns
  return(datl)
}
dat<-bind_rows(lapply(files,make_file))
names(dat)<-c("ID","DATE","CODE1")

#Phenotypes have been structured into levels to account for the inherent hierarchy. 
#Level 1 phenotypes can be defined using ICD codes. 
#Level 2 phenotypes require level 1 phenotypes to be defined. 
#Level 3 phenotypes require level 2 phenotypes to be defined. 


#For loop which iterates over each endpoint and extracts all instances where the code can be found in the longitudinal dataset. Extracts the earliest date for each participant. 
dataset <- data.frame(NULL)
for(i in c(1:nrow(endpoints))){
  endpoint<-endpoints[i,]
  print(endpoint$NAME)
  cases <- subset(dat, grepl(endpoint$CODE, dat$CODE1) & endpoint$CODE!="")
  cases <- cases %>% 
              group_by(ID) %>%
                arrange(DATE) %>%
                  slice(1L)
  exclusions <- subset(dat, grepl(endpoint$CODE_EXC, CODE1) & endpoint$CODE_EXC!="")
  exclusions <- exclusions[!duplicated(exclusions$ID),]
  cases <- subset(cases, !(ID %in% exclusions$ID))
  cases$ENDPOINT_LEVEL_1 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  dataset <- rbind(dataset, cases)
}

endpoints_level2 <- subset(endpoints, LEVEL==2) 
dataset$ENDPOINT_LEVEL_2 <- as.character(NA)

for(i in c(1:nrow(endpoints_level2))){
  endpoint <- endpoints_level2[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME)
  cases <- cases %>% 
              group_by(ID) %>%
                arrange(DATE) %>%
                  slice(1L)
  cases$ENDPOINT_LEVEL_2 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  
  exclusions <- subset(dataset, grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1))
  exclusions <- exclusions[!duplicated(exclusions$ID),]
  cases <- subset(cases, !(ID %in% exclusions$ID)) 
  dataset <- rbind(dataset, cases)
}

endpoints_level3 <- subset(endpoints, LEVEL==3) 
dataset$ENDPOINT_LEVEL_3 <- as.character(NA)

for(i in 1:nrow(endpoints_level3)){
  endpoint <- endpoints_level3[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME | grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_2) | ENDPOINT_LEVEL_2 == endpoint$NAME)
  cases <- cases %>% 
              group_by(ID) %>%
                  arrange(DATE) %>%
                      slice(1L)
  cases$ENDPOINT_LEVEL_3 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  exclusions <- subset(dataset, grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_2)) 
  exclusions <- exclusions[!duplicated(exclusions$ID),]
  cases <- subset(cases, !(ID %in% exclusions$ID)) 
  ####note all cancers had a problem due to level 3 typo in endpoint$EXCLUDE, should not be endpoints
  dataset <- rbind(dataset, cases)
}

#Save dataset
write.csv(dataset, 'endpointsDefinitionsHUNT.csv',row.names=FALSE)

##########################################################################################################################################################################################

#Aim for this section of code is to produce a wide table which will create an array of all phenotypes that belong to a person.
#No information is provided on earliest date for this table. 

endpoints <- fread('endpointsDefinitionsHUNT.csv', data.table=FALSE)
endpoints$ID<-as.character(endpoints$ID)
uniqueIDs <- as.data.frame(endpoints$ID[!duplicated(endpoints$ID)])
colnames(uniqueIDs) <- 'ID'

level1 <- subset(endpoints, LEVEL==1)
level1 <- select(level1, c("ID","DATE","CODE1","ENDPOINT_LEVEL_1", "LEVEL"))
for(i in unique(level1$ENDPOINT_LEVEL_1)){
  print(i)
  case <- subset(level1, i==ENDPOINT_LEVEL_1)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('ID',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case) #joins by ID
}

level2 <- subset(endpoints, LEVEL==2 & !is.na(ENDPOINT_LEVEL_2))
level2 <- select(level2, c("ID","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_2"))
for(i in unique(level2$ENDPOINT_LEVEL_2)){
  print(i)
  case <- subset(level2, i==ENDPOINT_LEVEL_2)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('ID',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

level3 <- subset(endpoints, LEVEL==3 & !is.na(ENDPOINT_LEVEL_3))
level3 <- select(level3, c("ID","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_3"))
for(i in unique(level3$ENDPOINT_LEVEL_3)){
  print(i)
  case <- subset(level3, i==ENDPOINT_LEVEL_3)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('ID',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

#Select endpoints of interest for flagship project
endpointswide <- select(uniqueIDs, contains(c("ID","C3_CANCER", "C3_COLORECTAL", "C3_BREAST", "T2D", "C3_PROSTATE", "I9_CHD", "I9_SAH", "C3_MELANOMA_SKIN", 
                                     "J10_ASTHMA", "I9_HEARTFAIL_NS", "I9_STR", "G6_AD_WIDE", "T1D", "I9_AF", "N14_CHRONKIDNEYDIS", "F5_DEPRESSIO", 
                                     "C3_BRONCHUS_LUNG", "RHEUMA_SEROPOS_OTH", "K11_IBD_STRICT", "I9_VTE", "I9_THAORTANEUR", "I9_ABAORTANEUR", 
                                     "COX_ARTHROSIS", "KNEE_ARTHROSIS", "M13_OSTEOPOROSIS", "AUD_SWEDISH", "E4_HYTHYNAS", "G6_SLEEPAPNO", "IPF", 
                                     "ILD", "GOUT", "H7_GLAUCOMA", "G6_EPLEPSY", "GE_STRICT", "FE_STRICT", "K11_APPENDACUT", "COVID",
                                     "C3_CANCER_DATE", "C3_COLORECTAL_DATE", "C3_BREAST_DATE", "T2D_DATE", "C3_PROSTATE_DATE", "I9_CHD_DATE", "I9_SAH_DATE", "C3_MELANOMA_SKIN_DATE", "J10_ASTHMA_DATE", 
                                     "I9_HEARTFAIL_NS_DATE", "I9_STR_DATE", "G6_AD_WIDE_DATE", "T1D_DATE", "I9_AF_DATE", "N14_CHRONKIDNEYDIS_DATE", "F5_DEPRESSIO_DATE", "C3_BRONCHUS_LUNG_DATE", 
                                     "RHEUMA_SEROPOS_OTH_DATE", "K11_IBD_STRICT_DATE", "I9_VTE_DATE", "I9_THAORTANEUR_DATE", "I9_ABAORTANEUR_DATE", "COX_ARTHROSIS_DATE", "KNEE_ARTHROSIS_DATE", 
                                     "M13_OSTEOPOROSIS_DATE", "AUD_SWEDISH_DATE", "E4_HYTHYNAS_DATE", "G6_SLEEPAPNO_DATE", "IPF_DATE", "ILD_DATE", "GOUT_DATE", "H7_GLAUCOMA_DATE", "G6_EPLEPSY_DATE", 
                                     "GE_STRICT_DATE", "FE_STRICT_DATE", "K11_APPENDACUT_DATE", "COVID_DATE")))
#NOTE: added contains because some of the phenotypes aren't working (AAA), need to figure out why
write.csv(uniqueIDs, 'endpointsWideFormatHUNT.csv',row.names=FALSE)
### note: editing seropositive in the definitions file to include M05.8 and .9 becuase not getting cases otherwise

#There will be a group who have received no ICD codes. 
#These participants should be included as controls also.
#Link these remaining people to the saved dataset above and include as controls. 
#Note: see below 

########## Make phenotype file 
#MISSING VALUES: missing values should be denoted with -
#NOT APPLICABLE FIELDS: if a field is not applicable (for example due to gender),
#this should be denoted with NA.
#``
######## https://docs.google.com/document/d/1GbZszpPeyf-hyb0V_YDx828YbM7woh8OBJhvzkEwo2g/edit

ew<-read.csv('endpointsWideFormatHUNT.csv')
ew$ID<-as.character(ew$ID)
#what IDs are not in uniqueIDs but are in the master file?
mdf<-fread(master_file) #70517
mdf$gid.current<-as.character(mdf$gid.current)
bdf<-read_sav(bridge_file) #79058
names(bdf)<-c("PID108485","PID105118")

#bring in HUNT data for BMI
hunt_file<-"/mnt/work/phenotypes/allin-phecode-2018_41492/kilde/hunt/2021-11-08_Data_delivery_1084851/2021-11-08_108485_Data.sav"
hunt<-read_sav(hunt_file)
hunt_file<-"/mnt/work/phenotypes/allin-phecode-2018_41492/kilde/hunt/2019 02 27 Datautlevering/2019-02-27_108485_Data.sav"
#do we have BMI for more than HUNT 3?

### impute DOB from birth year 
set.seed(1234)
hunt$DATE_OF_BIRTH<-impute_dob(as.Date(paste0(hunt$BirthYear+1,"01","01"),format='%Y%m%d'),age_years=0) #as.Date defaults to UTC
hunt$DATE_OF_BIRTH<-as.POSIXct(hunt$DATE_OF_BIRTH)

#identify start of follow up ie recruitment 
hunt<-hunt %>% filter(!(is.na(`PartAg@NT1BLQ1`) & is.na(`PartAg@NT2BLQ1`) & is.na(`PartAg@NT3BLQ1`))) #drop people in hospital records not in HUNT
hunt<-hunt %>% rowwise() %>% mutate(START_OF_FOLLOWUP=min(`PartAg@NT1BLQ1`,`PartAg@NT2BLQ1`,`PartAg@NT3BLQ1`,na.rm=TRUE))
hunt$START_OF_FOLLOWUP<-ceiling_date(ymd(hunt$DATE_OF_BIRTH) + dyears(hunt$START_OF_FOLLOWUP),unit="days")

#identify end of follow up 
#last linking date (i.e. the last date when the person was still known to be alive) or with the date of death.
#could use date that EHR was pulled
last_date<-dat %>% group_by(ID) %>% mutate(max=max(DATE)) %>% select(ID,max) %>% unique()
hunt<-hunt %>% left_join(last_date,by=c("PID@108485"="ID")) %>% rename(END_OF_FOLLOWUP=max) %>% select(-Sex)

#identify genotyped samples
fam<-fread(fam_file)

#smoking
hunt$SMOKING<-NA

#education 
hunt$EDUCATION_97<-NA
hunt$EDUCATION_11<-NA

#merge all
df<-left_join(ew, last_date,by="ID") %>% 
  full_join(bdf,by=c("ID"="PID108485")) %>% 
  left_join(mdf,by=c("PID105118"="gid.current")) %>% 
  left_join(hunt,by=c("ID"="PID@108485")) %>% 
  right_join(fam,by=c("IID"="V1"))

#delted RHEUM SERO POS doesn't exist,  "I9_THAORTANEUR","I9_ABAORTANEUR",  "E4_HYTHYNAS","G6_SLEEPAPNO",IPF"GE_STRICT","FE_STRICT","COVID"

###get COVID-19 ICD codes from COVID HGI data, not in these files from the hospital

#convert numeric dates back to date format
date_columns<-grep("_DATE",colnames(df))
df<-df %>% mutate_at(all_of(date_columns),~as.Date(as.POSIXct(.x,origin="1970-01-01")))

header<-c("ID","SEX","DATE_OF_BIRTH",	"PC1",	"PC2",	"PC3",	"PC4","PC5",	
          "PC6","PC7","PC8","PC9","PC10","ANCESTRY",
          "C3_CANCER",	"C3_COLORECTAL",	
          "C3_BREAST",	"T2D",	"C3_PROSTATE",	"I9_CHD",	"I9_SAH",	"C3_MELANOMA_SKIN",	"J10_ASTHMA",	"I9_HEARTFAIL_NS",	"I9_STR",	"G6_AD_WIDE",	
          "BMI",	"T1D",	"I9_AF",	"N14_CHRONKIDNEYDIS",	"COVID",	"F5_DEPRESSIO",	"C3_BRONCHUS_LUNG",	"RHEUMA_SEROPOS_OTH",	"K11_IBD_STRICT",	"I9_VTE",	
          "I9_THAORTANEUR",	"I9_ABAORTANEUR",	"COX_ARTHROSIS",	"KNEE_ARTHROSIS",	"M13_OSTEOPOROSIS",	"AUD_SWEDISH",	"E4_HYTHYNAS",	"E4_THYTOXGOITDIF",	"G6_SLEEPAPNO",	"IPF",	
          "ILD",	"GOUT",	"H7_GLAUCOMA",	"G6_EPLEPSY",	"GE_STRICT",	"FE_STRICT",	"K11_APPENDACUT",	"C3_CANCER_DATE",	"C3_COLORECTAL_DATE",	"C3_BREAST_DATE",	"T2D_DATE",	"C3_PROSTATE_DATE",	
          "I9_CHD_DATE",	"I9_SAH_DATE",	"C3_MELANOMA_SKIN_DATE",	"J10_ASTHMA_DATE",	"I9_HEARTFAIL_NS_DATE",	"I9_STR_DATE",	"G6_AD_WIDE_DATE",	"BMI_DATE",	"T1D_DATE",	"I9_AF_DATE",	"N14_CHRONKIDNEYDIS_DATE",	
          "COVID_DATE",	"F5_DEPRESSIO_DATE",	"C3_BRONCHUS_LUNG_DATE",	"RHEUMA_SEROPOS_OTH_DATE",	"K11_IBD_STRICT_DATE",	"I9_VTE_DATE",	"I9_THAORTANEUR_DATE",	"I9_ABAORTANEUR_DATE",	"COX_ARTHROSIS_DATE",	
          "KNEE_ARTHROSIS_DATE",	"M13_OSTEOPOROSIS_DATE",	"AUD_SWEDISH_DATE",	"E4_HYTHYNAS_DATE",	"E4_THYTOXGOITDIF_DATE",	"G6_SLEEPAPNO_DATE",	"IPF_DATE",	"ILD_DATE",	"GOUT_DATE",	"H7_GLAUCOMA_DATE",	"G6_EPLEPSY_DATE",	
          "GE_STRICT_DATE",	"FE_STRICT_DATE",	"K11_APPENDACUT_DATE",	
          "START_OF_FOLLOWUP",	"END_OF_FOLLOWUP","BATCH",	"SMOKING",	"EDUCATION_97",	"EDUCATION_11")	
names(df)[grep("^ID$",names(df))]<-"hospitalID" #need exact match
names(df)[grep("IID",names(df))]<-"ID" #matches genetics data (n=69715)
names(df)[grep("Sex",names(df))]<-"SEX" #2 and 1
names(df)[grep("Ancestry",names(df))]<-"ANCESTRY"
names(df)[grep("Bmi",names(df))]<-"BMI"
names(df)[grep("batch",names(df))]<-"BATCH"


#replace NA in the binary phenotype columns with 0 for controls
#can we assume that anyone who is genotyped but not in the EHR files is a control? or could they be missing from hospital records and should be NA?
binary<-c("C3_CANCER",	"C3_COLORECTAL",	
          "C3_BREAST",	"T2D",	"C3_PROSTATE",	"I9_CHD",	"I9_SAH",	"C3_MELANOMA_SKIN",	"J10_ASTHMA",	"I9_HEARTFAIL_NS",	"I9_STR",	"G6_AD_WIDE",	
          "T1D",	"I9_AF",	"N14_CHRONKIDNEYDIS",	"COVID",	"F5_DEPRESSIO",	"C3_BRONCHUS_LUNG",	"RHEUMA_SEROPOS_OTH",	"K11_IBD_STRICT",	"I9_VTE",	
          "I9_THAORTANEUR",	"I9_ABAORTANEUR",	"COX_ARTHROSIS",	"KNEE_ARTHROSIS",	"M13_OSTEOPOROSIS",	"AUD_SWEDISH",	"E4_HYTHYNAS",	"E4_THYTOXGOITDIF",	"G6_SLEEPAPNO",	"IPF",	
          "ILD",	"GOUT",	"H7_GLAUCOMA",	"G6_EPLEPSY",	"GE_STRICT",	"FE_STRICT",	"K11_APPENDACUT")
binary<-binary[binary %in% names(df)] #what actually exists
df2<-df %>% mutate_at(binary,~replace(.,is.na(.),0))

#make NA for males with breast cancer and females with prostate cancer, instead of controls
df2<-df2 %>% mutate(C3_BREAST=ifelse(SEX==1,NA,C3_BREAST))
df2<-df2 %>% mutate(C3_PROSTATE=ifelse(SEX==2,NA,C3_PROSTATE))

#make empty columns for missing phenotypes/variables 
columns<-header[!header %in% names(df2)]
empty<-data.frame(matrix(nrow=nrow(df2), ncol = length(columns))) 
names(empty)<-columns
df3<-cbind(df2,empty) %>% select(header) #subset just to headers of interest 

#write file
write.csv(df3,'endpointsPhenoFormatHUNT.csv',row.names=FALSE,quote=FALSE)

############## CALCULATE INTERQUARTILE RANGE FOR AGE AT ONSET FOR SELECT PHENOTYPES
p<-c("C3_BREAST","G6_EPLEPSY","GOUT","C3_PROSTATE","RHEUMA_SEROPOS_OTH","T1D","C3_CANCER","I9_AF","I9_CHD","T2D","I9_SAH","C3_MELANOMA_SKIN","J10_ASTHMA","F5_DEPRESSIO","C3_BRONCHUS_LUNG","COX_ARTHROSIS","KNEE_ARTHROSIS","K11_APPENDACUT","C3_COLORECTAL","ILD")

for (idx in 1:length(p)){
  print(p[idx])
  date_col=paste0(p[idx],"_DATE")
  if (p[idx] %in% colnames(df3)) {
    tmp<-df3 %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(as.POSIXct(get(date_col))-DATE_OF_BIRTH)/365.5)
   if (idx==1){
      age_df<-data.frame(p[idx],t(as.data.frame(quantile(tmp$age, probs = c(.25, .5, .75)))))
    } else{
      age_df<-rbind(age_df,data.frame(p[idx],t(as.data.frame(quantile(tmp$age, probs = c(.25, .5, .75))))))
    }
  }
}
names(age_df)<-c("trait","X25","X50","X75")
write.csv(format(age_df,digits=3),"age_quartiles.csv",row.names=FALSE,quote=FALSE)

##### other summary stats for phenotypes
p<-c("C3_CANCER","C3_COLORECTAL","C3_BREAST","T2D","C3_PROSTATE","I9_CHD","I9_SAH","C3_MELANOMA_SKIN","J10_ASTHMA","I9_HEARTFAIL_NS","I9_STR","G6_AD_WIDE","T1D","I9_AF","N14_CHRONKIDNEYDIS","COVID","F5_DEPRESSIO","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH",
     "K11_IBD_STRICT","I9_VTE","I9_THAORTANEUR","I9_ABAORTANEUR","COX_ARTHROSIS","KNEE_ARTHROSIS","M13_OSTEOPOROSIS","AUD_SWEDISH","E4_HYTHYNAS","E4_THYTOXGOITDIF","G6_SLEEPAPNO","IPF","ILD","GOUT","H7_GLAUCOMA","G6_EPLEPSY","GE_STRICT","FE_STRICT","K11_APPENDACUT")
for (idx in 1:length(p)){
  print(p[idx])
  if(sum(is.na(pull(df3,p[idx])))==nrow(df3)){ #if entire case/control designation is NA then the phenotype is missing 
    next
  } else {
    date_col=paste0(p[idx],"_DATE")
    # of cases and controls
    cases<-df3 %>% filter(get(p[idx])==1) %>% nrow()
    controls<-df3 %>% filter(get(p[idx])==0) %>% nrow()
    prev<-cases/(cases+controls)*100
    #Age distribution at time of recruitment/baseline (median, IQR)	
    tmp<-df3  %>% mutate(age=as.numeric((START_OF_FOLLOWUP-DATE_OF_BIRTH)/365.5))
    age_recruitment_median<-median(tmp$age,na.rm=TRUE)
    age_recruitment_IQR<-median(tmp$age,na.rm=TRUE)
    #Age of onset distribution (median, IQR)	ONLY CASES
    age_onset_median<-df3 %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(as.POSIXct(get(date_col))-DATE_OF_BIRTH)/365.5) %>% summarize(median(age))
    age_onset_IQR<-df3 %>% filter(get(p[idx])==1) %>% mutate(age=as.numeric(as.POSIXct(get(date_col))-DATE_OF_BIRTH)/365.5) %>% summarize(IQR(age))
    #Distribution of time of follow-up (median, IQR)	
    #why NA with followup dates?
    follow_up_median<-df3 %>% filter(!is.na(get(p[idx]))) %>% mutate(follow=as.numeric((END_OF_FOLLOWUP-START_OF_FOLLOWUP)/365.5)) %>% summarize(median(follow,na.rm=TRUE))
    follow_up_IQR<-df3 %>% filter(!is.na(get(p[idx]))) %>% mutate(follow=as.numeric((END_OF_FOLLOWUP-START_OF_FOLLOWUP)/365.5)) %>% summarize(IQR(follow,na.rm=TRUE))
    #correlations 
    age_corr<-cor.test(pull(df3,p[idx]),tmp$age)$estimate #has to be age of recruitment because age of onset wouldn't have controls
    age_cor_ci<-cor.test(pull(df3,p[idx]),tmp$age)$conf.int
    sex_corr<-cor.test(pull(df3,p[idx]),df3$SEX)$estimate
    sex_cor_ci<-cor.test(pull(df3,p[idx]),df3$SEX)$conf.int
    female_perc<-unlist(table(df3$SEX)/nrow(df3))[[2]]*100
  }
    if (idx==1){
      summary_stats_cases_df<-data.frame(p[idx],cases,controls,prev,age_recruitment_median,age_recruitment_IQR,age_onset_median,age_onset_IQR,follow_up_median,follow_up_IQR,age_corr,sex_corr,female_perc)
    } else{
      summary_stats_cases_df<-rbind(summary_stats_cases_df,data.frame(p[idx],cases,controls,prev,age_recruitment_median,age_recruitment_IQR,age_onset_median,age_onset_IQR,follow_up_median,follow_up_IQR,age_corr,sex_corr,female_perc))
}}

names(summary_stats_cases_df)<-c("trait","cases","controls","prevalence","age_recruitment_median","age_recruitment_IQR","age_onset_median","age_onset_IQR","follow_up_median","follow_up_IQR","age_corr","sex_corr","female_prev")
write.csv(format(summary_stats_cases_df,digits=3),"summary_stats_cases.csv",row.names=FALSE,quote=FALSE)

#################### files for the prspipeline
setwd("~/intervene")
df<-fread("endpointsWideFormatHUNT.csv")
config<-fread("/mnt/scratch/brooke/prspipe/prspipe/config/studies_for_methods_comparison.tsv")
phenos<-config$name
phenos<-unique(unlist(strsplit(phenos,",")))



"G6_ALZHEIMER"
"C3_BREAST"
"T2D"
"K11_IBD_STRICT"



