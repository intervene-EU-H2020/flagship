#Libraries
library(data.table)
library(tidyverse)
library(cmprsk)
library(lubridate)

#Read in Phenotype file
pheno <- fread("path/to/phenotype/file", data.table=FALSE)

#Death registry 
death <- fread(input = "/path/to/death_registry", data.table=FALSE)
death <- subset(death, EVENT_TYPE=='DEATH')
death <- unique(death$ID)

phenos <- c("C3_CANCER", "C3_COLORECTAL", "C3_BREAST", "T2D", "C3_PROSTATE",  "I9_CHD", "I9_SAH", "C3_MELANOMA_SKIN", "J10_ASTHMA", "T1D", "I9_AF", "F5_DEPRESSIO", "C3_BRONCHUS_LUNG", "RHEUMA_SEROPOS_OTH", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "ILD", "GOUT", "G6_EPLEPSY", "K11_APPENDACUT")

for(i in phenos){
  
pheno$death <- ifelse((pheno$ID %in% death & pheno[[i]]==0), 1, 0)

pheno$fstatus <- case_when(pheno[[i]] == 1 ~ 1,
                           pheno$death == 1 ~ 2,
                           TRUE ~ 0)

pheno$ftime<- ifelse(pheno[[i]] == 1, time_length(difftime(pheno[,paste0(i,"_DATE")], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP_ADJ, pheno$DATE_OF_BIRTH), 'years'))

cumulative_inc <- cuminc(ftime=pheno$ftime, fstatus=pheno$fstatus, cencode=0, na.action=na.omit)

print(i)
print(cumulative_inc)

}
