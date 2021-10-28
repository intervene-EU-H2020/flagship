
/apps/statistics2/R-4.0.0/bin/R

library(data.table)
library(tidyverse)

#Read in the information for each phenotype
endpoints <- fread("UKBB_definitions_demo_TEST.csv", data.table=FALSE)
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
dat <- fread("zcat ../tuomo/UKBB_longitudinal.txt", data.table=FALSE)
head(dat)

#Phenotypes have been structured into levels to account for the inherent hierarchy. Level 1 phenotypes can be defined using ICD codes. Level 2 phenotypes require level 1 phenotypes to be defined. Level 3 phenotypes require level 2 phenotypes to be defined. 

#For loop which iterates over each endpoint and extracts all instances where the code can be found in the longitudinal dataset. Extracts the earliest date for each participant. 
dataset <- data.frame(NULL)
for(i in c(1:nrow(endpoints))){
  endpoint<-endpoints[i,]
  print(endpoint$NAME)
  cases <- subset(dat, grepl(endpoint$CODE, CODE1) & endpoint$CODE!="")
  cases <- cases %>% 
              group_by(eid) %>%
                arrange(DATE) %>%
                  slice(1L)
  exclusions <- subset(dat, grepl(endpoint$CODE_EXC, CODE1) & endpoint$CODE_EXC!="")
  exclusions <- exclusions[!duplicated(exclusions$eid),]
  cases <- subset(cases, !(eid %in% exclusions$eid))
  cases$ENDPOINT_LEVEL_1 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  dataset <- rbind(dataset, cases)
}

endpoints_level2 <- subset(endpoints, LEVEL==2) 
dataset$ENDPOINT_LEVEL_2 <- NA 

for(i in c(1:nrow(endpoints_level2))){
  endpoint <- endpoints_level2[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME)
  cases <- cases %>% 
              group_by(eid) %>%
                arrange(DATE) %>%
                  slice(1L)
  cases$ENDPOINT_LEVEL_2 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  
  exclusions <- subset(dataset, grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1))
  exclusions <- exclusions[!duplicated(exclusions$eid),]
  cases <- subset(cases, !(eid %in% exclusions$eid)) 
  dataset <- rbind(dataset, cases)
}

endpoints_level3 <- subset(endpoints, LEVEL==3) 
dataset$ENDPOINT_LEVEL_3 <- NA 

for(i in 1:nrow(endpoints_level3)){
  endpoint <- endpoints_level3[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME | grepl(paste0("\\b",endpoint$INCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_2) | ENDPOINT_LEVEL_2 == endpoint$NAME)
  cases <- cases %>% 
              group_by(eid) %>%
                  arrange(DATE) %>%
                      slice(1L)
  cases$ENDPOINT_LEVEL_3 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  exclusions <- subset(dataset, grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_1) | grepl(paste0("\\b",endpoint$EXCLUDE,"\\b", sep=""), ENDPOINT_LEVEL_2)) 
  exclusions <- exclusions[!duplicated(exclusions$eid),]
  cases <- subset(cases, !(eid %in% exclusions$eid)) 
  
  dataset <- rbind(dataset, cases)
}

#Save dataset
write.csv(dataset, 'endpointsDefinitionsUKB.csv')

##########################################################################################################################################################################################

#Aim for this section of code is to produce a wide table which will create an array of all phenotypes that belong to a person. No information is provided on earliest date for this table. 

endpoints <- fread('endpointDefinitionsUKB.csv', data.table=FALSE)
endpoints <- endpoints[,-1]

uniqueIDs <- as.data.frame(endpoints$eid[!duplicated(endpoints$eid)])
colnames(uniqueIDs) <- 'eid'

level1 <- subset(endpoints, LEVEL==1)
level1 <- select(level1, c("eid","DATE","CODE1","ENDPOINT_LEVEL_1", "LEVEL"))
for(i in unique(level1$ENDPOINT_LEVEL_1)){
  print(i)
  case <- subset(level1, i==ENDPOINT_LEVEL_1)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('eid',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

level2 <- subset(endpoints, LEVEL==2 & !is.na(ENDPOINT_LEVEL_2))
level2 <- select(level2, c("eid","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_2"))
for(i in unique(level2$ENDPOINT_LEVEL_2)){
  print(i)
  case <- subset(level2, i==ENDPOINT_LEVEL_2)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('eid',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

level3 <- subset(endpoints, LEVEL==3 & !is.na(ENDPOINT_LEVEL_3))
level3 <- select(level3, c("eid","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_3"))
for(i in unique(level3$ENDPOINT_LEVEL_3)){
  print(i)
  case <- subset(level3, i==ENDPOINT_LEVEL_3)
  case[[paste0(i)]] <- 1
  case[[paste0(i,"_DATE",sep="")]] <- case$DATE
  case <- case[,c('eid',paste0(i),paste0(i,"_DATE",sep=""))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

#Select endpoints of interest for flagship project
endpointswide <- select(uniqueIDs, c("eid","C3_CANCER", "C3_COLORECTAL", "C3_BREAST", "T2D", "C3_PROSTATE", "I9_CHD", "I9_SAH", "C3_MELANOMA_SKIN", 
                                     "J10_ASTHMA", "I9_HEARTFAIL_NS", "I9_STR", "G6_AD_WIDE", "T1D", "I9_AF", "N14_CHRONKIDNEYDIS", "F5_DEPRESSIO", 
                                     "C3_BRONCHUS_LUNG", "RHEUMA_SEROPOS_OTH", "K11_IBD_STRICT", "I9_VTE", "I9_THAORTANEUR", "I9_ABAORTANEUR", 
                                     "COX_ARTHROSIS", "KNEE_ARTHROSIS", "M13_OSTEOPOROSIS", "AUD_SWEDISH", "E4_HYTHYNAS", "G6_SLEEPAPNO", "IPF", 
                                     "ILD", "GOUT", "H7_GLAUCOMA", "G6_EPLEPSY", "GE_STRICT", "FE_STRICT", "K11_APPENDACUT",
                                     "C3_CANCER_DATE", "C3_COLORECTAL_DATE", "C3_BREAST_DATE", "T2D_DATE", "C3_PROSTATE_DATE", "I9_CHD_DATE", "I9_SAH_DATE", "C3_MELANOMA_SKIN_DATE", "J10_ASTHMA_DATE", 
                                     "I9_HEARTFAIL_NS_DATE", "I9_STR_DATE", "G6_AD_WIDE_DATE", "T1D_DATE", "I9_AF_DATE", "N14_CHRONKIDNEYDIS_DATE", "F5_DEPRESSIO_DATE", "C3_BRONCHUS_LUNG_DATE", 
                                     "RHEUMA_SEROPOS_OTH_DATE", "K11_IBD_STRICT_DATE", "I9_VTE_DATE", "I9_THAORTANEUR_DATE", "I9_ABAORTANEUR_DATE", "COX_ARTHROSIS_DATE", "KNEE_ARTHROSIS_DATE", 
                                     "M13_OSTEOPOROSIS_DATE", "AUD_SWEDISH_DATE", "E4_HYTHYNAS_DATE", "G6_SLEEPAPNO_DATE", "IPF_DATE", "ILD_DATE", "GOUT_DATE", "H7_GLAUCOMA_DATE", "G6_EPLEPSY_DATE", 
                                     "GE_STRICT_DATE", "FE_STRICT_DATE", "K11_APPENDACUT_DATE"))

write.csv(uniqueIDs, 'endpointsWideFormat.csv')

#There will be a group who have received no ICD codes. These participants should be included as controls also. Link these remaining people to the saved dataset above and include as controls. 
