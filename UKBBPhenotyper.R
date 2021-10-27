
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
  case <- case[,c('eid',paste0(i))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

level2 <- subset(endpoints, LEVEL==2 & !is.na(ENDPOINT_LEVEL_2))
level2 <- select(level2, c("eid","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_2"))
for(i in unique(level2$ENDPOINT_LEVEL_2)){
  print(i)
  case <- subset(level2, i==ENDPOINT_LEVEL_2)
  case[[paste0(i)]] <- 1
  case <- case[,c('eid',paste0(i))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

level3 <- subset(endpoints, LEVEL==3 & !is.na(ENDPOINT_LEVEL_3))
level3 <- select(level3, c("eid","DATE","CODE1", "LEVEL", "ENDPOINT_LEVEL_3"))
for(i in unique(level3$ENDPOINT_LEVEL_3)){
  print(i)
  case <- subset(level3, i==ENDPOINT_LEVEL_3)
  case[[paste0(i)]] <- 1
  case <- case[,c('eid',paste0(i))]
  uniqueIDs <- left_join(uniqueIDs, case)
}

write.csv(uniqueIDs, 'endpointsWideFormat.csv')
