
/apps/statistics2/R-4.0.0/bin/R

library(data.table)
library(tidyverse)
library(doParallel)

print("Starting!!!!")

# Start the clock!
ptm <- proc.time()

endpoints <- fread("UKBB_definitions_demo_TEST.csv")
head(endpoints)

#Contains information on baseline age and age at end of follow up. 
#bl <- fread("zcat ../tuomo/UKBB_COV.txt.gz")
#head(bl)

#Transforms input into regex format for detecting ICD_10, ICD_9 and ATC codes.
endpoints$ICD_10 <- ifelse(!is.na(endpoints$COD_ICD_10), paste0("^10x(",endpoints$COD_ICD_10,")"), NA)
endpoints$ICD_10_EXC <- ifelse(!is.na(endpoints$COD_ICD_10_EXC), paste0("^10x(",endpoints$COD_ICD_10_EXC,")"), NA)
endpoints$ICD_9 <- ifelse(!is.na(endpoints$COD_ICD_9), paste0("^9x(",endpoints$COD_ICD_9,")"), NA)  
endpoints$ICD_9_EXC <- ifelse(!is.na(endpoints$COD_ICD_9_EXC), paste0("^9x(",endpoints$COD_ICD_9_EXC,")"), NA) 
endpoints$ATC <- ifelse(!is.na(endpoints$KELA_ATC), paste0("^ATCx(",endpoints$KELA_ATC,")"),NA) 
endpoints$INCLUDE <- ifelse(!is.na(endpoints$INCLUDE), paste0("(", endpoints$INCLUDE, ")"), NA)
endpoints$CONDITIONS <- ifelse(!is.na(endpoints$CONDITIONS), paste0("(", endpoints$CONDITIONS, ")"), NA)

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
dat <- fread("zcat ../tuomo/UKBB_longitudinal.txt")
head(dat)

# You must indicate .packages parameter.

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
  exclusions <- exclusions[!duplicated(exclusions$eid)]
  cases <- subset(cases, !(eid %in% exclusions$eid))
  cases$ENDPOINT_LEVEL_1 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  dataset <- rbind(dataset, cases)
}

endpoints_level2 <- subset(endpoints, LEVEL==2) 
for(i in 1:nrow(endpoints_level2)){
  endpoint <- endpoints_level2[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(endpoint$INCLUDE, ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME)
  cases <- cases %>% 
              group_by(eid) %>%
                arrange(DATE) %>%
                  slice(1L)
  cases$ENDPOINT_LEVEL_2 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  
  exclusions <- subset(dataset, grepl(endpoint$EXCLUDE, ENDPOINT_LEVEL_1))
  exclusions <- exclusions[!duplicated(exclusions$eid)]
  cases <- subset(cases, !(eid %in% exclusions$eid)) 
  
  dataset <- left_join(dataset, cases)
}

endpoints_level3 <- subset(endpoints, LEVEL==3) 
for(i in 1:nrow(endpoints_level3)){
  endpoint <- endpoints_level3[i,]
  print(endpoint$NAME)
  cases <- subset(dataset, grepl(endpoint$INCLUDE, ENDPOINT_LEVEL_1) | ENDPOINT_LEVEL_1 == endpoint$NAME | grepl(endpoint$INCLUDE, ENDPOINT_LEVEL_2) | ENDPOINT_LEVEL_2 == endpoint$NAME)
  cases <- cases %>% 
              group_by(eid) %>%
                  arrange(DATE) %>%
                      slice(1L)
  cases$ENDPOINT_LEVEL_3 <- endpoint$NAME
  cases$LEVEL <- endpoint$LEVEL
  
  exclusions <- subset(dataset, grepl(endpoint$EXCLUDE, ENDPOINT_LEVEL_1) | grepl(endpoint$EXCLUDE, ENDPOINT_LEVEL_1))
  exclusions <- exclusions[!duplicated(exclusions$eid)]
  cases <- subset(cases, !(eid %in% exclusions$eid)) 
  
  dataset <- left_join(dataset, cases)
}

#UP TO HERE!!

#Sex specific phenotypes - shouldn't remove any as they have the code. More of a santy check.  
dataset[sample(1:nrow(dataset), size = 20)]

dataset$value <- 1
dataset <- dataset[,c("eid", "ENDPOINT", "value")]
ukb_wide <- spread(dataset,ENDPOINT,value)

bleid <- bl[,1]
ukb_wide <- left_join(bleid, ukb_wide, by ="eid")
ukb_wide[is.na(ukb_wide)]<-0
ukb_wide <- left_join(bl,ukb_wide, by = "eid")

head(ukb_wide)

sapply(ukb_wide[,(ncol(bl)+1):ncol(ukb_wide)],sum)

# Stop the clock
proc.time() - ptm

print("ready")