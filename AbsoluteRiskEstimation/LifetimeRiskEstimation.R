# Please download the accompanying global burden of disease statistics from the github 

# Points to edit - Full Sample:
#   * Lines 41 and 42 - Subset to diseases included in the analysis of your biobank
#   * Lines 49 and 50 - Enter file path and subset location to the country of interest
#   * Lines 53 and 54 - Enter file path and subset location to the country of interest
#   * Lines 77 and 78 - Enter file path and subset location to the country of interest
#   * Lines 81 and 82 - Enter file path and subset location to the country of interest
#   * Lines 97 and 98 - Enter file path and subset location to the country of interest
#   * Lines 100 and 101 - Enter file path and subset location to the country of interest
#   * Line 142 - Enter file path
#   * Line 215 - Enter file path and biobank name
#   * Line 233 - Enter file path and biobank name

# Points to edit - Sex Stratified Samples:
#   * Lines 244 and 245 - Subset to diseases included in the analysis of your biobank
#   * Lines 256 and 259 - Enter file path and subset location to the country of interest
#   * Line 284 and 286 - Enter file path and subset location to the country of interest
#   * Line 302 and 306 - Enter file path and subset location to the country of interest
#   * Line 335 - Enter file path
#   * Line 408 - Enter file path and biobank name
#   * Line 425 - Enter file path and biobank name

# Points to edit - Age Stratified Samples:
#   * Line 437 and 438 - Subset to diseases included in the analysis of your biobank
#   * Line 444 and 445 - Enter file path and subset location to the country of interest
#   * Line 448 and 449 - Enter file path and subset location to the country of interest
#   * Line 472 and 473 - Enter file path and subset location to the country of interest
#   * Line 476 and 477 - Enter file path and subset location to the country of interest
#   * Line 492 and 493 - Enter file path and subset location to the country of interest
#   * Line 495 and 496 - Enter file path and subset location to the country of interest
#   * Line 537 _ Enter file path
#   * Line 609 - Enter file path and biobank name
#   * Line 627 - Enter file path and biobank name

#Libraries
library(data.table)
library(dplyr)
library(ggplot2)

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  
  #Read in GBD incidence data 
  incidence <- fread("/enter/file/path/for/GBD_Incidence.csv", data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread("/enter/file/path/for/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="ENTER_COUNTRY")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  
  population <- c()
  for(i in unique(incidence$age)){
    subby <- subset(incidence, age==i)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0
  
  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population)
  incidence <- incidence[,c("location","age","cause","metric","population","incidence")]
  
  prevalence <- fread("/enter/file/path/for/GBD_Prevalence.csv", data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread("/enter/file/path/for/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="ENTER_COUNTRY")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread("/enter/file/path/for/GBD_Mortality.csv", data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  
  bcpc_mortality <- fread("/enter/file/path/for/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="ENTER_COUNTRY")
  
  mortality <- rbind(mortality, bcpc_mortality)
  
  mortality <- subset(mortality, cause==gbd_phenos[j] | cause=="All causes")
  
  if(gbd_phenos[j]=="Breast cancer"){
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Female")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  } else if(gbd_phenos[j]=="Prostate cancer"){
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Male")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  } else {
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Both")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  }
  
  cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
  cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val")]
  colnames(cause_specific_mortality)[5] <- c("cause_specific_rate")
  
  mortality <- left_join(cause_specific_mortality, all_cause_mortality)
  mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
  mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
  
  mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
  mortality <- mortality[,c("location","age","cause","mortality_rate")]
  
  #Merge mortality data to incidence data
  incidence <- left_join(incidence, mortality)
  
  incidence <- subset(incidence, age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
  
  #######################################################################################################################################################################################
  
  #Hazard Ratios
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread("/enter/file/path/for/HazardRatios_FullSample", data.table = FALSE)
  hazrats <- hazrats[,-1]
  colnames(hazrats) <- c("phenotype", "prs", "group", "beta", "se", "pval", "HR", "CIpos", "CIneg")
  hazrats <- subset(hazrats, phenotype==hr_phenos[j])
  
  #Hazard Ratios
  hr01 <- hazrats[1,"HR"]
  hr02 <- hazrats[2,"HR"]
  hr03 <- hazrats[3,"HR"]
  hr04 <- hazrats[4,"HR"]
  hr05 <- hazrats[5,"HR"]
  hr07 <- hazrats[6,"HR"]
  hr08 <- hazrats[7,"HR"]
  hr09 <- hazrats[8,"HR"]
  hr10 <- hazrats[9,"HR"]
  hr11 <- hazrats[10,"HR"]
  
  #Proportions - 0.2 by definition of PRS group
  props01 <- 0.01
  props02 <- 0.04
  props03 <- 0.05
  props04 <- 0.1
  props05 <- 0.2
  props06 <- 0.2
  props07 <- 0.2
  props08 <- 0.1
  props09 <- 0.05
  props10 <- 0.04
  props11 <- 0.01
  
  #Estimate incidence attributable to different distributions of PRS 
  incidence$i6 <- (incidence$incidence*incidence$population) / ((props06 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population)) + (hr11 * (props11 * incidence$population))) 
  incidence$i6[is.na(incidence$i6)] <- 0
  incidence$i1 <- incidence$i6 * hr01
  incidence$i2 <- incidence$i6 * hr02
  incidence$i3 <- incidence$i6 * hr03
  incidence$i4 <- incidence$i6 * hr04
  incidence$i5 <- incidence$i6 * hr05
  incidence$i7 <- incidence$i6 * hr07
  incidence$i8 <- incidence$i6 * hr08
  incidence$i9 <- incidence$i6 * hr09
  incidence$i10 <- incidence$i6 * hr10
  incidence$i11 <- incidence$i6 * hr11
  
  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  for(i in 1:11){
    #Calculate hazard
    incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
    
    #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
    incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
    
    #Mortality and risk
    incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
    
    #Survival
    incidence[[paste0("survival",i)]] <- 1
    
    for(k in 2:nrow(incidence)){
      incidence[[paste0("survival",i)]][k] <- exp(-5*incidence[[paste0("mortandrisk",i)]][k-1])
    }
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
  write.csv(lifetimerisk, paste0("/enter/file/path/",hr_phenos[j],"_LifetimeRisk_ENTER_BIOBANK_NAME.csv"))
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age Range") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group') +
    scale_color_hue(labels = c("0-1%", "1-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0("/enter/file/path/",hr_phenos[j],"_LifetimeRisk_ENTER_BIOBANK_NAME.png"), height=10 , width=10)
  
}

################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################

#Stratified by sex

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "T2D", "ILD", "C3_BRONCHUS_LUNG")

for(j in 1:length(gbd_phenos)){
  
  for(k in c("Male", "Female")){
    
    print(gbd_phenos[j])
    
    print(k)
    
    #Incidence data - basic pre-processing
    incidence <- fread("/enter/file/path/for/Sex_Stratified_Incidence_GBD.csv", data.table=FALSE)
    
    #Incidence data to be replaced with that for males and females for prostate cancer and breast cancer respectively. Came from a separate dataset to reduce size of the full dataset. 
    incidence <- subset(incidence, sex==k & cause==gbd_phenos[j] & location=="ENTER COUNTRY")
    incidence <- incidence[,c("location","age","cause","metric","val")]
    
    incidence$val <- as.numeric(incidence$val)
    
    population <- c()
    for(i in unique(incidence$age)){
      subby <- subset(incidence, age==i)
      poppy <- subby$val[1]/(subby$val[2]/100000)
      population <- c(population, poppy)
    }
    
    population[is.na(population)] <- 0
    
    ##Subset to Rate only as no longer require absolute numbers.
    incidence <- subset(incidence, metric=="Rate")
    incidence <- cbind(incidence, population)
    
    incidence$val <- as.numeric(incidence$val)
    
    #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
    incidence$incidence <- incidence$val / 100000
    incidence <- incidence[,c("location","age","cause","metric","incidence","population")]
    
    #Prevalence - use to calculate hazard (incidence/(1-prevalence)) - The code is equivalent to that defined for incidence. 
    prevalence <- fread("/enter/file/path/for/Sex_Stratified_Prevalence_GBD.csv", data.table=FALSE)
    
    prevalence <- subset(prevalence, sex==k & cause==gbd_phenos[j] & location=="ENTER_COUNTRY")
    prevalence <- prevalence[,c("location","age","cause","metric","val")]
    
    ##Subset to Rate only as no longer require absolute number
    prevalence <- subset(prevalence, metric=="Rate")
    
    prevalence$val <- as.numeric(prevalence$val)
    
    #Divide by 100,000 to get prevalence in terms of probability. 
    prevalence$prevalence <- prevalence$val / 100000
    prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
    
    #Left join to incidence to calculate hazard 
    incidence <- left_join(incidence, prevalence)
    
    #Calculate age specific and disease specific mortality
    mortality <- fread("/enter/file/path/for/Sex_Stratified_Mortality_GBD", data.table=FALSE)
    
    mortality <- mortality[,c("location","sex","age","cause","metric","val")]
    
    mortality <- subset(mortality, sex==k & (cause==gbd_phenos[j] | cause=="All causes") & location=="ENTER_COUNTRY")
    
    #Separate the current dataset into all cause and cause specific mortality so that the table can be converted into a wide format. 
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate")
    all_cause_mortality <- all_cause_mortality[,c("location","age","val")]
    colnames(all_cause_mortality)[3] <- c("all_cause_rate")
    
    cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
    cause_specific_mortality <- cause_specific_mortality[,c("location","age","cause","val")]
    colnames(cause_specific_mortality)[4] <- c("cause_specific_rate")
    
    mortality <- left_join(cause_specific_mortality, all_cause_mortality)
    mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
    mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
    
    #Subtract cause specific mortality rate from all cause mortality rate and divide by 100,000 to put into probability. 
    mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
    mortality <- mortality[,c("location","age","cause","mortality_rate")]
    
    #Merge mortality data to incidence data
    incidence <- left_join(incidence, mortality)
    
    incidence <- subset(incidence, age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
    
    #######################################################################################################################################################################################
    
    #Hazard Ratios 
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(paste0("/enter/file/path/for/hazardratios"), data.table = FALSE)
    hazrats <- hazrats[,-1]
    colnames(hazrats) <- c("phenotype", "prs", "group", "beta", "se", "pval", "HR", "CIpos", "CIneg")
    hazrats <- subset(hazrats, phenotype==hr_phenos[j])
    
    #Hazard Ratios
    hr01 <- hazrats[1,"HR"]
    hr02 <- hazrats[2,"HR"]
    hr03 <- hazrats[3,"HR"]
    hr04 <- hazrats[4,"HR"]
    hr05 <- hazrats[5,"HR"]
    hr07 <- hazrats[6,"HR"]
    hr08 <- hazrats[7,"HR"]
    hr09 <- hazrats[8,"HR"]
    hr10 <- hazrats[9,"HR"]
    hr11 <- hazrats[10,"HR"]
    
    #Proportions - 0.2 by definition of PRS group
    props01 <- 0.01
    props02 <- 0.04
    props03 <- 0.05
    props04 <- 0.1
    props05 <- 0.2
    props06 <- 0.2
    props07 <- 0.2
    props08 <- 0.1
    props09 <- 0.05
    props10 <- 0.04
    props11 <- 0.01
    
    #Estimate incidence attributable to different distributions of PRS 
    incidence$i6 <- (incidence$incidence*incidence$population) / ((props06 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population)) + (hr11 * (props11 * incidence$population))) 
    incidence$i6[is.na(incidence$i6)] <- 0
    incidence$i1 <- incidence$i6 * hr01
    incidence$i2 <- incidence$i6 * hr02
    incidence$i3 <- incidence$i6 * hr03
    incidence$i4 <- incidence$i6 * hr04
    incidence$i5 <- incidence$i6 * hr05
    incidence$i7 <- incidence$i6 * hr07
    incidence$i8 <- incidence$i6 * hr08
    incidence$i9 <- incidence$i6 * hr09
    incidence$i10 <- incidence$i6 * hr10
    incidence$i11 <- incidence$i6 * hr11
    
    ###################################################
    
    lifetimerisk <- data.frame(NULL)
    for(i in 1:11){
      #Calculate hazard
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
      
      #Survival
      incidence[[paste0("survival",i)]] <- 1
    
      for(l in 2:nrow(incidence)){
        incidence[[paste0("survival",i)]][k] <- exp(-5*incidence[[paste0("mortandrisk",i)]][k-1])
      }
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
    
    #Plot all as well as overall lifetime risk
    lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
    
    write.csv(lifetimerisk, paste0("/enter/file/path/",hr_phenos[j],"_",k,"_LifetimeRisk_ENTER_BIOBANK_NAME.csv"))
    
    ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
      geom_point() +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group') +
      scale_color_hue(labels = c("0-1%", "1-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/enter/file/path/",hr_phenos[j],"_",k,"_LifetimeRisk_ENTER_BIOBANK_NAME.png"), height=10 , width=10)
    
  }
}

########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################

#Age stratification results

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Idiopathic epilepsy", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "G6_EPLEPSY", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "T2D", "ILD", "C3_BRONCHUS_LUNG")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  #Read in GBD incidence data 
  incidence <- fread("/enter/file/path/for/GBD_Incidence.csv", data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread("/enter/file/path/for/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="ENTER_COUNTRY")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  
  population <- c()
  for(i in unique(incidence$age)){
    subby <- subset(incidence, age==i)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0
  
  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population)
  incidence <- incidence[,c("location","age","cause","metric","population","incidence")]
  
  prevalence <- fread("/enter/file/path/for/GBD_Prevalence.csv", data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread("/enter/file/path/for/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="ENTER_COUNTRY")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread("/enter/file/path/for/GBD_Mortality.csv", data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="ENTER_COUNTRY")
  
  bcpc_mortality <- fread("/enter/file/path/for/GBD_Data/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="ENTER_COUNTRY")
  
  mortality <- rbind(mortality, bcpc_mortality)
  
  mortality <- subset(mortality, cause==gbd_phenos[j] | cause=="All causes")
  
  if(gbd_phenos[j]=="Breast cancer"){
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Female")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  } else if(gbd_phenos[j]=="Prostate cancer"){
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Male")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  } else {
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Both")
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
  }
  
  cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
  cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val")]
  colnames(cause_specific_mortality)[5] <- c("cause_specific_rate")
  
  mortality <- left_join(cause_specific_mortality, all_cause_mortality)
  mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
  mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
  
  mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
  mortality <- mortality[,c("location","age","cause","mortality_rate")]
  
  #Merge mortality data to incidence data
  incidence <- left_join(incidence, mortality)
  
  incidence <- subset(incidence, age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
  
  #######################################################################################################################################################################################
  
  #Hazard Ratios 
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread(paste0("/enter/file/path/",hr_phenos[j],"_Age_Specific_Hazards.csv"), data.table = FALSE)
  hazrats <- hazrats[,-1]
  colnames(hazrats) <- c("Age","Beta","HR","Group")
  
  #Hazard Ratios
  hr01 <- hazrats[hazrats$Group=="Group 1","HR"]
  hr02 <- hazrats[hazrats$Group=="Group 2","HR"]
  hr03 <- hazrats[hazrats$Group=="Group 3","HR"]
  hr04 <- hazrats[hazrats$Group=="Group 4","HR"]
  hr05 <- hazrats[hazrats$Group=="Group 5","HR"]
  hr07 <- hazrats[hazrats$Group=="Group 7","HR"]
  hr08 <- hazrats[hazrats$Group=="Group 8","HR"]
  hr09 <- hazrats[hazrats$Group=="Group 9","HR"]
  hr10 <- hazrats[hazrats$Group=="Group 10","HR"]
  hr11 <- hazrats[hazrats$Group=="Group 11","HR"]
  
  #Proportions of each PRS group.
  props01 <- 0.01
  props02 <- 0.04
  props03 <- 0.05
  props04 <- 0.1
  props05 <- 0.2
  props06 <- 0.2
  props07 <- 0.2
  props08 <- 0.1
  props09 <- 0.05
  props10 <- 0.04
  props11 <- 0.01
  
  #Estimate incidence attributable to different distributions of PRS 
  incidence$i6 <- (incidence$incidence*incidence$population) / ((props06 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population)) + (hr11 * (props11 * incidence$population))) 
  incidence$i6[is.na(incidence$i6)] <- 0
  incidence$i1 <- incidence$i6 * hr01
  incidence$i2 <- incidence$i6 * hr02
  incidence$i3 <- incidence$i6 * hr03
  incidence$i4 <- incidence$i6 * hr04
  incidence$i5 <- incidence$i6 * hr05
  incidence$i7 <- incidence$i6 * hr07
  incidence$i8 <- incidence$i6 * hr08
  incidence$i9 <- incidence$i6 * hr09
  incidence$i10 <- incidence$i6 * hr10
  incidence$i11 <- incidence$i6 * hr11
  
  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  for(i in 1:11){
    #Calculate hazard
    incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
    
    #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
    incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
    
    #Mortality and risk
    incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
    
    #Survival
    incidence[[paste0("survival",i)]] <- 1
    
    for(k in 2:nrow(incidence)){
      incidence[[paste0("survival",i)]][k] <- exp(-5*incidence[[paste0("mortandrisk",i)]][k-1])
    }
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
  write.csv(lifetimerisk, paste0("/enter/file/path/",hr_phenos[j],"_LifetimeRisk_AgeStratification_ENTER_BIOBANK_NAME.csv"))
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age Range") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group') +
    scale_color_hue(labels = c("0-1%", "1-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0("/enter/file/path/",hr_phenos[j],"_LifetimeRisk_AgeStratification_ENTER_BIOBANK_NAME.png"), height=10 , width=10)
  
}
