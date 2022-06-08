brewer.pal(n = 3, name = 'Dark2')

#Libraries
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4) 


read_sheet("https://docs.google.com/spreadsheets/d/13ZvAgt9oYVpPN_FJMoQKLHsH23HyDmyPo-jfPl2_6HY/edit#gid=896420200")

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")

### custom path to GBD data
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/"

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  #Read in GBD incidence data 
  incidence <- fread(paste0(path,"GBD_Incidence.csv"), data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  incidence$val <- as.numeric(incidence$val)
  incidence$upper <- as.numeric(incidence$upper)
  incidence$lower <- as.numeric(incidence$lower)
  
  bcpc_incidence <- fread(paste0(path,"BreastCancerProstateCancer_Incidence.csv"), data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  bcpc_incidence$upper <- as.numeric(bcpc_incidence$upper)
  bcpc_incidence$lower <- as.numeric(bcpc_incidence$lower)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  incidence$incidence_upper <- incidence$upper / 100000
  incidence$incidence_lower <- incidence$lower / 100000
  
  population <- c()
  
  for(i in unique(incidence$age)){
    subby <- subset(incidence, age==i)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0

  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population) 

  incidence <- incidence[,c("location","age","cause","metric","population","incidence","incidence_upper","incidence_lower")]
  
  prevalence <- fread(paste0(path,"GBD_Prevalence.csv"), data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  prevalence$val <- as.numeric(prevalence$val)
  prevalence$upper <- as.numeric(prevalence$upper)
  prevalence$lower <- as.numeric(prevalence$lower)
  
  bcpc_prevalence <- fread(paste0(path,"BreastCancerProstateCancer_Prevalence.csv"), data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  bcpc_prevalence$upper <- as.numeric(bcpc_prevalence$upper)
  bcpc_prevalence$lower <- as.numeric(bcpc_prevalence$lower)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence$prevalence_upper <- prevalence$upper / 100000
  prevalence$prevalence_lower <- prevalence$lower / 100000
  
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence","prevalence_upper","prevalence_lower")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread(paste0(path,"GBD_Mortality.csv"), data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  
  bcpc_mortality <- fread(paste0(path,"BreastCancerProstateCancer_Mortality.csv"), data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="Norway")
  
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
  hazrats <- fread("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/HazardRatios_FullSample_FinnGen.csv", data.table = FALSE)
  colnames(hazrats) <- c("phenotype", "prs", "group", "controls", "cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
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
  
  #Upper Confidence Intervals
  CIpos_hr01 <- hazrats[1,"CIpos"]
  CIpos_hr02 <- hazrats[2,"CIpos"]
  CIpos_hr03 <- hazrats[3,"CIpos"]
  CIpos_hr04 <- hazrats[4,"CIpos"]
  CIpos_hr05 <- hazrats[5,"CIpos"]
  CIpos_hr07 <- hazrats[6,"CIpos"]
  CIpos_hr08 <- hazrats[7,"CIpos"]
  CIpos_hr09 <- hazrats[8,"CIpos"]
  CIpos_hr10 <- hazrats[9,"CIpos"]
  CIpos_hr11 <- hazrats[10,"CIpos"]
  
  #Lower Confidence Intervals
  CIneg_hr01 <- hazrats[1,"CIneg"]
  CIneg_hr02 <- hazrats[2,"CIneg"]
  CIneg_hr03 <- hazrats[3,"CIneg"]
  CIneg_hr04 <- hazrats[4,"CIneg"]
  CIneg_hr05 <- hazrats[5,"CIneg"]
  CIneg_hr07 <- hazrats[6,"CIneg"]
  CIneg_hr08 <- hazrats[7,"CIneg"]
  CIneg_hr09 <- hazrats[8,"CIneg"]
  CIneg_hr10 <- hazrats[9,"CIneg"]
  CIneg_hr11 <- hazrats[10,"CIneg"]
  
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

  #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
  incidence$CIpos_i6 <- (incidence$incidence_upper*incidence$population) / ((props06 * incidence$population) + (CIpos_hr01 * (props01 * incidence$population)) + (CIpos_hr02 * (props02 * incidence$population)) + (CIpos_hr03 * (props03 * incidence$population)) + (CIpos_hr04 * (props04 * incidence$population)) + (CIpos_hr05 * (props05 * incidence$population)) + (CIpos_hr07 * (props07 * incidence$population)) + (CIpos_hr08 * (props08 * incidence$population)) + (CIpos_hr09 * (props09 * incidence$population)) + (CIpos_hr10 * (props10 * incidence$population)) + (CIpos_hr11 * (props11 * incidence$population))) 
  incidence$CIpos_i6[is.na(incidence$CIpos_i6)] <- 0
  incidence$CIpos_i1 <- incidence$CIpos_i6 * CIpos_hr01
  incidence$CIpos_i2 <- incidence$CIpos_i6 * CIpos_hr02
  incidence$CIpos_i3 <- incidence$CIpos_i6 * CIpos_hr03
  incidence$CIpos_i4 <- incidence$CIpos_i6 * CIpos_hr04
  incidence$CIpos_i5 <- incidence$CIpos_i6 * CIpos_hr05
  incidence$CIpos_i7 <- incidence$CIpos_i6 * CIpos_hr07
  incidence$CIpos_i8 <- incidence$CIpos_i6 * CIpos_hr08
  incidence$CIpos_i9 <- incidence$CIpos_i6 * CIpos_hr09
  incidence$CIpos_i10 <- incidence$CIpos_i6 * CIpos_hr10
  incidence$CIpos_i11 <- incidence$CIpos_i6 * CIpos_hr11
  
  #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
  incidence$CIneg_i6 <- (incidence$incidence_lower*incidence$population) / ((props06 * incidence$population) + (CIneg_hr01 * (props01 * incidence$population)) + (CIneg_hr02 * (props02 * incidence$population)) + (CIneg_hr03 * (props03 * incidence$population)) + (CIneg_hr04 * (props04 * incidence$population)) + (CIneg_hr05 * (props05 * incidence$population)) + (CIneg_hr07 * (props07 * incidence$population)) + (CIneg_hr08 * (props08 * incidence$population)) + (CIneg_hr09 * (props09 * incidence$population)) + (CIneg_hr10 * (props10 * incidence$population)) + (CIneg_hr11 * (props11 * incidence$population))) 
  incidence$CIneg_i6[is.na(incidence$CIneg_i6)] <- 0
  incidence$CIneg_i1 <- incidence$CIneg_i6 * CIneg_hr01
  incidence$CIneg_i2 <- incidence$CIneg_i6 * CIneg_hr02
  incidence$CIneg_i3 <- incidence$CIneg_i6 * CIneg_hr03
  incidence$CIneg_i4 <- incidence$CIneg_i6 * CIneg_hr04
  incidence$CIneg_i5 <- incidence$CIneg_i6 * CIneg_hr05
  incidence$CIneg_i7 <- incidence$CIneg_i6 * CIneg_hr07
  incidence$CIneg_i8 <- incidence$CIneg_i6 * CIneg_hr08
  incidence$CIneg_i9 <- incidence$CIneg_i6 * CIneg_hr09
  incidence$CIneg_i10 <- incidence$CIneg_i6 * CIneg_hr10
  incidence$CIneg_i11 <- incidence$CIneg_i6 * CIneg_hr11
  
  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  for(i in 1:11){
    #Calculate hazard
    incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
    
    #CIpos hazard
    incidence[[paste0("CIpos_hazard",i)]] <- incidence[[paste0("CIpos_i",i)]] / (1 - incidence$prevalence_upper)
    
    #CIneg hazard
    incidence[[paste0("CIneg_hazard",i)]] <- incidence[[paste0("CIneg_i",i)]] / (1 - incidence$prevalence_lower)
    
    #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
    incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
    
    #CIpos risk
    incidence[[paste0("CIpos_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIpos_hazard",i)]])
    
    #CIneg risk
    incidence[[paste0("CIneg_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIneg_hazard",i)]])
    
    #Calculate the probability of remaining disease free during the age interval - death as a competing interest not considered in this section. 
    #incidence[[paste0("cumulativerisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]])
    #incidence[[paste0("survival_noncomp",i)]] <- exp(-5*incidence[[paste0("cumulativerisk",i)]])
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    #incidence[[paste0("lifetimerisk_noncomp",i)]] <- cumsum(incidence[[paste0("survival_noncomp",i)]]*incidence[[paste0("risk",i)]])*100
    
    #Mortality and risk
    incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
    
    #CIpos mortandrisk
    incidence[[paste0("CIpos_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIpos_hazard",i)]] + incidence$mortality_rate)
    
    #CIneg mortandrisk
    incidence[[paste0("CIneg_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIneg_hazard",i)]] + incidence$mortality_rate)
    
    #Survival
    incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
    
    #CIpos survival
    incidence[[paste0("CIpos_survival",i)]] <- exp(-5*incidence[[paste0("CIpos_mortandrisk",i)]])
    
    #CIneg survival
    incidence[[paste0("CIneg_survival",i)]] <- exp(-5*incidence[[paste0("CIneg_mortandrisk",i)]])
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    #CIpos lifetime risk
    incidence[[paste0("CIpos_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIpos_survival",i)]]*incidence[[paste0("CIpos_risk",i)]])*100
    
    #CIneg lifetime risk
    incidence[[paste0("CIneg_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIneg_survival",i)]]*incidence[[paste0("CIneg_risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]], incidence[[paste0("CIpos_lifetimerisk",i)]], incidence[[paste0("CIneg_lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk", "LifetimeRisk_CIpos", "LifetimeRisk_CIneg")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group11","Group10","Group9","Group8","Group7","Group6","Group5","Group4","Group3","Group2","Group1"))
  
  write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_FinnGen.csv"))
  
  color <- brewer.pal(n = 12, name = 'Paired')
  color <- color[-11]
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age") + 
    ylab("Lifetime Risk (%)") + 
    theme_bw() +
    labs(color='PRS Strata') +
    scale_fill_manual(values=color, labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) + 
    scale_color_manual(values=color, labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) + 
    theme(title = element_text(size = 30),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 26),
          axis.title.x = element_text(size = 26),
          axis.text.x = element_text(size = 24, angle=-90, hjust=0),
          axis.title.y = element_text(size = 26),
          axis.text.y = element_text(size = 24))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_FinnGen.png"), height=10 , width=10, dpi=1200)
  
}

################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################

#Stratified by sex

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Diabetes mellitus type 2")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "T2D")

for(j in 1:length(gbd_phenos)){
  
  for(k in c("Male", "Female")){
    
    print(gbd_phenos[j])
    
    print(k)
    
    #Incidence data - basic pre-processing
    incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBD.csv", data.table=FALSE)
    
    #Incidence data to be replaced with that for males and females for prostate cancer and breast cancer respectively. Came from a separate dataset to reduce size of the full dataset. 
    incidence <- subset(incidence, sex==k & cause==gbd_phenos[j] & location=="Norway")
    incidence <- incidence[,c("location","age","cause","metric","val")]
    
    incidence$val <- as.numeric(incidence$val)
    
    population <- c()
    for(i in unique(incidence$age)){
      print(i)
      subby <- subset(incidence, age==i)
      print(subby)
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
    prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBD.csv", data.table=FALSE)
    
    prevalence <- subset(prevalence, sex==k & cause==gbd_phenos[j] & location=="Norway")
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
    mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBD.csv", data.table=FALSE)
    
    mortality <- mortality[,c("location","sex","age","cause","metric","val")]
    
    mortality <- subset(mortality, sex==k & (cause==gbd_phenos[j] | cause=="All causes") & location=="Norway")
    
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
    
    #Hazard Ratios - maddavat
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/HazardRatios_",k,"Sample_FinnGen.csv"), data.table = FALSE)
    colnames(hazrats) <- c("phenotype", "prs", "group", "controls", "cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
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
    
    #Upper Confidence Intervals
    CIpos_hr01 <- hazrats[1,"CIpos"]
    CIpos_hr02 <- hazrats[2,"CIpos"]
    CIpos_hr03 <- hazrats[3,"CIpos"]
    CIpos_hr04 <- hazrats[4,"CIpos"]
    CIpos_hr05 <- hazrats[5,"CIpos"]
    CIpos_hr07 <- hazrats[6,"CIpos"]
    CIpos_hr08 <- hazrats[7,"CIpos"]
    CIpos_hr09 <- hazrats[8,"CIpos"]
    CIpos_hr10 <- hazrats[9,"CIpos"]
    CIpos_hr11 <- hazrats[10,"CIpos"]
    
    #Lower Confidence Intervals
    CIneg_hr01 <- hazrats[1,"CIneg"]
    CIneg_hr02 <- hazrats[2,"CIneg"]
    CIneg_hr03 <- hazrats[3,"CIneg"]
    CIneg_hr04 <- hazrats[4,"CIneg"]
    CIneg_hr05 <- hazrats[5,"CIneg"]
    CIneg_hr07 <- hazrats[6,"CIneg"]
    CIneg_hr08 <- hazrats[7,"CIneg"]
    CIneg_hr09 <- hazrats[8,"CIneg"]
    CIneg_hr10 <- hazrats[9,"CIneg"]
    CIneg_hr11 <- hazrats[10,"CIneg"]
    
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
    
    #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
    incidence$CIpos_i6 <- (incidence$incidence*incidence$population) / ((props06 * incidence$population) + (CIpos_hr01 * (props01 * incidence$population)) + (CIpos_hr02 * (props02 * incidence$population)) + (CIpos_hr03 * (props03 * incidence$population)) + (CIpos_hr04 * (props04 * incidence$population)) + (CIpos_hr05 * (props05 * incidence$population)) + (CIpos_hr07 * (props07 * incidence$population)) + (CIpos_hr08 * (props08 * incidence$population)) + (CIpos_hr09 * (props09 * incidence$population)) + (CIpos_hr10 * (props10 * incidence$population)) + (CIpos_hr11 * (props11 * incidence$population))) 
    incidence$CIpos_i6[is.na(incidence$CIpos_i6)] <- 0
    incidence$CIpos_i1 <- incidence$CIpos_i6 * CIpos_hr01
    incidence$CIpos_i2 <- incidence$CIpos_i6 * CIpos_hr02
    incidence$CIpos_i3 <- incidence$CIpos_i6 * CIpos_hr03
    incidence$CIpos_i4 <- incidence$CIpos_i6 * CIpos_hr04
    incidence$CIpos_i5 <- incidence$CIpos_i6 * CIpos_hr05
    incidence$CIpos_i7 <- incidence$CIpos_i6 * CIpos_hr07
    incidence$CIpos_i8 <- incidence$CIpos_i6 * CIpos_hr08
    incidence$CIpos_i9 <- incidence$CIpos_i6 * CIpos_hr09
    incidence$CIpos_i10 <- incidence$CIpos_i6 * CIpos_hr10
    incidence$CIpos_i11 <- incidence$CIpos_i6 * CIpos_hr11
    
    #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
    incidence$CIneg_i6 <- (incidence$incidence*incidence$population) / ((props06 * incidence$population) + (CIneg_hr01 * (props01 * incidence$population)) + (CIneg_hr02 * (props02 * incidence$population)) + (CIneg_hr03 * (props03 * incidence$population)) + (CIneg_hr04 * (props04 * incidence$population)) + (CIneg_hr05 * (props05 * incidence$population)) + (CIneg_hr07 * (props07 * incidence$population)) + (CIneg_hr08 * (props08 * incidence$population)) + (CIneg_hr09 * (props09 * incidence$population)) + (CIneg_hr10 * (props10 * incidence$population)) + (CIneg_hr11 * (props11 * incidence$population))) 
    incidence$CIneg_i6[is.na(incidence$CIneg_i6)] <- 0
    incidence$CIneg_i1 <- incidence$CIneg_i6 * CIneg_hr01
    incidence$CIneg_i2 <- incidence$CIneg_i6 * CIneg_hr02
    incidence$CIneg_i3 <- incidence$CIneg_i6 * CIneg_hr03
    incidence$CIneg_i4 <- incidence$CIneg_i6 * CIneg_hr04
    incidence$CIneg_i5 <- incidence$CIneg_i6 * CIneg_hr05
    incidence$CIneg_i7 <- incidence$CIneg_i6 * CIneg_hr07
    incidence$CIneg_i8 <- incidence$CIneg_i6 * CIneg_hr08
    incidence$CIneg_i9 <- incidence$CIneg_i6 * CIneg_hr09
    incidence$CIneg_i10 <- incidence$CIneg_i6 * CIneg_hr10
    incidence$CIneg_i11 <- incidence$CIneg_i6 * CIneg_hr11
    
    ###################################################
    
    lifetimerisk <- data.frame(NULL)
    for(i in 1:11){
      #Calculate hazard
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
      
      #CIpos hazard
      incidence[[paste0("CIpos_hazard",i)]] <- incidence[[paste0("CIpos_i",i)]] / (1 - incidence$prevalence)
      
      #CIneg hazard
      incidence[[paste0("CIneg_hazard",i)]] <- incidence[[paste0("CIneg_i",i)]] / (1 - incidence$prevalence)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #CIpos risk
      incidence[[paste0("CIpos_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIpos_hazard",i)]])
      
      #CIneg risk
      incidence[[paste0("CIneg_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIneg_hazard",i)]])
      
      #Calculate the probability of remaining disease free during the age interval - death as a competing interest not considered in this section. 
      #incidence[[paste0("cumulativerisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]])
      #incidence[[paste0("survival_noncomp",i)]] <- exp(-5*incidence[[paste0("cumulativerisk",i)]])
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      #incidence[[paste0("lifetimerisk_noncomp",i)]] <- cumsum(incidence[[paste0("survival_noncomp",i)]]*incidence[[paste0("risk",i)]])*100
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
      
      #CIpos mortandrisk
      incidence[[paste0("CIpos_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIpos_hazard",i)]] + incidence$mortality_rate)
      
      #CIneg mortandrisk
      incidence[[paste0("CIneg_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIneg_hazard",i)]] + incidence$mortality_rate)
      
      #Survival
      incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
      
      #CIpos survival
      incidence[[paste0("CIpos_survival",i)]] <- exp(-5*incidence[[paste0("CIpos_mortandrisk",i)]])
      
      #CIneg survival
      incidence[[paste0("CIneg_survival",i)]] <- exp(-5*incidence[[paste0("CIneg_mortandrisk",i)]])
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      #CIpos lifetime risk
      incidence[[paste0("CIpos_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIpos_survival",i)]]*incidence[[paste0("CIpos_risk",i)]])*100
      
      #CIneg lifetime risk
      incidence[[paste0("CIneg_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIneg_survival",i)]]*incidence[[paste0("CIneg_risk",i)]])*100
      
      result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]], incidence[[paste0("CIpos_lifetimerisk",i)]], incidence[[paste0("CIneg_lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk", "LifetimeRisk_CIpos", "LifetimeRisk_CIneg")
    
    #Plot all as well as overall lifetime risk
    lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
    
    write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_FinnGen.csv"))
    
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
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_FinnGen.png"), height=10 , width=10, dpi=1200)
    
    #Considering confidence intervals
    riskwithintervals <- subset(lifetimerisk, Group=="Group1" | Group=="Group6" | Group=="Group11")
    
    ggplot(riskwithintervals, aes(Age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
      geom_point() +
      geom_ribbon(aes(ymin=LifetimeRisk_CIneg, ymax=LifetimeRisk_CIpos, fill=Group), alpha=0.2) +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group', fill='PRS Group') +
      scale_color_hue(labels = c("0-1%", "40-60%", "99-100%")) +
      scale_fill_hue(labels = c("0-1%", "40-60%", "99-100%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_ConfidenceIntervals_FinnGen.png"), height=10 , width=10, dpi=1200)
    
  }
}

########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################

#Rheumatoid Arthritis and T1D need to be run separately as only have 10 groups
gbd_phenos <- c("Rheumatoid arthritis", "Diabetes mellitus type 1")
hr_phenos <- c("RHEUMA_SEROPOS_OTH", "T1D")

for(j in 1:length(gbd_phenos)){
  
  for(k in c("Male", "Female")){
    
    print(gbd_phenos[j])
    
    print(k)
    
    #Incidence data - basic pre-processing
    incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBD.csv", data.table=FALSE)
    
    #Incidence data to be replaced with that for males and females for prostate cancer and breast cancer respectively. Came from a separate dataset to reduce size of the full dataset. 
    incidence <- subset(incidence, sex==k & cause==gbd_phenos[j] & location=="Norway")
    incidence <- incidence[,c("location","age","cause","metric","val")]
    
    incidence$val <- as.numeric(incidence$val)
    
    population <- c()
    for(i in unique(incidence$age)){
      print(i)
      subby <- subset(incidence, age==i)
      print(subby)
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
    prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBD.csv", data.table=FALSE)
    
    prevalence <- subset(prevalence, sex==k & cause==gbd_phenos[j] & location=="Norway")
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
    mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBD.csv", data.table=FALSE)
    
    mortality <- mortality[,c("location","sex","age","cause","metric","val")]
    
    mortality <- subset(mortality, sex==k & (cause==gbd_phenos[j] | cause=="All causes") & location=="Norway")
    
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
    
    #Hazard Ratios - maddavat
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/HazardRatios_",k,"Sample_FinnGen.csv"), data.table = FALSE)
    colnames(hazrats) <- c("phenotype", "prs", "group", "controls", "cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
    hazrats <- subset(hazrats, phenotype==hr_phenos[j])
    
    #Hazard Ratios
    hr01 <- hazrats[1,"HR"]
    hr02 <- hazrats[2,"HR"]
    hr03 <- hazrats[3,"HR"]
    hr04 <- hazrats[4,"HR"]
    hr06 <- hazrats[5,"HR"]
    hr07 <- hazrats[6,"HR"]
    hr08 <- hazrats[7,"HR"]
    hr09 <- hazrats[8,"HR"]
    hr10 <- hazrats[9,"HR"]
    
    #Upper Confidence Intervals
    CIpos_hr01 <- hazrats[1,"CIpos"]
    CIpos_hr02 <- hazrats[2,"CIpos"]
    CIpos_hr03 <- hazrats[3,"CIpos"]
    CIpos_hr04 <- hazrats[4,"CIpos"]
    CIpos_hr06 <- hazrats[5,"CIpos"]
    CIpos_hr07 <- hazrats[6,"CIpos"]
    CIpos_hr08 <- hazrats[7,"CIpos"]
    CIpos_hr09 <- hazrats[8,"CIpos"]
    CIpos_hr10 <- hazrats[9,"CIpos"]
    
    #Lower Confidence Intervals
    CIneg_hr01 <- hazrats[1,"CIneg"]
    CIneg_hr02 <- hazrats[2,"CIneg"]
    CIneg_hr03 <- hazrats[3,"CIneg"]
    CIneg_hr04 <- hazrats[4,"CIneg"]
    CIneg_hr06 <- hazrats[5,"CIneg"]
    CIneg_hr07 <- hazrats[6,"CIneg"]
    CIneg_hr08 <- hazrats[7,"CIneg"]
    CIneg_hr09 <- hazrats[8,"CIneg"]
    CIneg_hr10 <- hazrats[9,"CIneg"]

    #Proportions - 0.2 by definition of PRS group
    props01 <- 0.05
    props02 <- 0.05
    props03 <- 0.1
    props04 <- 0.2
    props05 <- 0.2
    props06 <- 0.2
    props07 <- 0.1
    props08 <- 0.05
    props9 <- 0.04
    props10 <- 0.01
    
    #Estimate incidence attributable to different distributions of PRS 
    incidence$i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population))) 
    incidence$i5[is.na(incidence$i5)] <- 0
    incidence$i1 <- incidence$i5 * hr01
    incidence$i2 <- incidence$i5 * hr02
    incidence$i3 <- incidence$i5 * hr03
    incidence$i4 <- incidence$i5 * hr04
    incidence$i6 <- incidence$i5 * hr06
    incidence$i7 <- incidence$i5 * hr07
    incidence$i8 <- incidence$i5 * hr08
    incidence$i9 <- incidence$i5 * hr09
    incidence$i10 <- incidence$i5 * hr10
    
    #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
    incidence$CIpos_i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (CIpos_hr01 * (props01 * incidence$population)) + (CIpos_hr02 * (props02 * incidence$population)) + (CIpos_hr03 * (props03 * incidence$population)) + (CIpos_hr04 * (props04 * incidence$population)) + (CIpos_hr06 * (props06 * incidence$population)) + (CIpos_hr07 * (props07 * incidence$population)) + (CIpos_hr08 * (props08 * incidence$population)) + (CIpos_hr09 * (props09 * incidence$population)) + (CIpos_hr10 * (props10 * incidence$population))) 
    incidence$CIpos_i5[is.na(incidence$CIpos_i5)] <- 0
    incidence$CIpos_i1 <- incidence$CIpos_i5 * CIpos_hr01
    incidence$CIpos_i2 <- incidence$CIpos_i5 * CIpos_hr02
    incidence$CIpos_i3 <- incidence$CIpos_i5 * CIpos_hr03
    incidence$CIpos_i4 <- incidence$CIpos_i5 * CIpos_hr04
    incidence$CIpos_i6 <- incidence$CIpos_i5 * CIpos_hr06
    incidence$CIpos_i7 <- incidence$CIpos_i5 * CIpos_hr07
    incidence$CIpos_i8 <- incidence$CIpos_i5 * CIpos_hr08
    incidence$CIpos_i9 <- incidence$CIpos_i5 * CIpos_hr09
    incidence$CIpos_i10 <- incidence$CIpos_i5 * CIpos_hr10
    
    #Estimate incidence attributable to different distributions of PRS - Upper confidence interval
    incidence$CIneg_i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (CIneg_hr01 * (props01 * incidence$population)) + (CIneg_hr02 * (props02 * incidence$population)) + (CIneg_hr03 * (props03 * incidence$population)) + (CIneg_hr04 * (props04 * incidence$population)) + (CIneg_hr06 * (props06 * incidence$population)) + (CIneg_hr07 * (props07 * incidence$population)) + (CIneg_hr08 * (props08 * incidence$population)) + (CIneg_hr09 * (props09 * incidence$population)) + (CIneg_hr10 * (props10 * incidence$population)) + (CIneg_hr11 * (props11 * incidence$population))) 
    incidence$CIneg_i5[is.na(incidence$CIneg_i5)] <- 0
    incidence$CIneg_i1 <- incidence$CIneg_i5 * CIneg_hr01
    incidence$CIneg_i2 <- incidence$CIneg_i5 * CIneg_hr02
    incidence$CIneg_i3 <- incidence$CIneg_i5 * CIneg_hr03
    incidence$CIneg_i4 <- incidence$CIneg_i5 * CIneg_hr04
    incidence$CIneg_i6 <- incidence$CIneg_i5 * CIneg_hr06
    incidence$CIneg_i7 <- incidence$CIneg_i5 * CIneg_hr07
    incidence$CIneg_i8 <- incidence$CIneg_i5 * CIneg_hr08
    incidence$CIneg_i9 <- incidence$CIneg_i5 * CIneg_hr09
    incidence$CIneg_i10 <- incidence$CIneg_i5 * CIneg_hr10
    
    ###################################################
    
    lifetimerisk <- data.frame(NULL)
    for(i in 1:10){
      #Calculate hazard
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
      
      #CIpos hazard
      incidence[[paste0("CIpos_hazard",i)]] <- incidence[[paste0("CIpos_i",i)]] / (1 - incidence$prevalence)
      
      #CIneg hazard
      incidence[[paste0("CIneg_hazard",i)]] <- incidence[[paste0("CIneg_i",i)]] / (1 - incidence$prevalence)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #CIpos risk
      incidence[[paste0("CIpos_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIpos_hazard",i)]])
      
      #CIneg risk
      incidence[[paste0("CIneg_risk",i)]] <- 1 - exp(-5*incidence[[paste0("CIneg_hazard",i)]])
      
      #Calculate the probability of remaining disease free during the age interval - death as a competing interest not considered in this section. 
      #incidence[[paste0("cumulativerisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]])
      #incidence[[paste0("survival_noncomp",i)]] <- exp(-5*incidence[[paste0("cumulativerisk",i)]])
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      #incidence[[paste0("lifetimerisk_noncomp",i)]] <- cumsum(incidence[[paste0("survival_noncomp",i)]]*incidence[[paste0("risk",i)]])*100
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
      
      #CIpos mortandrisk
      incidence[[paste0("CIpos_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIpos_hazard",i)]] + incidence$mortality_rate)
      
      #CIneg mortandrisk
      incidence[[paste0("CIneg_mortandrisk",i)]] <- cumsum(incidence[[paste0("CIneg_hazard",i)]] + incidence$mortality_rate)
      
      #Survival
      incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
      
      #CIpos survival
      incidence[[paste0("CIpos_survival",i)]] <- exp(-5*incidence[[paste0("CIpos_mortandrisk",i)]])
      
      #CIneg survival
      incidence[[paste0("CIneg_survival",i)]] <- exp(-5*incidence[[paste0("CIneg_mortandrisk",i)]])
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      #CIpos lifetime risk
      incidence[[paste0("CIpos_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIpos_survival",i)]]*incidence[[paste0("CIpos_risk",i)]])*100
      
      #CIneg lifetime risk
      incidence[[paste0("CIneg_lifetimerisk",i)]] <- cumsum(incidence[[paste0("CIneg_survival",i)]]*incidence[[paste0("CIneg_risk",i)]])*100
      
      result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]], incidence[[paste0("CIpos_lifetimerisk",i)]], incidence[[paste0("CIneg_lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk", "LifetimeRisk_CIpos", "LifetimeRisk_CIneg")
    
    #Plot all as well as overall lifetime risk
    lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10"))
    
    write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_FinnGen.csv"))
    
    #Not considering confidence intervals
    ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
      geom_point() +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group') +
      scale_color_hue(labels = c("0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_FinnGen.png"), height=10 , width=10, dpi=1200)
    
    #Considering confidence intervals
    riskwithintervals <- subset(lifetimerisk, Group=="Group1" | Group=="Group5" | Group=="Group10")
    
    ggplot(riskwithintervals, aes(Age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
      geom_point() +
      geom_ribbon(aes(ymin=LifetimeRisk_CIneg, ymax=LifetimeRisk_CIpos, fill=Group), alpha=0.2) +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group', fill='PRS Group') +
      scale_color_hue(labels = c("0-5%", "40-60%", "99-100%")) +
      scale_fill_hue(labels = c("0-5%", "40-60%", "99-100%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_",k,"_LifetimeRisk_ConfidenceIntervals_FinnGen.png"), height=10 , width=10, dpi=1200)
    
  }
}

################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################


#Age stratification results

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Idiopathic epilepsy", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Diabetes mellitus type 2")
hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "G6_EPLEPSY", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "T2D")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  #Read in GBD incidence data 
  incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Incidence.csv", data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  
  population <- c()
  for(i in unique(incidence$age)){
    #print(i)
    subby <- subset(incidence, age==i)
    #print(subby)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0
  
  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population)
  incidence <- incidence[,c("location","age","cause","metric","population","incidence")]
  
  prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Prevalence.csv", data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Mortality.csv", data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  
  bcpc_mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="Norway")
  
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
  
  #Hazard Ratios - maddavat
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/",hr_phenos[j],"_Age_Specific_Hazards.csv"), data.table = FALSE)
  hazrats <- hazrats[,-1]
  colnames(hazrats) <- c("Age","Beta","HR","Group")
  
  #Hazard Ratios
  hr01 <- hazrats[hazrats$Group=="< 1%","HR"]
  hr02 <- hazrats[hazrats$Group=="1-5%","HR"]
  hr03 <- hazrats[hazrats$Group=="5-10%","HR"]
  hr04 <- hazrats[hazrats$Group=="10-20%","HR"]
  hr05 <- hazrats[hazrats$Group=="20-40%","HR"]
  hr07 <- hazrats[hazrats$Group=="60-80%","HR"]
  hr08 <- hazrats[hazrats$Group=="80-90%","HR"]
  hr09 <- hazrats[hazrats$Group=="90-95%","HR"]
  hr10 <- hazrats[hazrats$Group=="95-99%","HR"]
  hr11 <- hazrats[hazrats$Group=="> 99%","HR"]
  
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
    incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
  write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.csv"))
  
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
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.png"), height=10 , width=10, dpi=1200)
  
}

########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################

#For the phenotypes with only 10 groupings. 

gbd_phenos <- c("Colon and rectum cancer", "Gout", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis")
hr_phenos <- c("C3_COLORECTAL", "GOUT", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  #Read in GBD incidence data 
  incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Incidence.csv", data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  
  population <- c()
  for(i in unique(incidence$age)){
    #print(i)
    subby <- subset(incidence, age==i)
    #print(subby)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0
  
  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population)
  incidence <- incidence[,c("location","age","cause","metric","population","incidence")]
  
  prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Prevalence.csv", data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Mortality.csv", data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  
  bcpc_mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="Norway")
  
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
  
  #Hazard Ratios - maddavat
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/",hr_phenos[j],"_Age_Specific_Hazards.csv"), data.table = FALSE)
  hazrats <- hazrats[,-1]
  colnames(hazrats) <- c("Age","Beta","HR","Group")
  
  #Hazard Ratios
  hr01 <- hazrats[hazrats$Group=="0-5%","HR"]
  hr02 <- hazrats[hazrats$Group=="5-10%","HR"]
  hr03 <- hazrats[hazrats$Group=="10-20%","HR"]
  hr04 <- hazrats[hazrats$Group=="20-40%","HR"]
  hr06 <- hazrats[hazrats$Group=="60-80%","HR"]
  hr07 <- hazrats[hazrats$Group=="80-90%","HR"]
  hr08 <- hazrats[hazrats$Group=="90-95%","HR"]
  hr09 <- hazrats[hazrats$Group=="95-99%","HR"]
  hr10 <- hazrats[hazrats$Group=="> 99%","HR"]
  
  #Proportions - 0.2 by definition of PRS group
  props01 <- 0.05
  props02 <- 0.05
  props03 <- 0.1
  props04 <- 0.2
  props05 <- 0.2
  props06 <- 0.2
  props07 <- 0.1
  props08 <- 0.05
  props09 <- 0.04
  props10 <- 0.01
  
  #Estimate incidence attributable to different distributions of PRS 
  incidence$i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population))) 
  incidence$i5[is.na(incidence$i5)] <- 0
  incidence$i1 <- incidence$i5 * hr01
  incidence$i2 <- incidence$i5 * hr02
  incidence$i3 <- incidence$i5 * hr03
  incidence$i4 <- incidence$i5 * hr04
  incidence$i6 <- incidence$i5 * hr06
  incidence$i7 <- incidence$i5 * hr07
  incidence$i8 <- incidence$i5 * hr08
  incidence$i9 <- incidence$i5 * hr09
  incidence$i10 <- incidence$i5 * hr10

  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  for(i in 1:10){
    #Calculate hazard
    incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
    
    #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
    incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
    
    #Mortality and risk
    incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
    
    #Survival
    incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10"))
  
  write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.csv"))
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age Range") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group') +
    scale_color_hue(labels = c("0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.png"), height=10 , width=10, dpi=1200)
  
}

########################################################################################################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################################################################################################

#For the phenotypes with only 8 groupings. 

gbd_phenos <- c("Diabetes mellitus type 1")
hr_phenos <- c("T1D")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  #Read in GBD incidence data 
  incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Incidence.csv", data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
  
  incidence <- rbind(incidence, bcpc_incidence)
  
  incidence <- subset(incidence, cause==gbd_phenos[j])
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  
  population <- c()
  for(i in unique(incidence$age)){
    #print(i)
    subby <- subset(incidence, age==i)
    #print(subby)
    poppy <- subby$val[1]/(subby$val[2]/100000)
    population <- c(population, poppy)
  }
  
  population[is.na(population)] <- 0
  
  incidence <- subset(incidence, metric=='Rate')
  incidence <- cbind(incidence, population)
  incidence <- incidence[,c("location","age","cause","metric","population","incidence")]
  
  prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Prevalence.csv", data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Mortality.csv", data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
  
  bcpc_mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="Norway")
  
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
  
  #Hazard Ratios - maddavat
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios/",hr_phenos[j],"_Age_Specific_Hazards.csv"), data.table = FALSE)
  hazrats <- hazrats[,-1]
  colnames(hazrats) <- c("Age","Beta","HR","Group")
  
  #Hazard Ratios
  hr01 <- hazrats[hazrats$Group=="0-20%","HR"]
  hr02 <- hazrats[hazrats$Group=="20-40%","HR"]
  hr04 <- hazrats[hazrats$Group=="60-80%","HR"]
  hr05 <- hazrats[hazrats$Group=="80-90%","HR"]
  hr06 <- hazrats[hazrats$Group=="90-95%","HR"]
  hr07 <- hazrats[hazrats$Group=="95-99%","HR"]
  hr08 <- hazrats[hazrats$Group=="> 99%","HR"]

  #Proportions - 0.2 by definition of PRS group
  props01 <- 0.2
  props02 <- 0.2
  props03 <- 0.2
  props04 <- 0.2
  props05 <- 0.1
  props06 <- 0.05
  props07 <- 0.04
  props08 <- 0.01
  
  #Estimate incidence attributable to different distributions of PRS 
  incidence$i3 <- (incidence$incidence*incidence$population) / ((props03 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population))) 
  incidence$i3[is.na(incidence$i3)] <- 0
  incidence$i1 <- incidence$i3 * hr01
  incidence$i2 <- incidence$i3 * hr02
  incidence$i4 <- incidence$i3 * hr04
  incidence$i5 <- incidence$i3 * hr05
  incidence$i6 <- incidence$i3 * hr06
  incidence$i7 <- incidence$i3 * hr07
  incidence$i8 <- incidence$i3 * hr08

  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  for(i in 1:8){
    #Calculate hazard
    incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
    
    #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
    incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
    
    #Mortality and risk
    incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
    
    #Survival
    incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
    
    #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
    incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
    
    result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
    lifetimerisk <- rbind(lifetimerisk, result)
  }
  
  colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8"))
  
  write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.csv"))
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age Range") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group') +
    scale_color_hue(labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen",hr_phenos[j],"_LifetimeRisk_AgeStratification_FinnGen.png"), height=10 , width=10, dpi=1200)
  
}
