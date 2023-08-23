library(data.table)
library(dplyr)
library(ggplot2)

#CHD across each country - bottom 20% vs top 5%

#Age and sex stratification results

gbd_phenos <- c("Rheumatoid arthritis")
hr_phenos <- c("RHEUMA_SEROPOS_OTH")
directory <- c("RheumatoidArthritis")
sexes <- c("Male", "Female")

countries <- c("Finland","Massachusetts","Norway","United Kingdom","Estonia")
biobank <- c("FinnGen", "PartnersBiobank", "HUNT", "UKBiobank", "EstonianBiobank")

for(j in 1:length(biobank)){
  
  print(biobank[j])
  
  for(k in 1:2){
    
    #Read in GBD incidence data 
    incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBDUpdated.csv", data.table=FALSE)
    incidence <- subset(incidence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
    incidence$age <- factor(incidence$age, levels=c("1-4 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years", "25-29 years", "30-34 years", "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "All ages"))
    incidence <- incidence[order(incidence$age),]
    incidence$val <- as.numeric(incidence$val)
    
    incidence <- subset(incidence, cause==gbd_phenos)
    
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
    
    prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBDUpdated.csv", data.table=FALSE)
    prevalence <- subset(prevalence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
    prevalence$val <- as.numeric(prevalence$val)
    
    prevalence <- subset(prevalence, cause==gbd_phenos)
    
    #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
    prevalence$prevalence <- prevalence$val / 100000
    prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
    
    #Left join to incidence to calculate hazard 
    incidence <- left_join(incidence, prevalence)
    
    #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
    mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBDUpdated.csv", data.table=FALSE)
    mortality <- subset(mortality, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
    
    mortality <- subset(mortality, cause==gbd_phenos | cause=="All causes")
    
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex==sexes[k])
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
    
    cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
    cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val")]
    colnames(cause_specific_mortality)[5] <- c("cause_specific_rate")
    
    mortality <- left_join(cause_specific_mortality, all_cause_mortality)
    mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
    mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
    
    mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
    mortality <- mortality[,c("location","sex","age","cause","mortality_rate")]
    
    #Merge mortality data to incidence data
    incidence <- left_join(incidence, mortality)
    
    incidence <- subset(incidence, age!="All ages" & age!="80-84 years" & age!="85-89 years" & age!="90-94 years")
    
    #######################################################################################################################################################################################
    
    #Hazard Ratios - maddavat
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[j],"/Age/",hr_phenos,"_Age_Specific_Hazards.csv"), data.table = FALSE)
    hazrats <- hazrats[,-1]
    colnames(hazrats) <- c("Age","Beta","HR","Group")
    
    #Hazard Ratios
    hr01 <- hazrats[hazrats$Group=="< 20%","HR"]
    hr02 <- hazrats[hazrats$Group=="20-40%","HR"]
    hr04 <- hazrats[hazrats$Group=="60-80%","HR"]
    hr05 <- hazrats[hazrats$Group=="80-90%","HR"]
    hr06 <- hazrats[hazrats$Group=="90-95%","HR"]
    hr07 <- hazrats[hazrats$Group=="> 95%","HR"]
    
    #Proportions - 0.2 by definition of PRS group
    props01 <- 0.2
    props02 <- 0.2
    props03 <- 0.2
    props04 <- 0.2
    props05 <- 0.1
    props06 <- 0.05
    props07 <- 0.05
    
    #Estimate incidence attributable to different distributions of PRS 
    incidence$i3 <- (incidence$incidence*incidence$population) / ((props03 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population))) 
    incidence$i3[is.na(incidence$i3)] <- 0
    incidence$i1 <- incidence$i3 * hr01
    incidence$i2 <- incidence$i3 * hr02
    incidence$i4 <- incidence$i3 * hr04
    incidence$i5 <- incidence$i3 * hr05
    incidence$i6 <- incidence$i3 * hr06
    incidence$i7 <- incidence$i3 * hr07
    
    
    ###################################################
    
    Groups <- c("< 20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "> 95%")
    lifetimerisk <- data.frame(NULL)
    for(i in 1:7){
      #Calculate hazard
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
      
      #Survival
      incidence[[paste0("survival",i)]] <- 1
      
      for(l in 2:nrow(incidence)){
        incidence[[paste0("survival",i)]][l] <- exp(-5*incidence[[paste0("mortandrisk",i)]][l-1])
      }
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      result <- data.frame(incidence$age, Groups[i], incidence[[paste0("lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
    
    #Plot all as well as overall lifetime risk
    lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"))
    lifetimerisk$age <- rep(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80),7)
    lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
    
    write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],".csv"))
    
    #Not considering confidence intervals
    ggplot(lifetimerisk, aes(age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
      geom_point() +
      xlab("Age") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Strata') +
      scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
      theme(legend.text = element_text(size = 24),
            legend.title = element_text(size = 28),
            axis.title.x = element_text(size = 28),
            axis.text.x = element_text(size = 24),
            axis.title.y = element_text(size = 28),
            axis.text.y = element_text(size = 24))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],".png"), height=6, width=8, dpi=300)
    
  }
}

###############################################################################################################################################################################################
###############################################################################################################################################################################################

library(data.table)
library(dplyr)
library(ggplot2)

#Supplementary analysis - Finngen extension to bottom and top 1%

gbd_phenos <- c("Rheumatoid arthritis")
hr_phenos <- c("RHEUMA_SEROPOS_OTH")
directory <- c("RheumatoidArthritis")
sexes <- c("Male", "Female")

countries <- c("Finland")
biobank <- c("FinnGen")

for(j in 1:length(biobank)){
  
  print(biobank[j])
  
  for(k in 1:2){
    
    print(sexes[k])
    
    #Read in GBD incidence data 
    incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBDUpdated.csv", data.table=FALSE)
    incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & sex==sexes[k] & location==countries[j])
    incidence$age <- factor(incidence$age, levels=c("1-4 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years", "25-29 years", "30-34 years", "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "All ages"))
    incidence <- incidence[order(incidence$age),]
    incidence$val <- as.numeric(incidence$val)
    
    incidence <- subset(incidence, cause==gbd_phenos)
    
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
    
    prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBDUpdated.csv", data.table=FALSE)
    prevalence <- subset(prevalence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
    prevalence$val <- as.numeric(prevalence$val)
    
    prevalence <- subset(prevalence, cause==gbd_phenos)
    
    #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
    prevalence$prevalence <- prevalence$val / 100000
    prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
    
    #Left join to incidence to calculate hazard 
    incidence <- left_join(incidence, prevalence)
    
    #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
    mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBDUpdated.csv", data.table=FALSE)
    mortality <- subset(mortality, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
    
    mortality <- subset(mortality, cause==gbd_phenos | cause=="All causes")
    
    all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex==sexes[k])
    all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val")]
    colnames(all_cause_mortality)[4] <- c("all_cause_rate")
    
    cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
    cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val")]
    colnames(cause_specific_mortality)[5] <- c("cause_specific_rate")
    
    mortality <- left_join(cause_specific_mortality, all_cause_mortality)
    mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
    mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
    
    mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
    mortality <- mortality[,c("location","sex","age","cause","mortality_rate")]
    
    #Merge mortality data to incidence data
    incidence <- left_join(incidence, mortality)
    
    incidence <- subset(incidence, age!="All ages" & age!="80-84 years" & age!="85-89 years" & age!="90-94 years")
    
    #######################################################################################################################################################################################
    
    #Hazard Ratios - maddavat
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[j],"/Age/",hr_phenos,"_Age_Specific_Hazards_OnePercentTails.csv"), data.table = FALSE)
    hazrats <- hazrats[,-1]
    colnames(hazrats) <- c("Age","Beta","HR","Group")
    
    #Hazard Ratios
    hr01 <- hazrats[hazrats$Group=="< 5%","HR"]
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
    
    Groups <- c("< 5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%")
    lifetimerisk <- data.frame(NULL)
    for(i in 1:10){
      #Calculate hazard
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate)
      
      #Survival
      incidence[[paste0("survival",i)]] <- 1
      
      for(l in 2:nrow(incidence)){
        incidence[[paste0("survival",i)]][l] <- exp(-5*incidence[[paste0("mortandrisk",i)]][l-1])
      }
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      result <- data.frame(incidence$age, Groups[i], incidence[[paste0("lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group","LifetimeRisk")
    
    #Plot all as well as overall lifetime risk
    lifetimerisk$Age <- factor(lifetimerisk$Age, levels=c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"))
    lifetimerisk$age <- rep(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80),10)
    lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "< 5%"))
    
    write.csv(lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],"_OnePercentTails.csv"))
    
    #Not considering confidence intervals
    ggplot(lifetimerisk, aes(age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
      geom_point() +
      xlab("Age") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Strata') +
      scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "0-5%")) +
      theme(legend.text = element_text(size = 24),
            legend.title = element_text(size = 28),
            axis.title.x = element_text(size = 28),
            axis.text.x = element_text(size = 24),
            axis.title.y = element_text(size = 28),
            axis.text.y = element_text(size = 24))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],"_OnePercentTails.png"), height=6, width=8, dpi=300)
    
  }
}