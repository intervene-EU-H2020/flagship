library(data.table)
library(dplyr)
library(ggplot2)

#Type 2 Diabetes across each country - bottom 20% vs top 5%

#Age and sex stratification results

gbd_phenos <- c("Rheumatoid arthritis")
hr_phenos <- c("RHEUMA_SEROPOS_OTH")
directory <- c("RheumatoidArthritis")
sexes <- c("Male", "Female")

#countries <- c("Finland","Massachusetts","Norway","United Kingdom","Estonia")
#biobank <- c("FinnGen", "PartnersBiobank", "HUNT", "UKBiobank", "EstonianBiobank")

countries <- c("Massachusetts")
biobank <- c("PartnersBiobank")

for(j in 1:length(biobank)){
  
  print(biobank[j])
  
  for(k in 1:2){
    
    lifetimerisks <- data.frame(Age=rep(c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"),7),
                                Group=c(rep("< 20%",16), rep("20-40%",16), rep("40-60%",16), rep("60-80%",16), rep("80-90%",16), rep("90-95%",16), rep("> 95%",16)))
    
    m <- 0 
    
    while(m < 1000){
      
      #Read in GBD incidence data 
      incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBDUpdated.csv", data.table=FALSE)
      incidence <- subset(incidence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      incidence$age <- factor(incidence$age, levels=c("1-4 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years", "25-29 years", "30-34 years", "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "All ages"))
      incidence <- incidence[order(incidence$age),]
      incidence$val <- as.numeric(incidence$val)
      incidence$upper <- as.numeric(incidence$upper)
      incidence$lower <- as.numeric(incidence$lower)
      
      incidence <- subset(incidence, cause==gbd_phenos)
      
      #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
      incidence$val <- incidence$val / 100000
      incidence$upper <- incidence$upper / 100000
      incidence$lower <- incidence$lower / 100000
      
      population <- c()
      for(i in unique(incidence$age)){
        subby <- subset(incidence, age==i)
        poppy <- subby$val[1]/(subby$val[2]/100000)
        population <- c(population, poppy)
      }
      
      population[is.na(population)] <- 0
      
      incidence <- subset(incidence, metric=='Rate')
      incidence <- cbind(incidence, population) 
      
      #Sample from the distribution to get an incidence value
      incidence$sd_pos <- (incidence$upper - incidence$val) / 1.96
      incidence$sd_neg <- (incidence$val - incidence$lower) / 1.96
      incidence$sd <- rowMeans(incidence[,c("sd_pos","sd_neg")], na.rm=T)
      
      incidence$incidence_sample <- with(incidence, rnorm(nrow(incidence), incidence$val, incidence$sd))
      incidence$incidence_sample[incidence$incidence_sample < 0] <- 0
      
      incidence <- incidence[,c("location","age","cause","metric","population","incidence_sample")]
      
      prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBDUpdated.csv", data.table=FALSE)
      prevalence <- subset(prevalence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      prevalence$val <- as.numeric(prevalence$val)
      prevalence$upper <- as.numeric(prevalence$upper)
      prevalence$lower <- as.numeric(prevalence$lower)
      
      prevalence <- subset(prevalence, cause==gbd_phenos)
      
      #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
      prevalence$val <- prevalence$val / 100000
      prevalence$upper <- prevalence$upper / 100000
      prevalence$lower <- prevalence$lower / 100000
      
      prevalence <- subset(prevalence, metric=='Rate')
      
      #Sample from the distribution to get a prevalence value
      prevalence$sd_pos <- (prevalence$upper - prevalence$val) / 1.96
      prevalence$sd_neg <- (prevalence$val - prevalence$lower) / 1.96
      prevalence$sd <- rowMeans(prevalence[,c("sd_pos","sd_neg")], na.rm=T)
      
      prevalence$prevalence_sample <- with(prevalence, rnorm(nrow(prevalence), prevalence$val, prevalence$sd))
      prevalence$prevalence_sample[prevalence$prevalence_sample < 0] <- 0
      
      prevalence <- prevalence[,c("location","age","cause","metric","prevalence_sample")]
      
      #Left join to incidence to calculate hazard 
      incidence <- left_join(incidence, prevalence)
      
      #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
      mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBDUpdated.csv", data.table=FALSE)
      mortality <- subset(mortality, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      
      mortality <- subset(mortality, cause==gbd_phenos | cause=="All causes")
      
      all_cause_mortality <- subset(mortality, sex==sexes[k] & cause=="All causes" & metric=="Rate")
      all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val","upper","lower")]
      
      #Sample from the distribution to get an all cause mortality value
      all_cause_mortality$sd_pos <- (all_cause_mortality$upper - all_cause_mortality$val) / 1.96
      all_cause_mortality$sd_neg <- (all_cause_mortality$val - all_cause_mortality$lower) / 1.96
      all_cause_mortality$sd <- rowMeans(all_cause_mortality[,c("sd_pos","sd_neg")], na.rm=T)
      
      all_cause_mortality$all_cause_rate_sample <- with(all_cause_mortality, rnorm(nrow(all_cause_mortality), all_cause_mortality$val, all_cause_mortality$sd))
      all_cause_mortality$all_cause_rate_sample[all_cause_mortality$all_cause_rate_sample < 0] <- 0
      
      cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
      cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val","upper","lower")]
      
      #Sample from the distribution to get an all cause mortality value
      cause_specific_mortality$sd_pos <- (cause_specific_mortality$upper - cause_specific_mortality$val) / 1.96
      cause_specific_mortality$sd_neg <- (cause_specific_mortality$val - cause_specific_mortality$lower) / 1.96
      cause_specific_mortality$sd <- rowMeans(cause_specific_mortality[,c("sd_pos","sd_neg")], na.rm=T)
      
      cause_specific_mortality$cause_specific_rate_sample <- with(cause_specific_mortality, rnorm(nrow(cause_specific_mortality), cause_specific_mortality$val, cause_specific_mortality$sd))
      cause_specific_mortality$cause_specific_rate_sample[cause_specific_mortality$cause_specific_rate_sample < 0] <- 0
      
      mortality <- left_join(cause_specific_mortality, all_cause_mortality, by=c("location","sex","age"))
      
      mortality$all_cause_rate_sample <- as.numeric(mortality$all_cause_rate_sample)
      mortality$cause_specific_rate_sample <- as.numeric(mortality$cause_specific_rate_sample)
      mortality$mortality_rate_sample <- (mortality$all_cause_rate_sample - mortality$cause_specific_rate_sample)/100000
      
      mortality <- mortality[,c("location","age","cause","mortality_rate_sample")]
      
      #Merge mortality data to incidence data
      incidence <- left_join(incidence, mortality)
      
      incidence <- subset(incidence, age!="All ages" & age!="80-84 years" & age!="85-89 years" & age!="90-94 years")
      
      #######################################################################################################################################################################################
      
      #Hazard Ratios - maddavat
      
      #Read in the hazard ratios and allocate to variables...
      hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[j],"/Age/",hr_phenos,"_Age_Specific_Hazards_Bootstrapped.csv"), data.table = FALSE)
      hazrats <- hazrats[,-1]
      hazrats$sdpos <- hazrats$CIpos - hazrats$Beta 
      hazrats$sdneg <- hazrats$Beta - hazrats$CIneg
      hazrats$SD <- rowMeans(hazrats[,c("sdpos","sdneg")], na.rm=T)
      
      #Hazard Ratios
      hr01 <- with(hazrats[hazrats$Group=="< 20%",], exp(rnorm(nrow(hazrats[hazrats$Group=="< 20%",]), hazrats[hazrats$Group=="< 20%","Beta"], hazrats[hazrats$Group=="< 20%","SD"])))
      hr02 <- with(hazrats[hazrats$Group=="20-40%",], exp(rnorm(nrow(hazrats[hazrats$Group=="20-40%",]), hazrats[hazrats$Group=="20-40%","Beta"], hazrats[hazrats$Group=="20-40%","SD"])))
      hr04 <- with(hazrats[hazrats$Group=="60-80%",], exp(rnorm(nrow(hazrats[hazrats$Group=="60-80%",]), hazrats[hazrats$Group=="60-80%","Beta"], hazrats[hazrats$Group=="60-80%","SD"])))
      hr05 <- with(hazrats[hazrats$Group=="80-90%",], exp(rnorm(nrow(hazrats[hazrats$Group=="80-90%",]), hazrats[hazrats$Group=="80-90%","Beta"], hazrats[hazrats$Group=="80-90%","SD"])))
      hr06 <- with(hazrats[hazrats$Group=="90-95%",], exp(rnorm(nrow(hazrats[hazrats$Group=="90-95%",]), hazrats[hazrats$Group=="90-95%","Beta"], hazrats[hazrats$Group=="90-95%","SD"])))
      hr07 <- with(hazrats[hazrats$Group=="> 95%",], exp(rnorm(nrow(hazrats[hazrats$Group=="> 95%",]), hazrats[hazrats$Group=="> 95%","Beta"], hazrats[hazrats$Group=="> 95%","SD"])))
      
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
        
        for(p in 2:nrow(incidence)){
          incidence[[paste0("survival",i)]][p] <- exp(-5*incidence[[paste0("mortandrisk",i)]][p-1])
        }
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
        
        result <- data.frame(incidence$age, Groups[i], incidence[[paste0("lifetimerisk",i)]])
        lifetimerisk <- rbind(lifetimerisk, result)
      }
      
      colnames(lifetimerisk) <- c("Age","Group",paste0("LifetimeRisk",m))
      
      m <- m+1
      
      lifetimerisks <- left_join(lifetimerisks, lifetimerisk)
      
    }
    
    lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2)])
    confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975))
    
    bootstrapped_lifetimerisk <- data.frame(Age=rep(c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"),7),
                                            Group=c(rep("< 20%",16), rep("20-40%",16), rep("40-60%",16), rep("60-80%",16), rep("80-90%",16), rep("90-95%",16), rep("> 95%",16)))
    
    bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
    bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]
    
    #Add in actual lifetime risks
    lifetimeriskactual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],".csv"), select=c("LifetimeRisk"), data.table=FALSE)
    bootstrapped_lifetimerisk <- cbind(bootstrapped_lifetimerisk, lifetimeriskactual)
    
    #Plot all as well as overall lifetime risk
    bootstrapped_lifetimerisk$Age <- factor(bootstrapped_lifetimerisk$Age, levels=c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"))
    bootstrapped_lifetimerisk$age <- rep(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80),7)
    bootstrapped_lifetimerisk$Group <- factor(bootstrapped_lifetimerisk$Group, levels=c(c("< 20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "> 95%")))
    
    write.csv(bootstrapped_lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_Bootstrapped_",biobank[j],".csv"))
    
    #Considering confidence intervals
    riskwithintervals <- subset(bootstrapped_lifetimerisk, Group=="< 20%" | Group=="40-60%" | Group=="> 95%")
    
    ggplot(riskwithintervals, aes(age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
      geom_point() +
      geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group', fill='PRS Group') +
      scale_color_hue(labels = c("< 20%", "40-60%", "> 95%")) +
      scale_fill_hue(labels = c("< 20%", "40-60%", "> 95%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_Bootstrapped_",biobank[j],".png"), height=10 , width=10)
    
  }
}

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

library(data.table)
library(dplyr)
library(ggplot2)

#Supplementary analysis - FinnGen top and bottom percentile

gbd_phenos <- c("Rheumatoid arthritis")
hr_phenos <- c("RHEUMA_SEROPOS_OTH")
directory <- c("RheumatoidArthritis")
sexes <- c("Male", "Female")

countries <- c("Finland")
biobank <- c("FinnGen")

for(j in 1:length(biobank)){
  
  print(biobank[j])
  
  for(k in 1:2){
    
    lifetimerisks <- data.frame(Age=rep(c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"),10),
                                Group=c(rep("< 5%",16), rep("5-10%",16), rep("10-20%",16), rep("20-40%",16), rep("40-60%",16), rep("60-80%",16), rep("80-90%",16), rep("90-95%",16), rep("95-99%",16), rep("> 99%",16)))
    
    m <- 0 
    
    while(m < 1000){
      
      #Read in GBD incidence data 
      incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBDUpdated.csv", data.table=FALSE)
      incidence <- subset(incidence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      incidence$age <- factor(incidence$age, levels=c("1-4 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years", "25-29 years", "30-34 years", "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "All ages"))
      incidence <- incidence[order(incidence$age),]
      incidence$val <- as.numeric(incidence$val)
      incidence$upper <- as.numeric(incidence$upper)
      incidence$lower <- as.numeric(incidence$lower)
      
      incidence <- subset(incidence, cause==gbd_phenos)
      
      #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
      incidence$val <- incidence$val / 100000
      incidence$upper <- incidence$upper / 100000
      incidence$lower <- incidence$lower / 100000
      
      population <- c()
      for(i in unique(incidence$age)){
        subby <- subset(incidence, age==i)
        poppy <- subby$val[1]/(subby$val[2]/100000)
        population <- c(population, poppy)
      }
      
      population[is.na(population)] <- 0
      
      incidence <- subset(incidence, metric=='Rate')
      incidence <- cbind(incidence, population) 
      
      #Sample from the distribution to get an incidence value
      incidence$sd_pos <- (incidence$upper - incidence$val) / 1.96
      incidence$sd_neg <- (incidence$val - incidence$lower) / 1.96
      incidence$sd <- rowMeans(incidence[,c("sd_pos","sd_neg")], na.rm=T)
      
      incidence$incidence_sample <- with(incidence, rnorm(nrow(incidence), incidence$val, incidence$sd))
      incidence$incidence_sample[incidence$incidence_sample < 0] <- 0
      
      incidence <- incidence[,c("location","age","cause","metric","population","incidence_sample")]
      
      prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBDUpdated.csv", data.table=FALSE)
      prevalence <- subset(prevalence, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      prevalence$val <- as.numeric(prevalence$val)
      prevalence$upper <- as.numeric(prevalence$upper)
      prevalence$lower <- as.numeric(prevalence$lower)
      
      prevalence <- subset(prevalence, cause==gbd_phenos)
      
      #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
      prevalence$val <- prevalence$val / 100000
      prevalence$upper <- prevalence$upper / 100000
      prevalence$lower <- prevalence$lower / 100000
      
      prevalence <- subset(prevalence, metric=='Rate')
      
      #Sample from the distribution to get a prevalence value
      prevalence$sd_pos <- (prevalence$upper - prevalence$val) / 1.96
      prevalence$sd_neg <- (prevalence$val - prevalence$lower) / 1.96
      prevalence$sd <- rowMeans(prevalence[,c("sd_pos","sd_neg")], na.rm=T)
      
      prevalence$prevalence_sample <- with(prevalence, rnorm(nrow(prevalence), prevalence$val, prevalence$sd))
      prevalence$prevalence_sample[prevalence$prevalence_sample < 0] <- 0
      
      prevalence <- prevalence[,c("location","age","cause","metric","prevalence_sample")]
      
      #Left join to incidence to calculate hazard 
      incidence <- left_join(incidence, prevalence)
      
      #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
      mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBDUpdated.csv", data.table=FALSE)
      mortality <- subset(mortality, sex==sexes[k] & cause!="Breast cancer" & cause!="Prostate cancer" & location==countries[j])
      
      mortality <- subset(mortality, cause==gbd_phenos | cause=="All causes")
      
      all_cause_mortality <- subset(mortality, sex==sexes[k] & cause=="All causes" & metric=="Rate")
      all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val","upper","lower")]
      
      #Sample from the distribution to get an all cause mortality value
      all_cause_mortality$sd_pos <- (all_cause_mortality$upper - all_cause_mortality$val) / 1.96
      all_cause_mortality$sd_neg <- (all_cause_mortality$val - all_cause_mortality$lower) / 1.96
      all_cause_mortality$sd <- rowMeans(all_cause_mortality[,c("sd_pos","sd_neg")], na.rm=T)
      
      all_cause_mortality$all_cause_rate_sample <- with(all_cause_mortality, rnorm(nrow(all_cause_mortality), all_cause_mortality$val, all_cause_mortality$sd))
      all_cause_mortality$all_cause_rate_sample[all_cause_mortality$all_cause_rate_sample < 0] <- 0
      
      cause_specific_mortality <- subset(mortality, cause!="All causes" & metric=="Rate")
      cause_specific_mortality <- cause_specific_mortality[,c("location","sex","age","cause","val","upper","lower")]
      
      #Sample from the distribution to get an all cause mortality value
      cause_specific_mortality$sd_pos <- (cause_specific_mortality$upper - cause_specific_mortality$val) / 1.96
      cause_specific_mortality$sd_neg <- (cause_specific_mortality$val - cause_specific_mortality$lower) / 1.96
      cause_specific_mortality$sd <- rowMeans(cause_specific_mortality[,c("sd_pos","sd_neg")], na.rm=T)
      
      cause_specific_mortality$cause_specific_rate_sample <- with(cause_specific_mortality, rnorm(nrow(cause_specific_mortality), cause_specific_mortality$val, cause_specific_mortality$sd))
      cause_specific_mortality$cause_specific_rate_sample[cause_specific_mortality$cause_specific_rate_sample < 0] <- 0
      
      mortality <- left_join(cause_specific_mortality, all_cause_mortality, by=c("location","sex","age"))
      
      mortality$all_cause_rate_sample <- as.numeric(mortality$all_cause_rate_sample)
      mortality$cause_specific_rate_sample <- as.numeric(mortality$cause_specific_rate_sample)
      mortality$mortality_rate_sample <- (mortality$all_cause_rate_sample - mortality$cause_specific_rate_sample)/100000
      
      mortality <- mortality[,c("location","age","cause","mortality_rate_sample")]
      
      #Merge mortality data to incidence data
      incidence <- left_join(incidence, mortality)
      
      incidence <- subset(incidence, age!="All ages" & age!="80-84 years" & age!="85-89 years" & age!="90-94 years")
      
      #######################################################################################################################################################################################
      
      #Hazard Ratios - maddavat
      
      #Read in the hazard ratios and allocate to variables...
      hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[j],"/Age/",hr_phenos,"_Age_Specific_Hazards_Bootstrapped_OnePercentTails.csv"), data.table = FALSE)
      hazrats <- hazrats[,-1]
      hazrats$sdpos <- hazrats$CIpos - hazrats$Beta 
      hazrats$sdneg <- hazrats$Beta - hazrats$CIneg
      hazrats$SD <- rowMeans(hazrats[,c("sdpos","sdneg")], na.rm=T)
      
      #Hazard Ratios
      hr01 <- with(hazrats[hazrats$Group=="< 5%",], exp(rnorm(nrow(hazrats[hazrats$Group=="< 5%",]), hazrats[hazrats$Group=="< 5%","Beta"], hazrats[hazrats$Group=="< 5%","SD"])))
      hr02 <- with(hazrats[hazrats$Group=="5-10%",], exp(rnorm(nrow(hazrats[hazrats$Group=="5-10%",]), hazrats[hazrats$Group=="5-10%","Beta"], hazrats[hazrats$Group=="5-10%","SD"])))
      hr03 <- with(hazrats[hazrats$Group=="10-20%",], exp(rnorm(nrow(hazrats[hazrats$Group=="10-20%",]), hazrats[hazrats$Group=="10-20%","Beta"], hazrats[hazrats$Group=="10-20%","SD"])))
      hr04 <- with(hazrats[hazrats$Group=="20-40%",], exp(rnorm(nrow(hazrats[hazrats$Group=="20-40%",]), hazrats[hazrats$Group=="20-40%","Beta"], hazrats[hazrats$Group=="20-40%","SD"])))
      hr06 <- with(hazrats[hazrats$Group=="60-80%",], exp(rnorm(nrow(hazrats[hazrats$Group=="60-80%",]), hazrats[hazrats$Group=="60-80%","Beta"], hazrats[hazrats$Group=="60-80%","SD"])))
      hr07 <- with(hazrats[hazrats$Group=="80-90%",], exp(rnorm(nrow(hazrats[hazrats$Group=="80-90%",]), hazrats[hazrats$Group=="80-90%","Beta"], hazrats[hazrats$Group=="80-90%","SD"])))
      hr08 <- with(hazrats[hazrats$Group=="90-95%",], exp(rnorm(nrow(hazrats[hazrats$Group=="90-95%",]), hazrats[hazrats$Group=="90-95%","Beta"], hazrats[hazrats$Group=="90-95%","SD"])))
      hr09 <- with(hazrats[hazrats$Group=="95-99%",], exp(rnorm(nrow(hazrats[hazrats$Group=="95-99%",]), hazrats[hazrats$Group=="95-99%","Beta"], hazrats[hazrats$Group=="95-99%","SD"])))
      hr10 <- with(hazrats[hazrats$Group=="> 99%",], exp(rnorm(nrow(hazrats[hazrats$Group=="> 99%",]), hazrats[hazrats$Group=="> 99%","Beta"], hazrats[hazrats$Group=="> 99%","SD"])))
      
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
        
        for(p in 2:nrow(incidence)){
          incidence[[paste0("survival",i)]][p] <- exp(-5*incidence[[paste0("mortandrisk",i)]][p-1])
        }
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
        
        result <- data.frame(incidence$age, Groups[i], incidence[[paste0("lifetimerisk",i)]])
        lifetimerisk <- rbind(lifetimerisk, result)
      }
      
      colnames(lifetimerisk) <- c("Age","Group",paste0("LifetimeRisk",m))
      
      m <- m+1
      
      lifetimerisks <- left_join(lifetimerisks, lifetimerisk)
      
    }
    
    lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2)])
    confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975))
    
    bootstrapped_lifetimerisk <- data.frame(Age=rep(c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"),10),
                                            Group=c(rep("< 5%",16), rep("5-10%",16), rep("10-20%",16), rep("20-40%",16), rep("40-60%",16), rep("60-80%",16), rep("80-90%",16), rep("90-95%",16), rep("95-99%",16), rep("> 99%",16)))
    
    bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
    bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]
    
    #Add in actual lifetime risks
    lifetimeriskactual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_",biobank[j],"_OnePercentTails.csv"), select=c("LifetimeRisk"), data.table=FALSE)
    bootstrapped_lifetimerisk <- cbind(bootstrapped_lifetimerisk, lifetimeriskactual)
    
    #Plot all as well as overall lifetime risk
    bootstrapped_lifetimerisk$Age <- factor(bootstrapped_lifetimerisk$Age, levels=c("1-4 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years","75-79 years"))
    bootstrapped_lifetimerisk$age <- rep(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80),10)
    bootstrapped_lifetimerisk$Group <- factor(bootstrapped_lifetimerisk$Group, levels=c(c("< 5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%")))
    
    write.csv(bootstrapped_lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_Bootstrapped_",biobank[j],"_OnePercentTails.csv"))
    
    #Considering confidence intervals
    riskwithintervals <- subset(bootstrapped_lifetimerisk, Group=="< 5%" | Group=="40-60%" | Group=="> 99%")
    
    ggplot(riskwithintervals, aes(age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
      geom_point() +
      geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
      xlab("Age Range") + 
      ylab("Cumulative Risk (%)") + 
      theme_bw() +
      labs(color='PRS Group', fill='PRS Group') +
      scale_color_hue(labels = c("< 5%", "40-60%", "> 99%")) +
      scale_fill_hue(labels = c("< 5%", "40-60%", "> 99%")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/",biobank[j],"/",directory,"/",hr_phenos,"_",sexes[k],"_LifetimeRisk_Bootstrapped_",biobank[j],"_OnePercentTails.png"), height=10 , width=10)
    
  }
}
