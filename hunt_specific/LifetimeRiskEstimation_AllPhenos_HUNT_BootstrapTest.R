#Libraries
library(data.table)
library(dplyr)
library(ggplot2)

gbd_phenos <- c("Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Diabetes mellitus type 1", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
hr_phenos <- c("J10_ASTHMA", "I9_AF", "C3_BREAST", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")

for(j in 1:length(gbd_phenos)){
  
  lifetimerisks <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                              Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))
  k <- 0 
  
  while(k < 5000){
    
    #Read in GBD incidence data 
    incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Incidence.csv", data.table=FALSE)
    incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
    incidence$val <- as.numeric(incidence$val)
    incidence$upper <- as.numeric(incidence$upper)
    incidence$lower <- as.numeric(incidence$lower)
    
    bcpc_incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
    bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
    bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
    bcpc_incidence$upper <- as.numeric(bcpc_incidence$upper)
    bcpc_incidence$lower <- as.numeric(bcpc_incidence$lower)
    
    incidence <- rbind(incidence, bcpc_incidence)
    
    incidence <- subset(incidence, cause==gbd_phenos[j])
    
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
    
    prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Prevalence.csv", data.table=FALSE)
    prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
    prevalence$val <- as.numeric(prevalence$val)
    prevalence$upper <- as.numeric(prevalence$upper)
    prevalence$lower <- as.numeric(prevalence$lower)
    
    bcpc_prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
    bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location=="Norway")
    bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
    bcpc_prevalence$upper <- as.numeric(bcpc_prevalence$upper)
    bcpc_prevalence$lower <- as.numeric(bcpc_prevalence$lower)
    
    prevalence <- rbind(prevalence, bcpc_prevalence)
    
    prevalence <- subset(prevalence, cause==gbd_phenos[j])
    
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
    mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/GBD_Mortality.csv", data.table=FALSE)
    mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
    
    bcpc_mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)
    bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location=="Norway")
    
    mortality <- rbind(mortality, bcpc_mortality)
    
    mortality <- subset(mortality, cause==gbd_phenos[j] | cause=="All causes")
    
    if(gbd_phenos[j]=="Breast cancer"){
      all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Female")
      all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val","upper","lower")]
    } else if(gbd_phenos[j]=="Prostate cancer"){
      all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Male")
      all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val","upper","lower")]
    } else {
      all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate" & sex=="Both")
      all_cause_mortality <- all_cause_mortality[,c("location","sex","age","val","upper","lower")]
    }
    
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
    
    incidence <- subset(incidence, age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
    
    #######################################################################################################################################################################################
    
    #Hazard Ratios 
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread("/Users/jermy/Documents/INTERVENE/Results/HUNT/HazardRatios_FullSample_HUNT.csv", data.table = FALSE)
    colnames(hazrats) <- c("phenotype", "prs", "group", "controls", "cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
    hazrats <- subset(hazrats, phenotype==hr_phenos[j])
    
    hazrats$beta <- log(hazrats$HR)
    hazrats$beta_pos <- log(hazrats$CIpos)
    
    hazrats$SD <- (hazrats[,"beta_pos"] - hazrats[,"beta"]) / 1.96
    
    #Sample from the hazard ratio distribution
    
    #Hazard Ratios
    hr01 <- exp(rnorm(1, mean=hazrats[1,"beta"], sd=hazrats[1,"SD"]))
    hr02 <- exp(rnorm(1, mean=hazrats[2,"beta"], sd=hazrats[2,"SD"]))
    hr03 <- exp(rnorm(1, mean=hazrats[3,"beta"], sd=hazrats[3,"SD"]))
    hr04 <- exp(rnorm(1, mean=hazrats[4,"beta"], sd=hazrats[4,"SD"]))
    hr05 <- exp(rnorm(1, mean=hazrats[5,"beta"], sd=hazrats[5,"SD"]))
    hr07 <- exp(rnorm(1, mean=hazrats[6,"beta"], sd=hazrats[6,"SD"]))
    hr08 <- exp(rnorm(1, mean=hazrats[7,"beta"], sd=hazrats[7,"SD"]))
    hr09 <- exp(rnorm(1, mean=hazrats[8,"beta"], sd=hazrats[8,"SD"]))
    hr10 <- exp(rnorm(1, mean=hazrats[9,"beta"], sd=hazrats[9,"SD"]))
    hr11 <- exp(rnorm(1, mean=hazrats[10,"beta"], sd=hazrats[10,"SD"]))
    
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
    incidence$i6 <- (incidence$incidence_sample*incidence$population) / ((props06 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population)) + (hr11 * (props11 * incidence$population))) 
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
      incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence_sample)
      
      #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
      incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
      
      #Mortality and risk
      incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate_sample)
      
      #Survival
      incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
      
      result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
      lifetimerisk <- rbind(lifetimerisk, result)
    }
    
    colnames(lifetimerisk) <- c("Age","Group",paste0("LifetimeRisk",k))
    
    k <- k+1
    
    lifetimerisks <- left_join(lifetimerisks, lifetimerisk)
    
  } 
  
  lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2)])
  confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975))
  
  bootstrapped_lifetimerisk <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                                          Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))
  
  bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
  bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]
  
  #Add in actual lifetime risks
  lifetimeriskactual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_LifetimeRisk_HUNT.csv"), select=c("LifetimeRisk"), data.table=FALSE)
  bootstrapped_lifetimerisk <- cbind(bootstrapped_lifetimerisk, lifetimeriskactual)
  
  #Plot all as well as overall lifetime risk
  bootstrapped_lifetimerisk$Age <- factor(bootstrapped_lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  bootstrapped_lifetimerisk$Group <- factor(bootstrapped_lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
  #write.csv(bootstrapped_lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_HUNT.csv"))
  
  #Considering confidence intervals
  riskwithintervals <- subset(bootstrapped_lifetimerisk, Group=="Group1" | Group=="Group6" | Group=="Group11")
  
  ggplot(riskwithintervals, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
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
  #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_HUNT.png"), height=10 , width=10, dpi=1200)
  
}

#########################################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################################################

gbd_phenos <- c()
hr_phenos <- c()

for(j in 1:length(gbd_phenos)){
  
  for(l in c("Male", "Female")){
    
    lifetimerisks <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                                Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))
    k <- 0 
    
    while(k < 5000){
      
      #Read in GBD incidence data 
      incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBD.csv", data.table=FALSE)
      incidence <- subset(incidence, sex==l & cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
      incidence$val <- as.numeric(incidence$val)
      incidence$upper <- as.numeric(incidence$upper)
      incidence$lower <- as.numeric(incidence$lower)
      
      incidence <- subset(incidence, cause==gbd_phenos[j])
      
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
      
      prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBD.csv", data.table=FALSE)
      prevalence <- subset(prevalence, sex==l & cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
      prevalence$val <- as.numeric(prevalence$val)
      prevalence$upper <- as.numeric(prevalence$upper)
      prevalence$lower <- as.numeric(prevalence$lower)
      
      prevalence <- subset(prevalence, cause==gbd_phenos[j])
      
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
      mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBD.csv", data.table=FALSE)
      mortality <- subset(mortality, sex==l & cause!="Breast cancer" & cause!="Prostate cancer" & location=="Norway")
      
      mortality <- subset(mortality, cause==gbd_phenos[j] | cause=="All causes")
      
      all_cause_mortality <- subset(mortality, cause=="All causes" & metric=="Rate")
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
      
      incidence <- subset(incidence, age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
      
      #######################################################################################################################################################################################
      
      #Hazard Ratios 
      
      #Read in the hazard ratios and allocate to variables...
      hazrats <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HUNT/HazardRatios/HazardRatios_",l,"Sample_HUNT.csv"), data.table = FALSE)
      colnames(hazrats) <- c("phenotype", "prs", "group", "controls", "cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
      hazrats <- subset(hazrats, phenotype==hr_phenos[j])
      
      hazrats$beta <- log(hazrats$HR)
      hazrats$beta_pos <- hazrats$beta + 1.96*hazrats$se
      
      hazrats$SD <- (hazrats[,"beta_pos"] - hazrats[,"beta"]) / 1.96
      
      #Sample from the hazard ratio distribution
      
      #Hazard Ratios
      hr01 <- exp(rnorm(1, mean=hazrats[1,"beta"], sd=hazrats[1,"SD"]))
      hr02 <- exp(rnorm(1, mean=hazrats[2,"beta"], sd=hazrats[2,"SD"]))
      hr03 <- exp(rnorm(1, mean=hazrats[3,"beta"], sd=hazrats[3,"SD"]))
      hr04 <- exp(rnorm(1, mean=hazrats[4,"beta"], sd=hazrats[4,"SD"]))
      hr05 <- exp(rnorm(1, mean=hazrats[5,"beta"], sd=hazrats[5,"SD"]))
      hr07 <- exp(rnorm(1, mean=hazrats[6,"beta"], sd=hazrats[6,"SD"]))
      hr08 <- exp(rnorm(1, mean=hazrats[7,"beta"], sd=hazrats[7,"SD"]))
      hr09 <- exp(rnorm(1, mean=hazrats[8,"beta"], sd=hazrats[8,"SD"]))
      hr10 <- exp(rnorm(1, mean=hazrats[9,"beta"], sd=hazrats[9,"SD"]))
      hr11 <- exp(rnorm(1, mean=hazrats[10,"beta"], sd=hazrats[10,"SD"]))
      
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
      incidence$i6 <- (incidence$incidence_sample*incidence$population) / ((props06 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr05 * (props05 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) + (hr10 * (props10 * incidence$population)) + (hr11 * (props11 * incidence$population))) 
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
        incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence_sample)
        
        #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
        incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
        
        #Mortality and risk
        incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate_sample)
        
        #Survival
        incidence[[paste0("survival",i)]] <- exp(-5*incidence[[paste0("mortandrisk",i)]])
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
        
        result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
        lifetimerisk <- rbind(lifetimerisk, result)
      }
      
      colnames(lifetimerisk) <- c("Age","Group",paste0("LifetimeRisk",k))
      
      k <- k+1
      
      lifetimerisks <- left_join(lifetimerisks, lifetimerisk)
      
    } 
    
    lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2)])
    confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975))
    
    bootstrapped_lifetimerisk <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                                            Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))
    
    bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
    bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]
    
    #Add in actual lifetime risks
    lifetimeriskactual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_",l,"_LifetimeRisk_HUNT.csv"), select=c("LifetimeRisk"), data.table=FALSE)
    bootstrapped_lifetimerisk <- cbind(bootstrapped_lifetimerisk, lifetimeriskactual)
    
    #Plot all as well as overall lifetime risk
    bootstrapped_lifetimerisk$Age <- factor(bootstrapped_lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    bootstrapped_lifetimerisk$Group <- factor(bootstrapped_lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
    
    write.csv(bootstrapped_lifetimerisk, paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_",l,"_LifetimeRisk_BootstrappedConfidenceIntervals_HUNT.csv"))
    
    #Considering confidence intervals
    riskwithintervals <- subset(bootstrapped_lifetimerisk, Group=="Group1" | Group=="Group6" | Group=="Group11")
    
    ggplot(riskwithintervals, aes(Age, LifetimeRisk, color=Group, group=Group)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
      geom_point() +
      geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
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
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",hr_phenos[j],"_",l,"LifetimeRisk_BootstrappedConfidenceIntervals_HUNT.png"), height=10 , width=10)
    
  }
}
