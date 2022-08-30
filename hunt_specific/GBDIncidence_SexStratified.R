#Libraries
library(data.table)
library(ggplot2)
library(dplyr)

for(k in c("Male", "Female")){
  
  print(k)
  
  #Incidence data - basic pre-processing
  incidence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Incidence_GBD.csv", data.table=FALSE)
  
  #Incidence data to be replaced with that for males and females for prostate cancer and breast cancer respectively. Came from a separate dataset to reduce size of the full dataset. 
  incidence <- subset(incidence, sex==k & cause!="Prostate cancer" & cause!="Breast cancer")
  incidence <- incidence[,c("location","age","cause","metric","val")]
  
  incidence$val <- as.numeric(incidence$val)
  
  ##Subset to Rate only as no longer require absolute numbers.
  incidence <- subset(incidence, metric=="Rate")
  
  incidence$val <- as.numeric(incidence$val)
  
  #For each cause plot incidence rates
  for(i in unique(incidence$cause)){
    disease <- subset(incidence, cause==i & age!="All Ages")
    disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79","80 to 84","85 to 89","90 to 94","95 plus","All Ages"))
    disease$location <- as.factor(disease$location)
    ggplot(disease, aes(age, val, color=location, group=location)) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 19), se = FALSE) +
      geom_point() +
      xlab("Age Range") + 
      ylab("New cases / 100000") + 
      theme_bw() +
      labs(color='Country') +
      scale_color_hue(labels = c("Estonia", "Finland", "Global", "Norway", "UK", "USA")) +
      theme(title = element_text(size = 22),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle=-90, hjust=0),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/LifetimeRisk/",i,"_",k,"_GBDRates.png"), height=10 , width=10)
  }
  
  #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
  incidence$incidence <- incidence$val / 100000
  incidence <- incidence[,c("location","age","cause","metric","incidence")]
  
  #Prevalence - use to calculate hazard (incidence/(1-prevalence)) - The code is equivalent to that defined for incidence. 
  prevalence <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Prevalence_GBD.csv", data.table=FALSE)
  
  prevalence <- subset(prevalence, sex==k & cause!="Prostate cancer" & cause!="Breast cancer")
  prevalence <- prevalence[,c("location","age","cause","metric","val")]
  
  ##Subset to Rate only as no longer require absolute number
  prevalence <- subset(prevalence, metric=="Rate")
  
  prevalence$val <- as.numeric(prevalence$val)
  
  #Divide by 100,000 to get prevalence in terms of probability. 
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- left_join(incidence, prevalence)
  
  #Calculate hazard
  incidence$hazard <- incidence$incidence / (1 - incidence$prevalence)
  
  #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
  incidence$risk <- 1 - exp(-5*incidence$hazard)
  
  #Subset to relevant columns only
  incidence <- incidence[,c("location","age","cause","incidence","hazard","risk")]
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  
  #Calculate age specific and disease specific mortality
  mortality <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GBD_Data/Sex_Stratified_GBD/Sex_Stratified_Mortality_GBD.csv", data.table=FALSE)
  
  mortality <- mortality[,c("location","sex","age","cause","metric","val")]
  
  mortality <- subset(mortality, sex==k & cause!="Prostate cancer" & cause!="Breast cancer")
  
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
  
  result <- c()
  
  #For each cause within each country:
  for(i in unique(incidence$location)){
    
    for(j in unique(incidence$cause)){
      
      subby <- subset(incidence, cause==j & location==i)
      
      #Calculate the probability of surviving disease free during the age interval.
      subby$mortandrisk <- cumsum(subby$hazard + subby$mortality_rate)
      subby$survival <- exp(-5*subby$mortandrisk)
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      subby$lifetimerisk <- cumsum(subby$survival*subby$risk)*100
      
      #Calculate the probability of remaining disease free during the age interval - death as a competing interest not considered in this section. 
      subby$cumulativerisk <- cumsum(subby$hazard)
      #subby$survival_noncomp <- exp(-5*subby$cumulativerisk)
      
      #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
      #subby$lifetimerisk_noncomp <- cumsum(subby$survival_noncomp*subby$risk)*100
      
      result <- rbind(result, subby)
    }
  }
  
  incidence <- result
  #Plot lifetime risk for each cause - stratified by country
  
  for(i in unique(incidence$cause)){
    disease <- subset(incidence, cause==i)
    disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    disease$location <- as.factor(disease$location)
    print(ggplot(disease, aes(age, lifetimerisk, color=location, group=location)) +
            stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) +
            geom_point() +
            xlab("Age Range") + 
            ylab("Lifetime Risk (%)") + 
            theme_bw() +
            labs(color='Country') +
            scale_color_hue(labels = c("Estonia", "Finland", "Global", "Norway", "UK", "USA")) +
            theme(title = element_text(size = 22),
                  legend.text = element_text(size = 16),
                  legend.title = element_text(size = 18),
                  axis.title.x = element_text(size = 18),
                  axis.text.x = element_text(size = 12, angle=-90, hjust=0),
                  axis.title.y = element_text(size = 18),
                  axis.text.y = element_text(size = 16)))
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/LifetimeRisk/",i,"_",k,"_LifetimeRisk.png"), height=10 , width=10)
  }
}


###############################################################################################################################################################################
###############################################################################################################################################################################