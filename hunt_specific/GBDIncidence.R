#Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(rcartocolor)
library(stringr)

##### put custom path for incidence, prevalence, mortality
output_dir<-"/mnt/work/workbench/bwolford/intervene/results/summary/"

#Incidence data - basic pre-processing
incidence <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/GBD_Incidence.csv",data.table=FALSE)

#Incidence data to be replaced with that for males and females for prostate cancer and breast cancer respectively. Came from a separate dataset to reduce size of the full dataset. 
incidence <- subset(incidence, cause!="Prostate cancer" & cause!="Breast cancer")
incidence <- incidence[,c("location","age","cause","metric","val")]

bc_pc_incidence <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/BreastCancerProstateCancer_Incidence.csv", data.table=FALSE)
bc_pc_incidence <- subset(bc_pc_incidence, (sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer"))
bc_pc_incidence <- bc_pc_incidence[,c("location","age","cause","metric","val")]

incidence <- rbind(incidence, bc_pc_incidence)

incidence$val <- as.numeric(incidence$val)

#For each country and age range calculate the incidence rate for stroke (intracerebral haemmorhage and ischemic stroke combined).
for(j in unique(incidence$location)){
  print(j)
  
  agebins <- c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")
  
  for(i in agebins){
    
    print(i)
    
    #Subset to cause and age 
    stroke <- subset(incidence, location==j & age==i & (cause=="Ischemic stroke" | cause=="Intracerebral hemorrhage"))
    stroke$val <- as.numeric(stroke$val)
    
    #Determine size of population for stroke (average of the population size for ischemic stroke and intracerebral hemorrhage)
    population <- c(stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Ischemic stroke","val"]/100000) , stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Intracerebral hemorrhage","val"]/100000))
    stroke_pop <- (population[1] + population[2]) / 2
    
    #The incidence of stroke is the sum of the incidence for two types of stroke.
    stroke_inc <- stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"] + stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]
    
    #Calculate the incidence rate according to the population and incidence
    stroke_perc <- (stroke_inc / stroke_pop)*100000
    
    result <- c(j,i,"Stroke (exc SAH)", "Rate", stroke_perc)
    
    #Join the estimate to the original incidence dataset
    incidence <- rbind(incidence, result)
    
  }
}

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
  #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/LifetimeRisk/",i,"_GBDRates.png"), height=10 , width=10)
}

#Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
incidence$incidence <- incidence$val / 100000
incidence <- incidence[,c("location","age","cause","metric","incidence")]

#Prevalence - use to calculate hazard (incidence/(1-prevalence)) - The code is equivalent to that defined for incidence. 
prevalence <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/GBD_Prevalence.csv", data.table=FALSE)
prevalence <- subset(prevalence, cause!="Prostate cancer" & cause!="Breast cancer")
prevalence <- prevalence[,c("location","age","cause","metric","val")]

bc_pc_prevalence <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/BreastCancerProstateCancer_Prevalence.csv", data.table=FALSE)
bc_pc_prevalence <- subset(bc_pc_prevalence, (sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer"))
bc_pc_prevalence <- bc_pc_prevalence[,c("location","age","cause","metric","val")]

prevalence <- rbind(prevalence, bc_pc_prevalence)

#Iterate over all countries/regions
for(j in unique(prevalence$location)){
  print(j)
  #Subset to cause name and age name
  agebins <- c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")
  
  for(i in agebins){
    
    print(i)
    
    stroke <- subset(prevalence, location==j & age==i & (cause=="Ischemic stroke" | cause=="Intracerebral hemorrhage"))
    stroke$val <- as.numeric(stroke$val)
    
    #Make population for stroke 
    population <- c(stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Ischemic stroke","val"]/100000) , stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Intracerebral hemorrhage","val"]/100000))
    percent <- c((stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"] / population)*100 , (stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"] / population)*100)
    stroke_pop <- (population[1] + population[2]) / 2
    stroke_inc <- stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"] + stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]
    print(stroke_pop)
    print(stroke_inc)
    stroke_perc <- (stroke_inc / stroke_pop)*100000
    result <- c(j,i,"Stroke (exc SAH)", "Rate", stroke_perc)
    prevalence <- rbind(prevalence, result)
    
  }
}

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
mortality <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/GBD_Mortality.csv", data.table=FALSE)
mortality <- mortality[,c("location","age","cause","metric","val")]

mortality <- subset(mortality, cause!="Prostate cancer" & cause!="Breast cancer")

#Iterate over all countries/regions to define mortality rate for stroke. 
for(j in unique(mortality$location)){
  print(j)
  #Subset to cause name and age name
  agebins <- c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")
  
  for(i in agebins){
    
    print(i)
    
    stroke <- subset(mortality, location==j & age==i & (cause=="Ischemic stroke" | cause=="Intracerebral hemorrhage"))
    stroke$val <- as.numeric(stroke$val)
    
    #Make population for stroke 
    population <- c(stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Ischemic stroke","val"]/100000) , stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]/(stroke[stroke$metric=="Rate" & stroke$cause=="Intracerebral hemorrhage","val"]/100000))
    percent <- c((stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"] / population)*100 , (stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"] / population)*100)
    stroke_pop <- (population[1] + population[2]) / 2
    stroke_inc <- stroke[stroke$metric=="Number" & stroke$cause=="Ischemic stroke","val"] + stroke[stroke$metric=="Number" & stroke$cause=="Intracerebral hemorrhage","val"]
    print(stroke_pop)
    print(stroke_inc)
    stroke_perc <- (stroke_inc / stroke_pop)*100000
    result <- c(j,i,"Stroke (exc SAH)", "Rate", stroke_perc)
    mortality <- rbind(mortality, result)
    
  }
}

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

#Perform the same as above but for sex specific diseases - breast cancer and prostate cancer 
bc_pc_mortality <- fread("/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/BreastCancerProstateCancer_Mortality.csv", data.table=FALSE)

bc_pc_mortality <- subset(bc_pc_mortality, (sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | cause=="All causes")
bc_pc_mortality <- bc_pc_mortality[,c("location","sex","age","cause","metric","val")]

all_cause_bc_pc_mortality <- subset(bc_pc_mortality, cause=="All causes" & metric=="Rate")
all_cause_bc_pc_mortality <- all_cause_bc_pc_mortality[,c("location","sex","age","val")]
colnames(all_cause_bc_pc_mortality)[4] <- c("all_cause_rate")

cause_specific_bc_pc_mortality <- subset(bc_pc_mortality, cause!="All causes" & metric=="Rate")
cause_specific_bc_pc_mortality <- cause_specific_bc_pc_mortality[,c("location","sex","age","cause","val")]
colnames(cause_specific_bc_pc_mortality)[5] <- c("cause_specific_rate")

bc_pc_mortality <- left_join(cause_specific_bc_pc_mortality, all_cause_bc_pc_mortality)
bc_pc_mortality$all_cause_rate <- as.numeric(bc_pc_mortality$all_cause_rate)
bc_pc_mortality$cause_specific_rate <- as.numeric(bc_pc_mortality$cause_specific_rate)

bc_pc_mortality$mortality_rate <- (bc_pc_mortality$all_cause_rate - bc_pc_mortality$cause_specific_rate)/100000
bc_pc_mortality <- bc_pc_mortality[,c("location","age","cause","mortality_rate")]

#Join breast cancer and prostate cancer to the original mortality dataset
mortality <- rbind(mortality, bc_pc_mortality)

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
#removed US, customized colors
for(i in unique(incidence$cause)){
  disease<-incidence %>% filter(cause==i,location!="United States of America")
  my_colors = carto_pal(n=length(unique(disease$location)), name="Safe")
  disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  disease$location <- as.factor(disease$location)
  
  pdf(file=paste0(output_dir,i,"_LifetimeRisk.pdf"), height=10 , width=10,useDingbats=TRUE)
  print(ggplot(disease, aes(age, lifetimerisk, color=location, group=location)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) +
    geom_point(size=5,alpha=0.7) +
    xlab("Age Range") + 
    ylab("Baseline Lifetime Risk (%)") + 
    theme_bw() +
    labs(color='Country',title=i) + 
    guides(color=guide_legend(reverse = TRUE)) +
    scale_color_manual(values=my_colors,labels = c("Estonia", "Finland", "Global", "Norway", "UK")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 18, angle=45, hjust=1),
          axis.title.y = element_text(size = 22),
          axis.text.y = element_text(size = 16)))
  dev.off()
}

prscols <- c("Asthma","AllCancers","Appendicitis", "Atrial_Fibrillation", "Breast_Cancer", "CHD","Colorectal_Cancer", "Epilepsy","Gout",
             "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer","Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", 
             "T1D","T2D", "ILD", "Lung_Cancer")

flagship_20<-c("Asthma","Breast cancer", "Diabetes mellitus type 1","Osteoarthritis hip","Atrial fibrillation and flutter","Subarachnoid hemorrhage","Rheumatoid arthritis",
"Gout","Major depressive disorder","Prostate cancer","Appendicitis","Ischemic heart disease","Colon and rectum cancer","Idiopathic epilepsy","Total cancers",
"Osteoarthritis knee","Malignant skin melanoma","Diabetes mellitus type 2","Interstitial lung disease and pulmonary sarcoidosis","Tracheal, bronchus, and lung cancer")

  disease<-incidence %>% filter(location!="United States of America") %>% 
    mutate(pretty_varname = as.factor(str_wrap(cause,width=30)))  %>% filter(cause %in% flagship_20)
  my_colors = carto_pal(n=length(unique(disease$location)), name="Safe")
  disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  disease$location <- as.factor(disease$location)
  
  pdf(file=paste0(output_dir,"LifetimeRisk_facet.pdf"),height=6,width=12)
  ggplot(disease, aes(age, lifetimerisk, color=location, group=location)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) + facet_wrap(~pretty_varname,nrow=4,scales="free_y") +
    geom_point(size=1,alpha=0.7) +
    xlab("Age Range") + 
    ylab("Baseline Lifetime Risk (%)") + 
    theme_bw() +
    labs(color='Country') +
    guides(color=guide_legend(reverse = TRUE)) +
    scale_color_manual(values=my_colors,labels = c("Estonia", "Finland", "Global", "Norway", "UK")) +
    theme(title = element_text(size = 12),
          strip.background = element_rect(color="black", fill="white"),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 8, angle=45, hjust=1),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  dev.off()

  select<-disease%>%filter(cause=="Gout"|cause=="Ischemic heart disease" |cause=="Diabetes mellitus type 2" | cause=="Prostate cancer")
                            
  pdf(file=paste0(output_dir,"LifetimeRisk_facet_select.pdf"),height=4,width=12)
  ggplot(select, aes(age, lifetimerisk, color=location, group=location)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) + facet_wrap(~pretty_varname,nrow=1,scales="free_y") +
    geom_point(size=3,alpha=0.7) +
    xlab("Age Range") + 
    ylab("Baseline Lifetime Risk (%)") + 
    theme_bw() +
    labs(color='Country') +
    guides(color=guide_legend(reverse = TRUE)) +
    scale_color_manual(values=my_colors,labels = c("Estonia", "Finland", "Global", "Norway", "UK")) +
    theme(title = element_text(size = 12),
          strip.background = element_rect(color="black", fill="white"),
          strip.text.x = element_text(size=14,margin = margin(.1, 0, .1, 0, "cm")),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 8, angle=45, hjust=1),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  dev.off()
  

# See how the results differ when not considering mortality as a competing interest. Plot lifetime risk for each cause - stratified by country 
#for(i in unique(incidence$cause)){
#  disease <- subset(incidence, cause==i)
#  disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
#  disease$location <- as.factor(disease$location)
#  ggplot(disease, aes(age, lifetimerisk_noncomp, color=location, group=location)) +
#          stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) +
#          geom_point() +
#          xlab("Age Range") + 
#          ylab("Lifetime Risk (%)") + 
#          theme_bw() +
#          labs(color='Country') +
#          scale_color_hue(labels = c("Estonia", "Finland", "Global", "Norway", "UK", "USA")) +
#          theme(title = element_text(size = 22),
#                legend.text = element_text(size = 16),
#                legend.title = element_text(size = 18),
#                axis.title.x = element_text(size = 18),
#                axis.text.x = element_text(size = 12, angle=-90, hjust=0),
#                axis.title.y = element_text(size = 18),
#                axis.text.y = element_text(size = 16))
#  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/LifetimeRisk/Archive/",i,"_LifetimeRisk_NonCompete.png"), height=10 , width=10)
#}
