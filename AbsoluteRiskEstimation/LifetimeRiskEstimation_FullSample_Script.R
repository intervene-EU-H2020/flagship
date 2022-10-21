
#Libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(optparse)

#options(warn = -1)

print(R.version)
set.seed(1234)

option_list <- list(
  make_option("--country", type="character", default=""),
  make_option("--path",type="character",default=""),
  make_option("--output_dir",type="character",default=""),
  make_option("--biobank",type="character",default=""),
  make_option("--k",type="numeric",default=""),
  make_option("--hr",type="character",default=""),
  make_option("--ancestry",type="character",default="EUR")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script creates a key for old reference genome coords to new reference genome coordinates after liftover\n")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

country<-opt$country
path<-opt$path #must have final backslash
full_HR_path<-opt$hr  #must go straight to csv
#   * Line 215 - Enter file path and biobank name
#   * Line 233 - Enter file path and biobank name
output_dir<-opt$output_dir
biobank<-opt$biobank
ancestry<-opt$ancestry
nk<-opt$k

gbd_phenos <- c("Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer", "Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
hr_phenos <- c("ILD", "C3_BRONCHUS_LUNG","C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")

################# Actual Lifetime Risk Estimates ###########
for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  
  #Read in GBD incidence data 
  incidence <- fread(paste0(path,"GBD_Incidence.csv"), data.table=FALSE)
  incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
  incidence$val <- as.numeric(incidence$val)
  
  bcpc_incidence <- fread(paste0(path,"BreastCancerProstateCancer_Incidence.csv"), data.table=FALSE)
  bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location==country)
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
  
  prevalence <- fread(paste0(path,"/GBD_Prevalence.csv"), data.table=FALSE)
  prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
  prevalence$val <- as.numeric(prevalence$val)
  
  bcpc_prevalence <- fread(paste0(path,"/BreastCancerProstateCancer_Prevalence.csv"), data.table=FALSE)
  bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location==country)
  bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
  
  prevalence <- rbind(prevalence, bcpc_prevalence)
  
  prevalence <- subset(prevalence, cause==gbd_phenos[j])
  
  #Divide prevalence rates by 100,000 to get the prevalence as a probability (note: prevalence rates are per year)
  prevalence$prevalence <- prevalence$val / 100000
  prevalence <- prevalence[,c("location","age","cause","metric","prevalence")]
  
  #Left join to incidence to calculate hazard 
  incidence <- suppressMessages(left_join(incidence, prevalence))
  
  #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
  mortality <- fread(paste0(path,"GBD_Mortality.csv"), data.table=FALSE)
  mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
  
  bcpc_mortality <- fread(paste0(path,"BreastCancerProstateCancer_Mortality.csv"), data.table=FALSE)
  bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location==country)
  
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
  
  mortality <- suppressMessages(left_join(cause_specific_mortality, all_cause_mortality))
  mortality$all_cause_rate <- as.numeric(mortality$all_cause_rate)
  mortality$cause_specific_rate <- as.numeric(mortality$cause_specific_rate)
  
  mortality$mortality_rate <- (mortality$all_cause_rate - mortality$cause_specific_rate)/100000
  mortality <- mortality[,c("location","age","cause","mortality_rate")]
  
  #Merge mortality data to incidence data
  incidence <- suppressMessages(left_join(incidence, mortality))
  
  incidence$age<-trimws(gsub("years","",gsub("-"," to ",incidence$age)))
  incidence <- subset(incidence, age!="All ages" & age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
  incidence<-incidence %>% arrange(factor(age,levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")))
  
  #######################################################################################################################################################################################
  
  #Hazard Ratios
  
  #Read in the hazard ratios and allocate to variables...
  hazrats <- fread(full_HR_path, data.table = FALSE)
  colnames(hazrats) <- c("phenotype", "prs", "group","controls","cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
  hazrats <- subset(hazrats, phenotype==hr_phenos[j])
  if(nrow(hazrats)<1){
    next
  }
  
  #use top and bottom 5% instead
  if (is.na(hazrats[hazrats$group=="< 1%",]$HR) | is.na(hazrats[hazrats$group=="> 99%",]$HR)){
    #Hazard Ratios
    hr01 <- hazrats[hazrats$group=="< 5%",]$HR
    hr02 <- hazrats[hazrats$group=="5-10%",]$HR
    hr03 <- hazrats[hazrats$group=="10-20%",]$HR
    hr04 <- hazrats[hazrats$group=="20-40%",]$HR
    hr06 <- hazrats[hazrats$group=="60-80%",]$HR
    hr07 <- hazrats[hazrats$group=="80-90%",]$HR
    hr08 <- hazrats[hazrats$group=="90-95%",]$HR
    hr09 <- hazrats[hazrats$group=="> 95%",]$HR
    
    #Proportions - 0.2 by definition of PRS group
    props01 <- 0.04
    props02 <- 0.05
    props03 <- 0.1
    props04 <- 0.1
    props05 <- 0.2
    props06 <- 0.2
    props07 <- 0.1
    props08 <- 0.05
    props09 <- 0.04
    
    #Estimate incidence attributable to different distributions of PRS, now group 5 is reference
    incidence$i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) )
    incidence$i5[is.na(incidence$i5)] <- 0
    incidence$i1 <- incidence$i5 * hr01
    incidence$i2 <- incidence$i5 * hr02
    incidence$i3 <- incidence$i5 * hr03
    incidence$i4 <- incidence$i5 * hr04
    incidence$i6 <- incidence$i5 * hr06
    incidence$i7 <- incidence$i5 * hr07
    incidence$i8 <- incidence$i5 * hr08
    incidence$i9 <- incidence$i5 * hr09
    
    groups<-9
  }else {
    #Hazard Ratios
    hr01 <- hazrats[hazrats$group=="< 1%",]$HR
    hr02 <- hazrats[hazrats$group=="1-5%",]$HR
    hr03 <- hazrats[hazrats$group=="5-10%",]$HR
    hr04 <- hazrats[hazrats$group=="10-20%",]$HR
    hr05 <- hazrats[hazrats$group=="20-40%",]$HR
    hr07 <- hazrats[hazrats$group=="60-80%",]$HR
    hr08 <- hazrats[hazrats$group=="80-90%",]$HR
    hr09 <- hazrats[hazrats$group=="90-95%",]$HR
    hr10 <- hazrats[hazrats$group=="95-99%",]$HR
    hr11 <- hazrats[hazrats$group=="> 99%",]$HR
    
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
    
    groups<-11
  }
  ###################################################
  
  lifetimerisk <- data.frame(NULL)
  
  for(i in 1:groups){
    print(i)
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
  #lifetimerisk$Age<-trimws(gsub("years","",gsub("-"," to ",lifetimerisk$Age)))
  #lifetimerisk <- subset(lifetimerisk, Age!="All ages" & Age!="80 to 84" & Age!="85 to 89" & Age!="90 to 94" & Age!="95 plus")
  
  
  #Plot all as well as overall lifetime risk
  lifetimerisk$Age <- factor(lifetimerisk$Age,levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  lifetimerisk$Group <- factor(lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
  write.csv(lifetimerisk, paste0(output_dir,hr_phenos[j],"_LifetimeRisk_",biobank,"_",ancestry,".csv"),row.names=FALSE)
  
  #Not considering confidence intervals
  ggplot(lifetimerisk, aes(Age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    xlab("Age Range") + 
    ylab("Absolute Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group') +
    scale_color_brewer(palette="RdYlBu",direction=-1,guide = guide_legend(reverse = TRUE), labels = c("0-1%", "1-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-90%", "90-95%", "95-99%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0(output_dir,hr_phenos[j],"_LifetimeRisk_",biobank,"_",ancestry,".png"), height=10 , width=10)
  
}


################ BOOOT STRAP ####################

print("bootstrap")

for(j in 1:length(gbd_phenos)){
  
  print(gbd_phenos[j])
  
  lifetimerisks <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                              Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))
  k <- 0 
  
  while(k < nk){
    
    #Read in GBD incidence data 
    incidence <- fread(paste0(path,"GBD_Incidence.csv"), data.table=FALSE)
    incidence <- subset(incidence, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
    incidence$val <- as.numeric(incidence$val)
    incidence$upper<-as.numeric(incidence$upper)
    incidence$lower<-as.numeric(incidence$lower)
    
    bcpc_incidence <- fread(paste0(path,"BreastCancerProstateCancer_Incidence.csv"), data.table=FALSE)
    bcpc_incidence <- subset(bcpc_incidence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location==country)
    bcpc_incidence$val <- as.numeric(bcpc_incidence$val)
    bcpc_incidence$upper<-as.numeric(bcpc_incidence$upper)
    bcpc_incidence$lower<-as.numeric(bcpc_incidence$lower)
    
    incidence <- rbind(incidence, bcpc_incidence)
    
    incidence <- subset(incidence, cause==gbd_phenos[j])
    
    #Divide incidence rates by 100,000 to get the incidence as a probability (note: incidence rates are per year)
    incidence$val <- incidence$val / 100000
    incidence$upper<-incidence$upper/100000
    incidence$lower<-incidence$lower/100000
    
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
    
    #### read in GBD prevalence data
    prevalence <- fread(paste0(path,"GBD_Prevalence.csv"), data.table=FALSE)
    prevalence <- subset(prevalence, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
    prevalence$val <- as.numeric(prevalence$val)
    prevalence$upper<-as.numeric(prevalence$upper)
    prevalence$lower<-as.numeric(prevalence$lower)
    
    bcpc_prevalence <- fread(paste0(path,"BreastCancerProstateCancer_Prevalence.csv"), data.table=FALSE)
    bcpc_prevalence <- subset(bcpc_prevalence, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer")) & location==country)
    bcpc_prevalence$val <- as.numeric(bcpc_prevalence$val)
    bcpc_prevalence$upper <-as.numeric(bcpc_prevalence$upper)
    bcpc_prevalence$lower <-as.numeric(bcpc_prevalence$lower)
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
    incidence <- suppressMessages(left_join(incidence, prevalence))
    
    #Use all cause and cause specific mortality incidence rates to calculate the competing risk of death during the age interval
    mortality <- fread(paste0(path,"GBD_Mortality.csv"), data.table=FALSE)
    mortality <- subset(mortality, cause!="Breast cancer" & cause!="Prostate cancer" & location==country)
    
    bcpc_mortality <- fread(paste0(path,"BreastCancerProstateCancer_Mortality.csv"), data.table=FALSE)
    bcpc_mortality <- subset(bcpc_mortality, ((sex=="Male" & cause=="Prostate cancer") | (sex=="Female" & cause=="Breast cancer") | ((sex=="Male" | sex=="Female") & cause=="All causes")) & location==country)
    
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
    incidence <- suppressMessages(left_join(incidence, mortality))
    
    incidence$age<-trimws(gsub("years","",gsub("-"," to ",incidence$age)))
    incidence <- subset(incidence, age!="All ages" & age!="All Ages" & age!="80 to 84" & age!="85 to 89" & age!="90 to 94" & age!="95 plus")
    incidence<-incidence %>% arrange(factor(age,levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")))
    
    #######################################################################################################################################################################################
    #Hazard Ratios
    
    #Read in the hazard ratios and allocate to variables...
    hazrats <- fread(full_HR_path, data.table = FALSE)
    colnames(hazrats) <- c("phenotype","prs", "group","controls","cases", "beta", "se", "pval", "HR", "CIpos", "CIneg")
    hazrats <- subset(hazrats, phenotype==hr_phenos[j])
    if(nrow(hazrats)<1) { #if hazard ratios don't exist for trait
     break
    }
    hazrats$beta <- log(hazrats$HR)
    hazrats$beta_pos <- log(hazrats$CIpos)
    hazrats$SD <- (hazrats[,"beta_pos"] - hazrats[,"beta"]) / 1.96
    if (sum(complete.cases(hazrats))<nrow(hazrats)) { #if there are NAs then skip trait
      print("test")
      next
    }
    if(is.na(hazrats[hazrats$group=="< 1%",]$HR) | is.na(hazrats[hazrats$group=="> 99%",]$HR) & is.finite(hazrats[hazrats$group=="< 5%",]$beta_pos)){
      #Hazard Ratios, #Sample from the hazard ratio distribution
      hr01 <- exp(rnorm(1, mean=hazrats[hazrats$group=="< 5%",]$beta, sd=hazrats[hazrats$group=="< 5%",]$SD))
      hr02 <- exp(rnorm(1, mean=hazrats[hazrats$group=="5-10%",]$beta, sd=hazrats[hazrats$group=="5-10%",]$SD))
      hr03 <- exp(rnorm(1, mean=hazrats[hazrats$group=="10-20%",]$beta, sd=hazrats[hazrats$group=="10-20%",]$SD))
      hr04 <- exp(rnorm(1, mean=hazrats[hazrats$group=="20-40%",]$beta, sd=hazrats[hazrats$group=="20-40%",]$SD))
      hr06 <- exp(rnorm(1, mean=hazrats[hazrats$group=="60-80%",]$beta, sd=hazrats[hazrats$group=="60-80%",]$SD))
      hr07 <- exp(rnorm(1, mean=hazrats[hazrats$group=="80-90%",]$beta, sd=hazrats[hazrats$group=="80-90%",]$SD))
      hr08 <- exp(rnorm(1, mean=hazrats[hazrats$group=="90-95%",]$beta, sd=hazrats[hazrats$group=="90-95%",]$SD))
      hr09 <- exp(rnorm(1, mean=hazrats[hazrats$group=="> 95%",]$beta, sd=hazrats[hazrats$group=="> 95%",]$SD))
  
      #Proportions - 0.2 by definition of PRS group
      props01 <- 0.04
      props02 <- 0.05
      props03 <- 0.1
      props04 <- 0.1
      props05 <- 0.2
      props06 <- 0.2
      props07 <- 0.1
      props08 <- 0.05
      props09 <- 0.04
      
      #Estimate incidence attributable to different distributions of PRS, now group 5 is reference
      incidence$i5 <- (incidence$incidence*incidence$population) / ((props05 * incidence$population) + (hr01 * (props01 * incidence$population)) + (hr02 * (props02 * incidence$population)) + (hr03 * (props03 * incidence$population)) + (hr04 * (props04 * incidence$population)) + (hr06 * (props06 * incidence$population)) + (hr07 * (props07 * incidence$population)) + (hr08 * (props08 * incidence$population)) + (hr09 * (props09 * incidence$population)) )
      incidence$i5[is.na(incidence$i5)] <- 0
      incidence$i1 <- incidence$i5 * hr01
      incidence$i2 <- incidence$i5 * hr02
      incidence$i3 <- incidence$i5 * hr03
      incidence$i4 <- incidence$i5 * hr04
      incidence$i6 <- incidence$i5 * hr06
      incidence$i7 <- incidence$i5 * hr07
      incidence$i8 <- incidence$i5 * hr08
      incidence$i9 <- incidence$i5 * hr09
      
      no_groups<-9
    } else {
      #Hazard Ratios
      hr01 <- exp(rnorm(1, mean=hazrats[hazrats$group=="< 1%",]$beta, sd=hazrats[hazrats$group=="< 1%",]$SD))
      hr02 <- exp(rnorm(1, mean=hazrats[hazrats$group=="1-5%",]$beta, sd=hazrats[hazrats$group=="1-5%",]$SD))
      hr03 <- exp(rnorm(1, mean=hazrats[hazrats$group=="5-10%",]$beta, sd=hazrats[hazrats$group=="5-10%",]$SD))
      hr04 <- exp(rnorm(1, mean=hazrats[hazrats$group=="10-20%",]$beta, sd=hazrats[hazrats$group=="10-20%",]$SD))
      hr05 <- exp(rnorm(1, mean=hazrats[hazrats$group=="20-40%",]$beta, sd=hazrats[hazrats$group=="20-40%",]$SD))
      hr07 <- exp(rnorm(1, mean=hazrats[hazrats$group=="60-80%",]$beta, sd=hazrats[hazrats$group=="60-80%",]$SD))
      hr08 <- exp(rnorm(1, mean=hazrats[hazrats$group=="80-90%",]$beta, sd=hazrats[hazrats$group=="80-90%",]$SD))
      hr09 <- exp(rnorm(1, mean=hazrats[hazrats$group=="90-95%",]$beta, sd=hazrats[hazrats$group=="90-95%",]$SD))
      hr10 <- exp(rnorm(1, mean=hazrats[hazrats$group=="95-99%",]$beta, sd=hazrats[hazrats$group=="95-99%",]$SD))
      hr11 <- exp(rnorm(1, mean=hazrats[hazrats$group=="> 99%",]$beta, sd= hazrats[hazrats$group=="> 99%",]$SD))
      
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
      
      no_groups<-11
    } 

    ###################################################
    if (exists("no_groups")){
      lifetimerisk <- data.frame(NULL)
      for(i in 1:no_groups){
        #Calculate hazard
        incidence[[paste0("hazard",i)]] <- incidence[[paste0("i",i)]] / (1 - incidence$prevalence_sample)
        
        #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
        incidence[[paste0("risk",i)]] <- 1 - exp(-5*incidence[[paste0("hazard",i)]])
        
        #Mortality and risk
        incidence[[paste0("mortandrisk",i)]] <- cumsum(incidence[[paste0("hazard",i)]] + incidence$mortality_rate_sample)
        
        #Survival
        incidence[[paste0("survival",i)]] <- 1
        
        for(r in 2:nrow(incidence)){
          incidence[[paste0("survival",i)]][r] <- exp(-5*incidence[[paste0("mortandrisk",i)]][r-1])
        }
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        incidence[[paste0("lifetimerisk",i)]] <- cumsum(incidence[[paste0("survival",i)]]*incidence[[paste0("risk",i)]])*100
        
        result <- data.frame(incidence$age, paste0("Group", i), incidence[[paste0("lifetimerisk",i)]])
        lifetimerisk <- rbind(lifetimerisk, result)
      }
      colnames(lifetimerisk) <- c("Age","Group",paste0("LifetimeRisk",k))
      k<-k+1
      lifetimerisks <- suppressMessages(left_join(lifetimerisks, lifetimerisk))
    }
  }
  
  lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2)])
  confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975), na.rm=TRUE)

  bootstrapped_lifetimerisk <- data.frame(Age=rep(c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"),11),
                                      Group=c(rep("Group1",16), rep("Group2",16), rep("Group3",16), rep("Group4",16), rep("Group5",16), rep("Group6",16), rep("Group7",16), rep("Group8",16), rep("Group9",16), rep("Group10",16), rep("Group11",16)))

  bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
  bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]

  #Add in actual lifetime risks (written earlier)
  if(file.exists(paste0(output_dir,hr_phenos[j],"_LifetimeRisk_",biobank,"_",ancestry,".csv"))){
    lifetimeriskactual <- fread(paste0(output_dir,hr_phenos[j],"_LifetimeRisk_",biobank,"_",ancestry,".csv"), select=c("LifetimeRisk"), data.table=FALSE)
    bootstrapped_lifetimerisk <- cbind(bootstrapped_lifetimerisk, lifetimeriskactual)

    #Plot all as well as overall lifetime risk
    bootstrapped_lifetimerisk$Age <- factor(bootstrapped_lifetimerisk$Age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
    bootstrapped_lifetimerisk$Group <- factor(bootstrapped_lifetimerisk$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))
  
    write.csv(bootstrapped_lifetimerisk, paste0(output_dir,hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_",biobank,"_",ancestry,".csv"),row.names=FALSE)
  
    #Considering confidence intervals
    riskwithintervals <- subset(bootstrapped_lifetimerisk, Group=="Group1" | Group=="Group6" | Group=="Group11")
  
    #colors<-c("light green","dark green","plum2","orchid4")
    #colors<-brewer.pal(3,"RdYlBu")
    colors<-c(brewer.pal(11,"RdYlBu")[11],"dark grey",brewer.pal(11,"RdYlBu")[1])
    ggplot(riskwithintervals, aes(Age, LifetimeRisk, fill=Group, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age Range") + 
    ylab("Absolute Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PRS Group', fill='PRS Group') +
    scale_color_manual(values=colors,guide = guide_legend(reverse = TRUE),labels = c("0-1%", "40-60%", "99-100%")) +
    scale_fill_manual(values=colors,guide = guide_legend(reverse = TRUE),labels = c("0-1%", "40-60%", "99-100%")) +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12, angle=-90, hjust=0),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16))
  ggsave(paste0(output_dir,hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_",biobank,"_",ancestry,".png"), height=10 , width=10, dpi=1200)
  }
}