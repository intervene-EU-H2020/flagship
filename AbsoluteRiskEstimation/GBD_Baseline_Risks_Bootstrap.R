#Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(rcartocolor)
library(stringr)

##### put custom path for incidence, prevalence, mortality
#output_dir<-"/mnt/work/workbench/bwolford/intervene/results/summary/"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/"
path<-"/mnt/scratch/brooke/flagship/AbsoluteRiskEstimation/"

gbd_phenos <- c("Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer", "Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
hr_phenos <- c("ILD", "C3_BRONCHUS_LUNG","C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
nk<-3

gbd_phenos<-c("I9_CHD","GOUT","C3_PROSTATE","T2D")
hr_phenos<-c("Ischemic heart disease","Gout","Prostate cancer","Diabetes mellitus type 2")
countries<-c("Massachusetts","Norway","Estonia","Finland","England","United Kingdom","Japan","Scotland","United States of America","Global")
#countries<-c("Global")
result<-c()
point_estimates<-c()
for(j in 1:length(gbd_phenos)){
  for (c in 1:length(countries)){
      lifetimerisks <- data.frame(age=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
      k <- 1 
      while(k < nk){
        country<-countries[c]
        print(country)
        print(gbd_phenos[j])
        
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
          poppy <- subby$val[1]/(subby$val[2]/100000) #number / (rate /100000)
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
        
        #save point estimates
        if(k==1){
          point<-incidence %>% select("location","age","cause","val","upper","lower") %>% rename("incidence"=val,"incidence_upper"=upper,"incidence_lower"=lower)
        }
        
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
        
        if(k==1){
          point<-prevalence %>% select("location","age","cause","val","upper","lower") %>% rename("prevalence"=val,"prevalence_upper"=upper,"prevalence_lower"=lower) %>% left_join(point)
        }
        
        prevalence <- prevalence[,c("location","age","cause","metric","prevalence_sample")]
        
        #Left join to incidence to calculate hazard 
        incidence <- suppressMessages(left_join(incidence, prevalence))
        
        #Calculate hazard
        incidence$hazard_sample <- incidence$incidence_sample / (1 - incidence$prevalence_sample)
        
        #Calculate probability of experiencing the endpoint within the age interval. hazard multiplied by 5 as that is the age interval and current probabilities are per year. 
        incidence$risk_sample <- 1 - exp(-5*incidence$hazard_sample)
        
        #Subset to relevant columns only
        incidence <- incidence[,c("location","age","cause","incidence_sample","prevalence_sample","hazard_sample","risk_sample")]
        
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
        
        if (k==1){
          point<-cause_specific_mortality %>% select("location","age","cause","val","upper","lower") %>% rename("cause_specific_mort"=val,"cause_specific_mort_upper"=upper,"cause_specific_mort_lower"=lower) %>% 
            left_join(point)
          point<-all_cause_mortality %>% select("location","age","val","upper","lower") %>% rename("all_cause_mort"=val,"all_cause_mort_upper"=upper,"all_cause_mort_lower"=lower) %>%
            left_join(point)
        }
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
        
        #Calculate the probability of surviving disease free during the age interval.
        incidence$mortandrisk_sample <- cumsum(incidence$hazard_sample + incidence$mortality_rate_sample)
        incidence$survival_sample <- exp(-5*incidence$mortandrisk_sample)
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        incidence$lifetimerisk_sample <- cumsum(incidence$survival_sample*incidence$risk_sample)*100
        
        #Calculate the probability of remaining disease free during the age interval - death as a competing interest not considered in this section. 
        incidence$cumulativerisk_sample <- cumsum(incidence$hazard_sample)
        #subby$survival_noncomp <- exp(-5*subby$cumulativerisk)
        
        #Calculate lifetime risk as the cumulative sum of the product of survival and risk.
        #subby$lifetimerisk_noncomp <- cumsum(subby$survival_noncomp*subby$risk)*100
        
        if(k==1){
          point_estimates<-rbind(point_estimates,point)
        }
        lifetimerisk<-incidence %>% select(location,age,cause,risk_sample) %>% rename(!!paste0("LifetimeRisk",k):=risk_sample)
        k<-k+1
        lifetimerisks <- suppressMessages(right_join(lifetimerisks, lifetimerisk))
      }
    lifetimerisk_percentile <- as.matrix(lifetimerisks[,-c(1,2,3)])
    confidenceintervals <- apply(lifetimerisk_percentile, 1, quantile, c(0.025, 0.975), na.rm=TRUE)
    bootstrapped_lifetimerisk <- point %>% filter(cause==gbd_phenos[j] & location==countries[c] & age!=c("80-84 years","85-89 years","90-94 years","All ages"))
  
    bootstrapped_lifetimerisk$CIneg <- confidenceintervals[1,]
    bootstrapped_lifetimerisk$CIpos <- confidenceintervals[2,]
    bootstrapped_lifetimerisk$location<-countries[c]
    bootstrapped_lifetimerisk$cause<-gbd_phenos[j]
    
    result<-rbind(result,bootstrapped_lifetimerisk)
  }
}
 
fwrite(result,paste0(output_dir,"GBD_LifetimeRisks.tab"),quote=FALSE,row.names=FALSE,sep="\t")


############################# PLOTTING #############################
print("plotting")
incidence<-fread(paste0(output_dir,"GBD_LifetimeRisks.tab"),fill=TRUE) %>% filter(location!="United States of America"&location!="Global") 
#Plot lifetime risk for each cause - stratified by country
#removed US, customized colors
for(i in unique(incidence$cause)){
  disease<-incidence %>% filter(cause==i)
  my_colors = carto_pal(n=length(unique(disease$location)), name="Safe")
  disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
  disease$location <- as.factor(disease$location)
  
  pdf(file=paste0(output_dir,i,"_LifetimeRisk.pdf"), height=10 , width=10,useDingbats=TRUE)
  print(ggplot(disease, aes(age, lifetimerisk, color=location, group=location)) +
          stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) +
          geom_point(size=5,alpha=0.7) +
          xlab("Age Range") + 
          ylab("Baseline Cumulative Risk (%)") + 
          theme_bw() +
          labs(color='Location',title=i) + 
          guides(color=guide_legend(reverse = TRUE)) +
          scale_color_manual(values=my_colors) +
          theme(title = element_text(size = 22),
                legend.position="bottom",
                legend.text = element_text(size = 28),
                legend.title = element_text(size = 28),
                axis.title.x = element_text(size = 28),
                axis.text.x = element_text(size = 18, angle=45, hjust=1),
                axis.title.y = element_text(size = 28),
                axis.text.y = element_text(size = 28)))
  dev.off()
}

prscols <- c("Asthma","AllCancers","Appendicitis", "Atrial_Fibrillation", "Breast_Cancer", "CHD","Colorectal_Cancer", "Epilepsy","Gout",
            "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer","Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", 
          "T1D","T2D", "ILD", "Lung_Cancer")

flagship_20<-c("Asthma","Breast cancer", "Diabetes mellitus type 1","Osteoarthritis hip","Atrial fibrillation and flutter","Subarachnoid hemorrhage","Rheumatoid arthritis",
            "Gout","Major depressive disorder","Prostate cancer","Appendicitis","Ischemic heart disease","Colon and rectum cancer","Idiopathic epilepsy","Total cancers",
           "Osteoarthritis knee","Malignant skin melanoma","Diabetes mellitus type 2","Interstitial lung disease and pulmonary sarcoidosis","Tracheal, bronchus, and lung cancer")

disease<-incidence %>% filter(location!="United States of America"&location!="Global") %>% 
  mutate(pretty_varname = as.factor(str_wrap(cause,width=30)))  %>% filter(cause %in% flagship_20)
my_colors = carto_pal(n=length(unique(disease$location)), name="Safe")
disease$age <- factor(disease$age, levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
disease$location <- as.factor(disease$location)

pdf(file=paste0(output_dir,"LifetimeRisk_facet.pdf"),height=6,width=12,useDingbats=TRUE)
ggplot(disease, aes(age, lifetimerisk, color=location, group=location)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) + facet_wrap(~pretty_varname,nrow=4,scales="free_y") +
  geom_point(size=1,alpha=0.7) +
  xlab("Age Range") + 
  ylab("Baseline Cumulative Risk (%)") + 
  theme_bw() +
  labs(color='Location') +
  guides(color=guide_legend(reverse = TRUE)) +
  scale_color_manual(values=my_colors) +
  theme(title = element_text(size = 12),legend.position="bottom",
        strip.background = element_rect(color="black", fill="white"),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 8, angle=45, hjust=1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))
dev.off()

#traits of interest
select<-disease%>%filter(cause=="Gout"|cause=="Ischemic heart disease" |cause=="Diabetes mellitus type 2" | cause=="Prostate cancer")

pdf(file=paste0(output_dir,"LifetimeRisk_facet_select.pdf"),height=4,width=12)
ggplot(select, aes(age, lifetimerisk, color=location, group=location)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 15), se = FALSE) + facet_wrap(~pretty_varname,nrow=1,scales="free_y") +
  geom_point(size=3,alpha=0.7) +
  xlab("Age Range") + 
  ylab("Baseline Cumulative Risk (%)") + 
  theme_bw() +
  labs(color='Location') +
  guides(color=guide_legend(reverse = TRUE)) +
  scale_color_manual(values=my_colors) +
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
    