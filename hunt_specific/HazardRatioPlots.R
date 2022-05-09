library(data.table)
library(ggplot2)
library(dplyr)

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

#UK Biobank - HRs per SDs

associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRsperSD_UKBiobank_FullSample.csv", data.table=FALSE)
associations <- associations[,-1]
colnames(associations) <- c("phenotype","prs","beta","se","pval","OR","Cipos","Cineg")
associations <- associations[,c("phenotype","prs","OR","Cipos","Cineg")]

associations$phenotype <- factor(associations$phenotype, levels=c("G6_EPLEPSY", "C3_BREAST", "GOUT", "RHEUMA_SEROPOS_OTH", "C3_PROSTATE", "T1D"))

phenotypelabels <- c("Epilepsy", "Breast Cancer", "Gout", "Rheumatoid Arthritis", "Prostate Cancer", "T1D")
  
ggplot(associations) +
  geom_point(aes(phenotype, OR)) +
  theme_bw() +
  geom_errorbar(aes(x=phenotype, ymin = Cineg, ymax = Cipos), size=0.5, width=0.5) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_x_discrete(labels= phenotypelabels) +
  ggtitle("HR Per SD") + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatiosPerSD_PRS_UKBiobank.png", height=10 , width=10)

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

#UK Biobank - HRs - Males vs Females vs Both - Top 1 Percent

male_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_MaleSample_top1percent.csv", data.table=FALSE)
male_sample_associations <- male_sample_associations[,-1]
male_sample_associations$Sex <- "Male"
colnames(male_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg","Sex")
male_sample_associations <- male_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

female_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FemaleSample_top1percent.csv", data.table=FALSE)
female_sample_associations <- female_sample_associations[,-1]
colnames(female_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg")
female_sample_associations$Sex <- "Female"
female_sample_associations <- female_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

associations <- rbind(male_sample_associations, female_sample_associations) 

associations$Sex <- as.factor(associations$Sex)
associations$phenotype <- as.factor(associations$phenotype)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5"),paste0(i,"_groupGroup 7"), 
                                                  paste0(i,"_groupGroup 8"), paste0(i,"_groupGroup 9"), paste0(i,"_groupGroup 10"), paste0(i,"_groupGroup 11")))
  
  phenotypelabels <- c("<1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", ">99%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=Sex, col=Sex), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, " - Sex Stratification")) + 
          scale_colour_manual(values=c("black","orange")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_MalevsFemale_UKBiobank.png"), height=10 , width=10)
}

####################################################################################################################################################################################################################

#UK Biobank - HRs - Males vs Females vs Both - Top 10 Percent

male_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_MaleSample_top10percent.csv", data.table=FALSE)
male_sample_associations <- male_sample_associations[,-1]
male_sample_associations$Sex <- "male"
colnames(male_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg","Sex")
male_sample_associations <- male_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

female_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FemaleSample_top10percent.csv", data.table=FALSE)
female_sample_associations <- female_sample_associations[,-1]
colnames(female_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg")
female_sample_associations$Sex <- "female"
female_sample_associations <- female_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

full_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FullSample_top10percent.csv", data.table=FALSE) 
full_sample_associations <- full_sample_associations[,-1]
colnames(full_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg")
full_sample_associations$Sex <- "both"
full_sample_associations <- full_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

associations <- rbind(male_sample_associations, female_sample_associations) %>%
  rbind(., full_sample_associations)

associations$Sex <- as.factor(associations$Sex)
associations$phenotype <- as.factor(associations$phenotype)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 5"), paste0(i,"_groupGroup 6"),paste0(i,"_groupGroup 7")))
  
  phenotypelabels <- c("0-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-100%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=Sex, col=Sex), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, " - Sex Stratification")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_UKBiobank_top10percent.png"), height=10 , width=10)
}

####################################################################################################################################################################################################################

#UK Biobank - HRs - Males vs Females vs Both - Top 20 Percent

male_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_MaleSample_top20percent.csv", data.table=FALSE)
male_sample_associations <- male_sample_associations[,-1]
male_sample_associations$Sex <- "male"
colnames(male_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg","Sex")
male_sample_associations <- male_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

female_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FemaleSample_top20percent.csv", data.table=FALSE)
female_sample_associations <- female_sample_associations[,-1]
colnames(female_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg")
female_sample_associations$Sex <- "female"
female_sample_associations <- female_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

full_sample_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FullSample_top20percent.csv", data.table=FALSE) 
full_sample_associations <- full_sample_associations[,-1]
colnames(full_sample_associations) <- c("phenotype","prs","group","beta","se","pval","OR","Cipos","Cineg")
full_sample_associations$Sex <- "both"
full_sample_associations <- full_sample_associations[,c("phenotype","prs","group","OR","Cipos","Cineg","Sex")]

associations <- rbind(male_sample_associations, female_sample_associations) %>%
  rbind(., full_sample_associations)

associations$Sex <- as.factor(associations$Sex)
associations$phenotype <- as.factor(associations$phenotype)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5")))
  
  phenotypelabels <- c("0-20%", "20-40%", "60-80%", "80-100%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=Sex, col=Sex), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          scale_colour_manual(values=c("black","orange","blue")) + 
          ggtitle(paste0(i, " - Sex Stratification")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_UKBiobank_top20percent.png"), height=10 , width=10)
}

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

#Age stratification - Top 20 percent

associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRSs_HR_UKBiobank_AgeStratifiedResults_top20percent.csv", data.table=FALSE)
associations$agecat <- paste0(associations$minage," - ", associations$maxage)
associations$agecat<- as.factor(associations$agecat)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5")))
  
  phenotypelabels <- c("0-20%", "20-40%", "60-80%", "80-100%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=agecat, col=agecat), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=agecat, col=agecat), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, " - Age Stratification")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_AgeStratification_top20percent_UKBiobank.png"), height=10 , width=10)
}

####################################################################################################################################################################################################################

#Age stratification - Top 10 percent

associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRSs_HR_UKBiobank_AgeStratifiedResults_top10percent.csv", data.table=FALSE)
associations$agecat <- paste0(associations$minage," - ", associations$maxage)
associations$agecat<- as.factor(associations$agecat)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 5"), paste0(i,"_groupGroup 6"), paste0(i,"_groupGroup 7")))
  
  phenotypelabels <- c("0-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-100%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=agecat, col=agecat), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=agecat, col=agecat), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, " - Age Stratification")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_AgeStratification_top10percent_UKBiobank.png"), height=10 , width=10)
}

####################################################################################################################################################################################################################

#Age stratification - Top 1 percent

associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRSs_HR_UKBiobank_AgeStratifiedResults_top1percent.csv", data.table=FALSE)
associations <- associations[,-1]
colnames(associations) <- c("phenotype","prs","minage","maxage","group","controls","cases","betas","se","pval","OR","Cipos","Cineg")
associations$agecat <- paste0(associations$minage," - ", associations$maxage)
associations$agecat<- as.factor(associations$agecat)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5"), paste0(i,"_groupGroup 7"), paste0(i,"_groupGroup 8"), paste0(i,"_groupGroup 9"), paste0(i,"_groupGroup 10"), paste0(i,"_groupGroup 11")))
  
  phenotypelabels <- c("<1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", ">99%")
  
  print(ggplot(disease) +
          geom_point(aes(group, OR, group=agecat, col=agecat), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=agecat, col=agecat), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          ggtitle(paste0(i, " - Age Stratification")) + 
          scale_x_discrete(labels= phenotypelabels) +
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_AgeStratification_top1percent_UKBiobank.png"), height=10 , width=10)
}

####################################################################################################################################################################################
####################################################################################################################################################################################

#UK Biobank vs FinnGen  - HRs - Top 1 Percent

ukb_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/UKBiobank/PRS_HRs_UKBiobank_FullSample_top1percent.csv", data.table=FALSE)
ukb_associations <- ukb_associations[,-1]
ukb_associations$Biobank <- "UKB"
colnames(ukb_associations) <- c("phenotype","prs","group","beta","se","pval","HR","Cipos","Cineg","Biobank")
ukb_associations <- ukb_associations[,c("phenotype","prs","group","HR","Cipos","Cineg","Biobank")]

finngen_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios_FullSample_FinnGen.csv", data.table=FALSE)
finngen_associations <- finngen_associations[,-1]
finngen_associations$Biobank <- "FinnGen"
colnames(finngen_associations) <- c("phenotype","prs","group","beta","se","pval","HR","Cipos","Cineg","Biobank")
phenos <- c("C3_BREAST", "G6_EPLEPSY", "GOUT", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D")
finngen_associations <- subset(finngen_associations, phenotype %in% phenos)
finngen_associations <- finngen_associations[,c("phenotype","prs","group","HR","Cipos","Cineg","Biobank")]

associations <- rbind(ukb_associations, finngen_associations) 

associations$Biobank <- as.factor(associations$Biobank)
associations$phenotype <- as.factor(associations$phenotype)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5"),paste0(i,"_groupGroup 7"), 
                                                  paste0(i,"_groupGroup 8"), paste0(i,"_groupGroup 9"), paste0(i,"_groupGroup 10"), paste0(i,"_groupGroup 11")))
  
  phenotypelabels <- c("<1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", ">99%")
  
  print(ggplot(disease) +
          geom_point(aes(group, HR, group=Biobank, col=Biobank), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=Biobank, col=Biobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          scale_colour_manual(values=c("black","orange")) + 
          ggtitle(paste0(i, "- UKBiobank vs FinnGen")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_FinngenvsUKBiobank.png"), height=10 , width=10)
}

####################################################################################################################################################################################
####################################################################################################################################################################################

#FinnGen - Male vs Female - HRs - Top 1 Percent

male_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios_MaleSample_FinnGen.csv", data.table=FALSE)
male_associations <- male_associations[,-1]
male_associations$Sex<- "Male"
colnames(male_associations) <- c("phenotype","prs","group","beta","se","pval","HR","Cipos","Cineg","Sex")
male_associations <- male_associations[,c("phenotype","prs","group","HR","Cipos","Cineg","Sex")]

female_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios_FemaleSample_FinnGen.csv", data.table=FALSE)
female_associations <- female_associations[,-1]
female_associations$Sex <- "Female"
colnames(female_associations) <- c("phenotype","prs","group","beta","se","pval","HR","Cipos","Cineg","Sex")
female_associations <- female_associations[,c("phenotype","prs","group","HR","Cipos","Cineg","Sex")]

associations <- rbind(male_associations, female_associations) 

associations$Sex <- as.factor(associations$Sex)
associations$phenotype <- as.factor(associations$phenotype)

for(i in unique(associations$prs)){
  print(i)
  disease <- subset(associations, prs==i)
  disease$group <- factor(disease$group, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5"),paste0(i,"_groupGroup 7"), 
                                                  paste0(i,"_groupGroup 8"), paste0(i,"_groupGroup 9"), paste0(i,"_groupGroup 10"), paste0(i,"_groupGroup 11")))
  
  phenotypelabels <- c("<1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", ">99%")
  
  print(ggplot(disease) +
          geom_point(aes(group, HR, group=Sex, col=Sex), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=group, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_colour_manual(values=c("black","orange")) + 
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, "- Male vs Female")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_Finngen_MalevsFemale.png"), height=10, width=10)
}

####################################################################################################################################################################################
####################################################################################################################################################################################

#FinnGen - Age Stratification

age_associations <- fread("/Users/jermy/Documents/INTERVENE/Results/FinnGen/HazardRatios_AgeStratifiedResults_FinnGen.csv", data.table=FALSE)
age_associations$AgeQuartile <- paste0(age_associations$MinAge," - ",age_associations$MaxAge)

age_associations$AgeQuartile <- as.factor(age_associations$AgeQuartile)
age_associations$Phenotype <- as.factor(age_associations$Phenotype)

age_associations$HR <- ifelse(age_associations$PRS=="AllCancers", 1/age_associations$HR, age_associations$HR)
age_associations$Cipos <- ifelse(age_associations$PRS=="AllCancers", 1/age_associations$Cipos, age_associations$Cipos)
age_associations$Cineg <- ifelse(age_associations$PRS=="AllCancers", 1/age_associations$Cineg, age_associations$Cineg)

for(i in unique(age_associations$PRS)){
  print(i)
  disease <- subset(age_associations, PRS==i)
  disease$PRS_GROUP <- factor(disease$PRS_GROUP, levels=c(paste0(i,"_groupGroup 1"), paste0(i,"_groupGroup 2"), paste0(i,"_groupGroup 3"), paste0(i,"_groupGroup 4"), paste0(i,"_groupGroup 5"),paste0(i,"_groupGroup 7"), 
                                                  paste0(i,"_groupGroup 8"), paste0(i,"_groupGroup 9"), paste0(i,"_groupGroup 10"), paste0(i,"_groupGroup 11")))
  
  phenotypelabels <- c("<1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", ">99%")
  
  print(ggplot(disease) +
          geom_point(aes(PRS_GROUP, HR, group=AgeQuartile, col=AgeQuartile), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=PRS_GROUP, ymin = Cineg, ymax = Cipos, group=AgeQuartile, col=AgeQuartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_x_discrete(labels= phenotypelabels) +
          ggtitle(paste0(i, "- Age Stratification")) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 16)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/Plots/HazardRatios_",i,"_PRS_Finngen_AgeStratified.png"), height=10, width=10)
}

