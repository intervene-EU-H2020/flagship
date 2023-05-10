library(ggpubr)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)
library(dplyr)
library(data.table)

uk_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/KneeOsteoarthritis/KNEE_ARTHROSIS_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)

directory <- c("AllCancers","Asthma","ColorectalCancer","HipOsteoarthritis","KneeOsteoarthritis","LungCancer","MajorDepression","Melanoma")
phenotype <- c("C3_CANCER","J10_ASTHMA","C3_COLORECTAL","COX_ARTHROSIS","KNEE_ARTHROSIS","C3_BRONCHUS_LUNG","F5_DEPRESSIO","C3_MELANOMA_SKIN")

#Read in all figures
for(i in 1:length(directory)){
  #MGB
  mgb_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_t2d_m$Sex <- "Male"
  mgb_t2d_m$Country <- "Massachusetts (USA)"
  
  mgb_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_t2d_f$Sex <- "Female"
  mgb_t2d_f$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_t2d_m$Sex <- "Male"
  hunt_t2d_m$Country <- "Norway"
  
  hunt_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_t2d_f$Sex <- "Female"
  hunt_t2d_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_t2d_m$Sex <- "Male"
  uk_t2d_m$Country <- "United Kingdom"
  
  uk_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_t2d_f$Sex <- "Female"
  uk_t2d_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_t2d_m$Sex <- "Male"
  finngen_t2d_m$Country <- "Finland"
  
  finngen_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_t2d_f$Sex <- "Female"
  finngen_t2d_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_t2d_m$Sex <- "Male"
  estbb_t2d_m$Country <- "Estonia"
  
  estbb_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_t2d_f$Sex <- "Female"
  estbb_t2d_f$Country <- "Estonia"
  
  #All
  all <- rbind(mgb_t2d_m, mgb_t2d_f) %>%
    rbind(., hunt_t2d_m) %>%
    rbind(., hunt_t2d_f) %>%
    rbind(., uk_t2d_m) %>%
    rbind(., uk_t2d_f) %>%
    rbind(., finngen_t2d_m) %>%
    rbind(., finngen_t2d_f) %>%
    rbind(., estbb_t2d_m) %>%
    rbind(., estbb_t2d_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","Massachusetts (USA)","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
 
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=9, width=6, dpi=300)
  
}

#################################################################################################################################################################################################################

directory <- c("Appendicitis")
phenotype <- c("K11_APPENDACUT")

#Read in all figures
for(i in 1:length(directory)){
  #MGB
  mgb_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_t2d_m$Sex <- "Male"
  mgb_t2d_m$Country <- "Massachusetts (USA)"
  
  mgb_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_t2d_f$Sex <- "Female"
  mgb_t2d_f$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_t2d_m$Sex <- "Male"
  hunt_t2d_m$Country <- "Norway"
  
  hunt_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_t2d_f$Sex <- "Female"
  hunt_t2d_f$Country <- "Norway"
  
  #FinnGen
  finngen_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_t2d_m$Sex <- "Male"
  finngen_t2d_m$Country <- "Finland"
  
  finngen_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_t2d_f$Sex <- "Female"
  finngen_t2d_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_t2d_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_t2d_m$Sex <- "Male"
  estbb_t2d_m$Country <- "Estonia"
  
  estbb_t2d_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_t2d_f$Sex <- "Female"
  estbb_t2d_f$Country <- "Estonia"
  
  #All
  all <- rbind(mgb_t2d_m, mgb_t2d_f) %>%
    rbind(., hunt_t2d_m) %>%
    rbind(., hunt_t2d_f) %>%
    rbind(., finngen_t2d_m) %>%
    rbind(., finngen_t2d_f) %>%
    rbind(., estbb_t2d_m) %>%
    rbind(., estbb_t2d_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","Massachusetts (USA)","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=9, width=6, dpi=300)
  
}

#################################################################################################################################################################################################################

#Gout

directory <- c("Gout")
phenotype <- c("GOUT")

for(i in 1){
  #MGB
  mgb_gout_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_gout_m$Sex <- "Male"
  mgb_gout_m$Country <- "Massachusetts (USA)"
  
  mgb_gout_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_gout_f$Sex <- "Female"
  mgb_gout_f$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_gout_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_gout_m$Sex <- "Male"
  hunt_gout_m$Country <- "Norway"
  
  hunt_gout_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_gout_f$Sex <- "Female"
  hunt_gout_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_gout_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_gout_m$Sex <- "Male"
  uk_gout_m$Country <- "United Kingdom"
  
  uk_gout_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_gout_f$Sex <- "Female"
  uk_gout_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_gout_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_gout_m$Sex <- "Male"
  finngen_gout_m$Country <- "Finland"
  
  finngen_gout_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_gout_f$Sex <- "Female"
  finngen_gout_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_gout_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_gout_m$Sex <- "Male"
  estbb_gout_m$Country <- "Estonia"
  
  estbb_gout_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_gout_f$Sex <- "Female"
  estbb_gout_f$Country <- "Estonia"
  
  #All
  all <- rbind(mgb_gout_m, mgb_gout_f) %>%
    rbind(., hunt_gout_m) %>%
    rbind(., hunt_gout_f) %>%
    rbind(., uk_gout_m) %>%
    rbind(., uk_gout_f) %>%
    rbind(., finngen_gout_m) %>%
    rbind(., finngen_gout_f) %>%
    rbind(., estbb_gout_m) %>%
    rbind(., estbb_gout_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","Massachusetts (USA)","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#Atrial Fibrillation

directory <- c("AtrialFibrillation")
phenotype <- c("I9_AF")

for(i in 1){
  
  #HUNT
  hunt_af_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_af_m$Sex <- "Male"
  hunt_af_m$Country <- "Norway"
  
  hunt_af_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_af_f$Sex <- "Female"
  hunt_af_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_af_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_af_m$Sex <- "Male"
  uk_af_m$Country <- "United Kingdom"
  
  uk_af_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_af_f$Sex <- "Female"
  uk_af_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_af_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_af_m$Sex <- "Male"
  finngen_af_m$Country <- "Finland"
  
  finngen_af_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_af_f$Sex <- "Female"
  finngen_af_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_af_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_af_m$Sex <- "Male"
  estbb_af_m$Country <- "Estonia"
  
  estbb_af_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_af_f$Sex <- "Female"
  estbb_af_f$Country <- "Estonia"
  
  #All
  all <- rbind(hunt_af_m,hunt_af_f) %>%
    rbind(., uk_af_m) %>%
    rbind(., uk_af_f) %>%
    rbind(., finngen_af_m) %>%
    rbind(., finngen_af_f) %>%
    rbind(., estbb_af_m) %>%
    rbind(., estbb_af_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#Coronary Heart Disease

directory <- c("CHD")
phenotype <- c("I9_CHD")

for(i in 1){
  
  #MGB
  mgb_chd_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_chd_m$Sex <- "Male"
  mgb_chd_m$Country <- "Massachusetts (USA)"
  
  mgb_chd_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_chd_f$Sex <- "Female"
  mgb_chd_f$Country <- "Massachusetts (USA)"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_chd_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_chd_m$Sex <- "Male"
  uk_chd_m$Country <- "United Kingdom"
  
  uk_chd_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_chd_f$Sex <- "Female"
  uk_chd_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_chd_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_chd_m$Sex <- "Male"
  finngen_chd_m$Country <- "Finland"
  
  finngen_chd_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_chd_f$Sex <- "Female"
  finngen_chd_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_chd_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_chd_m$Sex <- "Male"
  estbb_chd_m$Country <- "Estonia"
  
  estbb_chd_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_chd_f$Sex <- "Female"
  estbb_chd_f$Country <- "Estonia"
  
  #All
  all <- rbind(mgb_chd_m,mgb_chd_f) %>%
    rbind(., uk_chd_m) %>%
    rbind(., uk_chd_f) %>%
    rbind(., finngen_chd_m) %>%
    rbind(., finngen_chd_f) %>%
    rbind(., estbb_chd_m) %>%
    rbind(., estbb_chd_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Massachusetts (USA)","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#Epilepsy

directory <- c("Epilepsy")
phenotype <- c("G6_EPLEPSY")

for(i in 1){
  
  #HUNT
  hunt_epilepsy_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_epilepsy_m$Sex <- "Male"
  hunt_epilepsy_m$Country <- "Norway"
  
  hunt_epilepsy_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_epilepsy_f$Sex <- "Female"
  hunt_epilepsy_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_epilepsy_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_epilepsy_m$Sex <- "Male"
  uk_epilepsy_m$Country <- "United Kingdom"
  
  uk_epilepsy_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_epilepsy_f$Sex <- "Female"
  uk_epilepsy_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_epilepsy_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_epilepsy_m$Sex <- "Male"
  finngen_epilepsy_m$Country <- "Finland"
  
  finngen_epilepsy_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_epilepsy_f$Sex <- "Female"
  finngen_epilepsy_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_epilepsy_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_epilepsy_m$Sex <- "Male"
  estbb_epilepsy_m$Country <- "Estonia"
  
  estbb_epilepsy_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_epilepsy_f$Sex <- "Female"
  estbb_epilepsy_f$Country <- "Estonia"
  
  #All
  all <- rbind(hunt_epilepsy_m,hunt_epilepsy_f) %>%
    rbind(., uk_epilepsy_m) %>%
    rbind(., uk_epilepsy_f) %>%
    rbind(., finngen_epilepsy_m) %>%
    rbind(., finngen_epilepsy_f) %>%
    rbind(., estbb_epilepsy_m) %>%
    rbind(., estbb_epilepsy_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","United Kingdom"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#ILD
directory <- c("ILD")
phenotype <- c("ILD")

for(i in 1){
  
  #MGB
  mgb_ild_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_ild_m$Sex <- "Male"
  mgb_ild_m$Country <- "Massachusetts (USA)"
  
  mgb_ild_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_ild_f$Sex <- "Female"
  mgb_ild_f$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_ild_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_ild_m$Sex <- "Male"
  hunt_ild_m$Country <- "Norway"
  
  hunt_ild_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_ild_f$Sex <- "Female"
  hunt_ild_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_ild_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_ild_m$Sex <- "Male"
  uk_ild_m$Country <- "United Kingdom"
  
  uk_ild_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_ild_f$Sex <- "Female"
  uk_ild_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_ild_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_ild_m$Sex <- "Male"
  finngen_ild_m$Country <- "Finland"
  
  finngen_ild_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_ild_f$Sex <- "Female"
  finngen_ild_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_ild_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_ild_m$Sex <- "Male"
  estbb_ild_m$Country <- "Estonia"
  
  estbb_ild_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_ild_f$Sex <- "Female"
  estbb_ild_f$Country <- "Estonia"
  
  #All
  all <- rbind(hunt_ild_m,hunt_ild_f) %>%
    rbind(., uk_ild_m) %>%
    rbind(., uk_ild_f) %>%
    rbind(., finngen_ild_m) %>%
    rbind(., finngen_ild_f) %>%
    rbind(., estbb_ild_m) %>%
    rbind(., estbb_ild_f) %>%
    rbind(., mgb_ild_m) %>%
    rbind(., mgb_ild_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland","Estonia","Norway","United Kingdom", "Massachusetts (USA)"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#Rheumatoid Arthritis
directory <- c("RheumatoidArthritis")
phenotype <- c("RHEUMA_SEROPOS_OTH")

for(i in 1){
  
  #MGB
  mgb_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_rheart_m$Sex <- "Male"
  mgb_rheart_m$Country <- "Massachusetts (USA)"
  
  mgb_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_rheart_f$Sex <- "Female"
  mgb_rheart_f$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_rheart_m$Sex <- "Male"
  hunt_rheart_m$Country <- "Norway"
  
  hunt_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_rheart_f$Sex <- "Female"
  hunt_rheart_f$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_rheart_m$Sex <- "Male"
  uk_rheart_m$Country <- "United Kingdom"
  
  uk_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_rheart_f$Sex <- "Female"
  uk_rheart_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_rheart_m$Sex <- "Male"
  finngen_rheart_m$Country <- "Finland"
  
  finngen_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_rheart_f$Sex <- "Female"
  finngen_rheart_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_rheart_m$Sex <- "Male"
  estbb_rheart_m$Country <- "Estonia"
  
  estbb_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_rheart_f$Sex <- "Female"
  estbb_rheart_f$Country <- "Estonia"
  
  #All
  all <- rbind(hunt_rheart_m,hunt_rheart_f) %>%
    rbind(., uk_rheart_m) %>%
    rbind(., uk_rheart_f) %>%
    rbind(., finngen_rheart_m) %>%
    rbind(., finngen_rheart_f) %>%
    rbind(., estbb_rheart_m) %>%
    rbind(., estbb_rheart_f) %>%
    rbind(., mgb_rheart_m) %>%
    rbind(., mgb_rheart_f)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland", "Estonia", "Norway", "United Kingdom", "Massachusetts (USA)"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
}

#################################################################################################################################################################################################################

#Type 1 Diabetes
directory <- c("T1D")
phenotype <- c("T1D")

for(i in 1){
  
  #FinnGen
  finngen_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_rheart_m$Sex <- "Male"
  finngen_rheart_m$Country <- "Finland"
  
  finngen_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_rheart_f$Sex <- "Female"
  finngen_rheart_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_rheart_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_rheart_m$Sex <- "Male"
  estbb_rheart_m$Country <- "Estonia"
  
  estbb_rheart_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_rheart_f$Sex <- "Female"
  estbb_rheart_f$Country <- "Estonia"
  
  #All
  all <- rbind(finngen_rheart_m, finngen_rheart_f) %>%
    rbind(., estbb_rheart_m) %>%
    rbind(., estbb_rheart_f) 
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland", "Estonia"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_grid(Country ~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=7, width=7, dpi=300)
}

#################################################################################################################################################################################################################

#Prostate Cancer
directory <- c("ProstateCancer")
phenotype <- c("C3_PROSTATE")

for(i in 1){
  
  mgb_prostate_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_prostate_m$Sex <- "Male"
  mgb_prostate_m$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_prostate_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_prostate_m$Sex <- "Male"
  hunt_prostate_m$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_prostate_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_prostate_m$Sex <- "Male"
  uk_prostate_m$Country <- "United Kingdom"
  
  #FinnGen
  finngen_prostate_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_prostate_m$Sex <- "Male"
  finngen_prostate_m$Country <- "Finland"
  
  #Estonian Biobank
  estbb_prostate_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_prostate_m$Sex <- "Male"
  estbb_prostate_m$Country <- "Estonia"
  
  #All
  all <- rbind(mgb_prostate_m, uk_prostate_m) %>%
           rbind(., hunt_prostate_m) %>%
                rbind(., finngen_prostate_m) %>%
                  rbind(., estbb_prostate_m) 
  
  all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  all$Country <- factor(all$Country, levels=c("Finland", "Estonia", "Norway", "Massachusetts (USA)", "United Kingdom"))
  
  groupsy <- c("< 20%", "40-60%", "> 95%")
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
  
  ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap( ~ Country ) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Boostrap.png"), height=7, width=7, dpi=300)
}

#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################

