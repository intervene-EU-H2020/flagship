library(ggpubr)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)
library(dplyr)
library(data.table)

directory <- c("AllCancers","Appendicitis","Asthma","ColorectalCancer","HipOsteoarthritis","KneeOsteoarthritis","LungCancer","MajorDepression")
phenotype <- c("C3_CANCER","K11_APPENDACUT","J10_ASTHMA","C3_COLORECTAL","COX_ARTHROSIS","KNEE_ARTHROSIS","C3_BRONCHUS_LUNG","F5_DEPRESSIO")

#Read in all figures
for(i in 1:length(directory)){
  
  #FinnGen
  finngen_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_m$Sex <- "Male"
  finngen_m$Country <- "Finland"
  
  finngen_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_f$Sex <- "Female"
  finngen_f$Country <- "Finland"
  
  all <- rbind(finngen_m, finngen_f)
  #spline_int <- as.data.frame(spline(all$age, all$LifetimeRisk))
  
  all$Group <- factor(all$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "< 1%"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  suppfiga <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(formula = y ~ s(x, k = 14), method = "gam", se = FALSE) +
    geom_point() +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata') +
    scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  
  groupsy <- c("< 1%", "40-60%", "> 99%")
  all <- rbind(finngen_m, finngen_f)
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 99%", "40-60%", "< 1%"))
  all$Sex <- factor(all$Sex, levels=c("Female","Male"))
  
  suppfigb <- ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(formula = y ~ s(x, k = 14), method = "gam", se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles_wCI.png"), height=6, width=6, dpi=300)
  
  suppfig <- ggarrange(suppfiga, suppfigb,
                         nrow=2)
  print(suppfig)
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles.png"), height=8, width=6, dpi=300)
  
  
}

############################################################################################################################################################################################################################################################

directory <- c("AtrialFibrillation","CHD","Epilepsy","Gout","ILD","Melanoma","RheumatoidArthritis")
phenotype <- c("I9_AF","I9_CHD","G6_EPLEPSY","GOUT","ILD","C3_MELANOMA_SKIN","RHEUMA_SEROPOS_OTH")

#Read in all figures
for(i in 1:length(directory)){
  
  #FinnGen
  finngen_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_m$Sex <- "Male"
  finngen_m$Country <- "Finland"
  
  finngen_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_f$Sex <- "Female"
  finngen_f$Country <- "Finland"
  
  all <- rbind(finngen_m, finngen_f)
  
  all$Group <- factor(all$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "< 5%"))
  
  suppfiga <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(formula = y ~ s(x, k = 14), method = "gam", se = FALSE) +
    geom_point() +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata') +
    scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "0-5%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  
  groupsy <- c("< 5%", "40-60%", "> 99%")
  all <- rbind(finngen_m, finngen_f)
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 99%", "40-60%", "< 5%"))
  
  suppfigb <- ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(formula = y ~ s(x, k = 14), method = "gam", se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-5%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-5%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles_wCI.png"), height=6, width=6, dpi=300)
  
  suppfig <- ggarrange(suppfiga, suppfigb,
                       nrow=2)
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles.png"), height=8, width=6, dpi=300)
  
  
}

############################################################################################################################################################################################################################################################

directory <- c("ProstateCancer")
phenotype <- c("C3_PROSTATE")

#Read in all figures
for(i in 1:length(directory)){
  
  #FinnGen
  finngen_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_m$Sex <- "Male"
  finngen_m$Country <- "Finland"
  
  finngen_m$Group <- factor(finngen_m$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "< 5%"))
  
  suppfiga <-ggplot(finngen_m, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PGS Strata') +
    scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "0-5%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  
  groupsy <- c("< 5%", "40-60%", "> 99%")
  finngen_m <- subset(finngen_m, Group %in% groupsy)
  
  finngen_m$Group <- factor(finngen_m$Group, levels=c("> 99%", "40-60%", "< 5%"))
  
  suppfigb <- ggplot(finngen_m, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-5%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-5%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles_wCI.png"), height=6, width=6, dpi=300)
  
  suppfig <- ggarrange(suppfiga, suppfigb,
                       nrow=2)
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles.png"), height=8, width=6, dpi=300)
  
  
}

############################################################################################################################################################################################################################################################

directory <- c("T1D")
phenotype <- c("T1D")

#Read in all figures
for(i in 1:length(directory)){
  
  #FinnGen
  finngen_m <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_m$Sex <- "Male"
  finngen_m$Country <- "Finland"
  
  finngen_f <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
  finngen_f$Sex <- "Female"
  finngen_f$Country <- "Finland"
  
  all <- rbind(finngen_m, finngen_f)
  
  all$Group <- factor(all$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
  
  suppfiga <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata') +
    scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  
  groupsy <- c("< 20%", "40-60%", "> 99%")
  all <- rbind(finngen_m, finngen_f)
  all <- subset(all, Group %in% groupsy)
  
  all$Group <- factor(all$Group, levels=c("> 99%", "40-60%", "< 20%"))
  
  suppfigb <- ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +
    facet_wrap(~ Sex) + 
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles_wCI.png"), height=6, width=6, dpi=300)
  
  suppfig <- ggarrange(suppfiga, suppfigb,
                       nrow=2)
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_FinnGen_TopPercentiles.png"), height=8, width=6, dpi=300)
  
  
}
