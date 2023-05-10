library(data.table)
library(ggpubr)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)
library(dplyr)

#Read in all figures

#MGB
mgb_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_MaleLifetimeRisk_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_m$Sex <- "Male"
mgb_t2d_m$Country <- "Massachusetts (USA)"
mgb_t2d_m$Threshold <- 7.27714975

mgb_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_FemaleLifetimeRisk_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_f$Sex <- "Female"
mgb_t2d_f$Country <- "Massachusetts (USA)"
mgb_t2d_f$Threshold <- 5.40494655

#HUNT
hunt_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_MaleLifetimeRisk_HUNT.csv", data.table=FALSE)
hunt_t2d_m$Sex <- "Male"
hunt_t2d_m$Country <- "Norway"
hunt_t2d_m$Threshold <- 6.10085978

hunt_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_FemaleLifetimeRisk_HUNT.csv", data.table=FALSE)
hunt_t2d_f$Sex <- "Female"
hunt_t2d_f$Country <- "Norway"
hunt_t2d_f$Threshold <- 4.70159277

#Genomics England
ge_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenomicsEngland/T2D/T2D_MaleLifetimeRisk_GenomicsEngland.csv", data.table=FALSE)
ge_t2d_m$Sex <- "Male"
ge_t2d_m$Country <- "United Kingdom"
ge_t2d_m$Threshold <- 13.0264341

ge_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenomicsEngland/T2D/T2D_FemaleLifetimeRisk_GenomicsEngland.csv", data.table=FALSE)
ge_t2d_f$Sex <- "Female"
ge_t2d_f$Country <- "United Kingdom"
ge_t2d_f$Threshold <- 8.91818627

#Generation Scotland
gs_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenerationScotland/T2D/T2D_MaleLifetimeRisk_GenerationScotland.csv", data.table=FALSE)
gs_t2d_m$Sex <- "Male"
gs_t2d_m$Country <- "Scotland"
gs_t2d_m$Threshold <- 8.018056

gs_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenerationScotland/T2D/T2D_FemaleLifetimeRisk_GenerationScotland.csv", data.table=FALSE)
gs_t2d_f$Sex <- "Female"
gs_t2d_f$Country <- "Scotland"
gs_t2d_f$Threshold <- 5.66514572

#FinnGen
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_FinnGen.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_FinnGen.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

#Estonian Biobank
estbb_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_MaleLifetimeRisk_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_m$Sex <- "Male"
estbb_t2d_m$Country <- "Estonia"
estbb_t2d_m$Threshold <- 5.61893987

estbb_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_FemaleLifetimeRisk_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_f$Sex <- "Female"
estbb_t2d_f$Country <- "Estonia"
estbb_t2d_f$Threshold <- 4.73314501

#Biobank Japan
bbj_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/BiobankJapan/T2D/T2D_MaleLifetimeRisk_BiobankJapan.csv", data.table=FALSE)
bbj_t2d_m$Sex <- "Male"
bbj_t2d_m$Country <- "Japan"
bbj_t2d_m$Threshold <- 5.770162214

bbj_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/BiobankJapan/T2D/T2D_FemaleLifetimeRisk_BiobankJapan.csv", data.table=FALSE)
bbj_t2d_f$Sex <- "Female"
bbj_t2d_f$Country <- "Japan"
bbj_t2d_f$Threshold <- 3.440766689

#All
all <- rbind(mgb_t2d_m, mgb_t2d_f) %>%
        rbind(., hunt_t2d_m) %>%
          rbind(., hunt_t2d_f) %>%
            rbind(., ge_t2d_m) %>%
              rbind(., ge_t2d_f) %>%
                rbind(., gs_t2d_m) %>%
                  rbind(., gs_t2d_f) %>%
                    rbind(., finngen_t2d_m) %>%
                      rbind(., finngen_t2d_f) %>%
                        rbind(., estbb_t2d_m) %>%
                          rbind(., estbb_t2d_f) %>%
                            rbind(., bbj_t2d_m) %>%
                              rbind(., bbj_t2d_f) 

all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))

ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_grid(Country ~ Sex) + 
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/figure3.png"), height=10, width=6, dpi=300)

######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################

#Supplementary Figure 9

mgb_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_MaleLifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_m$Sex <- "Male"
mgb_t2d_m$Country <- "Massachusetts (USA)"
mgb_t2d_m$Threshold <- 7.27714975

mgb_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_f$Sex <- "Female"
mgb_t2d_f$Country <- "Massachusetts (USA)"
mgb_t2d_f$Threshold <- 5.40494655

#HUNT
hunt_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_MaleLifetimeRisk_Bootstrapped_HUNT.csv", data.table=FALSE)
hunt_t2d_m$Sex <- "Male"
hunt_t2d_m$Country <- "Norway"
hunt_t2d_m$Threshold <- 6.10085978

hunt_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_HUNT.csv", data.table=FALSE)
hunt_t2d_f$Sex <- "Female"
hunt_t2d_f$Country <- "Norway"
hunt_t2d_f$Threshold <- 4.70159277

#United Kingdom Meta Analysis
uk_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenomicsEngland/T2D/T2D_MaleLifetimeRisk_Bootstrapped_GenomicsEngland.csv", data.table=FALSE)
uk_t2d_m$Sex <- "Male"
uk_t2d_m$Country <- "United Kingdom"
uk_t2d_m$Threshold <- 13.0264341

uk_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenomicsEngland/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_GenomicsEngland.csv", data.table=FALSE)
uk_t2d_f$Sex <- "Female"
uk_t2d_f$Country <- "United Kingdom"
uk_t2d_f$Threshold <- 8.91818627

#FinnGen
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

#Estonian Biobank
estbb_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_MaleLifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_m$Sex <- "Male"
estbb_t2d_m$Country <- "Estonia"
estbb_t2d_m$Threshold <- 5.61893987

estbb_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_f$Sex <- "Female"
estbb_t2d_f$Country <- "Estonia"
estbb_t2d_f$Threshold <- 4.73314501

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

groupsy <- c("< 20%", "40-60%", "> 95%")
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))

ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  #geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
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
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure9.png"), height=10, width=6, dpi=300)

############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################

#Extension to top and bottom percentiles for FinnGen

#FinnGen
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

all <- rbind(finngen_t2d_m, finngen_t2d_f)

all$Group <- factor(all$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "< 1%"))

suppfig10a <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
                stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
                geom_point() +
                xlab("Age") + 
                ylab("Cumulative Risk (%)") + 
                geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
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
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure10a.png"), height=8, width=6, dpi=300)

#Bootstrapped
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

groupsy <- c("< 1%", "40-60%", "> 99%")
all <- rbind(finngen_t2d_m, finngen_t2d_f)
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 99%", "40-60%", "< 1%"))

suppfig10b <- ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
                stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
                geom_point() +
                geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
                xlab("Age") + 
                ylab("Cumulative Risk (%)") + 
                geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
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
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure10b.png"), height=6, width=6, dpi=300)

suppfig10 <- ggarrange(suppfig10a, suppfig10b,
                       nrow=2,
                       labels=c("a","b"))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure10.png"), height=8, width=6, dpi=300)

###################################################################################################################################################################################
###################################################################################################################################################################################

#Updated to include just finngen with and without confidence intervals as will make main figure for figure 4. 
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_FinnGen.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_FinnGen.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

all <- rbind(finngen_t2d_m, finngen_t2d_f)

all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))

suppfig10a <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Sex) + 
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/T2D_ScreeningFinnGenOnly_UpdatedFigure4a_withoutCI.png"), height=6, width=8, dpi=300)

###################################################################################################################################################################################

finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

groupsy <- c("< 20%", "40-60%", "> 95%")
all <- rbind(finngen_t2d_m, finngen_t2d_f)
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))

suppfig10b <- ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Sex) + 
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/T2D_ScreeningFinnGenOnly_UpdatedFigure4b_withCI.png"), dpi=300)

###################################################################################################################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################

#Test for inclusion of bootstrapped confidence intervals for all sections
finngen_t2d_m <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

all <- rbind(finngen_t2d_m, finngen_t2d_f)

all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))

suppfig10a <-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  theme_bw() +
  facet_wrap(~ Sex) + 
  labs(color='PGS Strata', fill="PGS Strata") +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  scale_fill_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/T2D_ScreeningFinnGenOnly_UpdatedFigure4_AllCILevels.png"), height=6, width=8, dpi=300)

