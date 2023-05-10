library(ggpubr)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)
library(dplyr)
library(data.table)

#Read in all figures

#MGB
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/BreastCancer/C3_BREAST_Female_LifetimeRisk_PartnersBiobank.csv", data.table=FALSE)
mgb$Country <- "Massachusetts (USA)"
mgb$Threshold <- 1.85913380706649

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/BreastCancer/C3_BREAST_Female_LifetimeRisk_HUNT.csv", data.table=FALSE)
hunt$Country <- "Norway"
hunt$Threshold <- 1.47119459974407

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenomicsEngland/BreastCancer/C3_BREAST_Female_LifetimeRisk_GenomicsEngland.csv", data.table=FALSE)
ge$Country <- "UK (Genomics England)"
ge$Threshold <- 2.05112936231737

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/GenerationScotland/BreastCancer/C3_BREAST_Female_LifetimeRisk_GenerationScotland.csv", data.table=FALSE)
gs$Country <- "Scotland"
gs$Threshold <- 2.38637060091801

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_FinnGen.csv", data.table=FALSE)
finngen$Country <- "Finland"
finngen$Threshold <- 1.70678589956649

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/BreastCancer/C3_BREAST_Female_LifetimeRisk_EstonianBiobank.csv", data.table=FALSE)
estbb$Country <- "Estonia"
estbb$Threshold <- 1.48670771704085

#UK Biobank
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKBiobank/BreastCancer/C3_BREAST_Female_LifetimeRisk_UKBiobank.csv", data.table=FALSE)
ukb$Country <- "UK (UK Biobank)"
ukb$Threshold <- 2.05112936231737

#South Asians
sas <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/SouthAsian/BreastCancer/C3_BREAST_LifetimeRisk_SouthAsian_UK.csv", data.table=FALSE)
sas$Country <- "UK (SAS)"
sas$Threshold <- 2.05112936231737

#All
all <- rbind(mgb) %>%
        rbind(., hunt) %>%
          rbind(., ge) %>%
            rbind(., gs) %>%
              rbind(., finngen) %>%
                rbind(., estbb) %>%
                  rbind(., ukb) %>%
                    rbind(., sas)

all$Group <- factor(all$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))
all$Country <- factor(all$Country, levels=c("Finland", "Estonia", "Norway", "Massachusetts (USA)", "UK (UK Biobank)", "UK (Genomics England)", "UK (SAS)", "Scotland"))

ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Country) + 
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/figure4.png"), height=7, width=7, dpi=300)

######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################

#Supplementary Figure 9

#MGB
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb$Country <- "Massachusetts (USA)"
mgb$Threshold <- 1.85913380706649

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_HUNT.csv", data.table=FALSE)
hunt$Country <- "Norway"
hunt$Threshold <- 1.47119459974407

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen$Country <- "Finland"
finngen$Threshold <- 1.70678589956649

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb$Country <- "Estonia"
estbb$Threshold <- 1.48670771704085

#United Kingdom meta analysis
uk <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv", data.table=FALSE)
uk$Country <- "United Kingdom"
uk$Threshold <- 2.05112936231737


#All
all <- rbind(mgb, hunt) %>%
         rbind(., finngen) %>%
           rbind(., estbb) %>%
             rbind(., uk) 

groupsy <- c("< 20%", "40-60%", "> 95%")
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 95%", "40-60%", "< 20%"))
all$Country <- factor(all$Country, levels=c("Finland", "Estonia", "Norway", "Massachusetts (USA)", "United Kingdom"))

print(ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  #stat_smooth(formula = y ~ s(x, k = 4), method = "gam", se = FALSE) + 
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  #geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Country) + 
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure11_BreCanLREwithCI.png"), height=8, width=7, dpi=300)

############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################

#Extension to top and bottom percentiles for FinnGen

#FinnGen
finngen_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_f$Country <- "Finland"
finngen_f$Threshold <- 1.70678589956649

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "< 1%"))

suppfiga <-ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

groupsy <- c("< 1%", "40-60%", "> 99%")
finngen_f <- subset(finngen_f, Group %in% groupsy)

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "40-60%", "< 1%"))

suppfigb <- ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

suppfig10 <- ggarrange(suppfiga, suppfigb,
                       nrow=2,
                       labels=c("a","b"))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure12_BreCanOnePercentTails.png"), height=8, width=6, dpi=300)

###################################################################################################################################################################################
###################################################################################################################################################################################
#FinnGen
finngen_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv", data.table=FALSE)
finngen_f$Country <- "Finland"
finngen_f$Threshold <- 1.70678589956649

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "< 1%"))

suppfiga <-ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

groupsy <- c("< 1%", "40-60%", "> 99%")
finngen_f <- subset(finngen_f, Group %in% groupsy)

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "40-60%", "< 1%"))

suppfigb <- ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/BreCanOnePercentTailswCI.png"), height=6, width=6, dpi=300)

suppfig10 <- ggarrange(suppfiga, suppfigb,
                       nrow=2,
                       labels=c("a","b"))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure12_BreCanOnePercentTails.png"), height=6, width=6, dpi=300)

###################################################################################################################################################################################
#Figure 4

#FinnGen
finngen_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_f$Country <- "Finland"
finngen_f$Threshold <- 1.70678589956649

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))

suppfiga <-ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata') +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/BreastCancer_Finland_WithoutCI_UpdatedFigure4c.png"), height=6, width=8, dpi=300)

groupsy <- c("< 20%", "40-60%", "> 95%")
finngen_f <- subset(finngen_f, Group %in% groupsy)

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 95%", "40-60%", "< 20%"))

suppfigb <- ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/BreastCancer_Finland_WithCI_UpdatedFigure4d.png"), height=6, width=8, dpi=300)

###################################################################################################################################################################################
###################################################################################################################################################################################

#FinnGen
finngen_f <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_f$Country <- "Finland"
finngen_f$Threshold <- 1.70678589956649

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 95%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "< 20%"))

suppfiga <-ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  scale_fill_hue(labels = c("95-100%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "0-20%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/BreastCancer_Finland_WithoutCI_UpdatedFigure4AllCILevels.png"), height=6, width=8, dpi=300)

