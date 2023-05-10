#T2D

#Biobank Specific

#MGB
mgb_t2d_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_MaleLifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_m_bs$Sex <- "Male"
mgb_t2d_m_bs$Country <- "Massachusetts (USA)"
mgb_t2d_m_bs$Threshold <- 7.27714975
mgb_t2d_m_bs$Test <- "Biobank Specific" 

mgb_t2d_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_t2d_f_bs$Sex <- "Female"
mgb_t2d_f_bs$Country <- "Massachusetts (USA)"
mgb_t2d_f_bs$Threshold <- 5.40494655
mgb_t2d_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
mgb_t2d_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_MaleLifetimeRisk_MetaAnalysis_Bootstrapped_Massachusetts.csv", data.table=FALSE)
mgb_t2d_m_ma$Sex <- "Male"
mgb_t2d_m_ma$Country <- "Massachusetts (USA)"
mgb_t2d_m_ma$Threshold <- 7.27714975
mgb_t2d_m_ma$Test <- "Meta-Analysis" 

mgb_t2d_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_FemaleLifetimeRisk_MetaAnalysis_Bootstrapped_Massachusetts.csv", data.table=FALSE)
mgb_t2d_f_ma$Sex <- "Female"
mgb_t2d_f_ma$Country <- "Massachusetts (USA)"
mgb_t2d_f_ma$Threshold <- 5.40494655
mgb_t2d_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#HUNT
hunt_t2d_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_MaleLifetimeRisk_Bootstrapped_HUNT.csv", data.table=FALSE)
hunt_t2d_m_bs$Sex <- "Male"
hunt_t2d_m_bs$Country <- "Norway"
hunt_t2d_m_bs$Threshold <- 6.10085978
hunt_t2d_m_bs$Test <- "Biobank Specific" 

hunt_t2d_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/HUNT/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_HUNT.csv", data.table=FALSE)
hunt_t2d_f_bs$Sex <- "Female"
hunt_t2d_f_bs$Country <- "Norway"
hunt_t2d_f_bs$Threshold <- 4.70159277
hunt_t2d_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
hunt_t2d_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_MaleLifetimeRisk_MetaAnalysis_Bootstrapped_Norway.csv", data.table=FALSE)
hunt_t2d_m_ma$Sex <- "Male"
hunt_t2d_m_ma$Country <- "Norway"
hunt_t2d_m_ma$Threshold <- 6.10085978
hunt_t2d_m_ma$Test <- "Meta-Analysis" 

hunt_t2d_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_FemaleLifetimeRisk_MetaAnalysis_Bootstrapped_Norway.csv", data.table=FALSE)
hunt_t2d_f_ma$Sex <- "Female"
hunt_t2d_f_ma$Country <- "Norway"
hunt_t2d_f_ma$Threshold <- 4.70159277
hunt_t2d_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#Genomics England
uk_t2d_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/T2D/T2D_MaleLifetimeRisk_Bootstrapped_UKMetaAnalysis.csv", data.table=FALSE)
uk_t2d_m_bs$Sex <- "Male"
uk_t2d_m_bs$Country <- "United Kingdom"
uk_t2d_m_bs$Threshold <- 13.0264341
uk_t2d_m_bs$Test <- "Biobank Specific" 

uk_t2d_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_UKMetaAnalysis.csv", data.table=FALSE)
uk_t2d_f_bs$Sex <- "Female"
uk_t2d_f_bs$Country <- "United Kingdom"
uk_t2d_f_bs$Threshold <- 8.91818627
uk_t2d_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
uk_t2d_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_MaleLifetimeRisk_MetaAnalysis_Bootstrapped_United Kingdom.csv", data.table=FALSE)
uk_t2d_m_ma$Sex <- "Male"
uk_t2d_m_ma$Country <- "United Kingdom"
uk_t2d_m_ma$Threshold <- 13.0264341
uk_t2d_m_ma$Test <- "Meta-Analysis" 

uk_t2d_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_FemaleLifetimeRisk_MetaAnalysis_Bootstrapped_United Kingdom.csv", data.table=FALSE)
uk_t2d_f_ma$Sex <- "Female"
uk_t2d_f_ma$Country <- "United Kingdom"
uk_t2d_f_ma$Threshold <- 8.91818627
uk_t2d_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#FinnGen
finngen_t2d_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_m_bs$Sex <- "Male"
finngen_t2d_m_bs$Country <- "Finland"
finngen_t2d_m_bs$Threshold <- 7.380580738985
finngen_t2d_m_bs$Test <- "Biobank Specific" 

finngen_t2d_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_t2d_f_bs$Sex <- "Female"
finngen_t2d_f_bs$Country <- "Finland"
finngen_t2d_f_bs$Threshold <- 6.40599016524789
finngen_t2d_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
finngen_t2d_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_MaleLifetimeRisk_MetaAnalysis_Bootstrapped_Finland.csv", data.table=FALSE)
finngen_t2d_m_ma$Sex <- "Male"
finngen_t2d_m_ma$Country <- "Finland"
finngen_t2d_m_ma$Threshold <- 7.380580738985
finngen_t2d_m_ma$Test <- "Meta-Analysis" 

finngen_t2d_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_FemaleLifetimeRisk_MetaAnalysis_Bootstrapped_Finland.csv", data.table=FALSE)
finngen_t2d_f_ma$Sex <- "Female"
finngen_t2d_f_ma$Country <- "Finland"
finngen_t2d_f_ma$Threshold <- 6.40599016524789
finngen_t2d_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#Estonian Biobank
estbb_t2d_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_MaleLifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_m_bs$Sex <- "Male"
estbb_t2d_m_bs$Country <- "Estonia"
estbb_t2d_m_bs$Threshold <- 5.61893987
estbb_t2d_m_bs$Test <- "Biobank Specific" 

estbb_t2d_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_t2d_f_bs$Sex <- "Female"
estbb_t2d_f_bs$Country <- "Estonia"
estbb_t2d_f_bs$Threshold <- 4.73314501
estbb_t2d_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
estbb_t2d_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_MaleLifetimeRisk_MetaAnalysis_Bootstrapped_Estonia.csv", data.table=FALSE)
estbb_t2d_m_ma$Sex <- "Male"
estbb_t2d_m_ma$Country <- "Estonia"
estbb_t2d_m_ma$Threshold <- 5.61893987
estbb_t2d_m_ma$Test <- "Meta-Analysis" 

estbb_t2d_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/T2D/T2D_FemaleLifetimeRisk_MetaAnalysis_Bootstrapped_Estonia.csv", data.table=FALSE)
estbb_t2d_f_ma$Sex <- "Female"
estbb_t2d_f_ma$Country <- "Estonia"
estbb_t2d_f_ma$Threshold <- 4.73314501
estbb_t2d_f_ma$Test <- "Meta-Analysis" 

#All
all <- rbind(mgb_t2d_m_bs, mgb_t2d_f_bs) %>%
         rbind(., mgb_t2d_m_ma) %>%
          rbind(., mgb_t2d_f_ma) %>%
            rbind(., hunt_t2d_m_bs) %>%
              rbind(., hunt_t2d_m_ma) %>%
                rbind(., hunt_t2d_f_bs) %>%
                  rbind(., hunt_t2d_f_ma) %>%
                    rbind(., uk_t2d_m_bs) %>%
                      rbind(., uk_t2d_m_ma) %>%
                        rbind(., uk_t2d_f_bs) %>%
                          rbind(., uk_t2d_f_ma) %>%
                            rbind(., finngen_t2d_m_bs) %>%
                              rbind(., finngen_t2d_m_ma) %>%
                                rbind(., finngen_t2d_f_bs) %>%
                                  rbind(., finngen_t2d_f_ma) %>%
                                    rbind(., estbb_t2d_m_bs) %>%
                                      rbind(., estbb_t2d_m_ma) %>%
                                        rbind(., estbb_t2d_f_bs) %>%
                                          rbind(., estbb_t2d_f_ma)

groupsy <- c("< 20%", "> 95%")
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 95%", "< 20%"))
all$Test <- factor(all$Test, levels=c("Biobank Specific", "Meta-Analysis"))

ggplot(all, aes(age, LifetimeRisk, color=interaction(Group,Test), group=interaction(Group,Test))) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=interaction(Group,Test)), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_grid(Country ~ Sex) + 
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("lightblue", "red", "darkblue", "darkred"), labels=c("Biobank Specific (95-100%)", "Biobank Specific (0-20%)", "Meta-Analysis (95-100%)", "Meta-Analysis (0-20%)")) +
  scale_fill_manual(values=c("lightblue", "red", "darkblue", "darkred"), labels=c("Biobank Specific (95-100%)", "Biobank Specific (0-20%)", "Meta-Analysis (95-100%)", "Meta-Analysis (0-20%)")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure_20a_T2D_MetaAnalysisvsBiobankSpecific.png"), height=8, width=10, dpi=300)

###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################

#CHD

#Biobank Specific

#MGB
mgb_CHD_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/CHD/I9_CHD_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_CHD_m_bs$Sex <- "Male"
mgb_CHD_m_bs$Country <- "Massachusetts (USA)"
mgb_CHD_m_bs$Test <- "Biobank Specific" 

mgb_CHD_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/PartnersBiobank/CHD/I9_CHD_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv", data.table=FALSE)
mgb_CHD_f_bs$Sex <- "Female"
mgb_CHD_f_bs$Country <- "Massachusetts (USA)"
mgb_CHD_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
mgb_CHD_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_MaleLifetimeRisk_Bootstrapped_MetaAnalysis_Massachusetts.csv", data.table=FALSE)
mgb_CHD_m_ma$Sex <- "Male"
mgb_CHD_m_ma$Country <- "Massachusetts (USA)"
mgb_CHD_m_ma$Test <- "Meta-Analysis" 

mgb_CHD_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_FemaleLifetimeRisk_Bootstrapped_MetaAnalysis_Massachusetts.csv", data.table=FALSE)
mgb_CHD_f_ma$Sex <- "Female"
mgb_CHD_f_ma$Country <- "Massachusetts (USA)"
mgb_CHD_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#Genomics England
uk_CHD_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/CHD/I9_CHD_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv", data.table=FALSE)
uk_CHD_m_bs$Sex <- "Male"
uk_CHD_m_bs$Country <- "United Kingdom"
uk_CHD_m_bs$Test <- "Biobank Specific" 

uk_CHD_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/UKMetaAnalysis/CHD/I9_CHD_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv", data.table=FALSE)
uk_CHD_f_bs$Sex <- "Female"
uk_CHD_f_bs$Country <- "United Kingdom"
uk_CHD_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
uk_CHD_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_MaleLifetimeRisk_Bootstrapped_MetaAnalysis_United Kingdom.csv", data.table=FALSE)
uk_CHD_m_ma$Sex <- "Male"
uk_CHD_m_ma$Country <- "United Kingdom"
uk_CHD_m_ma$Test <- "Meta-Analysis" 

uk_CHD_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_FemaleLifetimeRisk_Bootstrapped_MetaAnalysis_United Kingdom.csv", data.table=FALSE)
uk_CHD_f_ma$Sex <- "Female"
uk_CHD_f_ma$Country <- "United Kingdom"
uk_CHD_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#FinnGen
finngen_CHD_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/CHD/I9_CHD_Male_LifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_CHD_m_bs$Sex <- "Male"
finngen_CHD_m_bs$Country <- "Finland"
finngen_CHD_m_bs$Test <- "Biobank Specific" 

finngen_CHD_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/FinnGen/CHD/I9_CHD_Female_LifetimeRisk_Bootstrapped_FinnGen.csv", data.table=FALSE)
finngen_CHD_f_bs$Sex <- "Female"
finngen_CHD_f_bs$Country <- "Finland"
finngen_CHD_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
finngen_CHD_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_MaleLifetimeRisk_Bootstrapped_MetaAnalysis_Finland.csv", data.table=FALSE)
finngen_CHD_m_ma$Sex <- "Male"
finngen_CHD_m_ma$Country <- "Finland"
finngen_CHD_m_ma$Test <- "Meta-Analysis" 

finngen_CHD_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_FemaleLifetimeRisk_Bootstrapped_MetaAnalysis_Finland.csv", data.table=FALSE)
finngen_CHD_f_ma$Sex <- "Female"
finngen_CHD_f_ma$Country <- "Finland"
finngen_CHD_f_ma$Test <- "Meta-Analysis" 

#Biobank Specific

#Estonian Biobank
estbb_CHD_m_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/CHD/I9_CHD_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_CHD_m_bs$Sex <- "Male"
estbb_CHD_m_bs$Country <- "Estonia"
estbb_CHD_m_bs$Test <- "Biobank Specific" 

estbb_CHD_f_bs <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/EstonianBiobank/CHD/I9_CHD_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv", data.table=FALSE)
estbb_CHD_f_bs$Sex <- "Female"
estbb_CHD_f_bs$Country <- "Estonia"
estbb_CHD_f_bs$Test <- "Biobank Specific" 

#Meta Analysis
estbb_CHD_m_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_MaleLifetimeRisk_Bootstrapped_MetaAnalysis_Estonia.csv", data.table=FALSE)
estbb_CHD_m_ma$Sex <- "Male"
estbb_CHD_m_ma$Country <- "Estonia"
estbb_CHD_m_ma$Test <- "Meta-Analysis" 

estbb_CHD_f_ma <- fread("/Users/jermy/Documents/INTERVENE/Results/GBD_Incidence/MetaAnalysis/CHD/CHD_FemaleLifetimeRisk_Bootstrapped_MetaAnalysis_Estonia.csv", data.table=FALSE)
estbb_CHD_f_ma$Sex <- "Female"
estbb_CHD_f_ma$Country <- "Estonia"
estbb_CHD_f_ma$Test <- "Meta-Analysis" 

#All
all <- rbind(mgb_CHD_m_bs, mgb_CHD_f_bs) %>%
        rbind(., mgb_CHD_m_ma) %>%
          rbind(., mgb_CHD_f_ma) %>%
            rbind(., uk_CHD_m_bs) %>%
              rbind(., uk_CHD_m_ma) %>%
                rbind(., uk_CHD_f_bs) %>%
                  rbind(., uk_CHD_f_ma) %>%
                    rbind(., finngen_CHD_m_bs) %>%
                      rbind(., finngen_CHD_m_ma) %>%
                        rbind(., finngen_CHD_f_bs) %>%
                          rbind(., finngen_CHD_f_ma) %>%
                            rbind(., estbb_CHD_m_bs) %>%
                              rbind(., estbb_CHD_m_ma) %>%
                                rbind(., estbb_CHD_f_bs) %>%
                                  rbind(., estbb_CHD_f_ma)

groupsy <- c("< 20%", "> 95%")
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 95%", "< 20%"))
all$Test <- factor(all$Test, levels=c("Biobank Specific", "Meta-Analysis"))

ggplot(all, aes(age, LifetimeRisk, color=interaction(Group,Test), group=interaction(Group,Test))) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=interaction(Group,Test)), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  theme_bw() +
  facet_grid(Country ~ Sex) + 
  labs(color='PGS Strata', fill='PGS Strata') +
  scale_color_manual(values=c("lightblue", "red", "darkblue", "darkred"), labels=c("Biobank Specific (95-100%)", "Biobank Specific (0-20%)", "Meta-Analysis (95-100%)", "Meta-Analysis (0-20%)")) +
  scale_fill_manual(values=c("lightblue", "red", "darkblue", "darkred"), labels=c("Biobank Specific (95-100%)", "Biobank Specific (0-20%)", "Meta-Analysis (95-100%)", "Meta-Analysis (0-20%)")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure_20b_CHD_MetaAnalysisvsBiobankSpecific.png"), height=8, width=10, dpi=300)

