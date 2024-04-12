inp<-"~/Documents/MyStuff/NTNU/Research/INTERVENE/INTERVENE_flagship_zips_from_bradley_computer/"
out<-"~/Documents/MyStuff/NTNU/Research/INTERVENE/INTERVENE_flagship_zips_from_bradley_computer/Figures/"

library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)
library(ggtext) 


################# FIGURE 1a ################

# Main analysis
meta <- fread(paste0(inp,"MetaAnalysis/SexStratified/HRperSD/SexStratifiedMetaAnalysisFEandRE.csv"), data.table=FALSE)
meta_eur <- subset(meta, Test=="Fixed Effect" & Ancestry=="EUR" & Phenotype!="ILD")
meta_eur2b <- meta_eur[,c("Ancestry","Sex","Phenotype","Beta","SE","Pval","QHet","HetPval","HR","Cipos","Cineg")]
#write.csv(meta_eur, "/Users/jermy/Documents/INTERVENE/Write-up/Supplementary Tables/SupplementaryTable7b.csv")

#Join with Prostate Cancer and Breast Cancer to create one panel
meta_full <- fread(paste0(inp,"MetaAnalysis/FullSample/HRperSD/metaanalysisFEandRE.csv"), data.table=FALSE)
meta_full <- subset(meta_full, Test=="Fixed Effect" & Ancestry=="EUR")
bcpc <- subset(meta_full, Test=="Fixed Effect" & Ancestry=="EUR" & (Phenotype=="C3_BREAST" | Phenotype=="C3_PROSTATE"))
bcpc[bcpc$Phenotype=="C3_PROSTATE", "Sex"] <- "Males" 
bcpc[bcpc$Phenotype=="C3_BREAST", "Sex"] <- "Females" 
bcpc <- bcpc[,colnames(meta_eur2b)]

meta_eur2b <- rbind(meta_eur2b, bcpc)
orderPheno <- subset(meta_full, Ancestry=="EUR")
meta_eur2b$Phenotype <- factor(meta_eur2b$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))

#mark what is significantly different by sex
meta_eur2b<-meta_eur2b %>% mutate(signif=case_when(Sex=="Females" & (Phenotype=="COX_ARTHROSIS" | Phenotype=="I9_CHD" | Phenotype=="T2D" | Phenotype=="GOUT" |  Phenotype=="J10_ASTHMA")~"*"))



#Plot
pdf(file=paste0(out,"Figure1a.pdf"),useDingbats = FALSE,height=8,width=8)
ggplot(meta_eur2b,aes(Phenotype, HR, group=Sex, col=Sex,label=signif)) +
  geom_point(position=position_dodge(width=0.7)) +
  geom_text(position=position_dodge(width=0.7),vjust=-0.5,color="black") +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  ylab("Meta-Analyzed Hazard Ratio per Standard Deviation (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#D55E00", "#56B4E9")) +
  scale_x_discrete(labels=c("Appendicitis  
                            <span style='font-size:12pt'>*Total Cases=48,239*</span>",
                            "Epilepsy  
                            <span style='font-size:12pt'>*Total Cases=27,882*</span>",
                            "All Cancers  
                            <span style='font-size:12pt'>*Total Cases=129,563*</span>",
                            "Major Depression  
                            <span style='font-size:12pt'>*Total Cases=119,419*</span>",
                            "Lung Cancer  
                            <span style='font-size:12pt'>*Total Cases=8,703*</span>",
                            "Skin Melanoma  
                            <span style='font-size:12pt'>*Total Cases=11,473*</span>",
                            "Knee Osteoarthritis  
                            <span style='font-size:12pt'>*Total Cases=100,324*</span>",
                            "Hip Osteoarthritis  
                            <span style='font-size:12pt'>*Total Cases=55,383*</span>",
                            "Coronary Heart Disease  
                            <span style='font-size:12pt'>*Total Cases=64,261*</span>",
                            "Asthma  
                            <span style='font-size:12pt'>*Total Cases=85,543*</span>",
                            "Colorectal Cancer  
                            <span style='font-size:12pt'>*Total Cases=13,009*</span>",
                            "Atrial Fibrillation  
                            <span style='font-size:12pt'>*Total Cases=54,434*</span>",
                            "Breast Cancer  
                            <span style='font-size:12pt'>*Total Cases=42,843*</span>",
                            "Rheumatoid Arthritis  
                            <span style='font-size:12pt'>*Total Cases=14,497*</span>",
                            "Gout  
                            <span style='font-size:12pt'>*Total Cases=30,503*</span>",
                            "Type 2 Diabetes  
                            <span style='font-size:12pt'>*Total Cases=80,821*</span>",
                            "Prostate Cancer  
                            <span style='font-size:12pt'>*Total Cases=32,876*</span>",
                            "Type 1 Diabetes  
                            <span style='font-size:12pt'>*Total Cases=6,719*</span>")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_markdown(size = 14, hjust=0.5)) +
  coord_flip()
dev.off()


####### Figure 1b #########

# Main analysis
meta <- fread(paste0(inp,"/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv"), data.table=FALSE)
meta_eur <- subset(meta, Ancestry=="EUR" & Phenotype!="ILD")
#write.csv(meta_eur, "/Users/jermy/Documents/INTERVENE/Write-up/Supplementary Tables/SupplementaryTable7c.csv")
meta_eur$Phenotype <- factor(meta_eur$Phenotype, levels=c("K11_APPENDACUT","G6_EPLEPSY","C3_CANCER", "F5_DEPRESSIO", "C3_BRONCHUS_LUNG", "C3_MELANOMA_SKIN", "KNEE_ARTHROSIS", "COX_ARTHROSIS",
                                                          "I9_CHD", "J10_ASTHMA", "C3_COLORECTAL", "I9_AF", "C3_BREAST","RHEUMA_SEROPOS_OTH", "GOUT", "T2D", "C3_PROSTATE", "T1D"))
meta_eur$Quartile <- factor(meta_eur$Quartile, levels=c(1,2,3,4))

meta_eur2c <- subset(meta_eur, Test=="Fixed Effect")
meta_eur2c <- meta_eur2c[,-1]

pdf(paste0(out,"Figure1b.pdf"),useDingbats = FALSE,height=8,width=8)
ggplot(meta_eur2c) +
  geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Meta-Analyzed Hazard Ratio per Standard Deviation (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  #guides(color = guide_legend(reverse = TRUE)) + 
  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", "2nd Quartile", "3rd Quartile", "Oldest")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Major Depression","Lung Cancer","Skin Melanoma","Knee Osteoarthritis","Hip Osteoarthritis","Coronary Heart Disease","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
dev.off()

############# Figure 2 ############

#Plot Hazard Ratios similar to original figure but with hazard ratios percentiles top and bottom and a log scale
metas <- fread(paste0(inp,"MetaAnalysis/HRPercentiles/MetaAnalysisPercentilesAgeandSexStratified.csv"), data.table=FALSE)

#More extensive subsetting required here...
#metas <- subset(metas, Phenotype %in% c("J10_ASTHMA","I9_AF","T2D") | (Phenotype=="I9_CHD" & Sex=="male") | (Phenotype=="GOUT" & Sex=="male") | (Phenotype=="KNEE_ARTHROSIS" & Sex=="female"))
metas <- subset(metas, Phenotype=="T2D"|Phenotype=="I9_CHD")

#Read in prostate cancer for the main plot
prostate <- fread(paste0(inp,"MetaAnalysis/HRPercentiles/MetaAnalysisPercentilesAgeStratified.csv"), data.table=FALSE)
prostate <- subset(prostate, Phenotype=="C3_PROSTATE" | Phenotype=="C3_BREAST")
prostate$Sex <- case_when(prostate$Phenotype=="C3_PROSTATE" ~ "male",
                          prostate$Phenotype=="C3_BREAST" ~ "female")
prostate <- prostate[,names(metas)]

metas <- rbind(metas,prostate)
#metas$Phenotype <- factor(metas$Phenotype, levels=c())
metas$Quartile <- factor(metas$Quartile, levels=c(1,2,3,4))
metas$Group <- factor(metas$Group, levels=c("< 20%","20-40%","40-60%","60-80%","80-90%","90-95%","> 95%"))

#Main Plot 
pdf(file=paste0(out,"Figure2.pdf"),height=13,width=13,useDingbats = FALSE)
ggplot(metas) +
  geom_point(aes(interaction(Sex,Phenotype), HR, col=Group, shape=Quartile, group=Quartile), position=position_dodge(width=0.7), size=3) +
  theme_bw() +
  geom_errorbar(aes(x=interaction(Sex,Phenotype), ymin = Cineg, ymax = Cipos, col=Group, group=Quartile), position=position_dodge(width=0.7), size=0.75, width=0.125) +
  ylab("Meta-Analyzed Hazard Ratio (95% CI)") +
  xlab("") +
  geom_vline(xintercept = seq(2.5,(length(unique(metas$Phenotype))*2)+0.5,by=2), col='black') + 
  geom_hline(yintercept = 1.0) +
  scale_shape_manual(name="Age Quartiles", values=c(15,16,17,8), labels=c("Youngest", "2nd Quartile", "3rd Quartile", "Oldest")) +
  labs(color="PGS Strata") +
  scale_x_discrete(labels=c("Breast Cancer  
                            <span style='font-size:12pt'>(Females)</span>",
                            "Prostate Cancer  
                            <span style='font-size:12pt'>(Males)</span>",
                            "Coronary Heart Disease  
                            <span style='font-size:12pt'>(Females)</span>",
                            "Coronary Heart Disease  
                            <span style='font-size:12pt'>(Males)</span>",
                            "Type 2 Diabetes  
                            <span style='font-size:12pt'>(Females)</span>",
                            "Type 2 Diabetes  
                            <span style='font-size:12pt'>(Males)</span>")) + 
  scale_y_continuous(trans="log2") +             
  theme(title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_markdown(size = 16, hjust=0.5)) +
  coord_flip()
dev.off()

#scale_color_manual(values=rev(hcl.colors(n=6,palette="Red-Blue")))

######## Figure 3a #######
### Figure 3 code from SupplementaryFigure12_LRE_Bootstrap_AllPhenos.R ####

#Prostate Cancer
directory <- c("ProstateCancer")
phenotype <- c("C3_PROSTATE")

for(i in 1){
  
  mgb_prostate_m <- fread(paste0(inp,"/GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_prostate_m$Sex <- "Male"
  mgb_prostate_m$Country <- "Massachusetts (USA)"
  
  #HUNT
  hunt_prostate_m <- fread(paste0(inp,"GBD_Incidence/HUNT/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_HUNT.csv"), data.table=FALSE)
  hunt_prostate_m$Sex <- "Male"
  hunt_prostate_m$Country <- "Norway"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_prostate_m <- fread(paste0(inp,"GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_prostate_m$Sex <- "Male"
  uk_prostate_m$Country <- "United Kingdom"
  
  #FinnGen
  finngen_prostate_m <- fread(paste0(inp,"GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_prostate_m$Sex <- "Male"
  finngen_prostate_m$Country <- "Finland"
  
  #Estonian Biobank
  estbb_prostate_m <- fread(paste0(inp,"GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
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
  
  pdf(file=paste0(out,"Figure_3a.pdf"),useDingbats=FALSE, height=7,width=7)
  print(ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
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
          axis.text.y = element_text(size = 12)))
  #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Boostrap.png"), height=7, width=7, dpi=300)
  dev.off()
}

######## Figure 3b #######

#Coronary Heart Disease

directory <- c("CHD")
phenotype <- c("I9_CHD")

for(i in 1){
  
  #MGB
  mgb_chd_m <- fread(paste0(inp,"GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_chd_m$Sex <- "Male"
  mgb_chd_m$Country <- "Massachusetts (USA)"
  
  mgb_chd_f <- fread(paste0(inp,"GBD_Incidence/PartnersBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_PartnersBiobank.csv"), data.table=FALSE)
  mgb_chd_f$Sex <- "Female"
  mgb_chd_f$Country <- "Massachusetts (USA)"
  
  #UK Meta Analysis - UK Biobank, Genomics England and Generation Scotland
  uk_chd_m <- fread(paste0(inp,"GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_chd_m$Sex <- "Male"
  uk_chd_m$Country <- "United Kingdom"
  
  uk_chd_f <- fread(paste0(inp,"GBD_Incidence/UKMetaAnalysis/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_UKMetaAnalysis.csv"), data.table=FALSE)
  uk_chd_f$Sex <- "Female"
  uk_chd_f$Country <- "United Kingdom"
  
  #FinnGen
  finngen_chd_m <- fread(paste0(inp,"GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_chd_m$Sex <- "Male"
  finngen_chd_m$Country <- "Finland"
  
  finngen_chd_f <- fread(paste0(inp,"GBD_Incidence/FinnGen/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_FinnGen.csv"), data.table=FALSE)
  finngen_chd_f$Sex <- "Female"
  finngen_chd_f$Country <- "Finland"
  
  #Estonian Biobank
  estbb_chd_m <- fread(paste0(inp,"GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Male_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
  estbb_chd_m$Sex <- "Male"
  estbb_chd_m$Country <- "Estonia"
  
  estbb_chd_f <- fread(paste0(inp,"GBD_Incidence/EstonianBiobank/",directory[i],"/",phenotype[i],"_Female_LifetimeRisk_Bootstrapped_EstonianBiobank.csv"), data.table=FALSE)
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
  
  pdf(file=paste0(out,"Figure_3b.pdf"),useDingbats=FALSE, height=10,width=6)
  print(ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
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
          axis.text.y = element_text(size = 12)))
  #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Bootstrap.png"), height=10, width=6, dpi=300)
  dev.off()
}


########## Figure 4a ############
##T2D
#code from Figure3_T2D_LRE.R

#Bootstrapped
finngen_t2d_m <- fread(paste0(inp,"GBD_Incidence/FinnGen/T2D/T2D_MaleLifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
finngen_t2d_m$Sex <- "Male"
finngen_t2d_m$Country <- "Finland"
finngen_t2d_m$Threshold <- 7.380580738985

finngen_t2d_f <- fread(paste0(inp,"GBD_Incidence/FinnGen/T2D/T2D_FemaleLifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
finngen_t2d_f$Sex <- "Female"
finngen_t2d_f$Country <- "Finland"
finngen_t2d_f$Threshold <- 6.40599016524789

groupsy <- c("< 1%", "40-60%", "> 99%")
all <- rbind(finngen_t2d_m, finngen_t2d_f)
all <- subset(all, Group %in% groupsy)

all$Group <- factor(all$Group, levels=c("> 99%", "40-60%", "< 1%"))

pdf(file=paste0(out,"Figure4a.pdf"),height=6,width=6, useDingbats=FALSE)
ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Sex) + 
  labs(color='PGS Strata', fill='PGS Strata',title="Type 2 Diabetes") +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
dev.off()
#ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFigure10b.png"), height=6, width=6, dpi=300)

###### Figure 4b #####
##breast cancer
###code from Figure4_BreCan_LRE.R

finngen_f <- fread(paste0(inp,"GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
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
  labs(color='PGS Strata',title="Breast Cancer") +
  scale_color_hue(labels = c("99-100%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
groupsy <- c("< 1%", "40-60%", "> 99%")
finngen_f <- subset(finngen_f, Group %in% groupsy)

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "40-60%", "< 1%"))

pdf(file=paste0(out,"Figure4b.pdf"),height=6,width=6, useDingbats=FALSE)
ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  labs(color='PGS Strata', fill='PGS Strata',title="Breast Cancer") +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
dev.off()
#ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/BreCanOnePercentTailswCI.png"), height=6, width=6, dpi=300)

