inp<-"~/Documents/MyStuff/NTNU/Research/INTERVENE/INTERVENE_flagship_zips_from_bradley_computer/"
out<-"~/Documents/MyStuff/NTNU/Research/INTERVENE/INTERVENE_flagship_zips_from_bradley_computer/Figures/"

library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)
library(ggtext) 
library(patchwork)
library(lemon) #reposition_legend()


################# FIGURE 1 ################

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
a<-ggplot(meta_eur2b,aes(Phenotype, HR, group=Sex, col=Sex,label=signif)) +
  geom_point(position=position_dodge(width=0.7),size=2) +
  geom_text(position=position_dodge(width=0.7),vjust=-0.5,color="black") +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Sex, col=Sex), size=0.7, width=0.4, position=position_dodge(width=0.7)) +
  ylab("Meta-Analyzed Hazard Ratio\nper Standard Deviation (95% CI)") +
  xlab("") + 
  guides(color = guide_legend(nrow = 2))+
  geom_hline(yintercept = 1.0) +
  scale_color_manual(name="Sex",values=c("#D55E00", "#56B4E9")) +
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
        legend.position="bottom",
        legend.justification = "left",
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_markdown(size = 14, hjust=0.5)) +
  labs(tag="a") +
  coord_flip()


#1b

# Main analysis
meta <- fread(paste0(inp,"/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv"), data.table=FALSE)
meta_eur <- subset(meta, Ancestry=="EUR" & Phenotype!="ILD")
#write.csv(meta_eur, "/Users/jermy/Documents/INTERVENE/Write-up/Supplementary Tables/SupplementaryTable7c.csv")
meta_eur$Phenotype <- factor(meta_eur$Phenotype, levels=c("K11_APPENDACUT","G6_EPLEPSY","C3_CANCER", "F5_DEPRESSIO", "C3_BRONCHUS_LUNG", "C3_MELANOMA_SKIN", "KNEE_ARTHROSIS", "COX_ARTHROSIS",
                                                          "I9_CHD", "J10_ASTHMA", "C3_COLORECTAL", "I9_AF", "C3_BREAST","RHEUMA_SEROPOS_OTH", "GOUT", "T2D", "C3_PROSTATE", "T1D"))
meta_eur$Quartile <- factor(meta_eur$Quartile, levels=c(1,2,3,4))

meta_eur2c <- subset(meta_eur, Test=="Fixed Effect")
meta_eur2c <- meta_eur2c[,-1]


b<-ggplot(meta_eur2c) +
  geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7),size=2) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.7, width=0.7, position=position_dodge(width=0.7)) +
  ylab("Meta-Analyzed Hazard Ratio\nper Standard Deviation (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  #guides(color = guide_legend(reverse = TRUE)) + 
  scale_color_manual(name="Age\nQuartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", "2nd Quartile", "3rd Quartile", "Oldest")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Major Depression","Lung Cancer","Skin Melanoma","Knee Osteoarthritis","Hip Osteoarthritis","Coronary Heart Disease","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  guides(color = guide_legend(nrow = 4)) +
  theme(title = element_text(size = 18),
        legend.position="bottom",
        legend.justification = "left",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_blank()) +
  coord_flip() + labs(tag="b")

pdf(file=paste0(out,"Figure1.pdf"),useDingbats = FALSE,height=11.7,width=8.3)
a+b + plot_layout(axis_titles="collect")
dev.off()


######################### FIGURE 3 version 1 #####################


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
  

  a<-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
          stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
          geom_point() +
          geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
          xlab("Age (Years)") + 
          ylab("Cumulative Risk (%)") + 
          facet_wrap( ~ Country ) +
          theme_bw() +  labs(tag="a") +
          labs(color='PGS Strata', fill='PGS Strata') +
          scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
          scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
          theme(legend.text = element_text(size = 12),
                legend.position.inside=c(0,1),
                strip.text=element_text(size=14),
                strip.background =element_rect(fill="white"),
                legend.title = element_text(size = 14),
                axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                axis.title.y = element_text(size = 14),
                axis.text.y = element_text(size = 12))
  
  #lemon::reposition_legend(a,'center',panel='panel-3-2')
   #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Boostrap.png"), height=7, width=7, dpi=300)

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
  
  
  b<-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
          stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
          geom_point() +
          geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
          xlab("Age (Years)") + 
          ylab("Cumulative Risk (%)") + 
          theme_bw() +  labs(tag="b") +
          facet_grid(Country ~ Sex, labeller = label_wrap_gen(width=8,multi_line = TRUE)) +
          labs(color='PGS Strata', fill='PGS Strata') +
          scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%"),guide="none") +
          scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%"),guide="none") +
          theme(legend.position="none",
                strip.text=element_text(size=13),
                strip.background =element_rect(fill="white"),
                axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                axis.title.y = element_text(size = 14),
                axis.text.y = element_text(size = 12)) 
  
}


pdf(file=paste0(out,"Figure3.pdf"),useDingbats=FALSE,height=11.7,width=8.3)
a/b/guide_area() + plot_layout(guides="collect",heights=c(0.35,0.6,0.05)) & theme(legend.position = 'bottom')
dev.off()

##### 

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
  
  
  a<-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age (Years)") + 
    ylab("Cumulative Risk (%)") + 
    facet_wrap( ~ Country,ncol=1,strip.position="right") +
    theme_bw() +  labs(tag="a") +
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%")) +
    theme(legend.text = element_text(size = 12),
          legend.position.inside=c(0,1),
          strip.text=element_text(size=14),
          strip.background =element_rect(fill="white"),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
  
  #lemon::reposition_legend(a,'center',panel='panel-3-2')
  #ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/",phenotype[i],"_LRE_Boostrap.png"), height=7, width=7, dpi=300)
  
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
  
  
  b<-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
    geom_point() +
    geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
    xlab("Age (Years)") + 
    ylab("Cumulative Risk (%)") + 
    theme_bw() +  labs(tag="b") +
    facet_grid(Country ~ Sex) +
    labs(color='PGS Strata', fill='PGS Strata') +
    scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%"),guide="none") +
    scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("95-100%", "40-60%", "0-20%"),guide="none") +
    theme(legend.position="none",
          strip.text=element_text(size=13),
          strip.background =element_rect(fill="white"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12)) 
  
}


pdf(file=paste0(out,"Figure3.pdf"),useDingbats=FALSE,height=11.7,width=8.3)
(a + b )/guide_area() + plot_layout(guides="collect",heights=c(0.95,0.05)) & theme(legend.position = 'bottom')
dev.off()



######################### FIGURE 4 #####################


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

a<-ggplot(all, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age (Years)") + 
  ylab("Cumulative Risk (%)") + 
  geom_hline(data=all, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  theme_bw() +
  facet_wrap(~ Sex) + labs(tag="a") +
  labs(color='PGS Strata', fill='PGS Strata',title="Type 2 Diabetes") +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  theme(legend.text = element_text(size = 12),
        legend.position="none",
        strip.text=element_text(size=14),
        strip.background =element_rect(fill="white"),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))



###### Figure 4b #####
##breast cancer

finngen_f <- fread(paste0(inp,"GBD_Incidence/FinnGen/BreastCancer/C3_BREAST_Female_LifetimeRisk_Bootstrapped_FinnGen_OnePercentTails.csv"), data.table=FALSE)
finngen_f$Country <- "Finland"
finngen_f$Threshold <- 1.70678589956649

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "95-99%", "90-95%", "80-90%", "60-80%", "40-60%", "20-40%", "10-20%", "5-10%", "1-5%", "< 1%"))

finngen_f <- subset(finngen_f, Group %in% groupsy)

finngen_f$Group <- factor(finngen_f$Group, levels=c("> 99%", "40-60%", "< 1%"))


b<-ggplot(finngen_f, aes(age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 8), se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age (Years)") + 
  ylab("Cumulative Risk (%)") + 
  guides(colo=guide_legend(ncol=2)) +
  geom_hline(data=finngen_f, aes(yintercept=Threshold), linetype="dashed", color="red", size=1) + 
  labs(color='PGS Strata', fill='PGS Strata',title="Breast Cancer") +
  labs(tag="b") +
  scale_color_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) +
  scale_fill_manual(values=c("#88CCEE", "#000000", "#CC6677"), labels = c("99-100%", "40-60%", "0-1%")) + theme_bw() +
  theme(legend.text = element_text(size = 12),
        legend.position="bottom",
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))


pdf(file=paste0(out,"Figure4.pdf"),useDingbats = FALSE,height=11.7,width=8.3)
a+b + plot_layout(heights = c(7,7,1))
dev.off()




