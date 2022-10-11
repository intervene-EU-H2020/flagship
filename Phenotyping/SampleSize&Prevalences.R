library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
#library(ggpubr)
library(rcartocolor)

#color <- brewer.pal(n = 9, name = 'Paired')
#color <- color[-1]

#Read in full sample hazard ratios per standard deviation
bbj<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/Biobank_Japan_HazardRatios/HRperSD_BBJ.csv")
ukb<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/UKB_HazardRatios/PRS_HRsperSD_UKBiobank_AllAncestries.csv")
finngen<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/FinnGen_HazardRatios/HRperSD_FinnGen.csv")
gnh<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/GNH_HazardRatios/HRperSD_GNH.csv")
estbb<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/HRperSD_EstBB.csv")
gs<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/GenerationScotland_HazardRatios/HRperSD_GS.csv")
ge<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/GenomicsEngland_HazardRatios/HRperSD_GenomicsEngland.csv")
mgb_eur<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/MGB_HazardRatios/HRperSD_MGBB_AFR.csv")
mgb_afr<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/MGB_HazardRatios/HRperSD_MGBB_AFR.csv")
hunt<-fread("/mnt/work/workbench/bwolford/intervene/GoogleDrive/HUNT_HazardRatios/HRperSD_HUNT.csv")


output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/"
#Biobank Japan
#bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HRperSD_BBJ.csv", data.table=FALSE)
dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
bbj <- subset(bbj, !(Phenotype %in% dropbbj))
bbj$Biobank <- "Biobank Japan"
bbj$Ancestry <- "EAS"

#UKB
#ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/PRS_HRsperSD_UKBiobank_AllAncestries.csv", data.table=FALSE)
ukb$Biobank <- "UK Biobank"
ukb <- subset(ukb, (Phenotype %in% bbj$Phenotype & Ancestry=="EAS") | Ancestry=="EUR" | Ancestry=="SAS" | Ancestry=="AFR")
#ukb <- ukb[,names(bbj)]

#FinnGen
#finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HRperSD_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"
#finngen <- finngen[,names(bbj)]

#Genes&Health
#gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HRperSD_GNH.csv", data.table=FALSE)
gnh$Biobank <- "Genes & Health"
gnh$Ancestry <- "SAS"
#gnh <- gnh[,names(bbj)]

#Estonian Biobank
#estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HRperSD_EstBB.csv", data.table=FALSE)
estbb$Biobank <- "Estonian Biobank"
estbb$Ancestry <- "EUR"
#estbb <- estbb[,names(bbj)]

#Generation Scotland
#gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HRperSD_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs$Biobank <- "Generation Scotland"
gs$Ancestry <- "EUR"
#gs <- gs[,names(bbj)]

#Mass General Brigham
#mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_MGBB_EUR.csv", data.table=FALSE)
dropmgb <- c("I9_AF")
mgb_eur <- subset(mgb_eur, !(Phenotype %in% dropmgb))
mgb_eur$Biobank <- "Mass General Brigham"
mgb_eur$Ancestry <- "EUR"
#mgb_eur <- mgb_eur[,names(bbj)]

mgb_afr<- subset(mgb_afr, !(Phenotype %in% dropmgb))
mgb_afr$Biobank <- "Mass General Brigham"
mgb_afr$Ancestry <- "AFR"
#mgb_afr <- mgb_afr[,names(bbj)]

#Genomics England
#ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HRperSD_GenomicsEngland.csv", data.table=FALSE)
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
#ge <- ge[,names(bbj)]

#HUNT
drophunt <- c("I9_CHD","T1D")
hunt$Biobank <- "HUNT"
hunt$Ancestry <- "EUR"


#Combine into one dataset
all <- rbind(ukb, finngen,fill=TRUE) %>%
        rbind(., gnh,fill=TRUE) %>% 
          rbind(., estbb,fill=TRUE) %>%
            rbind(., gs,fill=TRUE) %>%
              rbind(., mgb_eur,fill=TRUE) %>% 
              rbind(., mgb_afr,fill=TRUE) %>% 
                rbind(., bbj,fill=TRUE) %>% 
                  rbind(., ge,fill=TRUE) %>% rbind(.,hunt,fill=TRUE)

full <- subset(all, Sample=="Full Sample")

##################################################################################################################################################################################
##################################################################################################################################################################################

#For Sample Size graphs (Supplementary Figure 2)
sample <- full[,c("Phenotype","Cases","Ancestry","Biobank")]
sample$Cases <- sample$Cases/1000

summies <- full %>% 
  group_by(Phenotype) %>%
  summarize(count=sum(Cases))

max<- full %>%
  group_by(Phenotype,Biobank) %>% summarize(total_case=sum(Cases),total_control=sum(Controls)) %>% group_by(Biobank) %>%
  summarize(total=max(total_case+total_control))


#group_by(Category) %>% 
 # summarise(Frequency = sum(Frequency))

summies <- summies[order(summies$count),]
summies$Phenotype <- as.factor(summies$Phenotype) 

sample$Phenotype <- factor(sample$Phenotype, levels=summies$Phenotype)
sample$Ancestry <- factor(sample$Ancestry, levels=c("AFR","SAS","EAS","EUR"))
sample$Biobank <- factor(sample$Biobank, levels=c("Generation Scotland","HUNT","Genes & Health","Genomics England","Mass General Brigham","UK Biobank","Biobank Japan","Estonian Biobank","FinnGen"))
color<-brewer.pal("Dark2",n=length(unique(sample$Biobank)))
caseno <- ggplot(data=sample, aes(x=Phenotype, y=Cases, fill=Biobank)) + 
            geom_bar(stat="identity") + 
              theme_bw() + 
                scale_fill_manual(values=rev(color)) + 
                  xlab("") + 
                    ylab("Number of Cases (per 1000)") + 
                      theme(legend.text = element_text(size = 28),
                            legend.title = element_blank(),
                            legend.spacing.y = unit(2.0, 'cm'),
                            axis.title.x = element_text(size = 24),
                            axis.text.x = element_text(size = 20),
                            axis.title.y = element_text(size = 20),
                            axis.text.y = element_text(size = 26)) +
                      guides(fill = guide_legend(reverse=TRUE, byrow = TRUE)) + 
                          scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Lung Cancer","Skin Melanoma","Rheumatoid Arthritis","Colorectal Cancer","Epilepsy","Gout","Prostate Cancer","Breast Cancer","Appendicitis","Hip Osteoarthritis","Atrial Fibrillation","CHD","Asthma","Knee Osteoarthritis","Major Depression","Type 2 Diabetes","All Cancers")) + 
                            coord_flip()
ggsave(filename=paste0(output_dir,"prevalences.png"), plot = caseno, height=10, width = 16, dpi=300)

ancno <- ggplot(data=sample, aes(x=Phenotype, y=Cases, fill=Ancestry)) + 
          geom_bar(stat="identity") + 
            theme_bw() + 
              scale_fill_manual(values=brewer.pal("Dark2",n=4)) + 
                xlab("") + 
                  ylab("Number of Cases (per 1000)") + 
  theme( axis.title.x = element_text(size = 18),
         axis.text.x = element_text(size = 14),
         axis.title.y = element_text(size = 18),
         axis.text.y = element_text(size = 18)) +
                    scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Lung Cancer","Skin Melanoma","Rheumatoid Arthritis","Colorectal Cancer","Epilepsy","Gout","Prostate Cancer","Atrial Fibrillation","Breast Cancer","Appendicitis","Hip Osteoarthritis","CHD","Asthma","Knee Osteoarthritis","Type 2 Diabetes","Major Depression","All Cancers")) + 
                      coord_flip()
ggsave(filename=paste0(output_dir,"ancestry.png"), ancno, height=10, width = 16, dpi=300)

eur <- ggplot(data=sample[sample$Ancestry=="EUR",], aes(x=Phenotype, y=Cases, fill=Ancestry)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values="#E7298A") + 
  xlab("") + 
  ylab("Number of Cases (per 1000)") + 
  theme( axis.title.x = element_text(size = 18),
         axis.text.x = element_text(size = 14),
         axis.title.y = element_text(size = 18),
         axis.text.y = element_text(size = 18)) +
  scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Lung Cancer","Skin Melanoma","Rheumatoid Arthritis","Colorectal Cancer","Epilepsy","Gout","Prostate Cancer","Atrial Fibrillation","Breast Cancer","Appendicitis","Hip Osteoarthritis","CHD","Asthma","Knee Osteoarthritis","Type 2 Diabetes","Major Depression","All Cancers")) + 
  coord_flip()
ggsave(filename=paste0(output_dir,"ancestry_EUR.png"), eur, height=10, width = 16, dpi=300)


##################################################################################################################################################################################
##################################################################################################################################################################################

color <- brewer.pal(n = 10, name = 'Paired')

#For prevalence variation - supplementary figure 3
prev <- full[,c("Phenotype","Controls","Cases","Ancestry","Biobank")]
prev$Prevalence <- (prev$Cases/prev$Controls)*100

highprev <- prev %>% 
              group_by(Phenotype) %>%
                summarize(meanprev=mean(Prevalence))
highprev <- highprev[order(highprev$meanprev),]

highprev$Phenotype <- as.factor(highprev$Phenotype) 

prev$Phenotype <- factor(prev$Phenotype, levels=c(highprev$Phenotype))

#Change Biobank to split UKB into three different factors given differences in ancestry also
prev$BiobankUpd <- case_when(prev$Biobank=="UK Biobank" & prev$Ancestry=="EAS" ~ "UK Biobank (EAS)",
                             prev$Biobank=="UK Biobank" & prev$Ancestry=="EUR" ~ "UK Biobank (EUR)",
                             prev$Biobank=="UK Biobank" & prev$Ancestry=="SAS" ~ "UK Biobank (SAS)",
                             TRUE ~ prev$Biobank)

prev$BiobankUpd <- as.factor(prev$BiobankUpd)

#Absolute prevalence differences
prevplot <- ggplot(data=prev, aes(x=Phenotype, y=Prevalence)) +
                geom_point(aes(color=BiobankUpd), size=1) + 
                  geom_boxplot(alpha=0.2) +
                    theme_bw() + 
                      scale_color_manual(values=c(color)) + 
                        labs(color="Biobank") + 
                          scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Rheumatoid Arthritis","Epilepsy","Lung Cancer","Colorectal Cancer","Gout","Skin Melanoma","Prostate Cancer","Appendicitis","Atrial Fibrillation","Hip Osteoarthritis","Breast Cancer","CHD","Knee Osteoarthritis","Major Depression","Asthma","Type 2 Diabetes","All Cancers")) + 
                            xlab("") +
                              ylab("Prevalence (%)") + 
                                coord_flip()

#Relative differences
prev <- left_join(prev, highprev)
prev$RelPrev <- (prev$Prevalence/prev$meanprev)*100

relprevplot <- ggplot(data=prev, aes(x=Phenotype, y=RelPrev)) +
              geom_point(aes(color=BiobankUpd), size=1) + 
                geom_boxplot(alpha=0.2) +
                  theme_bw() + 
                    scale_color_manual(values=c(color)) + 
                      labs(color="Biobank") + 
                        geom_hline(yintercept=100, linetype="dashed", col="red") + 
                          scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Rheumatoid Arthritis","Epilepsy","Lung Cancer","Colorectal Cancer","Gout","Skin Melanoma","Prostate Cancer","Appendicitis","Atrial Fibrillation","Hip Osteoarthritis","Breast Cancer","CHD","Knee Osteoarthritis","Major Depression","Asthma","Type 2 Diabetes","All Cancers")) + 
                            xlab("") +
                              ylab("Relative Difference to Mean Prevalence (%)") + 
                                coord_flip()

#supp_fig_3 <- ggarrange(prevplot, relprevplot, labels=c("A","B"), align = "hv", nrow = 2, ncol = 1)
ggsave(filename="/Users/jermy/Documents/INTERVENE/Write-up/supplementary_figure_3.png", plot = supp_fig_3, height=10, width = 10, dpi=300)





