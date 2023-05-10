library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)

#Read in age stratified hazard ratios

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HRperSD_AgeandSexStratified_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HRperSD_AgeandSexStratified_UKB_AllAncestries.csv", data.table=FALSE)
ukb <- ukb[,names(finngen)]

#Genes&Health
gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HRperSD_AgeandSexStratified_GNH.csv", data.table=FALSE)
gnh$Biobank <- "Genes & Health"
gnh$Ancestry <- "SAS"

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HRperSD_AgeandSexStratified_BBJ.csv", data.table=FALSE)
dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
bbj <- subset(bbj, !(Phenotype %in% dropbbj))

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HRperSD_AgeandSexStratified_EstBB.csv", data.table=FALSE)

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HRperSD_AgeandSexStratified_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))

#Mass General Brigham - European
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeandSexStratified_MGBB_EUR.csv", data.table=FALSE)
mgb$Beta <- ifelse(mgb$Phenotype=="C3_CANCER", mgb$Beta*-1, mgb$Beta)
mgb$HR <- exp(mgb$Beta)
mgb$Cipos <- exp(mgb$Beta + 1.96*mgb$SE)
mgb$Cineg <- exp(mgb$Beta - 1.96*mgb$SE)
dropmgb <- c("I9_AF")
mgb <- subset(mgb, !(Phenotype %in% dropmgb))
mgb$Biobank <- "Mass General Brigham"
mgb$Ancestry <- "EUR"
mgb <- mgb[,names(ukb)]

#Mass General Brigham - African
#mgb_afr <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeandSexStratified_MGBB_AFR.csv", data.table=FALSE)
#dropmgb <- c("I9_AF")
#mgb_afr <- subset(mgb_afr, !(Phenotype %in% dropmgb))
#mgb_afr$Biobank <- "Mass General Brigham"
#mgb_afr$Ancestry <- "AFR"
#mgb_afr <- mgb_afr[,names(ukb)]

#Mass General Brigham
#mgb <- rbind(mgb_eur, mgb_afr)

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HRperSD_AgeandSexStratified_Genomics_England.csv", data.table=FALSE)
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
ge <- ge[,names(finngen)]

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/HUNT/HRperSD_AgeandSexStratified_HUNT.csv", data.table=FALSE)
drophunt <- c("I9_CHD")
hunt <- subset(hunt, !(Phenotype %in% drophunt))
hunt$Biobank <- "HUNT"
hunt$Ancestry <- "EUR"
hunt <- hunt[,names(ukb)]

#Combine into one dataset
all <- rbind(finngen, gnh) %>%
        rbind(., estbb) %>% 
          rbind(., gs) %>%
            rbind(., mgb) %>%
              rbind(., bbj) %>%
                rbind(., ukb) %>%
                  rbind(., ge) %>%
                    rbind(., hunt)

all$Quartile <- factor(all$Quartile, levels=c(1,2,3,4))

#Make sure all quartiles across both sexes have been tested before including in the meta-analysis
q <- c()
for(i in unique(all$Phenotype)){
  for(k in unique(all$Biobank)){
    for(j in unique(all$Ancestry)){
      pheno <- subset(all, Phenotype==i & Biobank==k & Ancestry==j)
      print(pheno)
      
      if(any(is.na(pheno[,c(9:14)])) | dim(pheno)[1] < 8){
        next
      }
      
      q <- rbind(q,pheno)
    }
  }
}

all <- as.data.frame(q)

metaresults <- c()
for(k in c(1,2,3,4)){
  
  for(l in c("male","female")){
    
    quartile <- subset(all, Quartile==k & Sex==l)
    
    for(j in c("AFR","EUR","SAS","EAS")){
      
      for(i in unique(quartile$Phenotype)){
        print(i)
        #Test with T2D phenotype and EUR ancestry
        if(j=="EAS"){
          disease <- subset(quartile, Phenotype==i & Phenotype %in% bbj$Phenotype & Ancestry==j & !(is.na(Beta)))
        } else{
          disease <- subset(quartile, Phenotype==i & Ancestry==j & !(is.na(Beta)))
        }
        
        if(dim(disease)[1]==0){
          next
        }
        
        #Meta analysis should be done at the beta level and stratified by ancestry 
        metaFE <- rma(yi=Beta, sei=SE, data=disease, method="FE")
        metaRE <- rma(yi=Beta, sei=SE, data=disease)
        
        #png(file=paste0("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/ForestPlots/",l,"_Quartile_",k,"_",i,"_",j,"_Forest.png"), 
        #    width = 7,
        #    height    = 7,
        #    units     = "in",
        #    res       = 300)
        #forest(meta, slab=disease$Biobank)
        #dev.off()
        
        results <- matrix(c("Fixed Effect", "Random Effect", j, j, k, k, l, l, i, i, metaFE$b, metaRE$b, metaFE$se, metaRE$se, metaFE$pval, metaRE$pval, metaFE$QE, metaRE$QE, metaFE$QEp, metaRE$QEp), nrow = 2, ncol = 10)
        
        metaresults <- rbind(metaresults, results)
        
      }
      
    }
  }
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Test","Ancestry","Quartile","Sex","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(6:10),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/AgeandSexStratifiedMetaAnalysisFEandRE.csv")

######################################################################################################################################

alld <- all[,c("Ancestry","Quartile","Sex","Phenotype","Beta","SE","Biobank")]
diffs <- metaresults[,c("Ancestry","Sex","Quartile","Phenotype","Beta","SE")]
diffs$Biobank <- "All"
diffs <- rbind(diffs, alld)

differences <- c()
for(j in c("EUR","SAS","EAS")){
  for(l in c("male","female")){
    for(k in c(1,2,3,4)){
      for(i in unique(diffs$Phenotype)){
        
        disease <- subset(diffs, Phenotype==i & Ancestry==j & Quartile==k & Sex==l)
        
        if(dim(disease)[1]==0){
          next
        }
        
        disease$BetaDiff <- ifelse(!(is.na(disease$Beta)), disease$Beta - disease$Beta[disease$Biobank=="All"], NA)
        disease$SEDiff <- ifelse(!(is.na(disease$Beta)), sqrt(disease$SE**2 + disease$SE[disease$Biobank=="All"]**2), NA)
        disease$ZDiff <- ifelse(!(is.na(disease$Beta)), disease$BetaDiff/disease$SEDiff, NA)
        disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
        differences <- rbind(differences, disease)
        
      }
    }
  }
}

fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/AgeandSex_BiobankVariation.csv")

######################################################################################################################################

# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/AgeandSexStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
meta_eur <- subset(meta, Ancestry=="EUR" & Phenotype!="ILD")
#write.csv(meta_eur, "/Users/jermy/Documents/INTERVENE/Write-up/Supplementary Tables/SupplementaryTable7d.csv")

orderPheno <- subset(meta_eur, Quartile==1 & Sex=="female" & Test=="Fixed Effect" & Phenotype!="ILD")
meta_eur$Phenotype <- factor(meta_eur$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))
meta_eur$Quartile <- factor(meta_eur$Quartile, levels=c(4,3,2,1))
meta_eur$HR <- as.numeric(meta_eur$HR)

#Subset to phenotypes which are going to be age and sex stratified
phenotypes <- c("I9_CHD","GOUT","J10_ASTHMA","T2D","I9_AF","KNEE_ARTHROSIS")
meta_eur_sub <- subset(meta_eur, Phenotype %in% phenotypes & Test=="Fixed Effect") 

#Plot
ggplot(meta_eur_sub) +
        geom_point(aes(Phenotype, HR, group=interaction(Quartile,Sex), col=Quartile, shape=Sex), position=position_dodge(width=0.7), size=3) +
        theme_bw() +
        geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=interaction(Quartile,Sex), col=Quartile), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
        ylab("Meta-Analyzed Hazard Ratio per Standard Deviation (95% CI)") +
        xlab("") +
        geom_hline(yintercept = 1.0) +
        scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", " |", "V", "Oldest")) +
        #guides(color = guide_legend(reverse = TRUE)) + 
        scale_shape_discrete(labels=c("Females","Males")) + 
        scale_x_discrete(labels=c("Knee Osteoarthritis","Coronary Heart Disease","Asthma","Gout","Atrial Fibrillation","Type 2 Diabetes")) + 
        theme(title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 14)) +
        coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/Figure_2d.png", height=8, width=9, dpi=300)

######################################################################################################################################
######################################################################################################################################

#Significant differences between meta-analysis estimates?

diffs <- meta[,c("Ancestry","Sex","Quartile","Phenotype","Beta","SE")]

differences <- c()
for(j in unique(diffs$Ancestry)){
  for(k in unique(diffs$Sex)){
    for(i in unique(diffs$Phenotype)){
      disease <- subset(diffs, Phenotype==i & Ancestry==j & Sex==k & (Quartile==1 | Quartile==4))
      if(dim(disease)[1]==0){
        next
      }
      
      disease$delta <- disease$Beta[1] - disease$Beta[2]
      disease$se_diff <- sqrt((disease$SE[1]**2) + (disease$SE[2]**2))
      disease$z <- disease$delta/disease$se_diff
      disease$p <- 2*pnorm(abs(disease$z), lower.tail=FALSE)
      differences <- rbind(differences, disease)
    }
  }
}
fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/AgeandSex_HRDifferences_MA.csv")

######################################################################################################################################

#Related to this but looking for heterogeneity statistics stratified by sex to determine if age-specific effects exist
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/AgeandSexStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
meta <- subset(meta, Test=="Fixed Effect")

metaresults <- c()
for(k in unique(meta$Ancestry)){
  for(j in unique(meta$Sex)){
    for(i in unique(meta$Phenotype)){
      disease <- subset(meta, Phenotype==i & Sex==j & Ancestry==k)
      
      if(nrow(disease)==0){
        next
      }
      
      metaFE <- rma(yi=Beta, sei=SE, data=disease, method="FE")
      results <- c(k, i, j, metaFE$b, metaFE$se, metaFE$pval, metaFE$QE, metaFE$QEp)
      
      metaresults <- rbind(metaresults, results)
    }
  }
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Ancestry","Phenotype","Sex","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratified/HeterogeneityinAgeandSexEffects_ModelSelection_AllAncestries.csv")

######################################################################################################################################
######################################################################################################################################

#Create forest plots for each pheno - facet_wrap
eur <- subset(all, Ancestry=="EUR" & Phenotype!="ILD")
eur <- eur[,c("Phenotype","Sex","Quartile","HR","Cipos","Cineg","Biobank")]

meta_eur$Biobank <- case_when(meta_eur$Test=="Fixed Effect" ~ "Meta Analysis FE",
                              meta_eur$Test=="Random Effect" ~ "Meta Analysis RE")

meta_eur <- meta_eur[,c("Phenotype","Sex","Quartile","HR","Cipos","Cineg","Biobank")]
eur <- rbind(eur, meta_eur)
eur$Biobank <- factor(eur$Biobank, levels=c("Meta Analysis RE", "Meta Analysis FE", "FinnGen", "Estonian Biobank", "UK Biobank", "HUNT", "Mass General Brigham","Genomics England", "Generation Scotland"))
eur$Quartile <- factor(eur$Quartile, levels=c(4,3,2,1))

eur_males <- subset(eur, Sex=="male")
eur_males$Cipos <- ifelse(eur_males$Cipos > 15, 15, eur_males$Cipos)

phenotypes <- c(
  `C3_BRONCHUS_LUNG` = "Lung Cancer",
  `C3_CANCER` = "All Cancers",
  `C3_COLORECTAL` = "Colorectal Cancer",
  `C3_MELANOMA_SKIN` = "Skin Melanoma",
  `COX_ARTHROSIS` = "Hip Osteoarthritis",
  `F5_DEPRESSIO` = "Major Depression",
  `G6_EPLEPSY` = "Epilepsy",
  `GOUT` = "Gout",
  `I9_AF` = "Atrial Fibrillation",
  `I9_CHD` = "Coronary Heart Disease",
  #`ILD` = "Interstitial Lung Disease",
  `J10_ASTHMA` = "Asthma",
  `K11_APPENDACUT` = "Appendicitis",
  `KNEE_ARTHROSIS` = "Knee Osteoarthritis",
  `RHEUMA_SEROPOS_OTH` = "Rheumatoid Arthritis",
  `T1D` = "Type 1 Diabetes",
  `T2D` = "Type 2 Diabetes"
)

ggplot(eur_males) +
  geom_point(aes(Biobank, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  xlab(c("")) + 
  ylab(c("Hazard Ratio per Standard Deviation (95% CI)")) + 
  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", " |", "V", "Oldest")) +
  #guides(color = guide_legend(reverse = TRUE)) + 
  geom_errorbar(aes(x=Biobank, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  facet_wrap(~ Phenotype, scales="free_x", labeller = as_labeller(phenotypes))
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_5d_noILD.png",  height=10, width=10)

eur_females <- subset(eur, Sex=="female")
eur_females$Cipos <- ifelse(eur_females$Cipos > 15, 15, eur_females$Cipos)

ggplot(eur_females) +
  geom_point(aes(Biobank, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  xlab(c("")) + 
  ylab(c("Hazard Ratio per Standard Deviation (95% CI)")) + 
  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", " |", "V", "Oldest")) +
  #guides(color = guide_legend(reverse = TRUE)) + 
  geom_errorbar(aes(x=Biobank, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  facet_wrap(~ Phenotype, scales="free_x", labeller = as_labeller(phenotypes))
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_5e_noILD.png",  height=10, width=10)

######################################################################################################################################
######################################################################################################################################

#As above but including biobank japan and facet_grid sex and quartile
#Create forest plots for each pheno - facet_wrap
#Read in all figures

#MGB
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HR_AgeandSexStratified_MGBB_EUR.csv", data.table=FALSE)
mgb_t2d <- subset(mgb, Phenotype=="T2D" & Group=="> 95%")
mgb_t2d$Quartile <- rep(c(1,2,3,4),2)
mgb_t2d$Biobank <- "Mass General Brigham"

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/HUNT/HR_AgeandSexStratifiedHUNT.csv", data.table=FALSE)
hunt_t2d <- subset(hunt, Phenotype=="T2D" & Group=="> 95%")
hunt_t2d$Quartile <- rep(c(1,2,3,4),2)
hunt_t2d$Biobank <- "HUNT"

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HR_AgeandSexStratified_GenomicsEngland.csv", data.table=FALSE)
ge_t2d <- subset(ge, Phenotype=="T2D" & Group=="> 95%")
ge_t2d$Quartile <- rep(c(1,2,3,4),2)
ge_t2d$Biobank <- "Genomics England"

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HR_AgeandSexStratified_GS.csv", data.table=FALSE)
gs_t2d <- subset(gs, Phenotype=="T2D" & Group=="> 95%")
gs_t2d$Quartile <- rep(c(1,2,3,4),2)
gs_t2d$Biobank <- "Generation Scotland"

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HR_AgeandSex_Stratified_FinnGen.csv", data.table=FALSE)
finngen_t2d <- subset(finngen, Phenotype=="T2D" & Group=="> 95%")
finngen_t2d$Quartile <- rep(c(1,2,3,4),2)
finngen_t2d$Biobank <- "FinnGen"

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_AgeandSexStratified_EstBB.csv", data.table=FALSE)
estbb_t2d <- subset(estbb, Phenotype=="T2D" & Group=="> 95%")
estbb_t2d$Quartile <- rep(c(1,2,3,4),2)
estbb_t2d$Biobank <- "Estonian Biobank"

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HR_AgeandSexStratified_BBJ.csv", data.table=FALSE)
colnames(bbj) <- c("Phenotype","PRS","Sex","MinAge","MaxAge","MedianAAO","Group","Controls","Cases","Beta","SE","Pval","HR","Cipos","Cineg")
bbj_t2d <- subset(bbj, Phenotype=="T2D" & Group=="> 95%")
bbj_t2d$Quartile <- rep(c(1,2,3,4),2)
bbj_t2d$Biobank <- "Biobank Japan"

#All
comparison <- rbind(mgb_t2d, hunt_t2d) %>%
                rbind(., ge_t2d) %>%
                  rbind(., gs_t2d) %>%
                    rbind(., finngen_t2d) %>%
                      rbind(., estbb_t2d) %>%
                        rbind(., bbj_t2d) 

comparison <- comparison[,c("Phenotype","Sex","Quartile","HR","Cipos","Cineg","Biobank")]

comparison$Biobank <- factor(comparison$Biobank, levels=c("Generation Scotland","Genomics England", "Mass General Brigham","HUNT", "Estonian Biobank", "FinnGen" ,"Biobank Japan"))
comparison$Ages <- as.factor(case_when(comparison$Quartile==1 ~ "0 - 54.36",
                             comparison$Quartile==2 ~ "54.36 - 63.04",
                             comparison$Quartile==3 ~ "63.04 - 71.13",
                             comparison$Quartile==4 ~ "71.13 - 80"))

scaleFUN <- function(x) sprintf("%.2f", x)

ggplot(comparison) +
  geom_point(aes(Biobank, HR), position=position_dodge(width=0.7)) +
  theme_bw() +
  expand_limits(y=1) + 
  xlab(c("")) + 
  ylab(c("Hazard Ratio (95% CI)")) + 
  scale_y_continuous(labels=scaleFUN) + 
  geom_errorbar(aes(x=Biobank, ymin = Cineg, ymax = Cipos), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)) + 
  facet_grid(Ages ~ Sex)
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_19a.png",  height=7, width=7, dpi=300)

######################################################################################################################################
######################################################################################################################################
#Prostate Cancer

#MGB
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HR_AgeStratified_MGBB_EUR.csv", data.table=FALSE)
mgb_pro <- subset(mgb, Phenotype=="C3_PROSTATE" & Group=="> 95%")
mgb_pro$Quartile <- c(1,2,3,4)
mgb_pro$Biobank <- "Mass General Brigham"

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/HUNT/HR_AgeStratified_HUNT.csv", data.table=FALSE)
colnames(hunt) <- colnames(mgb)
hunt_pro <- subset(hunt, Phenotype=="C3_PROSTATE" & Group=="> 95%")
hunt_pro$Quartile <- c(1,2,3,4)
hunt_pro$Biobank <- "HUNT"

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HR_AgeStratified_GenomicsEngland.csv", data.table=FALSE)
ge_pro <- subset(ge, Phenotype=="C3_PROSTATE" & Group=="> 95%")
ge_pro$Quartile <- c(1,2,3,4)
ge_pro$Biobank <- "Genomics England"

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HR_AgeStratified_GS.csv", data.table=FALSE)
colnames(gs) <- colnames(mgb)
gs_pro <- subset(gs, Phenotype=="C3_PROSTATE" & Group=="> 95%")
gs_pro$Quartile <- c(1,2,3,4)
gs_pro$Biobank <- "Generation Scotland"

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HR_AgeStratified_FinnGen.csv", data.table=FALSE)
finngen_pro <- subset(finngen, Phenotype=="C3_PROSTATE" & Group=="> 95%")
finngen_pro$Quartile <- c(1,2,3,4)
finngen_pro$Biobank <- "FinnGen"

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_AgeStratified_EstBB.csv", data.table=FALSE)
colnames(estbb) <- colnames(mgb)
estbb_pro <- subset(estbb, Phenotype=="C3_PROSTATE" & Group=="> 95%")
estbb_pro$Quartile <- c(1,2,3,4)
estbb_pro$Biobank <- "Estonian Biobank"

#UK Biobank
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_AgeStratified_EstBB.csv", data.table=FALSE)
colnames(ukb) <- colnames(mgb)
ukb_pro <- subset(ukb, Phenotype=="C3_PROSTATE" & Group=="> 95%")
ukb_pro$Quartile <- c(1,2,3,4)
ukb_pro$Biobank <- "UK Biobank"

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HR_AgeStratified_BBJ.csv", data.table=FALSE)
colnames(bbj) <- colnames(mgb)
bbj_pro <- subset(bbj, Phenotype=="C3_PROSTATE" & Group=="> 95%")
bbj_pro$Quartile <- c(1,2,3,4)
bbj_pro$Biobank <- "Biobank Japan"

#All
comparison <- rbind(mgb_pro, hunt_pro) %>%
  rbind(., ge_pro) %>%
  rbind(., finngen_pro) %>%
  rbind(., estbb_pro) %>%
  rbind(., ukb_pro) %>%
  rbind(., bbj_pro) 

comparison <- comparison[,c("Phenotype","Quartile","HR","Cipos","Cineg","Biobank")]

comparison$Biobank <- factor(comparison$Biobank, levels=c("Generation Scotland","Genomics England", "Mass General Brigham","HUNT", "Estonian Biobank", "UK Biobank", "FinnGen", "Biobank Japan"))
comparison$Ages <- as.factor(case_when(comparison$Quartile==1 ~ "0 - 62.60",
                                       comparison$Quartile==2 ~ "62.60 - 68.25",
                                       comparison$Quartile==3 ~ "68.25 - 73.89",
                                       comparison$Quartile==4 ~ "73.89 - 80"))

scaleFUN <- function(x) sprintf("%.2f", x)

ggplot(comparison) +
  geom_point(aes(Biobank, HR), position=position_dodge(width=0.7)) +
  theme_bw() +
  expand_limits(y=1) + 
  xlab(c("")) + 
  ylab(c("Hazard Ratio (95% CI)")) + 
  scale_y_continuous(labels=scaleFUN) + 
  geom_errorbar(aes(x=Biobank, ymin = Cineg, ymax = Cipos), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)) + 
  facet_wrap(~ Ages)
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_19b.png",  height=7, width=7, dpi=300)

######################################################################################################################################
######################################################################################################################################
# Leave one out meta-analysis - Only really relevant for European meta-analysis - Each mark can show the effect on the effect size

metaresults <- c()
european <- subset(all, Ancestry=="EUR")

#Make sure that both males and females have been tested before meta-analysing
eur <- c()
for(i in unique(european$Phenotype)){
  for(k in unique(european$Biobank)){
    pheno <- subset(european, Phenotype==i & Biobank==k)
    print(pheno)
      
    if(any(is.na(pheno[,c(9:14)])) | dim(pheno)[1]<8){
      next
    }
      
    eur <- rbind(eur,pheno)
  }
}

european <- as.data.frame(eur)

#Also restrict to phenotypes which have results for both male and female sample
for(k in c(1,2,3,4)){
  for(l in c("male","female")){
    europeanAge <- subset(european, Quartile==k & Sex==l)
    for(i in unique(europeanAge$Biobank)){
      
      loo <- subset(europeanAge, Biobank!=i)
      
      for(j in unique(loo$Phenotype)){
        print(j)
        
        #Check to see if the biobank has been tested
        if(dim(subset(europeanAge, Phenotype==j & Biobank==i))[1]==0){
          next
        }
        
        disease <- subset(loo, Phenotype==j & !(is.na(Beta)))
        
        if(dim(disease)[1]==0){
          next
        }
        
        #Meta analysis should be done at the beta level and stratified by ancestry 
        meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
        
        forest(meta, slab=disease$Biobank)
        grid.text(paste0(i," ",k," - EUR Ancestry"), .5, .9, gp=gpar(cex=2))
        
        metaresults <- rbind(metaresults, c(j, k, l, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
        
      }
    }  
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","Quartile","Sex","LeftOutBiobank","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(5:9),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_AgeandSexStratified.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_AgeandSexStratified.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeandSexStratifiedMetaAnalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
orderPheno <- subset(meta, Quartile==1 & LeftOutBiobank=="All" & Sex=="female")
meta$Phenotype <- factor(meta$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))
meta$Quartile <- factor(meta$Quartile, levels=c(1,2,3,4))

#Plot leave one out results 
for(k in c(1,2,3,4)){
  for(l in c("male","female")){
    age <- subset(meta, Quartile==k & Sex==l)
    
    print(ggplot(age) +
      geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
      theme_bw() +
      geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
      ylab("Hazard Ratio (95% CI)") +
      xlab("") +
      geom_hline(yintercept = 1.0) +
      scale_x_discrete(labels=c("Appendicitis","ILD","Lung Cancer","All Cancers","Major Depression","Epilepsy","Skin Melanoma","Hip Osteoarthritis","Knee Osteoarthritis","CHD","Colorectal Cancer","Asthma","Gout","Atrial Fibrillation","Type 2 Diabetes","Type 1 Diabetes")) + 
      scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
      theme(title = element_text(size = 18),
            legend.text = element_text(size = 14),
            legend.title = element_blank(),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 14)) +
      coord_flip())
    ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_LOOMetaAnalysisEUR_",l,"_Quartile",k,".png"), height=8 , width=10)
    
  }
}

age <- subset(meta, Quartile==1 & Sex=="female" & Phenotype=="C3_CANCER")
