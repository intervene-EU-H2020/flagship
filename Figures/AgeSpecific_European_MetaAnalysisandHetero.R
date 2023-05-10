library(data.table)
library(dplyr)
library(metafor)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

#Read in age stratified hazard ratios

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HRperSD_AgeStratified_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HRperSD_AgeStratified_UKB_AllAncestries.csv", data.table=FALSE)
ukb <- ukb[,names(finngen)]

#Genes&Health
gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HRperSD_AgeStratified_GNH.csv", data.table=FALSE)
gnh$Biobank <- "Genes & Health"
gnh$Ancestry <- "SAS"
gnh <- gnh[,names(finngen)]

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HRperSD_AgeStratified_BBJ.csv", data.table=FALSE)
bbj$Biobank <- "Biobank Japan"
bbj$Ancestry <- "EAS"
bbj <- bbj[,names(finngen)]

dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
bbj <- subset(bbj, !(Phenotype %in% dropbbj))

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HRperSD_AgeStratified_EstBB.csv", data.table=FALSE)
estbb <- estbb[,names(finngen)]

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HRperSD_AgeStratified_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs <- gs[,names(finngen)]

#Mass General Brigham - European
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeStratified_MGBB_EUR.csv", data.table=FALSE)
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
#mgb_afr <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeStratified_MGBB_AFR.csv", data.table=FALSE)
#dropmgb <- c("I9_AF")
#mgb_afr <- subset(mgb_afr, !(Phenotype %in% dropmgb))
#mgb_afr$Biobank <- "Mass General Brigham"
#mgb_afr$Ancestry <- "AFR"
#mgb_afr <- mgb_afr[,names(ukb)]

#Mass General Brigham
#mgb <- rbind(mgb_eur, mgb_afr)

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HRperSD_AgeStratified_GenomicsEngland.csv", data.table=FALSE)
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
ge<- ge[,names(finngen)]

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/HUNT/HRperSD_AgeStratified_HUNT.csv", data.table=FALSE)
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

#Make sure all quartiles have been tested before including in the meta-analysis
q <- c()
for(i in unique(all$Phenotype)){
  for(k in unique(all$Biobank)){
    for(j in unique(all$Ancestry)){
      pheno <- subset(all, Phenotype==i & Biobank==k & Ancestry==j)
      print(pheno)
        
      if(any(is.na(pheno[,c(8:13)])) | dim(pheno)[1]<4){
        next
      }
        
      q <- rbind(q,pheno)
    }
  }
}

all <- as.data.frame(q)

metaresults <- c()
for(k in c(1,2,3,4)){
  
  quartile <- subset(all, Quartile==k)
  
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
      
      #png(file=paste0("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/ForestPlots/Quartile_",k,"_",i,"_",j,"_Forest.png"), width = 7,
      #    height    = 7,
      #    units     = "in",
      #    res       = 300)
      #forest(meta, slab=disease$Biobank)
      #dev.off()
      
      results <- results <- matrix(c("Fixed Effect", "Random Effect", j, j, k, k, i, i, metaFE$b, metaRE$b, metaFE$se, metaRE$se, metaFE$pval, metaRE$pval, metaFE$QE, metaRE$QE, metaFE$QEp, metaRE$QEp), nrow = 2, ncol = 9)
      
      metaresults <- rbind(metaresults, results)
      
    }
    
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Test","Ancestry","Quartile","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(5:9),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv")

######################################################################################################################################

#alld <- all[,c("Ancestry","Quartile","Phenotype","Beta","SE","Biobank")]
#diffs <- metaresults[,c("Ancestry","Quartile","Phenotype","Beta","SE")]
#diffs$Biobank <- "All"
#diffs <- rbind(diffs, alld)

#differences <- c()
#for(j in c("EUR","SAS","EAS")){
#  for(k in c(1,2,3,4)){
#    for(i in unique(diffs$Phenotype)){
      
#      disease <- subset(diffs, Phenotype==i & Ancestry==j & Quartile==k)
      
#      if(dim(disease)[1]==0){
#        next
#      }
      
#      disease$BetaDiff <- ifelse(!(is.na(disease$Beta)), disease$Beta - disease$Beta[disease$Biobank=="All"], NA)
#      disease$SEDiff <- ifelse(!(is.na(disease$Beta)), sqrt(disease$SE**2 + disease$SE[disease$Biobank=="All"]**2), NA)
#      disease$ZDiff <- ifelse(!(is.na(disease$Beta)), disease$BetaDiff/disease$SEDiff, NA)
#      disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
#      differences <- rbind(differences, disease)
      
#    }
#  }
#}

#fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/Age_BiobankVariation.csv")

######################################################################################################################################

# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
meta_eur <- subset(meta, Ancestry=="EUR" & Phenotype!="ILD")
#write.csv(meta_eur, "/Users/jermy/Documents/INTERVENE/Write-up/Supplementary Tables/SupplementaryTable7c.csv")
meta_eur$Phenotype <- factor(meta_eur$Phenotype, levels=c("K11_APPENDACUT","G6_EPLEPSY","C3_CANCER", "F5_DEPRESSIO", "C3_BRONCHUS_LUNG", "C3_MELANOMA_SKIN", "KNEE_ARTHROSIS", "COX_ARTHROSIS",
                                                          "I9_CHD", "J10_ASTHMA", "C3_COLORECTAL", "I9_AF", "C3_BREAST","RHEUMA_SEROPOS_OTH", "GOUT", "T2D", "C3_PROSTATE", "T1D"))
meta_eur$Quartile <- factor(meta_eur$Quartile, levels=c(1,2,3,4))

meta_eur2c <- subset(meta_eur, Test=="Fixed Effect")
meta_eur2c <- meta_eur2c[,-1]

ggplot(meta_eur2c) +
  geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Meta-Analyzed Hazard Ratio per Standard Deviation (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  #guides(color = guide_legend(reverse = TRUE)) + 
  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", " |", "V", "Oldest")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Major Depression","Lung Cancer","Skin Melanoma","Knee Osteoarthritis","Hip Osteoarthritis","Coronary Heart Disease","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/Figure2c_noILD.png", height=8, width=8)

######################################################################################################################################

# Supplementary analysis - EAS and SAS significant age effects
#meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
#meta_sas <- subset(meta, Ancestry=="SAS" & (Phenotype=="T2D" | Phenotype=="T1D" | Phenotype=="I9_CHD") & Test=="Fixed Effect")
#meta_sas$Quartile <- factor(meta_sas$Quartile, levels=c(4,3,2,1))

#sas <- ggplot(meta_sas) +
#          geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
#          theme_bw() +
#          geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
#          ylab("Meta-Analyzed Hazard Ratio per Standard Deviation (95% CI)") +
#          xlab("") +
#          geom_hline(yintercept = 1.0) +
#          scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Oldest", ">", " --", "Youngest")) +
#          guides(color = guide_legend(reverse = TRUE)) + 
#          scale_x_discrete(labels=c("Coronary Heart Disease", "Type 1 Diabetes", "Type 2 Diabetes")) + 
#          theme(title = element_text(size = 18),
#                legend.text = element_text(size = 14),
#                legend.title = element_text(size = 18),
#                axis.title.x = element_text(size = 18),
#                axis.text.x = element_text(size = 14),
#                axis.title.y = element_text(size = 18),
#                axis.text.y = element_text(size = 14)) +
#          coord_flip()

#meta_eas <- subset(meta, Ancestry=="EAS" & (Phenotype=="T2D" | Phenotype=="C3_PROSTATE") & Test=="Fixed Effect")
#meta_eas$Quartile <- factor(meta_eas$Quartile, levels=c(4,3,2,1))

#eas <- ggplot(meta_eas) +
#  geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
#  theme_bw() +
#  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
#  ylab("") +
#  xlab("") +
#  geom_hline(yintercept = 1.0) +
#  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Oldest", ">", " --", "Youngest")) +
#  scale_x_discrete(labels=c("Prostate Cancer","Type 2 Diabetes")) + 
#  guides(color = guide_legend(reverse = TRUE)) + 
#  theme(title = element_text(size = 18),
#        legend.text = element_text(size = 14),
#        legend.title = element_text(size = 18),
#        axis.title.x = element_text(size = 18),
#        axis.text.x = element_text(size = 14),
#        axis.title.y = element_text(size = 18),
#        axis.text.y = element_text(size = 14)) +
#  coord_flip()

#suppfig18 <- ggarrange(eas, sas, labels=c("East Asian","South Asian"), common.legend=TRUE, align = "hv", nrow = 2, ncol = 1)
#ggsave(filename="/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/SupplementaryFig18_EAS_SAS_AgeStrat.png", plot = suppfig18, height=10, width=10,dpi=300)

######################################################################################################################################

#Create forest plots for each pheno - facet_wrap
eur <- subset(all, Ancestry=="EUR" & Phenotype!="ILD")
eur <- eur[,c("Phenotype","Quartile","HR","Cipos","Cineg","Biobank")]

meta_eur$Biobank <- case_when(meta_eur$Test=="Fixed Effect" ~ "Meta Analysis FE",
                              meta_eur$Test=="Random Effect" ~ "Meta Analysis RE")

meta_eur <- meta_eur[,c("Phenotype","Quartile","HR","Cipos","Cineg","Biobank")]
eur <- rbind(eur, meta_eur)
eur$Biobank <- factor(eur$Biobank, levels=c("Meta Analysis RE", "Meta Analysis FE", "FinnGen", "Estonian Biobank", "UK Biobank", "HUNT", "Mass General Brigham","Genomics England", "Generation Scotland"))
eur$Quartile <- factor(eur$Quartile, levels=c(4,3,2,1))

scaleFUN <- function(x) sprintf("%.2f", x)

phenotypes <- c(
  `C3_BREAST` = "Breast Cancer",
  `C3_BRONCHUS_LUNG` = "Lung Cancer",
  `C3_CANCER` = "All Cancers",
  `C3_COLORECTAL` = "Colorectal Cancer",
  `C3_MELANOMA_SKIN` = "Skin Melanoma",
  `C3_PROSTATE` = "Prostate Cancer",
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

ggplot(eur) +
  geom_point(aes(Biobank, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  xlab(c("")) + 
  ylab(c("Hazard Ratio per Standard Deviation (95% CI)")) + 
  scale_color_manual(name="Age Quartiles", values=c("#DDCC77", "#117733", "#332288" ,"#AA4499"), labels=c("Youngest", " |", "V", "Oldest")) +
  geom_errorbar(aes(x=Biobank, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  #guides(color = guide_legend(reverse = TRUE)) + 
  facet_wrap(~ Phenotype, scales="free_x", labeller = as_labeller(phenotypes)) +
  theme(axis.text.x = element_text(size = 8))
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_5c_noILD.png",  height=10, width=10)

######################################################################################################################################

#Difference in estimates for each biobank relative to the meta-analytic estimates?




######################################################################################################################################

ukbaframr <- subset(all, Ancestry=="AFR" | Ancestry=="AMR")
orderPheno <- subset(ukbaframr, Ancestry=="AFR" & Quartile==1)
ukbaframr$Ancestry <- factor(ukbaframr$Ancestry, levels=c("AMR","AFR"))
ukbaframr$Phenotype <- factor(ukbaframr$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))
ukbaframr$Quartile <- as.factor(ukbaframr$Quartile)

ancestry <- c("AFR","AMR")
europelabels <- list(c("ILD","Knee Osteoarthritis","CHD","All Cancers","Prostate Cancer","Hip Osteoarthritis","Type 2 Diabetes","Atrial Fibrillation","Breast Cancer","Gout"),
                     c("Knee Osteoarthritis","CHD","All Cancers","Type 2 Diabetes"))

for(i in 1:2){
  
  anc <- subset(ukbaframr, Ancestry==ancestry[i])
  
  #Plot
  print(ggplot(anc) +
          geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
          theme_bw() +
          geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
          ylab("Hazard Ratio (95% CI)") +
          xlab("") +
          geom_hline(yintercept = 1.0) +
          scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#000000")) +
          scale_x_discrete(labels=europelabels[i]) + 
          theme(title = element_text(size = 18),
                legend.text = element_text(size = 14),
                legend.title = element_blank(),
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 14),
                axis.title.y = element_text(size = 18),
                axis.text.y = element_text(size = 14)) +
          coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_AgeStratified_",ancestry[i],".png"), height=8, width=8)
  
}

######################################################################################################################################
######################################################################################################################################

#Significant differences between meta-analysis estimates? 

diffs <- meta[,c("Ancestry","Quartile","Phenotype","Beta","SE")]

differences <- c()
for(j in unique(diffs$Ancestry)){
  for(i in unique(diffs$Phenotype)){
    disease <- subset(diffs, Ancestry==j & Phenotype==i & (Quartile==1 | Quartile==4))
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
fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/Age_HRDifferences_MA.csv")

######################################################################################################################################

#Related to this but looking for heterogeneity statistics
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
meta <- subset(meta, Test=="Fixed Effect")

metaresults <- c()
for(j in unique(meta$Ancestry)){
  for(i in unique(meta$Phenotype)){
    disease <- subset(meta, Phenotype==i & Ancestry==j)
    
    if(nrow(disease)==0){
      next
    }
    
    metaFE <- rma(yi=Beta, sei=SE, data=disease, method="FE")
    results <- c(j, i, metaFE$b, metaFE$se, metaFE$pval, metaFE$QE, metaFE$QEp)
    
    metaresults <- rbind(metaresults, results)
  }
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Ancestry","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(3:7),as.numeric)

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/HeterogeneityinAgeEffects_ModelSelection_AllAncestries.csv")

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
    
    if(any(is.na(pheno[,c(8:12)])) | dim(pheno)[1]<4){
      next
    }
    
    eur <- rbind(eur,pheno)
    
  }
}

european <- as.data.frame(eur)

#Also restrict to phenotypes which have results for both male and female sample

for(k in c(1,2,3,4)){
  
  europeanAge <- subset(european, Quartile==k)
  
  for(i in unique(europeanAge$Biobank)){
    
    loo <- subset(europeanAge, Biobank!=i)
    
    for(j in unique(loo$Phenotype)){
      print(j)
      
      #Check to see if the biobank has been tested
      if(dim(subset(europeanAge, Phenotype==j & Biobank==i))[1]==0){
        next
      }
      
      #Test with T2D phenotype and EUR ancestry
      disease <- subset(loo, Phenotype==j & !(is.na(Beta)))
      
      if(dim(disease)[1]==0){
        next
      }
      
      #Meta analysis should be done at the beta level and stratified by ancestry 
      meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
      
      forest(meta, slab=disease$Biobank)
      grid.text(paste0(i," ",k," - EUR Ancestry"), .5, .9, gp=gpar(cex=2))
      
      metaresults <- rbind(metaresults, c(j, k, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
      
    }
    
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","Quartile","LeftOutBiobank","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_AgeStratified.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_AgeStratified.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratifiedMetaAnalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
orderPheno <- subset(meta, Quartile==1 & LeftOutBiobank=="All")
meta$Phenotype <- factor(meta$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))
meta$Quartile <- factor(meta$Quartile, levels=c(1,2,3,4))

#Plot leave one out results 
for(k in c(1,2,3,4)){
  
  age <- subset(meta, Quartile==k)
  
  print(ggplot(age) +
    geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
    theme_bw() +
    geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
    ylab("Hazard Ratio (95% CI)") +
    xlab("") +
    geom_hline(yintercept = 1.0) +
    scale_x_discrete(labels=c("ILD","Appendicitis","All Cancers","Lung Cancer","Epilepsy","Major Depression","Skin Melanoma","Knee Osteoarthritis","Hip Osteoarthritis","CHD","Colorectal Cancer","Asthma","Breast Cancer","Gout","Atrial Fibrillation","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
    scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
    theme(title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 14)) +
    coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_LOOMetaAnalysisEUR_Quartile",k,".png"), height=8 , width=10)
  
}