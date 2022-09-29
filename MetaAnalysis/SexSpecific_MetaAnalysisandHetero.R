library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)

#Read in full sample hazard ratios per standard deviation

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/PRS_HRsperSD_UKBiobank_AllAncestries.csv", data.table=FALSE)
ukb$Biobank <- "UK Biobank"

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/HRperSD_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"
finngen <- finngen[,names(ukb)]

#Genes&Health
gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HRperSD_GNH.csv", data.table=FALSE)
gnh$Biobank <- "Genes & Health"
gnh$Ancestry <- "SAS"
gnh <- gnh[,names(ukb)]

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HRperSD_BBJ.csv", data.table=FALSE)
dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
bbj <- subset(bbj, !(Phenotype %in% dropbbj))
bbj$Biobank <- "Biobank Japan"
bbj$Ancestry <- "EAS"
bbj <- bbj[,names(ukb)]

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HRperSD_EstBB.csv", data.table=FALSE)
estbb$Biobank <- "Estonian Biobank"
estbb$Ancestry <- "EUR"
estbb <- estbb[,names(ukb)]

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HRperSD_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs$Biobank <- "Generation Scotland"
gs$Ancestry <- "EUR"
gs <- gs[,names(ukb)]

#Mass General Brigham
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_MGBB.csv", data.table=FALSE)
dropmgb <- c("I9_AF")
mgb <- subset(mgb, !(Phenotype %in% dropmgb))
mgb$Biobank <- "Mass General Brigham"
mgb$Ancestry <- "EUR"
mgb <- mgb[,names(ukb)]

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HRperSD_GenomicsEngland.csv", data.table=FALSE)
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
ge <- ge[,names(ukb)]

#Combine into one dataset
all <- rbind(ukb, finngen) %>%
        rbind(., gnh) %>% 
          rbind(., estbb) %>%
            rbind(., gs) %>%
              rbind(., mgb) %>% 
                rbind(., bbj) %>% 
                  rbind(., ge)

both <- subset(all, Sample=="Male Sample" | Sample=="Female Sample")
both$Sample[both$Sample=="Male Sample"] <- "Males"
both$Sample[both$Sample=="Female Sample"] <- "Females"

b <- c()
for(i in unique(both$Phenotype)){
  for(k in unique(both$Biobank)){
    for(j in unique(both$Ancestry)){
      pheno <- subset(both, Phenotype==i & Biobank==k & Ancestry==j)
      print(pheno)
      
      if(any(is.na(pheno[,c(7:12)]))){
        next
      }
      
      b <- rbind(b,pheno)
    }
  }
}

both <- as.data.frame(b)

metaresults <- c()
for(k in c("Males","Females")){
  
  sex <- subset(both, Sample==k)
  
  for(j in c("EUR","SAS","EAS")){
    
    for(i in unique(sex$Phenotype)){
      print(i)
 
      if(j=="EAS"){
        disease <- subset(sex, Phenotype==i & Phenotype %in% bbj$Phenotype & Ancestry==j & !(is.na(Beta)))
      } else{
        disease <- subset(sex, Phenotype==i & Ancestry==j & !(is.na(Beta)))
      }
      
      if(dim(disease)[1]==0){
        next
      }
      
      #Meta analysis should be done at the beta level and stratified by ancestry 
      meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
      
      png(file=paste0("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexStratified/ForestPlots/",k,"_",i,"_",j,"_Forest.png"), width = 7,
          height    = 7,
          units     = "in",
          res       = 300)
      forest(meta, slab=disease$Biobank)
      dev.off()
      
      metaresults <- rbind(metaresults, c(j, k, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
      
    }
    
  }

}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Ancestry","Sex","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexStratifiedMetaAnalysisFE.csv")

######################################################################################################################################

# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexStratifiedMetaAnalysisFE.csv", data.table=FALSE)
meta$Ancestry <- factor(meta$Ancestry, levels=c("EUR","SAS","EAS"))
orderPheno <- subset(meta, Ancestry=="EUR" & Sex=="Females")
meta$Phenotype <- factor(meta$Phenotype, levels=c(meta[order(orderPheno$HR),"Phenotype"]))

#Plot
figure2b <- ggplot(meta) +
  geom_point(aes(Phenotype, HR, group=interaction(Sex,Ancestry), col=Ancestry, shape=Sex), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=interaction(Sex,Ancestry), col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#D55E00", "#56B4E9", "#009E73")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","ILD","Major Depression","Lung Cancer","Skin Melanoma","Hip Osteoarthritis","Knee Osteoarthritis","CHD","Asthma","Colorectal Cancer","Gout","Atrial Fibrillation","Rheumatoid Arthritis","Type 2 Diabetes","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/Figure2b.png", height=8, width=8)

######################################################################################################################################
######################################################################################################################################

#Significant differences between meta-analysis estimates and individual biobanks

both <- both[,c("Ancestry","Sample","Phenotype","Beta","SE","Biobank")]
colnames(both)[2] <- "Sex"
diffs <- meta[,c("Ancestry","Sex","Phenotype","Beta","SE")]
diffs$Biobank <- "All"
diffs <- rbind(diffs, both)

differences <- c()
for(j in c("EUR","SAS","EAS")){
  for(k in c("Males","Females")){
    for(i in unique(diffs$Phenotype)){
      
      disease <- subset(diffs, Phenotype==i & Ancestry==j & Sex==k)
      
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

fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/Sex_BiobankVariation.csv")

######################################################################################################################################
######################################################################################################################################

#Significant differences between meta-analysis estimates across sexes?

diffs <- meta[,c("Ancestry","Sex","Phenotype","Beta","SE")]

differences <- c()
for(j in c("EUR","SAS","EAS")){
  for(i in unique(meta$Phenotype)){
    disease <- subset(diffs, Phenotype==i & Ancestry==j)
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

fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/Sex_HRDifferences_MA.csv")

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

#Supplementary figures for amr and afr

ukbaframr <- subset(both, (Ancestry=="AFR" | Ancestry=="AMR") & Biobank=="UK Biobank" & !(is.na(HR)))
afronly <- subset(ukbaframr, Ancestry=="AFR" & Sample=="Females")
ukbaframr$Phenotype <- factor(ukbaframr$Phenotype, afronly[order(afronly$HR),"Phenotype"])

#Supplementary Plot - AFR and AMR
ggplot(ukbaframr) +
  geom_point(aes(Phenotype, HR, group=interaction(Sample,Ancestry), col=Ancestry, shape=Sample), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=interaction(Sample,Ancestry), col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#88CCEE", "#CC6677")) +
  scale_x_discrete(labels=c("Appendicitis","ILD","Epilepsy","CHD","All Cancers","Knee Osteoarthritis","Hip Osteoarthritis","Gout","Major Depression","Lung Cancer","Type 1 Diabetes","Colorectal Cancer","Atrial Fibrillation","Type 2 Diabetes","Asthma")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigureAFR&AMRAncestriesHRperSDSexStratified.png", height=8 , width=8)

######################################################################################################################################
# Leave one out meta-analysis - Only really relevant for European meta-analysis - Each mark can show the effect on the effect size
metaresults <- c()
european <- subset(both, Ancestry=="EUR")

#Make sure that both males and females have been tested before meta-analysing
eur <- c()
for(i in unique(european$Phenotype)){
  for(k in unique(european$Biobank)){
    pheno <- subset(european, Phenotype==i & Biobank==k)
    print(pheno)
    
    if(any(is.na(pheno[,c(7:12)]))){
      next
    }
    
    eur <- rbind(eur,pheno)
  
  }
}

european <- as.data.frame(eur)

#Also restrict to phenotypes which have results for both male and female sample

for(k in c("Males","Females")){
  
  europeanSex <- subset(european, Sample==k)
  
  for(i in unique(europeanSex$Biobank)){
    
    loo <- subset(europeanSex, Biobank!=i)
    
    for(j in unique(loo$Phenotype)){
      print(j)
      
      #Check to see if the biobank has been tested
      if(dim(subset(european, Phenotype==j & Biobank==i))[1]==0){
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
colnames(metaresults) <- c("Phenotype","Sex","LeftOutBiobank","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_SexStratified.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE_SexStratified.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexStratifiedMetaAnalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
orderPheno <- subset(meta, Sex=="Females" & LeftOutBiobank=="All")
meta$Phenotype <- factor(meta$Phenotype, levels=c(orderPheno[order(orderPheno$HR),"Phenotype"]))

#Plot leave one out results - do separately for males and females otherwise gets too confusing...
for(k in c("Males","Females")){

  sex <- subset(meta, Sex==k)
  
  print(ggplot(sex) +
    geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
    theme_bw() +
    geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
    ylab("Hazard Ratio (95% CI)") +
    xlab("") +
    geom_hline(yintercept = 1.0) +
    scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","ILD","Lung Cancer","Major Depression","Skin Melanoma","Hip Osteoarthritis","Knee Osteoarthritis","CHD","Asthma","Colorectal Cancer","Gout","Atrial Fibrillation","Rheumatoid Arthritis","Type 2 Diabetes","Type 1 Diabetes")) + 
    scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
    theme(title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 14)) +
    coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_LOOMetaAnalysisEUR_",k,".png"), height=8 , width=10)

}

