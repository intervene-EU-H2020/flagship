library(data.table)
library(metafor)
library(dplyr)
library(grid)
library(gridExtra)

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

#Meta analyse Interaction

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/HR_SexInteraction_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HR_SexInteraction_UKB_AllAncestries.csv", data.table=FALSE)

#Genes&Health
gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HR_SexInteraction_GNH.csv", data.table=FALSE)

#Biobank Japan
bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HR_SexInteraction_BBJ.csv", data.table=FALSE)

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_SexInteraction_EstBB.csv", data.table=FALSE)

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HR_SexInteraction_GS.csv", data.table=FALSE)

#Mass General Brigham
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HR_SexInteraction_MGBB.csv", data.table=FALSE)

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HR_SexInteraction_GenomicsEngland.csv", data.table=FALSE)

#Combine into one dataset
all <- rbind(finngen, gnh) %>%
        rbind(., estbb) %>% 
          rbind(., gs) %>%
            rbind(., mgb) %>%
              rbind(., bbj) %>%
                rbind(., ukb) %>%
                  rbind(., ge)


all$Beta <- ifelse(all$Biobank=="Genomics England" & (all$Test=="Sex_Main" | all$Test=="PRS*Sex"), all$Beta*-1, all$Beta)

#Meta analysis
metaresults <- c()

for(k in unique(all$Test)){
  for(j in c("EUR","SAS","EAS")){
    
    for(i in unique(all$Phenotype)){
      print(i)
      #Test with T2D phenotype and EUR ancestry
      if(j=="EAS"){
        disease <- subset(all, Test==k & Phenotype==i & Phenotype %in% bbj$Phenotype & Ancestry==j & !(is.na(Beta)))
      } else{
        disease <- subset(all, Test==k & Phenotype==i & Ancestry==j & !(is.na(Beta)))
      }
      
      if(dim(disease)[1]==0){
        next
      }
      
      #Meta analysis should be done at the beta level and stratified by ancestry 
      meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
      
      forest(meta, slab=disease$Biobank)
      grid.text(paste0(i," - ",j," Ancestry"), .5, .9, gp=gpar(cex=2))
      
      metaresults <- rbind(metaresults, c(k, j, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
      
    }
    
  }
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Test","Ancestry","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexInteractionMetaAnalysisFE.csv")

######################################################################################################################################
# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexInteractionMetaAnalysisFE.csv", data.table=FALSE)
meta$Ancestry <- factor(meta$Ancestry, levels=c("EUR","SAS","EAS"))
orderPheno <- subset(meta, Ancestry=="EUR")
meta$Phenotype <- factor(meta$Phenotype, levels=c(meta[order(orderPheno$HR),"Phenotype"]))

#Plot
ggplot(meta) +
  geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77")) +
  scale_x_discrete(labels=c("Type 2 Diabetes","Lung Cancer","Epilepsy","Type 1 Diabetes","All Cancers","Rheumatoid Arthritis","Knee Osteoarthritis","Major Depression","Appendicitis","Atrial Fibrillation","Asthma","Hip Osteoarthritis","Skin Melanoma","Colorectal Cancer","Gout","CHD","ILD")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_SexInteraction.png", height=8, width=8)

######################################################################################################################################
#Supplementary figures for amr and afr

ukbaframr <- subset(all, (Ancestry=="AFR" | Ancestry=="AMR") & Biobank=="UK Biobank" & !(is.na(HR)))
afronly <- subset(ukbaframr, Ancestry=="AFR")
ukbaframr$Phenotype <- factor(ukbaframr$Phenotype, afronly[order(afronly$HR),"Phenotype"])

#Supplementary Plot - AFR and AMR
ggplot(ukbaframr) +
  geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#88CCEE", "#CC6677")) +
  scale_x_discrete(labels=c("ILD","Type 1 Diabetes","Major Depression","Lung Cancer","Atrial Fibrillation","Knee Osteoarthritis","Epilepsy","Type 2 Diabetes","Colorectal Cancer","All Cancers","Appendicitis","Asthma","CHD","Gout","Hip Osteoarthritis")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigureAFR&AMRAncestriesHRperSDSexInteraction.png", height=8 , width=8)

######################################################################################################################################

#Leave one out analysis

metaresults <- c()
european <- subset(all, Ancestry=="EUR")

for(i in unique(european$Biobank)){
  
  loo <- subset(european, Biobank!=i)
  
  #Repeat meta-analysis and save down the results but allude to the biobank
  
  #European ancestry meta-analysis
  
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
    grid.text(paste0(i," - EUR Ancestry"), .5, .9, gp=gpar(cex=2))
    
    metaresults <- rbind(metaresults, c(j, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
    
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","LeftOutBiobank","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(3:7),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexInteraction_LOOmetaanalysisFE.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexInteraction_LOOmetaanalysisFE.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/SexInteractionmetaanalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
meta$Phenotype <- factor(meta$Phenotype, levels=c(allmeta[order(allmeta$HR),"Phenotype"]))

#Plot leave one out results
ggplot(meta) +
  geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
  scale_x_discrete(labels=c("Type 2 Diabetes","Lung Cancer","Epilepsy","Type 1 Diabetes","All Cancers","Rheumatoid Arthritis","Knee Osteoarthritis","Major Depression","Appendicitis","Atrial Fibrillation","Asthma","Hip Osteoarthritis","Skin Melanoma","Colorectal Cancer","Gout","CHD","ILD")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_SexInteraction_LOOMetaAnalysisEUR.png", height=8 , width=10)
