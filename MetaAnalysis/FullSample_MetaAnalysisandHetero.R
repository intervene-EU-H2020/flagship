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
                    rbind(.,ge)

full <- subset(all, Sample=="Full Sample")

metaresults <- c()
for(j in c("EUR","SAS","EAS")){

  for(i in unique(full$Phenotype)){
    print(i)
    #Test with T2D phenotype and EUR ancestry
    if(j=="EAS"){
      disease <- subset(full, Phenotype==i & Phenotype %in% bbj$Phenotype & Ancestry==j & !(is.na(Beta)))
    } else{
      disease <- subset(full, Phenotype==i & Ancestry==j & !(is.na(Beta)))
    }
    
    if(dim(disease)[1]==0){
      next
    }

    #Meta analysis should be done at the beta level and stratified by ancestry 
    meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
    
    png(file=paste0("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/FullSample/HRperSD/ForestPlots/",i,"_",j,"_Forest.png"), width = 7,
        height    = 7,
        units     = "in",
        res       = 300)
    forest(meta, slab=disease$Biobank)
    dev.off()
    
    #Save the forest plot so can compare with the table
    
    metaresults <- rbind(metaresults, c(j, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
    
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Ancestry","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(3:7),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/metaanalysisFE.csv")

######################################################################################################################################
#Combine the meta analysed results with each biobank and then compare - reduce to those with differences and save file only
metaresults$Biobank <- "All"
metadiff <- metaresults[,c("Ancestry","Phenotype","Beta","SE","Pval","Biobank")]
fulldiff <- full[,names(metadiff)]
diffs <- rbind(fulldiff, metadiff)

diff <- c()
for(i in c("EUR","SAS","EAS")){
  for(j in unique(diffs$Phenotype)){
    disease <- subset(diffs, Phenotype==j & Ancestry==i)
    disease$BetaDiff <- ifelse(!(is.na(disease$Beta)), disease$Beta - disease$Beta[disease$Biobank=="All"], NA)
    disease$SEDiff <- ifelse(!(is.na(disease$Beta)), sqrt(disease$SE**2 + disease$SE[disease$Biobank=="All"]**2), NA)
    disease$ZDiff <- ifelse(!(is.na(disease$Beta)), disease$BetaDiff/disease$SEDiff, NA)
    disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
    diff <- rbind(diff, disease)
  }
}

diffs <- as.data.frame(diff)

fwrite(diffs, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/FullSampleIndBiobankComparisonwMAEffectSize.csv")

######################################################################################################################################
# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/metaanalysisFE.csv", data.table=FALSE)
meta$Ancestry <- factor(meta$Ancestry, levels=c("EUR","SAS","EAS"))
orderPheno <- subset(meta, Ancestry=="EUR")
meta$Phenotype <- factor(meta$Phenotype, levels=c(meta[order(orderPheno$HR),"Phenotype"]))

#Main Plot - EUR, SAS and EAS
figure2a <- ggplot(meta) +
              geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
              theme_bw() +
              geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
              ylab("Hazard Ratio (95% CI)") +
              xlab("") +
              geom_hline(yintercept = 1.0) +
              scale_color_manual(values=c("#D55E00", "#56B4E9", "#009E73")) +
              scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
              theme(title = element_text(size = 18),
                    legend.text = element_text(size = 14),
                    legend.title = element_blank(),
                    axis.title.x = element_text(size = 18),
                    axis.text.x = element_text(size = 14),
                    axis.title.y = element_text(size = 18),
                    axis.text.y = element_text(size = 14)) +
              coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/Figure2a.png", height=8 , width=8)

######################################################################################################################################

ukbaframr <- subset(ukb, (Ancestry=="AFR" | Ancestry=="AMR") & Sample=="Full Sample" & !(is.na(HR)))
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
  scale_x_discrete(labels=c("ILD","Appendicitis","Epilepsy","CHD","All Cancers","Lung Cancer","Type 1 Diabetes","Major Depression","Knee Osteoarthritis","Rheumatoid Arthritis","Prostate Cancer","Hip Osteoarthritis","Colorectal Cancer","Atrial Fibrillation","Type 2 Diabetes","Breast Cancer","Gout","Asthma")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigureAFR&AMRAncestriesHRperSDFullSample.png", height=8 , width=8)

######################################################################################################################################

# Leave one out meta-analysis - Only really relevant for European meta-analysis - Each mark can show the effect on the effect size
metaresults <- c()
european <- subset(full, Ancestry=="EUR")

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
  
fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/metaanalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
meta$Phenotype <- factor(meta$Phenotype, levels=c(allmeta[order(allmeta$HR),"Phenotype"]))

metaupdate <- c()
#Test to see if the individual biobank association is significantly different to that of the meta-analysed association
for(i in unique(meta$Phenotype)){
  disease <- subset(meta, Phenotype==i)
  disease$BetaDiff <- disease$Beta - disease$Beta[disease$LeftOutBiobank=="All"]
  disease$SEDiff <- sqrt(disease$SE**2 + disease$SE[disease$LeftOutBiobank=="All"]**2)
  disease$ZDiff <- disease$BetaDiff/disease$SEDiff
  disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
  metaupdate <- rbind(metaupdate, disease)
}

meta <- as.data.frame(metaupdate)

#Plot leave one out results
ggplot(meta) +
  geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_LOOMetaAnalysisEUR.png", height=8 , width=10)

################################################################################################################################################################################################################################################################################################################################
#Drop the All results and resave the file
meta <- subset(meta, Phenotype!="All")
fwrite(meta, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv")

################################################################################################################################################################################################################################################################################################################################
