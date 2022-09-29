library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)

#Read in age stratified hazard ratios

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/HRperSD_AgeStratified_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HRperSD_AgeStratified_UKB_AllAncestries.csv", data.table=FALSE)
ukb <- ukb[,names(finngen)]

#Genes&Health
gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HRperSD_AgeStratified_GNH.csv", data.table=FALSE)
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

#Mass General Brigham
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeStratified_MGBB.csv", data.table=FALSE)

dropmgb <- c("I9_AF")
mgb <- subset(mgb, !(Phenotype %in% dropmgb))
mgb <- mgb[,names(finngen)]

#Genomics England
#ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HRperSD_AgeStratified_MGBB.csv", data.table=FALSE)
#mgb <- subset(mgb, !(Phenotype %in% dropmgb))
#mgb <- mgb[,names(finngen)]

#Combine into one dataset
all <- rbind(finngen, gnh) %>%
        rbind(., estbb) %>% 
          rbind(., gs) %>%
            rbind(., mgb) %>%
              rbind(., bbj) %>%
                rbind(., ukb)

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
  
  for(j in c("EUR","SAS","EAS")){
    
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
      meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
      
      png(file=paste0("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/ForestPlots/Quartile_",k,"_",i,"_",j,"_Forest.png"), width = 7,
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
colnames(metaresults) <- c("Ancestry","Quartile","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(4:8),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

#fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratifiedMetaAnalysisFE.csv")

######################################################################################################################################

all <- all[,c("Ancestry","Quartile","Phenotype","Beta","SE","Biobank")]
diffs <- metaresults[,c("Ancestry","Quartile","Phenotype","Beta","SE")]
diffs$Biobank <- "All"
diffs <- rbind(diffs, all)

differences <- c()
for(j in c("EUR","SAS","EAS")){
  for(k in c(1,2,3,4)){
    for(i in unique(diffs$Phenotype)){
      
      disease <- subset(diffs, Phenotype==i & Ancestry==j & Quartile==k)
      
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

fwrite(differences, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/Age_BiobankVariation.csv")


######################################################################################################################################

# Main analysis
meta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratifiedMetaAnalysisFE.csv", data.table=FALSE)
meta$Ancestry <- factor(meta$Ancestry, levels=c("EUR","SAS","EAS"))
orderPheno <- subset(meta, Ancestry=="EUR" & Quartile==1)
meta$Phenotype <- factor(meta$Phenotype, levels=c(meta[order(orderPheno$HR),"Phenotype"]))
meta$Quartile <- as.factor(meta$Quartile)

ancestry <- c("EUR","SAS","EAS")
europelabels <- list(c("ILD","Appendicitis","All Cancers","Lung Cancer","Epilepsy","Major Depression","Hip Osteoarthritis","Knee Osteoarthritis","Skin Melanoma","CHD","Colorectal Cancer","Asthma","Breast Cancer","Gout","Atrial Fibrillation","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes"),
                     c("ILD","Appendicitis","All Cancers","Lung Cancer","Epilepsy","Major Depression","Hip Osteoarthritis","Knee Osteoarthritis","CHD","Asthma","Breast Cancer","Gout","Atrial Fibrillation","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes"),
                     c("All Cancers","Epilepsy","CHD","Colorectal Cancer","Asthma","Breast Cancer","Type 2 Diabetes","Prostate Cancer")
                     )

for(i in 1:3){
  
  anc <- subset(meta, Ancestry==ancestry[i])

  #Plot
  print(ggplot(anc) +
    geom_point(aes(Phenotype, HR, group=Quartile, col=Quartile), position=position_dodge(width=0.7)) +
    theme_bw() +
    geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Quartile, col=Quartile), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
    ylab("Hazard Ratio (95% CI)") +
    xlab("") +
    geom_hline(yintercept = 1.0) +
    scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#000000"), labels=c("1 (youngest)", "2", "3", "4 (oldest)")) +
    scale_x_discrete(labels=europelabels[i]) + 
    theme(title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 14)) +
    coord_flip())
  ggsave(paste0("/Users/jermy/Documents/INTERVENE/Write-up/Figure3_",ancestry[i],".png"), height=8, width=8)
  
}


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
    scale_x_discrete(labels=c("ILD","Appendicitis","All Cancers","Lung Cancer","Epilepsy","Major Depression","Hip Osteoarthritis","Knee Osteoarthritis","Skin Melanoma","CHD","Colorectal Cancer","Asthma","Breast Cancer","Gout","Atrial Fibrillation","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
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
