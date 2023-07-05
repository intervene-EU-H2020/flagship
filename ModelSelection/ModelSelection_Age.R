#Focus on european to begin with
library(data.table)
library(ggplot2)
library(metafor)

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HRperSD_AgeStratified_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HRperSD_AgeStratified_UKB_AllAncestries.csv", data.table=FALSE)
ukb <- ukb[,names(finngen)]

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HRperSD_AgeStratified_EstBB.csv", data.table=FALSE)
estbb <- estbb[,names(finngen)]

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HRperSD_AgeStratified_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs <- gs[,names(finngen)]

#Mass General Brigham
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
hunt <- hunt[,names(finngen)]

#Combine
all <- rbind(finngen, ukb) %>%
        rbind(., mgb) %>%
          rbind(., gs) %>%
            rbind(., estbb) %>%
              rbind(., ge) %>%
                rbind(., hunt)

all$Quartile <- factor(all$Quartile, levels=c(1,2,3,4))

#Make sure all quartiles have been tested before including in the meta-analysis
q <- c()
for(i in unique(all$Phenotype)){
  for(k in unique(all$Biobank)){
      pheno <- subset(all, Phenotype==i & Biobank==k)
      print(pheno)
      
      if(any(is.na(pheno[,c(8:13)])) | dim(pheno)[1]<4){
        next
      }
      
      q <- rbind(q,pheno)
  }
}

all <- as.data.frame(q)

#Meta-analyse MedianAAO for each quartile - add to meta analysis results - only meta-analyse if none of the quartiles are NA
metaresults <- c()
for(k in c(1,2,3,4)){
  
  quartile <- subset(all, Quartile==k)
    
  for(i in unique(quartile$Phenotype)){
    disease <- subset(quartile, Phenotype==i & !(is.na(Beta)))
    
    if(dim(disease)[1]==0){
      next
    }
    
    #Meta analysis should be done at the beta level and stratified by ancestry 
    meta <- rma(yi=MedianAAO, sei=SE, data=disease, method="FE")
    
    metaresults <- rbind(metaresults, c(i, k, meta$b))
    
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","Quartile","MedianAAO")
metaresults$Quartile <- as.numeric(metaresults$Quartile)
metaresults$MedianAAO <- as.numeric(metaresults$MedianAAO)

#Check makes sense - it does - great
#test <- sum((finngen[1,"MedianAAO"]/(finngen[1,"SE"]**2)), (estbb[1,"MedianAAO"]/(estbb[1,"SE"]**2)), (mgb[1,"MedianAAO"]/(mgb[1,"SE"]**2)), (gs[1,"MedianAAO"]/(gs[1,"SE"]**2))) / sum((1/(finngen[1,"SE"]**2)), (1/(estbb[1,"SE"]**2)), (1/(mgbb[1,"SE"]**2)), (1/(gs[1,"SE"]**2)))

#Read in the dataset for ages. 
age <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeStratifiedMetaAnalysisFEandRE.csv", data.table=FALSE)
age <- subset(age, Test=="Fixed Effect" & Ancestry=="EUR")
age <- left_join(age, metaresults)
age$BCipos <- log(age$Cipos)
age$BCineg <- log(age$Cineg)

results <- c()
hazards <- c()
for(i in unique(age$Phenotype)){
  
  disease <- subset(age, Phenotype==i) 
  
  #Make Betas and relevant confidence intervals from Hazard Ratios
  disease$BCipos <- log(disease$Cipos)
  disease$BCineg <- log(disease$Cineg)
  
  #Fit weighted linear regression - taking weights as inverse of the variance
  
  ##Create weights
  disease$weights <- 1/(disease$SE**2)
  
  weighted_fit <- lm(Beta~MedianAAO,data=disease,weights=weights)
  
  disease$predicted <- predict(weighted_fit, disease)
  
  ##If the betas are before the first age midpoint or after the last then keep the first predicted beta
  AgeIntervals <- data.frame(MedianAAO = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 47.89, 52.5, 57.5, 58.96, 62.5, 66.67, 67.5, 72.5, 74.51, 77.5))
  pred <- predict(weighted_fit, newdata=AgeIntervals)
  
  hazard <- cbind(AgeIntervals, pred)
  colnames(hazard) <- c("Age", "Beta")
  hazard$Phenotype <- i
  
  hazard$Beta <- ifelse(hazard$Age < disease$MedianAAO[1], disease$predicted[1], hazard$Beta) 
  hazard$Beta <- ifelse(hazard$Age > disease$MedianAAO[4], disease$predicted[4], hazard$Beta) 
  hazards <- rbind(hazards, hazard)
  
  results <- rbind(results, c(i, summary(weighted_fit)$coeff["MedianAAO","Estimate"], summary(weighted_fit)$coeff["MedianAAO","Std. Error"], summary(weighted_fit)$coeff["MedianAAO","Pr(>|t|)"]))
}

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
  `I9_CHD` = "CHD",
  #`ILD` = "ILD",
  `J10_ASTHMA` = "Asthma",
  `K11_APPENDACUT` = "Appendicitis",
  `KNEE_ARTHROSIS` = "Knee Osteoarthritis",
  `RHEUMA_SEROPOS_OTH` = "Rheumatoid Arthritis",
  `T1D` = "Type 1 Diabetes",
  `T2D` = "Type 2 Diabetes"
)

age <- subset(age, Phenotype!="ILD")
hazards <- subset(hazards, Phenotype!="ILD")

#Create facet_wrapped version of the age specific effects.
ggplot() +
  geom_point(data=age, aes(MedianAAO, Beta)) +
  theme_bw() +
  geom_errorbar(data=age, aes(x=MedianAAO, ymin = BCineg, ymax = BCipos), width=0.5) +
  ylab("log(HR) (95% CI)") +
  xlab("Age") +
  geom_hline(yintercept = 0) +
  geom_line(data=hazards, color='red', aes(x=Age, y=Beta)) +
  theme(title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) + 
  facet_wrap(~ Phenotype, labeller = as_labeller(phenotypes))
#ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/Supplementary_Figure_6_noILD.png",  height=8, width=10)

results <- as.data.frame(results)
colnames(results) <- c("Phenotype","Beta", "SE", "Pval")

fwrite(results, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/AgeStratified/AgeSpecificEffects_Test.csv")
