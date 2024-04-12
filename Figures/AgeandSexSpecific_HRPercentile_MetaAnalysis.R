library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)
library(ggtext) 

#UKB
ukb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/UKBiobank/HR_AgeandSexStratified_UKB_AllAncestries.csv", data.table=FALSE)
ukb$Biobank <- "UK Biobank"
colnames(ukb) <- tolower(colnames(ukb))

#FinnGen
finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/R10_HR_AgeandSex_Stratified_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"
colnames(finngen) <- tolower(colnames(finngen))
finngen <- finngen[,names(ukb)]

#Estonian Biobank
estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_AgeandSexStratified_EstBB.csv", data.table=FALSE)
estbb$Biobank <- "Estonian Biobank"
estbb$Ancestry <- "EUR"
colnames(estbb) <- tolower(colnames(estbb))
estbb <- estbb[,names(ukb)]

#Generation Scotland
gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HR_AgeandSexStratified_GS.csv", data.table=FALSE)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs$Biobank <- "Generation Scotland"
gs$Ancestry <- "EUR"
colnames(gs) <- tolower(colnames(gs))
gs <- gs[,names(ukb)]

#Mass General Brigham - European
mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HR_AgeandSexStratified_MGBB_EUR.csv", data.table=FALSE)
mgb$Beta <- ifelse(mgb$Phenotype=="C3_CANCER", mgb$Beta*-1, mgb$Beta)
mgb$HR <- exp(mgb$Beta)
mgb$Cipos <- exp(mgb$Beta + 1.96*mgb$SE)
mgb$Cineg <- exp(mgb$Beta - 1.96*mgb$SE)
dropmgb <- c("I9_AF")
mgb <- subset(mgb, !(Phenotype %in% dropmgb))
mgb$Biobank <- "Mass General Brigham"
mgb$Ancestry <- "EUR"
colnames(mgb) <- tolower(colnames(mgb))
mgb <- mgb[,names(ukb)]

#Genomics England
ge <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenomicsEngland/HR_AgeandSexStratified_GenomicsEngland.csv", data.table=FALSE)
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
colnames(ge) <- tolower(colnames(ge))
ge <- ge[,names(ukb)]

#HUNT
hunt <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/HUNT/HR_AgeandSexStratifiedHUNT.csv", data.table=FALSE)
drophunt <- c("I9_CHD")
hunt <- subset(hunt, !(Phenotype %in% drophunt))
hunt$Biobank <- "HUNT"
hunt$Ancestry <- "EUR"
colnames(hunt) <- tolower(colnames(hunt))
hunt <- hunt[,names(ukb)]

#Combine into one dataset
all <- rbind(ukb, finngen) %>%
        rbind(., estbb) %>%
          rbind(., gs) %>%
            rbind(., mgb) %>% 
              rbind(.,ge) %>%
                rbind(., hunt)

all$cipos[all$cipos==Inf] <- NA

percentile <- c("< 20%","20-40%","60-80%","80-90%","90-95%","> 95%")
q <- c()
for(i in c("J10_ASTHMA","I9_AF","I9_CHD","GOUT","F5_DEPRESSIO","KNEE_ARTHROSIS","T2D")){
  for(j in unique(all$biobank)){
    for(k in unique(all$sex)){
      disease <- subset(all, group %in% percentile & phenotype==i & ancestry=="EUR" & biobank==j & sex==k & group %in% percentile)
      
      if(any(is.na(disease[,c(11:16)])) | dim(disease)[1] < 24){
        next
      }
      
      q <- rbind(q,disease)
    }
  }
}

for(i in c("J10_ASTHMA","I9_AF","I9_CHD","GOUT","KNEE_ARTHROSIS","T2D")){

  disease <- subset(q, phenotype==i & ancestry=="EUR" & group %in% percentile)
  print(i)
  print(dim(disease)[1]/48)
  print(table(disease$biobank))
  
}

#Meta-analyse Q by quartile
#percentile <- c("< 20%","> 95%")

percentile <- c("< 20%","20-40%","60-80%","80-90%","90-95%","> 95%")

metaresults <- c()

for(i in unique(q$phenotype)){
  
  print(i)
  
  for(j in unique(q$quartile)){
    for(k in percentile){
      for(l in unique(q$sex)){
      
        print(j)
        disease <- subset(q, phenotype==i & ancestry=="EUR" & quartile==j & group==k & sex==l)
        
        #Meta analysis should be done at the beta level and stratified by ancestry 
        metaFE <- rma(yi=beta, sei=se, data=disease, method="FE")
        
        results <- c(i, j, k, l, metaFE$b, metaFE$se, metaFE$pval, metaFE$QE, metaFE$QEp)
        
        metaresults <- rbind(metaresults, results)
      }
    }
  }
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","Quartile","Group","Sex","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(5:9),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/HRPercentiles/MetaAnalysisPercentilesAgeandSexStratified.csv")

#Plot Hazard Ratios similar to original figure but with hazard ratios percentiles top and bottom and a log scale
metas <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/HRPercentiles/MetaAnalysisPercentilesAgeandSexStratified.csv", data.table=FALSE)

#More extensive subsetting required here...
#metas <- subset(metas, Phenotype %in% c("J10_ASTHMA","I9_AF","T2D") | (Phenotype=="I9_CHD" & Sex=="male") | (Phenotype=="GOUT" & Sex=="male") | (Phenotype=="KNEE_ARTHROSIS" & Sex=="female"))
metas <- subset(metas, Phenotype=="T2D"|Phenotype=="I9_CHD")

#Read in prostate cancer for the main plot
prostate <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/HRPercentiles/MetaAnalysisPercentilesAgeStratified.csv", data.table=FALSE)
prostate <- subset(prostate, Phenotype=="C3_PROSTATE" | Phenotype=="C3_BREAST")
prostate$Sex <- case_when(prostate$Phenotype=="C3_PROSTATE" ~ "male",
                          prostate$Phenotype=="C3_BREAST" ~ "female")
prostate <- prostate[,names(metas)]

metas <- rbind(metas,prostate)
#metas$Phenotype <- factor(metas$Phenotype, levels=c())
metas$Quartile <- factor(metas$Quartile, levels=c(1,2,3,4))
metas$Group <- factor(metas$Group, levels=c("< 20%","20-40%","40-60%","60-80%","80-90%","90-95%","> 95%"))

#Main Plot 
ggplot(metas) +
  geom_point(aes(interaction(Sex,Phenotype), HR, col=Group, shape=Quartile, group=Quartile), position=position_dodge(width=0.7), size=3) +
  theme_bw() +
  geom_errorbar(aes(x=interaction(Sex,Phenotype), ymin = Cineg, ymax = Cipos, col=Group, group=Quartile), position=position_dodge(width=0.7), size=0.75, width=0.125) +
  ylab("Meta-Analyzed Hazard Ratio (95% CI)") +
  xlab("") +
  geom_vline(xintercept = seq(2.5,(length(unique(metas$Phenotype))*2)+0.5,by=2), col='black') + 
  geom_hline(yintercept = 1.0) +
  scale_shape_manual(name="Age Quartiles", values=c(15,16,17,8), labels=c("Youngest", " |", "V", "Oldest")) +
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
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/HR_AgeandSexStratified_Percentiles_AllPercentiles.png", height=13 , width=13, dpi=300)

#################################################################################################################################################################################
#################################################################################################################################################################################

#Forest plot of top 5% and bottom 20% for CHD
#Create forest plots for each pheno - facet_wrap
chd <- subset(q, phenotype=="I9_CHD" & group %in% c("< 20%","90-95%","> 95%") & sex=="male")

chd$biobank <- as.factor(chd$biobank)
chd$quartile <- factor(chd$quartile, levels=c(4,3,2,1))
chd$group <- factor(chd$group, levels=c("< 20%","90-95%","> 95%"))

ggplot(chd) +
  geom_point(aes(biobank, hr, group=interaction(group,quartile), col=group, shape=quartile), position=position_dodge(width=0.7)) +
  theme_bw() +
  xlab(c("")) + 
  ylab(c("Hazard Ratio per Standard Deviation (95% CI)")) + 
  scale_shape_manual(name="Age Quartiles", values=c(15,16,17,8), labels=c("Oldest", "V", " |", "Youngest")) +
  scale_color_manual(name="PGS Strata", values=c("#117733", "#332288" ,"#AA4499")) +
  geom_errorbar(aes(x=biobank, ymin = cineg, ymax = cipos, group=interaction(group,quartile), col=group), size=0.5, width=0.25, position=position_dodge(width=0.7)) +
  coord_flip() + 
  scale_y_continuous(trans="log2") +    
  guides(color = guide_legend(reverse = TRUE)) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) 
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigures/CHD_Percentile_ForestPlot.png", dpi=300, height=8, width=8)
