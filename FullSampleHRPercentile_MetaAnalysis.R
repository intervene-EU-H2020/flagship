#When doing the meta-analysis check all the forest plots to understand which is driving the analysis. 
library(data.table)
library(ggplot2)
library(dplyr)
library(metafor)

#Meta-analysis of top 1% only - only keep when phenotype has a result for all percentiles - if not drop

finngen <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/HR_FullSample_FinnGen.csv", data.table=FALSE)
finngen$Biobank <- "FinnGen"

gs <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/GenerationScotland/HR_FullSampleGS.csv", data.table=FALSE)
gs$Biobank <- "Generation Scotland"
colnames(gs) <- colnames(finngen)
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT","F5_DEPRESSIO","C3_MELANOMA_SKIN","C3_PROSTATE","RHEUMA_SEROPOS_OTH","T1D","C3_BRONCHUS_LUNG")
gs <- subset(gs, !(Phenotype %in% dropgs))

estbb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/EstonianBiobank/HR_FullSample_EstBB.csv", data.table=FALSE)
estbb <- estbb[,-1]
estbb$Biobank <- "Estonian Biobank"
colnames(estbb) <- colnames(finngen)

mgb <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/PartnersBiobank/HR_FullSample_MGBB.csv", data.table=FALSE)
mgb$Biobank <- "Mass General Brigham"
dropmgb <- c("I9_AF","T1D")
mgb <- subset(mgb, !(Phenotype %in% dropmgb))

#gnh <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/Genes_and_Health/HR_FullSample_GNH.csv", data.table=FALSE)
#gnh <- gnh[,-1]
#gnh$Biobank <- "Genes & Health"
#gnh$Ancestry <- "SAS"
#colnames(gnh) <- colnames(finngen)

#bbj <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/BiobankJapan/HR_FullSampleBBJ.csv", data.table=FALSE)
#bbj <- bbj[,-1]
#bbj$Biobank <- "Biobank Japan"
#bbj$Ancestry <- "EAS"
#colnames(bbj) <- colnames(finngen)
#dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
#bbj <- subset(bbj, !(Phenotype %in% dropbbj))

all <- rbind(finngen, gs) %>%
        rbind(., estbb) %>%
          rbind(., mgb)

#######################################################################################################################################################

#Want to perform meta-analysis on each percentile in turn
percentile <- c("< 1%","1-5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","95-99%","> 99%","< 5%","> 95%","< 10%","> 90%","< 20%","> 80%")

metaresults <- c()
 
for(i in unique(all$Phenotype)){
    
  for(k in percentile){
    
    disease <- subset(all, Phenotype==i & Group==k)
    
    print(i)
    if(dim(disease)[1]==0){
      next
    }
    
    #Meta analysis should be done at the beta level and stratified by ancestry 
    meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
    
    #Save the forest plot so can compare with the table
    
    metaresults <- rbind(metaresults, c(i, k, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
  
  }  

}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","Group","Beta","SE","Pval","Heterogeneity","HetPval")

fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/FullSample/HRPercentiles/MetaAnalysis.csv")


