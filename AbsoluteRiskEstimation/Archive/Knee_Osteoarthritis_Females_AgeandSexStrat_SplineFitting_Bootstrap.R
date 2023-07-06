library(data.table)
library(ggplot2)
library(dplyr)

#biobank <- c("FinnGen", "PartnersBiobank", "HUNT", "GenomicsEngland", "GenerationScotland", "EstonianBiobank", "UKMetaAnalysis")
#filename <- c("/R10_HR_AgeandSex_Stratified_FinnGen.csv", "/HR_AgeandSexStratified_MGBB_EUR.csv", "/HR_AgeandSexStratifiedHUNT.csv","/HR_AgeandSexStratified_GenomicsEngland.csv", "/HR_AgeandSexStratified_GS.csv", "/HR_AgeandSexStratified_EstBB.csv", "/UKMetaAnalysisPercentilesAgeandSexStratified.csv")

biobank <- c("PartnersBiobank")
filename <- c("/HR_AgeandSexStratified_MGBB_EUR.csv")

for(l in 1:length(biobank)){
  
  print(biobank[l])
  
  #Read in the dataset for ages. 
  age <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],filename[l]), data.table=FALSE)
  if(biobank[l]=="UKMetaAnalysis"){
    colnames(age) <- c("Phenotype", "Quartile", "Group", "Sex", "MedianAAO", "Beta", "SE", "Pval", "QHet", "HetPval", "HR", "CIpos", "CIneg")
  } else{
    colnames(age) <- c("Phenotype", "PRS", "Sex", "MinAge", "MaxAge", "MedianAAO", "Group", "Controls", "Cases", "Beta", "SE", "Pval", "HR", "CIpos", "CIneg", "Quartile")
  }
  
  phenocols <- c("KNEE_ARTHROSIS")
  groups <- c("< 20%", "20-40%", "60-80%", "80-90%", "90-95%", "> 95%")
  
  #Subset to one disorder - something that definitely shows an age effect - T2D
  
  for(j in 1:length(phenocols)){
    for(k in c("female")){
      
      print(phenocols[j])
      
      disease <- subset(age, Phenotype==phenocols[j] & Sex==k & Group %in% groups) 
      
      #Plot 
      disease$Group <- factor(disease$Group, levels=c("< 20%", "20-40%", "60-80%", "80-90%", "90-95%","> 95%"))
      disease$agecat <- paste0(disease$MinAge," - ", disease$MaxAge)
      disease$agecat<- as.factor(disease$agecat)
      
      #Make Betas and relevant confidence intervals from Hazard Ratios
      disease$Beta <- log(disease$HR)
      disease$BCipos <- log(disease$CIpos)
      disease$BCineg <- log(disease$CIneg)
      
      b <- 0 
      
      print(k)
      
      bootstrapped_results <- data.frame(Age = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
      
      while(b < 1000){
        
        hazards <- c()
        
        for(i in unique(disease$Group)){
          
          print(i)
          
          #Subset to top percentile (Group 11)
          disease_group <- subset(disease, Group==i)
          
          #Fit weighted linear regression - taking weights as inverse of the variance
          
          ##Create weights
          disease_group[[paste0("Beta_",b)]] <- with(disease_group, rnorm(nrow(disease_group), disease_group$Beta, disease_group$SE))
          
          disease_group$weights <- 1/(disease_group$SE**2)
          
          weighted_fit <- lm(as.formula(paste0("Beta_",b," ~ MedianAAO")), data=disease_group, weights=weights)
          
          # save predictions of the model in the new data frame 
          # together with variable you want to plot against
          disease_group$predicted <- predict(weighted_fit, disease_group)
          
          #Construct betas and hazard ratios from predicted trend
          AgeIntervals <- data.frame(MedianAAO = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
          pred <- predict(weighted_fit, newdata=AgeIntervals)
          
          hazard <- cbind(AgeIntervals, pred)
          
          colnames(hazard) <- c("Age", paste0("Beta_",b))
          
          hazard[[paste0("Beta_",b)]] <- ifelse(hazard$Age < disease_group$MedianAAO[1], disease_group$predicted[1], hazard[[paste0("Beta_",b)]]) 
          hazard[[paste0("Beta_",b)]] <- ifelse(hazard$Age > disease_group$MedianAAO[4], disease_group$predicted[4], hazard[[paste0("Beta_",b)]]) 
          hazard$Group <- paste0(i)
          
          hazards <- rbind(hazards, hazard)
          
        }
        
        bootstrapped_results <- inner_join(hazards, bootstrapped_results)
        
        b <- b+1
        
      }
      
      bootstrapped_betas <- as.matrix(bootstrapped_results[,-c(1,3)])
      confidenceintervals <- apply(bootstrapped_betas, 1, quantile, c(0.025, 0.975))
      
      bootstrapped_results <- bootstrapped_results[,c(1,3)]
      
      bootstrapped_results$CIneg <- confidenceintervals[1,]
      bootstrapped_results$CIpos <- confidenceintervals[2,]
      
      #Add in actual lifetime risks
      actual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/AgeandSex/",phenocols[j],"_",k,"_Age_Specific_Hazards.csv"), select=c("Age","Beta","Group"), data.table=FALSE)
      bootstrapped_results <- inner_join(bootstrapped_results, actual)
      
      write.csv(bootstrapped_results, paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/AgeandSex/",phenocols[j],"_",k,"_Age_Specific_Hazards_Bootstrapped.csv"))
      
    }
    
  }
}


##########################################################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################################################

#Supplementary analysis - Finngen extension to bottom and top 1%

library(data.table)
library(ggplot2)

biobank <- c("FinnGen")
filename <- c("/R10_HR_AgeandSex_Stratified_FinnGen.csv")

for(l in 1:length(biobank)){
  
  print(biobank[l])
  
  #Read in the dataset for ages. 
  age <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],filename[l]), data.table=FALSE)
  colnames(age) <- c("Phenotype", "PRS", "Sex", "MinAge", "MaxAge", "MedianAAO", "Group", "Controls", "Cases", "Beta", "SE", "Pval", "HR", "CIpos", "CIneg")
  
  phenocols <- c("KNEE_ARTHROSIS")
  groups <- c("< 1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%")
  
  #Subset to one disorder - something that definitely shows an age effect - T2D
  
  for(j in 1:length(phenocols)){
    for(k in c("female")){
      
      print(phenocols[j])
      
      disease <- subset(age, Phenotype==phenocols[j] & Sex==k & Group %in% groups) 
      
      #Plot 
      disease$Group <- factor(disease$Group, levels=c("< 1%", "1-5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%"))
      disease$agecat <- paste0(disease$MinAge," - ", disease$MaxAge)
      disease$agecat<- as.factor(disease$agecat)
      
      #Make Betas and relevant confidence intervals from Hazard Ratios
      disease$Beta <- log(disease$HR)
      disease$BCipos <- log(disease$CIpos)
      disease$BCineg <- log(disease$CIneg)
      
      b <- 0 
      
      print(k)
      
      bootstrapped_results <- data.frame(Age = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
      
      while(b < 1000){
        
        hazards <- c()
        
        for(i in unique(disease$Group)){
          
          print(i)
          
          #Subset to top percentile (Group 11)
          disease_group <- subset(disease, Group==i)
          
          #Fit weighted linear regression - taking weights as inverse of the variance
          
          ##Create weights
          disease_group[[paste0("Beta_",b)]] <- with(disease_group, rnorm(nrow(disease_group), disease_group$Beta, disease_group$SE))
          
          disease_group$weights <- 1/(disease_group$SE**2)
          
          weighted_fit <- lm(as.formula(paste0("Beta_",b," ~ MedianAAO")), data=disease_group, weights=weights)
          
          # save predictions of the model in the new data frame 
          # together with variable you want to plot against
          disease_group$predicted <- predict(weighted_fit, disease_group)
          
          #Construct betas and hazard ratios from predicted trend
          AgeIntervals <- data.frame(MedianAAO = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
          pred <- predict(weighted_fit, newdata=AgeIntervals)
          
          hazard <- cbind(AgeIntervals, pred)
          
          colnames(hazard) <- c("Age", paste0("Beta_",b))
          
          hazard[[paste0("Beta_",b)]] <- ifelse(hazard$Age < disease_group$MedianAAO[1], disease_group$predicted[1], hazard[[paste0("Beta_",b)]]) 
          hazard[[paste0("Beta_",b)]] <- ifelse(hazard$Age > disease_group$MedianAAO[4], disease_group$predicted[4], hazard[[paste0("Beta_",b)]]) 
          hazard$Group <- paste0(i)
          
          hazards <- rbind(hazards, hazard)
          
        }
        
        bootstrapped_results <- inner_join(bootstrapped_results, hazards)
        
        b <- b+1
        
      }
      
      bootstrapped_betas <- as.matrix(bootstrapped_results[,-c(1,3)])
      confidenceintervals <- apply(bootstrapped_betas, 1, quantile, c(0.025, 0.975))
      
      bootstrapped_results <- bootstrapped_results[,c(1,3)]
      
      bootstrapped_results$CIneg <- confidenceintervals[1,]
      bootstrapped_results$CIpos <- confidenceintervals[2,]
      
      #Add in actual lifetime risks
      actual <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/AgeandSex/",phenocols[j],"_",k,"_Age_Specific_Hazards_OnePercentTails.csv"), select=c("Age","Beta","Group"), data.table=FALSE)
      bootstrapped_results <- inner_join(bootstrapped_results, actual)
      
      write.csv(bootstrapped_results, paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/AgeandSex/",phenocols[j],"_",k,"_Age_Specific_Hazards_Bootstrapped_OnePercentTails.csv"))
      
    }
    
  }
}
