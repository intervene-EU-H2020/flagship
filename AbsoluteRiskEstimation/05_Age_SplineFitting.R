library(data.table)
library(ggplot2)

biobank <- c("FinnGen", "PartnersBiobank", "HUNT", "EstonianBiobank", "UKBiobank")
filename <- c("/R10_HR_AgeStratified_FinnGen.csv", "/HR_AgeStratified_MGBB_EUR.csv", "/HR_AgeStratified_HUNT.csv","/HR_AgeStratified_EstBB.csv", "/HR_AgeStratified_UKB_AllAncestries.csv")

for(l in 1:length(biobank)){
  
  print(biobank[l])
  
  #Read in the dataset for ages. 
  age <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],filename[l]), data.table=FALSE)
  
  if(biobank[l]=="UKBiobank"){
    age <- age[,-c(1,2)]} else{
      age <- age
    }
  
  colnames(age) <- c("Phenotype", "PRS", "MinAge", "MaxAge", "MedianAAO", "Group", "Controls", "Cases", "Beta", "SE", "Pval", "HR", "CIpos", "CIneg")
  
  phenocols <- c("RHEUMA_SEROPOS_OTH")
  groups <- c("< 20%", "20-40%", "60-80%", "80-90%", "90-95%", "> 95%")
  
  #Subset to one disorder - something that definitely shows an age effect - T2D
  
  for(j in 1:length(phenocols)){
      
    print(phenocols[j])
    
    disease <- subset(age, Phenotype==phenocols[j] & Group %in% groups) 
    
    #Plot 
    disease$Group <- factor(disease$Group, levels=c("< 20%", "20-40%", "60-80%", "80-90%", "90-95%","> 95%"))
    disease$agecat <- paste0(disease$MinAge," - ", disease$MaxAge)
    disease$agecat<- as.factor(disease$agecat)
    
    #Make Betas and relevant confidence intervals from Hazard Ratios
    disease$Beta <- log(disease$HR)
    disease$BCipos <- log(disease$CIpos)
    disease$BCineg <- log(disease$CIneg)
    
    hazards <- c()
    
    for(i in unique(disease$Group)){
      
      print(i)
      
      #Subset to top percentile (Group 11)
      disease_group <- subset(disease, Group==i)
      
      #Fit weighted linear regression - taking weights as inverse of the variance
      
      ##Create weights
      disease_group$weights <- 1/(disease_group$SE**2)
      
      weighted_fit <- lm(Beta~MedianAAO,data=disease_group,weights=weights)
      
      # save predictions of the model in the new data frame 
      # together with variable you want to plot against
      disease_group$predicted <- predict(weighted_fit, disease_group)
      
      print(ggplot(disease_group) +
              geom_point(aes(MedianAAO, Beta), position=position_dodge(width=0.7)) +
              theme_bw() +
              geom_errorbar(aes(x=MedianAAO, ymin = BCineg, ymax = BCipos), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
              ylab("log(HR) (95% CI)") +
              xlab("Age") +
              geom_hline(yintercept = 0) +
              geom_line(color='red', aes(x=MedianAAO, y=predicted)) +
              ggtitle(paste0(phenocols[j]," - Group ",i)) +
              theme(title = element_text(size = 18),
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size = 18),
                    axis.title.x = element_text(size = 18),
                    axis.text.x = element_text(size = 16),
                    axis.title.y = element_text(size = 18),
                    axis.text.y = element_text(size = 16)))
      
      #Construct betas and hazard ratios from predicted trend
      
      ##If the betas are before the first age midpoint or after the last then keep the first predicted beta
      AgeIntervals <- data.frame(MedianAAO = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
      pred <- predict(weighted_fit, newdata=AgeIntervals)
      
      hazard <- cbind(AgeIntervals, pred)
      colnames(hazard) <- c("Age", "Beta")
      
      hazard$Beta <- ifelse(hazard$Age < disease_group$MedianAAO[1], disease_group$predicted[1], hazard$Beta) 
      hazard$Beta <- ifelse(hazard$Age > disease_group$MedianAAO[4], disease_group$predicted[4], hazard$Beta) 
      hazard$HR <- exp(hazard$Beta)
      hazard$Group <- paste0(i)
      
      hazards <- rbind(hazards, hazard)
      
    }
    
    write.csv(hazards, paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/Age/",phenocols[j],"_Age_Specific_Hazards.csv"))
    
  }
  
}

##########################################################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################################################

#Supplementary analysis - Finngen extension to bottom and top 1%

library(data.table)
library(ggplot2)

biobank <- c("FinnGen")
filename <- c("/R10_HR_AgeStratified_FinnGen.csv")

for(l in 1:length(biobank)){
  
  print(biobank[l])
  
  #Read in the dataset for ages. 
  age <- fread(paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],filename[l]), data.table=FALSE)
  colnames(age) <- c("Phenotype", "PRS", "MinAge", "MaxAge", "MedianAAO", "Group", "Controls", "Cases", "Beta", "SE", "Pval", "HR", "CIpos", "CIneg")
  
  phenocols <- c("RHEUMA_SEROPOS_OTH")
  groups <- c("< 5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%")
  
  #Subset to one disorder - something that definitely shows an age effect - T2D
  
  for(j in 1:length(phenocols)){
      
    print(phenocols[j])
    
    disease <- subset(age, Phenotype==phenocols[j] & Group %in% groups) 
    
    #Plot 
    disease$Group <- factor(disease$Group, levels=c("< 5%", "5-10%", "10-20%", "20-40%", "60-80%", "80-90%", "90-95%", "95-99%", "> 99%"))
    disease$agecat <- paste0(disease$MinAge," - ", disease$MaxAge)
    disease$agecat<- as.factor(disease$agecat)
    
    #Make Betas and relevant confidence intervals from Hazard Ratios
    disease$Beta <- log(disease$HR)
    disease$BCipos <- log(disease$CIpos)
    disease$BCineg <- log(disease$CIneg)
    
    hazards <- c()
    
    for(i in unique(disease$Group)){
      
      print(i)
      
      #Subset to top percentile (Group 11)
      disease_group <- subset(disease, Group==i)
      
      #Fit weighted linear regression - taking weights as inverse of the variance
      
      ##Create weights
      disease_group$weights <- 1/(disease_group$SE**2)
      
      weighted_fit <- lm(Beta~MedianAAO,data=disease_group,weights=weights)
      
      # save predictions of the model in the new data frame 
      # together with variable you want to plot against
      disease_group$predicted <- predict(weighted_fit, disease_group)
      
      print(ggplot(disease_group) +
              geom_point(aes(MedianAAO, Beta), position=position_dodge(width=0.7)) +
              theme_bw() +
              geom_errorbar(aes(x=MedianAAO, ymin = BCineg, ymax = BCipos), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
              ylab("log(HR) (95% CI)") +
              xlab("Age") +
              geom_hline(yintercept = 0) +
              geom_line(color='red', aes(x=MedianAAO, y=predicted)) +
              ggtitle(paste0(phenocols[j]," - Group ",i)) +
              theme(title = element_text(size = 18),
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size = 18),
                    axis.title.x = element_text(size = 18),
                    axis.text.x = element_text(size = 16),
                    axis.title.y = element_text(size = 18),
                    axis.text.y = element_text(size = 16)))
      
      #Construct betas and hazard ratios from predicted trend
      
      ##If the betas are before the first age midpoint or after the last then keep the first predicted beta
      AgeIntervals <- data.frame(MedianAAO = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5))
      pred <- predict(weighted_fit, newdata=AgeIntervals)
      
      hazard <- cbind(AgeIntervals, pred)
      colnames(hazard) <- c("Age", "Beta")
      
      hazard$Beta <- ifelse(hazard$Age < disease_group$MedianAAO[1], disease_group$predicted[1], hazard$Beta) 
      hazard$Beta <- ifelse(hazard$Age > disease_group$MedianAAO[4], disease_group$predicted[4], hazard$Beta) 
      hazard$HR <- exp(hazard$Beta)
      hazard$Group <- paste0(i)
      
      hazards <- rbind(hazards, hazard)
      
    }
    
    write.csv(hazards, paste0("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/",biobank[l],"/Age/",phenocols[j],"_Age_Specific_Hazards_OnePercentTails.csv"))
    
  }
  
}
