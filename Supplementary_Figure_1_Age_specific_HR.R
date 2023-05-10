library(data.table)
library(ggplot2)

#Read in the dataset for ages. 
age <- fread("/Users/jermy/Documents/INTERVENE/Results/HazardRatios/FinnGen/HRperSD_AgeStratified_FinnGen.csv", data.table=FALSE)
colnames(age) <- c("Phenotype", "PRS", "MinAge", "MaxAge", "MedianAAO", "Controls", "Cases", "Beta", "SE", "Pval", "HR", "CIpos", "CIneg")

disease <- subset(age, Phenotype=="T2D") 
  
#Plot 
disease$agecat <- paste0(disease$MinAge," - ", disease$MaxAge)
disease$agecat<- as.factor(disease$agecat)

#Make Betas and relevant confidence intervals from Hazard Ratios
disease$BCipos <- log(disease$CIpos)
disease$BCineg <- log(disease$CIneg)

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

hazard$Beta <- ifelse(hazard$Age < disease$MedianAAO[1], disease$predicted[1], hazard$Beta) 
hazard$Beta <- ifelse(hazard$Age > disease$MedianAAO[4], disease$predicted[4], hazard$Beta) 

ggplot() +
  geom_point(data=disease, aes(MedianAAO, Beta)) +
  theme_bw() +
  geom_errorbar(data=disease, aes(x=MedianAAO, ymin = BCineg, ymax = BCipos), width=0.5) +
  ylab("log(HR) (95% CI)") +
  xlab("Age") +
  geom_hline(yintercept = 0) +
  geom_line(data=hazard, color='red', aes(x=Age, y=Beta)) +
            ggtitle("Type 2 Diabetes - Age Interval Associations") +
            theme(title = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 16))
#ggsave("/Users/jermy/Documents/INTERVENE/Write-up/Figures/Supplementary_Figure_1_AgeSpecificHRs.png", height=7, width=7)
    
