#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

results <- c()

for(i in 1:length(phenocols)){
  
  print(phenocols[i])
  print(prscols[i])
  
  #Read in phenotype file
  pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"SEX","END_OF_FOLLOWUP"), data.table=FALSE)
  
  pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
  
  #Read in PRS scores
  PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)
  
  #Subset columns to the IDs and score only. Note: columns FID or IID may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c("#FID","IID","SCORE1_SUM")]
    
  #Rename ID column to the name of the ID column in the phenotype file
  colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))
  
  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS)
  
  pheno <- subset(pheno, !(is.na(pheno[[paste0(phenocols[i])]]) | is.na(pheno[[paste0(prscols[i],"_prs")]])))
  
  #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
  #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
  #Feel free to subset using your own code: only provided as a reminder.
  pheno <- subset(pheno, ANCESTRY=='EUR')
  
  pheno[[paste0(prscols[i],"_prs")]] <- scale(pheno[[paste0(prscols[i],"_prs")]])
  
  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
  
  #Adjust to censor at age 80
  pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
  
  #Perform survival analysis
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  
  if(phenocols[i] != "C3_BREAST" & phenocols[i] != "C3_PROSTATE"){
    
    males <- subset(pheno, SEX=="male")
    females <- subset(pheno, SEX=="female")
    
    malesurvival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=males, na.action=na.exclude)
    femalesurvival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=females, na.action=na.exclude)
    
    #Extract hazard ratios, betas, standard errors and p-vals
    phenotype <- rep(phenocols[i],3)
    prs <- rep(prscols[i],3)
    test <- c("Full Sample", "Male Sample", "Female Sample")
    controls <- c(table(pheno[[phenocols[i]]])[1], table(males[[phenocols[i]]])[1], table(females[[phenocols[i]]])[1])  
    cases <- c(table(pheno[[phenocols[i]]])[2], table(males[[phenocols[i]]])[2], table(females[[phenocols[i]]])[2])  
    betas <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"coef"],summary(malesurvival)$coefficients[paste0(prscols[i],"_prs"),"coef"],summary(femalesurvival)$coefficients[paste0(prscols[i],"_prs"),"coef"])
    std_errs <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"se(coef)"],summary(malesurvival)$coefficients[paste0(prscols[i],"_prs"),"se(coef)"],summary(femalesurvival)$coefficients[paste0(prscols[i],"_prs"),"se(coef)"])
    pvals <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"],summary(malesurvival)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"],summary(femalesurvival)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"])
    OR <- exp(betas)
    CIpos <- exp(betas+1.96*std_errs)
    CIneg <- exp(betas-1.96*std_errs)
    result <- matrix(c(phenotype, prs, test, controls, cases, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=3, ncol=11)
    results <- rbind(results, result)
  } else {
    phenotype <- phenocols[i]
    controls <- table(pheno[[phenocols[i]]])[1]
    cases <- table(pheno[[phenocols[i]]])[2]
    prs <- prscols[i]
    test <- "Full Sample"
    betas <- summary(survival)$coefficients[paste0(prscols[i],"_prs"),"coef"]
    std_errs <- summary(survival)$coefficients[paste0(prscols[i],"_prs"),"se(coef)"]
    pvals <- summary(survival)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
    OR <- exp(betas)
    CIpos <- exp(betas+1.96*std_errs)
    CIneg <- exp(betas-1.96*std_errs)
    result <- c(phenotype, prs, test, controls, cases, betas, std_errs, pvals, OR, CIpos, CIneg)
    results <- rbind(results, result)
  }
  
}

write.csv(results, "file/path/to/output/HRperSD_[ENTER_BIOBANK_NAME].csv")
