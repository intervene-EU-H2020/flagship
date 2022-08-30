#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

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
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_prs + SEX + ", prscols[i], "_prs:SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  print(summary(survival))
  
  #Extract hazard ratios, betas, standard errors and p-vals
  phenotype <- rep(phenocols[i], 3)
  prs <- rep(prscols[i], 3)
  test <- c("PRS_Main","Sex_Main","PRS*Sex")
  betas <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"coef"],summary(survival)$coefficients["SEXmale","coef"],summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"coef"])
  std_errs <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"se(coef)"],summary(survival)$coefficients["SEXmale","se(coef)"],summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"se(coef)"])
  pvals <- c(summary(survival)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"],summary(survival)$coefficients["SEXmale","Pr(>|z|)"],summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"Pr(>|z|)"])
  OR <- exp(betas)
  CIpos <- exp(betas+1.96*std_errs)
  CIneg <- exp(betas-1.96*std_errs)
  result <- matrix(c(phenotype, prs, test, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=3, ncol=9)
  results <- rbind(results, result)
  
  
}

write.csv(results, "file/path/to/output/HR_SexInteraction_[ENTER_BIOBANK_NAME].csv")
