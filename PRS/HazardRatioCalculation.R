#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTHER", "I9_SAH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D","T2D", "ILD", "Lung_Cancer)

results <- c()

for(i in 1:length(phenocols)){
  
  print(phenocols[i])
  print(prscols[i])
  
  #Read in phenotype file
  pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP","BATCH","COHORT"), data.table=FALSE)
    
  pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
 
  #Read in PRS scores
  PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)

  #Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c(1,2,5)]

  #Rename ID column to the name of the ID column in the 
  colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))

  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS)

  pheno <- subset(pheno, !is.na(pheno[[paste0(prscols[i],"_prs")]]))
 
  #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
  #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
  #Feel free to subset using your own code: only provided as a reminder.
  pheno <- subset(pheno, ANCESTRY=='EUR')
  
  #Assign PRS into percentiles
  q <- quantile(pheno[[paste0(prscols[i],"_prs")]], probs=c(0,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99,1))

  pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[[paste0(prscols[i],"_prs")]], q, include.lowest=TRUE,
                                              labels=paste("Group",1:11))
  
  #Make all necessary variables factors
  pheno$BATCH <- as.factor(pheno$BATCH)
  pheno$COHORT <- as.factor(pheno$COHORT)
  pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
  pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref="Group 6")
  
  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))

  #Adjust to censor at age 80
  pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
  
  #Perform survival analysis
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT")), data=pheno, na.action=na.exclude)

  #Extract hazard ratios, betas, standard errors and p-vals
  phenotype <- rep(phenocols[i],10)
  prs <- rep(prscols[i],10)
  group <- c(paste0(prscols[i],"_groupGroup ",c(1:5:7:11)))
  betas <- summary(survival)$coefficients[group,"coef"]
  std_errs <- summary(survival)$coefficients[group,"se(coef)"]
  pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+1.96*std_errs)
  CIneg <- exp(betas-1.96*std_errs)
  result <- matrix(c(phenotype, prs, group, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=9)
  results <- rbind(results, result)
  
}

write.csv(results, "file/path/to/output_FullSample.csv")

###########################################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################################

#Sex specific HRs

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "RHEUMA_SEROPOS_OTHER", "I9_SAH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D", "T2D", "ILD", "Lung_Cancer")

maleresults <- c()
femaleresults <- c()

for(i in 1:length(phenocols)){
  
  print(phenocols[i])
  print(prscols[i])
  
  #Read in phenotype file
  pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP","BATCH","COHORT"), data.table=FALSE)
  
  pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
 
  #Read in PRS scores
  PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)
  
  #Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c(1,2,5)]
  
  #Rename ID column to the name of the ID column in the 
  colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))
  
  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS)
  
  pheno <- subset(pheno, !is.na(pheno[[paste0(prscols[i],"_prs")]]))
 
  #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
  #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
  #Feel free to subset using your own code: only provided as a reminder.
  pheno <- subset(pheno, ANCESTRY=='EUR')
  
  #Assign PRS into percentiles
  q <- quantile(pheno[[paste0(prscols[i],"_prs")]], probs=c(0,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99,1))
  
  pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[[paste0(prscols[i],"_prs")]], q, include.lowest=TRUE,
                                              labels=paste("Group",1:11))
  
  #Make all necessary variables factors
  pheno$BATCH <- as.factor(pheno$BATCH)
  pheno$COHORT <- as.factor(pheno$COHORT)
  pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
  pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref="Group 6")
  
  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
  
  #Adjust to censor at age 80
  pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
  
  males <- subset(pheno, SEX="male")
  
  #Perform survival analysis
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT")), data=males, na.action=na.exclude)
  
  #Extract hazard ratios, betas, standard errors and p-vals
  phenotype <- rep(phenocols[i],10)
  prs <- rep(prscols[i],10)
  group <- c(paste0(prscols[i],"_groupGroup ",c(1:5,7:11)))
  betas <- summary(survival)$coefficients[group,"coef"]
  std_errs <- summary(survival)$coefficients[group,"se(coef)"]
  pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+1.96*std_errs)
  CIneg <- exp(betas-1.96*std_errs)
  maleresult <- matrix(c(phenotype, prs, group, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=9)
  maleresults <- rbind(maleresults, maleresult)
  
  females <- subset(pheno, SEX="female")
  
  #Perform survival analysis
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT")), data=females, na.action=na.exclude)
  
  #Extract hazard ratios, betas, standard errors and p-vals
  phenotype <- rep(phenocols[i],10)
  prs <- rep(prscols[i],10)
  group <- c(paste0(prscols[i],"_groupGroup ",c(1:5,7:11)))
  betas <- summary(survival)$coefficients[group,"coef"]
  std_errs <- summary(survival)$coefficients[group,"se(coef)"]
  pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+1.96*std_errs)
  CIneg <- exp(betas-1.96*std_errs)
  femaleresult <- matrix(c(phenotype, prs, group, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=9)
  femaleresults <- rbind(femaleresults, femaleresult)
  
}

write.csv(maleresults, "file/path/to/output_MaleSample.csv")
write.csv(femaleresults, "file/path/to/output_FemaleSample.csv")

###########################################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################################

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "RHEUMA_SEROPOS_OTHER", "I9_SAH", "T1D", "T2D")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D","T2D")

results <- c()

for(i in 1:length(phenocols)){
  
  print(phenocols[i])
  print(prscols[i])
  
  #Read in phenotype file
  pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP","BATCH","COHORT"), data.table=FALSE)
  
  pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
 
  #Read in PRS scores
  PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)
  
  #Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c(1,2,5)]
  
  #Rename ID column to the name of the ID column in the 
  colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))
  
  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS)
  
  pheno <- subset(pheno, !is.na(pheno[[paste0(prscols[i],"_prs")]]))
 
  #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
  #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
  #Feel free to subset using your own code: only provided as a reminder.
  pheno <- subset(pheno, ANCESTRY=='EUR')
  
  pheno[[paste0(prscols[i],"_prs")]] <- scale(pheno[[paste0(prscols[i],"_prs")]])

  #Make all necessary variables factors
  pheno$BATCH <- as.factor(pheno$BATCH)
  pheno$COHORT <- as.factor(pheno$COHORT)
 
  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
  
  #Adjust to censor at age 80
  pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
  
  #Perform survival analysis
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_prs + SEX + ", prscols[i], "_prs:SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT")), data=pheno, na.action=na.exclude)
  
  #Extract hazard ratios, betas, standard errors and p-vals
  phenotype <- rep(phenocols[i],10)
  prs <- rep(prscols[i],10)
  betas <- summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"coef"]
  std_errs <- summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"se(coef)"]
  pvals <- summary(survival)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+1.96*std_errs)
  CIneg <- exp(betas-1.96*std_errs)
  result <- c(phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
  results <- rbind(results, result)
  
}

write.csv(results, "file/path/to/output_SexInteraction.csv")

###########################################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################################

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTHER", "I9_SAH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D", "T2D", "ILD", "Lung_Cancer")

#TO BE COMPLETED - will be based on mean quartiles from biobanks to be used in the lifetime risk estimation
ages <- list(c(c(),c(),c(),c()), #C3_CANCER
             c(c(),c(),c(),c()), #K11_APPENDACUT
             c(c(),c(),c(),c()), #J10_ASTHMA
             c(c(),c(),c(),c()), #I9_AF
             c(c(),c(),c(),c()), #C3_BREAST
             c(c(),c(),c(),c()), #I9_CHD
             c(c(),c(),c(),c()), #C3_COLORECTAL
             c(c(),c(),c(),c()), #G6_EPLEPSY
             c(c(),c(),c(),c()), #GOUT
             c(c(),c(),c(),c()), #COX_ARTHROSIS
             c(c(),c(),c(),c()), #KNEE_ARTHROSIS
             c(c(),c(),c(),c()), #F5_DEPRESSIO
             c(c(),c(),c(),c()), #C3_MELANOMA_SKIN
             c(c(),c(),c(),c()), #C3_PROSTATE
             c(c(),c(),c(),c()), #RHEUMA_SEROPOS_OTHER
             c(c(),c(),c(),c()), #I9_SAH
             c(c(),c(),c(),c()), #T1D
             c(c(),c(),c(),c()), #T2D
             c(c(),c(),c(),c()), #ILD
             c(c(),c(),c(),c())) #C3_BRONCHUS_LUNG

results <- c()

for(i in 1:length(phenocols)){
  agestrat <- ages[i]
  for(j in c(1,3,5,7)){
  
    print(phenocols[i])
    print(prscols[i])
    print(agestrat[[1]][[j]])
    print(agestrat[[1]][[j+1]])

    #Read in phenotype file
    pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP","BATCH","COHORT"), data.table=FALSE)

    pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")

    #Read in PRS scores
    PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)

    #Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
    PRS <- PRS[,c(1,2,5)]

    #Rename ID column to the name of the ID column in the 
    colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))

    #left_join to the phenotype file
    pheno <- left_join(pheno, PRS)

    pheno <- subset(pheno, !is.na(pheno[[paste0(prscols[i],"_prs")]]))

    #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
    #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
    #Feel free to subset using your own code: only provided as a reminder.
    pheno <- subset(pheno, ANCESTRY=='EUR')

    #Assign PRS into percentiles
    q <- quantile(pheno[[paste0(prscols[i],"_prs")]], probs=c(0,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99,1))

    pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[[paste0(prscols[i],"_prs")]], q, include.lowest=TRUE,
                                                labels=paste("Group",1:11))

    #Make all necessary variables factors
    pheno$BATCH <- as.factor(pheno$BATCH)
    pheno$COHORT <- as.factor(pheno$COHORT)
    pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
    pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref="Group 6")

    #Specify age as either the Age at Onset or End of Follow-up (if not a case)
    pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))

    sample <- subset(pheno, is.na(AGE) | (pheno[[phenocols[i]]] == 1 & AGE > agestrat[[1]][[j]] & AGE <= agestrat[[1]][[j+1]]) | (pheno[[phenocols[i]]] == 0 & AGE > agestrat[[1]][[j]]))

    #Adjust to censor at upper bound of age group for controls
    sample$AGE <- ifelse(sample$AGE > agestrat[[1]][[j+1]], agestrat[[1]][[j+1]], sample$AGE)

    #Adjust to censor at age 80 as that is when we finish lifetime risk estimation
    sample$AGE <- ifelse(sample$AGE > 80, 80, sample$AGE)

    #Subtract the lower bound of age so that each age group only counts follow-up time once
    sample$AGE <- sample$AGE - agestrat[[1]][[j]]

    #Perform survival analysis
    survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT")), data=sample, na.action=na.exclude)

    controls <- table(sample[[paste0(prscols[i],"_group")]], sample[[phenocols[i]]])[2:11,1]
    cases <- if(sum(nrow(sample[sample[[paste0(phenocols[i])]]==0,])) == length(sample[[paste0(phenocols[i])]])){ 
      rep(0,10)} else {table(sample[[paste0(prscols[i],"_group")]], sample[[paste0(phenocols[i])]])[2:11,2]}

    #Extract hazard ratios, betas, standard errors and p-vals
    phenotype <- rep(phenocols[i],10)
    prs <- rep(prscols[i],10)
    minage <- rep(agestrat[[1]][[j]], 10)
    maxage <- rep(agestrat[[1]][[j+1]], 10)
    group <- c(paste0(prscols[i],"_groupGroup ",c(1:5:,7:11)))
    betas <- summary(survival)$coefficients[group,"coef"]
    std_errs <- summary(survival)$coefficients[group,"se(coef)"]
    pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
    OR <- exp(betas)
    CIpos <- exp(betas+1.96*std_errs)
    CIneg <- exp(betas-1.96*std_errs)
    result <- matrix(c(phenotype, prs, minage, maxage, controls, cases, group, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=13)
    results <- rbind(results, result)
  }
}

write.csv(results, "file/path/to/output_AgeStratifiedResults.csv")




