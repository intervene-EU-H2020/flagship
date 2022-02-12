#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)

#Read in phenotype file
pheno <- fread(input="path/to/pheno_file", data.table=FALSE)

phenotypes <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "BMI", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Educational_Attainment", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lifespan", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "Pain", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "Stroke", "Subarachnoid_Haemmorhage", "T2D", "Thyroid_Stimulating_Hormone")

for(i in phenotypes){
#Read in PRS scores
PRS <- fread(input=paste0("path/to/PRS/",i,"_PRS.sscore"), data.table=FALSE)

#Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
PRS <- PRS[,c(1,2,5)]

#Rename ID column to the name of the ID column in the 
colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(i,"_prs"))

#left_join to the phenotype file
pheno <- left_join(pheno, PRS)
}

#Read in file which contains the principal components computed on the subset of individuals who are of european ancestry
pcs <- fread("path/to/european/principal/components", data.table=FALSE)

#Subset to first 10 PCs - change subset if this doesnt work!!
pcs <- pcs[c(1:11)]

pheno <- left_join(pheno, pcs)

#Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
#As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
#Feel free to subset using your own code: only provided as a reminder.
pheno <- subset(pheno, ANCESTRY=='EUR')

#Standardise PRS now they are subset to european ancestry participants.
for(i in phenotypes){
  pheno[[paste0(i,"_prs")]] <- scale(pheno[[paste0(i,"_prs")]])
}

#Perform regressions
phenocols <- c("AUD_SWEDISH", "G6_AD_WIDE", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "N14_CHRONKIDNEYDIS", "G6_EPLEPSY", "FE_STRICT", "GE_STRICT", "GOUT", "I9_HEARTFAIL_NS", "COX_ARTHROSIS", "ILD", "K11_IBD_STRICT", "KNEE_ARTHROSIS", "C3_BRONCHUS_LUNG", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "M13_OSTEOPOROSIS", "H7_GLAUCOMA", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "G6_SLEEPAPNO", "I9_STR", "I9_SAH", "T2D", "E4_HYTHYNAS", "E4_THYTOXGOITDIF")
prscols <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "Stroke", "Subarachnoid_Haemmorhage", "T2D", "Thyroid_Stimulating_Hormone", "Thyroid_Stimulating_Hormone") #Lifepsan, Pain and maybe BMI and educational attainment to be assessed by all phenotypes so do separately

#Loop through each phenotype and corresponding prs to perform regression,
results <- c()
for(i in 1:29){
  print(phenocols[i])
  print(prscols[i])
  regression <- glm(as.formula(paste(phenocols[i], " ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)
  phenotype <- phenocols[i]
  prs <- prscols[i]
  betas <- summary(regression)$coefficients[prscols[paste0(i,"_prs")],"Estimate"]
  std_errs <- summary(regression)$coefficients[prscols[paste0(i,"_prs")],"Std. Error"]
  pvals <- summary(regression)$coefficients[prscols[paste0(i,"_prs")],"Pr(>|z|)"]
  
  # Bradley, please see, i think correct is this
  #betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Estimate"]
  #std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Std. Error"]
  #pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
 
  OR <- exp(betas)
  CIpos <- exp(betas+(1.96*std_errs))
  CIneg <- exp(betas-(1.96*std_errs))
  result <- c(phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
  results <- rbind(results, result)
}

#For the lifespan, pain, educational attainment and BMI PRS, run regressions for all phenotypes
broadriskPRS <- c("BMI", "Educational_Attainment", "Lifespan", "Pain")
for(i in 1:29){
  for(j in broadriskPRS){
    print(phenocols[i])
    print(prscols[i])
    regression <- glm(as.formula(paste(phenocols[i], " ~ ", broadriskPRS, "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)
      #for Bradley, i think this instead
    #regression <- glm(as.formula(paste(phenocols[i], " ~ ", j, "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)

    phenotype <- phenocols[i]
    prs <- prscols[i]
    betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs")],"Estimate"]
    std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs")],"Std. Error"]
    pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs")],"Pr(>|z|)"]
    #for Bradley, i propose this
    #betas <- summary(regression)$coefficients[paste0(j,"_prs"),"Estimate"]
    #std_errs <- summary(regression)$coefficients[paste0(j,"_prs"),"Std. Error"]
    #pvals <- summary(regression)$coefficients[paste0(j,"_prs"),"Pr(>|z|)"]

    OR <- exp(betas)
    CIpos <- exp(betas+(1.96*std_errs))
    CIneg <- exp(betas-(1.96*std_errs))
    result <- c(phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
    results <- rbind(results, result)
  }
}

#Aggregate results
write.csv(results, "prs_associations_ENTER_BIOBANK_NAME.csv")

#Sex stratified analysis 

#Perform regressions
phenocols <- c("AUD_SWEDISH", "G6_AD_WIDE", "J10_ASTHMA", "I9_AF", "I9_CHD", "N14_CHRONKIDNEYDIS", "G6_EPLEPSY", "FE_STRICT", "GE_STRICT", "GOUT", "I9_HEARTFAIL_NS", "COX_ARTHROSIS", "ILD", "K11_IBD_STRICT", "KNEE_ARTHROSIS", "C3_BRONCHUS_LUNG", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "M13_OSTEOPOROSIS", "H7_GLAUCOMA", "RHEUMA_SEROPOS_OTH", "G6_SLEEPAPNO", "I9_STR", "I9_SAH", "T2D", "E4_HYTHYNAS", "E4_THYTOXGOITDIF")
prscols <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "CHD", "Chronic_Kidney_Disease", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "POAG", "Rheumatoid_Arthritis", "Sleep_Apnoea", "Stroke", "Subarachnoid_Haemmorhage", "T2D", "Thyroid_Stimulating_Hormone", "Thyroid_Stimulating_Hormone") #Lifepsan, Pain and maybe BMI and educational attainment to be assessed by all phenotypes so do separately

results <- c()
#Test for sex*prs interaction 
for(i in 1:length(phenocols)){
  print(phenocols[i])
  print(prscols[i])
  regression <- glm(as.formula(paste(phenocols[i], " ~ ", prscols[i], "_prs + SEX + ", prscols[i], "_prs:SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)
  phenotype <- phenocols[i]
  prs <- prscols[i]
  betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"Estimate"]
  std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"Std. Error"]
  pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs:SEXmale"),"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+(1.96*std_errs))
  CIneg <- exp(betas-(1.96*std_errs))
  result <- c('male', phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
  results <- rbind(results, result)
}

write.csv(results, "sex_interaction_ORs_ENTER_BIOBANK_NAME.csv")

#Male only regressions
male_results <- c()
males <- subset(pheno, SEX=='male')
for(i in 1:length(phenocols)){
  print(phenocols[i])
  print(prscols[i])
  regression <- glm(as.formula(paste(phenocols[i], " ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=males, na.action=na.exclude)
  phenotype <- phenocols[i]
  prs <- prscols[i]
  betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Estimate"]
  std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Std. Error"]
  pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+(1.96*std_errs))
  CIneg <- exp(betas-(1.96*std_errs))
  male_result <- c('male', phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
  male_results <- rbind(male_results, male_result)
}

#Female only regressions
female_results <- c()
females <- subset(pheno, SEX=='female')
for(i in 1:length(phenocols)){
  print(phenocols[i])
  print(prscols[i])
  regression <- glm(as.formula(paste(phenocols[i], " ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=females, na.action=na.exclude)
  phenotype <- phenocols[i]
  prs <- prscols[i]
  betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Estimate"]
  std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Std. Error"]
  pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
  OR <- exp(betas)
  CIpos <- exp(betas+(1.96*std_errs))
  CIneg <- exp(betas-(1.96*std_errs))
  female_result <- c('female', phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
  female_results <- rbind(female_results, female_result)
}

results <- rbind(male_results, female_results)
write.csv(results, "sex_specific_ORs_ENTER_BIOBANK_NAME.csv")

#Age stratified results

#Calculate age at diagnosis and age at end of follow up
phenocols <- c("AUD_SWEDISH", "G6_AD_WIDE", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "N14_CHRONKIDNEYDIS", "G6_EPLEPSY", "FE_STRICT", "GE_STRICT", "GOUT", "I9_HEARTFAIL_NS", "COX_ARTHROSIS", "ILD", "K11_IBD_STRICT", "KNEE_ARTHROSIS", "C3_BRONCHUS_LUNG", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "M13_OSTEOPOROSIS", "H7_GLAUCOMA", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "G6_SLEEPAPNO", "I9_STR", "I9_SAH", "T2D", "E4_HYTHYNAS", "E4_THYTOXGOITDIF")
prscols <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "Stroke", "Subarachnoid_Haemmorhage", "T2D", "Thyroid_Stimulating_Hormone", "Thyroid_Stimulating_Hormone") #Lifepsan, Pain and maybe BMI and educational attainment to be assessed by all phenotypes so do separately

#For diagnosis
for(i in phenocols){
  pheno[,paste0(i,"_AGE")] <- time_length(difftime(pheno[,paste0(i,"_DATE")], pheno[,"DATE_OF_BIRTH"]), 'years')
  #For Bradley - for me it was necessary to define it as Date and define format structure
  #pheno[,paste0(i,"_AGE")] <-time_length(difftime(as.Date(pheno[,paste0(i,"_DATE")],format="%Y-%m-%d"), as.Date(pheno[,"DATE_OF_BIRTH"],format="%Y-%m-%d")), 'years')

} 

#For end of follow up
pheno$AGE_FU_END <- time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years')
#For bradley, again as.date conversion
#pheno$AGE_FU_END <- time_length(difftime(as.Date(pheno$END_OF_FOLLOWUP,format="%Y-%m-%d"), as.Date(pheno$DATE_OF_BIRTH,format="%Y-%m-%d")), 'years')


results <- c()

for(i in 1:length(phenocols)){
  #Age bins 0-10. 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  for(j in list(c(0,10),c(10,20),c(20,30),c(30,40),c(40,50),c(50,60),c(60,70),c(70,80))){
    print(phenocols[i])
    print(prscols[i])
    
    #Subset to those who have follow-up from the baseline age and do not have the disease prior.
    sample <- subset(pheno, AGE_FU_END >= j[2] & (is.na(pheno[[paste0(phenocols[i],"_AGE")]]) | pheno[[paste0(phenocols[i],"_AGE")]] >= j[1]))
    sample <- sample[,c("FINNGENID", paste0(prscols[i],"_prs"), paste0(phenocols[i],"_AGE"), "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "BATCH", "COHORT", "AGE_FU_END")]
    
    #Assign case control status  
    sample$disease_status <- case_when(is.na(sample[,paste0(phenocols[i],"_AGE")]) | sample[,paste0(phenocols[i],"_AGE")] > j[2] ~ 0,
                                       sample[,paste0(phenocols[i],"_AGE")] <= j[2] ~ 1,
                                       TRUE ~ NA_real_)
    
    #Perform regression
    regression <- glm(as.formula(paste("disease_status ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=females, na.action=na.exclude)
    #For Bradley:
    # regression <- glm(as.formula(paste("disease_status ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BATCH + COHORT", sep="")), family=binomial(link='logit'), data=sample, na.action=na.exclude)

    phenotype <- phenocols[i]
    prs <- prscols[i]
    betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Estimate"]
    std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Std. Error"]
    pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
    OR <- exp(betas)
    CIpos <- exp(betas+(1.96*std_errs))
    CIneg <- exp(betas-(1.96*std_errs))
    result <- c(prs, phenotype, j[1], j[2], table(sample$disease_status)[1], table(sample$disease_status)[2], betas, std_errs, pvals, OR, CIpos, CIneg)
    results <- rbind(results, result)
  }
}

colnames(results) <- c("prs", "phenotype", "min_age", "max_age", "controls", "cases", "betas", "std_errs", "pvals", "OR", "CIpos", "CIneg")

write.csv(results, "age_stratified_ORs_ENTER_BIOBANK_NAME.csv")

