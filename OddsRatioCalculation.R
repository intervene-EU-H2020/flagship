#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(fmsb)

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

#Read in batch and assessment centre for inclusion as covariates (if you have them).
covariates <- fread("path/to/covariate/file", data.table=FALSE)

#Change column names to batch and assessment centre. - Note you will need to identify the appropriate index.
colnames(covariates)[c('ENTER INDICES')] <- c("batch", "assessment_centre")

covariates <- full_join(pcs, covariates)

pheno <- left_join(pheno, covariates)

#Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
#As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
#Feel free to subset using your own code: only provided as a reminder.
pheno <- subset(pheno, ANCESTRY=='EUR')

#Subset to unrelated individuals - this will be specific to the biobank as was the case for ancestries. 

#Standardise PRS now they are subset to european ancestry participants.
for(i in phenotypes){
  pheno[[paste0(i,"_prs")]] <- scale(pheno[[paste0(i,"_prs"]])
}

#Perform regressions
phenocols <- c("AUD_SWEDISH", "G6_AD_WIDE", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "N14_CHRONKIDNEYDIS", "G6_EPLEPSY", "FE_STRICT", "GE_STRICT", "GOUT", "I9_HEARTFAIL_NS", "COX_ARTHROSIS", "ILD", "K11_IBD_STRICT", "KNEE_ARTHROSIS", "C3_BRONCHUS_LUNG", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "M13_OSTEOPOROSIS", "H7_GLAUCOMA", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "G6_SLEEPAPNO", "I9_STR", "I9_SAH", "T2D", "E4_HYTHYNAS", "E4_THYTOXGOITDIF")
prscols <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "Stroke", "Subarachnoid_Haemmorhage", "T2D", "Thyroid_Stimulating_Hormone", "Thyroid_Stimulating_Hormone") #Lifepsan, Pain and maybe BMI and educational attainment to be assessed by all phenotypes so do separately

#Loop through each phenotype and corresponding prs to perform regression,
results <- c()
for(i in 1:29){
  print(phenocols[i])
  print(prscols[i])
  regression <- glm(as.formula(paste(phenocols[i], " ~ ", prscols[i], "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + batch + assessment_centre", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)
  phenotype <- phenocols[i]
  prs <- prscols[i]
  betas <- summary(regression)$coefficients[paste(prscols[i],"_prs"),"Estimate"]
  std_errs <- summary(regression)$coefficients[paste(prscols[i],"_prs"),"Std. Error"]
  pvals <- summary(regression)$coefficients[paste(prscols[i],"_prs"),"Pr(>|z|)"]
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
    regression <- glm(as.formula(paste(phenocols[i], " ~ ", broadriskPRS, "_prs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + batch + assessment_centre", sep="")), family=binomial(link='logit'), data=pheno, na.action=na.exclude)
    phenotype <- phenocols[i]
    prs <- prscols[i]
    betas <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Estimate"]
    std_errs <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Std. Error"]
    pvals <- summary(regression)$coefficients[paste0(prscols[i],"_prs"),"Pr(>|z|)"]
    OR <- exp(betas)
    CIpos <- exp(betas+(1.96*std_errs))
    CIneg <- exp(betas-(1.96*std_errs))
    result <- c(phenotype, prs, betas, std_errs, pvals, OR, CIpos, CIneg)
    results <- rbind(results, result)
  }
}

#Aggregate results
write.csv(results, "prs_associations_ENTER_BIOBANK_NAME.csv")
