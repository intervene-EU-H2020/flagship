#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

percentiles <- list(c(0,0.01,0.05,0.1,0.2,0.4), #1%
                    c(0,0.05,0.1,0.2,0.4), #5%
                    c(0,0.1,0.2,0.4), #10%
                    c(0,0.2,0.4) #20%
)

maleresults <- c()
femaleresults <- c()

for(i in 1:length(phenocols)){
  for(p in percentiles){
    
    print(phenocols[i])
    print(prscols[i])
    
    #Read in phenotype file
    pheno <- fread(input="path/to/pheno_file", select=c("ID","DATE_OF_BIRTH","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP"), data.table=FALSE)
    
    pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
    
    #Read in PRS scores
    PRS <- fread(input=paste0("path/to/PRS/",prscols[i],"_PRS.sscore"), data.table=FALSE)
    
    #Subset columns to the IDs and score only. Note: columns 1 or 2 may be redundant and can be removed if necessary. Kept in to avoid bugs.
    PRS <- PRS[,c(1,2,5)]
    
    #Rename ID column to the name of the ID column in the 
    colnames(PRS) <- c("ENTER_ID", "ENTER_ID", paste0(prscols[i],"_prs"))
    
    #left_join to the phenotype file
    pheno <- left_join(pheno, PRS)
    
    pheno <- subset(pheno, !is.na(pheno[[paste0(phenocols[i])]]) | !is.na(pheno[[paste0(prscols[i],"_prs")]]))
    
    #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
    #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
    #Feel free to subset using your own code: only provided as a reminder.
    pheno <- subset(pheno, ANCESTRY=='EUR')
    
    #Assign PRS into percentiles
    q <- quantile(pheno[[paste0(prscols[i],"_prs")]], probs=c(p,rev(1-p)))
    
    pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[[paste0(prscols[i],"_prs")]], q, include.lowest=TRUE,
                                                labels=paste("Group",1:(2*length(p)-1)))
    
    #Make all necessary variables factors
    pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
    pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref=paste("Group",length(p)))
    
    #Specify age as either the Age at Onset or End of Follow-up (if not a case)
    pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
    
    #Adjust to censor at age 80
    pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
    pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
    
    males <- subset(pheno, SEX=="male")
    controls <- table(males[[paste0(prscols[i],"_group")]], males[[paste0(phenocols[i])]])[2:(2*length(p)-1),1]
    cases <- if(sum(nrow(males[males[[paste0(phenocols[i])]]==0,]))==length(males[[paste0(phenocols[i])]])){
      rep(0,(2*length(p)-2))} else {table(males[[paste0(prscols[i],"_group")]], males[[paste0(phenocols[i])]])[2:(2*length(p)-1),2]}
    
    #Perform survival analysis
    survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ARRAY + ASSESSMENT_CENTRE")), data=males, na.action=na.exclude)
    
    #Extract hazard ratios, betas, standard errors and p-vals - in the first instance extract all results, for the latter just take the 
    if(p[2] == 0.01){
      phenotype <- rep(phenocols[i],(2*length(p)-2))
      prs <- rep(prscols[i],(2*length(p)-2))
      group <- c(paste0(prscols[i],"_groupGroup ",c(1:(length(p)-1),(length(p)+1):(2*length(p)-1))))
      betas <- summary(survival)$coefficients[group,"coef"]
      std_errs <- summary(survival)$coefficients[group,"se(coef)"]
      pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
      groups <- c("< 1%","1-5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","95-99%", "> 99%")
      OR <- exp(betas)
      CIpos <- exp(betas+1.96*std_errs)
      CIneg <- exp(betas-1.96*std_errs)
      maleresult <- matrix(c(phenotype, prs, groups, controls, cases, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=11)
      maleresults <- rbind(maleresults, maleresult)
    } else {
      phenotype <- rep(phenocols[i],2)
      prs <- rep(prscols[i],2)
      group <- c(paste0(prscols[i],"_groupGroup 1"), paste0(prscols[i],"_groupGroup ", (2*length(p)-1)))
      betas <- summary(survival)$coefficients[group,"coef"]
      std_errs <- summary(survival)$coefficients[group,"se(coef)"]
      pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
      groups <- c(paste0("< ",(p[2]*100),"%"),paste0("> ",((1-p[2])*100),"%"))
      OR <- exp(betas)
      CIpos <- exp(betas+1.96*std_errs)
      CIneg <- exp(betas-1.96*std_errs)
      maleresult <- matrix(c(phenotype, prs, groups, controls[c(1,length(controls))], cases[c(1,length(cases))], betas, std_errs, pvals, OR, CIpos, CIneg), nrow=2, ncol=11)
      maleresults <- rbind(maleresults, maleresult)
    }
    
    females <- subset(pheno, SEX=="female")
    controls <- table(females[[paste0(prscols[i],"_group")]], females[[paste0(phenocols[i])]])[2:(2*length(p)-1),1]
    cases <- if(sum(nrow(females[females[[paste0(phenocols[i])]]==0,]))==length(females[[paste0(phenocols[i])]])){
      rep(0,(2*length(p)-2))} else {table(females[[paste0(prscols[i],"_group")]], females[[paste0(phenocols[i])]])[2:(2*length(p)-1),2]}
    
    #Perform survival analysis
    survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ARRAY + ASSESSMENT_CENTRE")), data=females, na.action=na.exclude)
    
    #Extract hazard ratios, betas, standard errors and p-vals - in the first instance extract all results, for the latter just take the 
    if(p[2] == 0.01){
      phenotype <- rep(phenocols[i],(2*length(p)-2))
      prs <- rep(prscols[i],(2*length(p)-2))
      group <- c(paste0(prscols[i],"_groupGroup ",c(1:(length(p)-1),(length(p)+1):(2*length(p)-1))))
      betas <- summary(survival)$coefficients[group,"coef"]
      std_errs <- summary(survival)$coefficients[group,"se(coef)"]
      pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
      groups <- c("< 1%","1-5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","95-99%", "> 99%")
      OR <- exp(betas)
      CIpos <- exp(betas+1.96*std_errs)
      CIneg <- exp(betas-1.96*std_errs)
      femaleresult <- matrix(c(phenotype, prs, groups, controls, cases, betas, std_errs, pvals, OR, CIpos, CIneg), nrow=10, ncol=11)
      femaleresults <- rbind(femaleresults, femaleresult)
    } else {
      phenotype <- rep(phenocols[i],2)
      prs <- rep(prscols[i],2)
      group <- c(paste0(prscols[i],"_groupGroup 1"), paste0(prscols[i],"_groupGroup ", (2*length(p)-1)))
      betas <- summary(survival)$coefficients[group,"coef"]
      std_errs <- summary(survival)$coefficients[group,"se(coef)"]
      pvals <- summary(survival)$coefficients[group,"Pr(>|z|)"]
      groups <- c(paste0("< ",(p[2]*100),"%"),paste0("> ",((1-p[2])*100),"%"))
      OR <- exp(betas)
      CIpos <- exp(betas+1.96*std_errs)
      CIneg <- exp(betas-1.96*std_errs)
      femaleresult <- matrix(c(phenotype, prs, groups, controls[c(1,length(controls))], cases[c(1,length(cases))], betas, std_errs, pvals, OR, CIpos, CIneg), nrow=2, ncol=11)
      femaleresults <- rbind(femaleresults, femaleresult)
    }
  }
}

write.csv(maleresults, "file/path/to/output/HR_MaleSample_[ENTER_BIOBANK_NAME].csv")
write.csv(femaleresults, "file/path/to/output/HR_FemaleSample_[ENTER_BIOBANK_NAME].csv")


