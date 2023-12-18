
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

print(R.version)

option_list <- list(
  make_option("--phenofile", type="character", default=""), 
  make_option("--prs", type="character", default=""),
  make_option("--output_dir",type="character",default=""),
  make_option("--biobank",type="character",default="")
)

### example Cindex.R --phenofile </path/to/file.txt> --prs </path/tpfile.txt> --output_dr /path/to/out --biobank "string"

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script calculates c-index for several models")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

out<-data.frame(model=NA,biobank=NA,cindex=NA,cindex_SE=NA,pheno=NA)

p<-c(0,0.01,0.05,0.1,0.2,0.4)

#assumes fixed column names
pheno <- fread(input=opt$phenofile, data.table=FALSE)

for(i in 1:length(phenocols)){        
  
  #you may need to customize this line depending on the format of your date variable
  pheno<-pheno %>% mutate_at(paste(sep="_",phenocols[i],"DATE"),as.Date,origin="1970-01-01")
  
  #Read in PRS scores 
  PRS <- fread(input=paste0(opt$prs), data.table=FALSE)

  #some people may have SCORE1_SUM
  #Subset columns to the IDs and score only. Note: columns FID or IID may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c("#FID","IID","SCORE1_AVG")]

  #Rename ID column to the name of the ID column in the phenotype file
  colnames(PRS) <- c("FID", "ID", "PRS")

  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS) #joins by "ID" column
  
  pheno <- subset(pheno, !(is.na(pheno[[paste0(phenocols[i])]]) | is.na(pheno[["PRS"]])))

  #Assign PRS into percentiles
  q <- quantile(pheno[["PRS"]], probs=c(p,rev(1-p)))

  pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[["PRS"]], q, include.lowest=TRUE,
                                            labels=paste("Group",1:(2*length(p)-1)))


  #Make all necessary variables factors
  pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
  pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref=paste("Group",length(p)))

  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))

  #Adjust to censor at age 80
  pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)

  #set birth year (also af function of how your date of birth is coded)
  pheno$BY<-year(pheno$DATE_OF_BIRTH)

  #Perform survival analysis


  ###null model
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="null",pheno=phenocols[i])

  
  ###group and raw PGS 
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group",pheno=phenocols[i])
  
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS",pheno=phenocols[i])
  
  ###group and raw PGS now add sex
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + SEX+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group_and_sex",pheno=phenocols[i])
  
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ PRS + SEX+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_and_sex",pheno=phenocols[i])
  
  
  ### stick with group and add birth year as surrogate for age, first drop the PRS grouping
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~  + BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="birthyear",pheno=phenocols[i])
  
  
  #then add grouping back
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group_and_birthyear",pheno=phenocols[i])

  
  #now birth year and sex only
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ SEX+ BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="birthyear_and_sex",pheno=phenocols[i])
  
  
  #now full model
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + SEX+ BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group_and_birthyear_and_sex",pheno=phenocols[i])
  
}
out$biobank<-out$biobank
write.table(out,paste0(opt$output_dir,"/cindex_output.tab"),quote=FALSE,sep="\t",row.names=FALSE)
