
library(optparse)
library(data.table)
library(dplyr)
library(lubridate)
library(survival)

print(R.version)

option_list <- list(
  make_option("--phenofile", type="character", default="",help="Provide the tab delimited file with set column names that has phenotypes and dates of diagnosis"), 
  make_option("--prs_path", type="character", default="",help="Provide the path to the .sscore files. Assums files end in _PRS.sscore and has header SCORE1_AVG"),
  make_option("--output_dir",type="character",default=".",help="path to where the output file 'cindex_output.tab' will be written."),
  make_option("--biobank",type="character",default="",help="String to describe your biobank")
)

### example Cindex.R --phenofile </path/to/file.txt> --prs </path/to/scores/> --output_dr /path/to/out --biobank "string"

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script calculates c-index for several models")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

out<-data.frame(model=NA,biobank=NA,cindex=NA,cindex_SE=NA,pheno=NA)

p<-c(0,0.01,0.05,0.1,0.2,0.4)

#assumes fixed column names
pheno <- fread(input=opt$phenofile, data.table=FALSE)

for(i in 1:length(phenocols)){        
  
  print(phenocols[i])
  
  #you may need to customize this line depending on the format of your date variable
  pheno<-pheno %>% mutate_at(paste(sep="_",phenocols[i],"DATE"),as.Date,origin="1970-01-01")
  
  #Read in PRS scores 
  PRS <- fread(input=paste0(opt$prs_path,prscols[i],"_PRS.sscore"), data.table=FALSE)

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


####breast and prostate cancer run separately because only in males/females
phenocols<-c("C3_BREAST","C3_PROSTATE")
prscols<-c("Breast_Cancer","Prostate_Cancer")
#sex_string<-c("FEMALE","MALE") #assumes this is the coding for SEX column
sex_string<-c(2,1) #comment this out and choose FEMALE, MALE if you use that convention

for(l in 1:length(phenocols)){     
  pheno <- fread(input=opt$phenofile, data.table=FALSE)
  #you may need to customize this line depending on the format of your date variable, this line is also filtering by sex
  pheno<-pheno %>% filter(SEX==!!sex_string[l]) %>% mutate_at(paste(sep="_",phenocols[l],"DATE"),as.Date,origin="1970-01-01")

  #Read in PRS scores 
  PRS <- fread(input=paste0(opt$prs_path,prscols[l],"_PRS.sscore"), data.table=FALSE)
  PRS <- fread(input=paste0(  "/home/bwolford/scratch/brooke/scores/",prscols[l],"_PRS.sscore"), data.table=FALSE)
  

  #some people may have SCORE1_SUM
  #Subset columns to the IDs and score only. Note: columns FID or IID may be redundant and can be removed if necessary. Kept in to avoid bugs.
  PRS <- PRS[,c("#FID","IID","SCORE1_AVG")]

  #Rename ID column to the name of the ID column in the phenotype file
  colnames(PRS) <- c("FID", "ID", "PRS")

  #left_join to the phenotype file
  pheno <- left_join(pheno, PRS) #joins by "ID" column

  pheno <- subset(pheno, !(is.na(pheno[[paste0(phenocols[l])]]) | is.na(pheno[["PRS"]])))

  #Assign PRS into percentiles
  q <- quantile(pheno[["PRS"]], probs=c(p,rev(1-p)))

  pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[["PRS"]], q, include.lowest=TRUE,
                                            labels=paste("Group",1:(2*length(p)-1)))
  
  
  
  #Make all necessary variables factors
  pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
  pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref=paste("Group",length(p)))
  
  #Specify age as either the Age at Onset or End of Follow-up (if not a case)
  pheno$AGE <- ifelse(pheno[[phenocols[l]]]==1, time_length(difftime(pheno[[paste0(phenocols[l],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
  
  #Adjust to censor at age 80
  pheno[[paste0(phenocols[l])]] <- ifelse(pheno[[paste0(phenocols[l])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[l])]])
  pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)
  
  #set birth year (also af function of how your date of birth is coded)
  pheno$BY<-year(pheno$DATE_OF_BIRTH)
  
  
  ##null model
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[l],") ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="null",pheno=phenocols[l])
  
  
  ###group and raw PGS 
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[l],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group",pheno=phenocols[l])
  
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[l],") ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS",pheno=phenocols[l])
  
  ## stick with group and add birth year as surrogate for age, first drop the PRS grouping
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[l],") ~  + BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="birthyear",pheno=phenocols[l])
  
  
  #then add grouping back
  survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[l],") ~ ",prscols[i],"_group + BY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
  out <- out %>% add_row(cindex=concordance(survival)$concordance, cindex_SE=sqrt(concordance(survival)$var), model="PRS_percentile_group_and_birthyear",pheno=phenocols[l])
  

}

#write the output 
out$biobank<-out$biobank
write.table(out,paste0(opt$output_dir,"cindex_output.tab"),quote=FALSE,sep="\t",row.names=FALSE)
