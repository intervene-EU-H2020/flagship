library(data.table)
library(dplyr)

#Mean is not associative so need to get back to the original sum to estimate. 

#Iterate over phenotype
pheno <- c("AllCancers", "Appendicitis", "Colorectal_Cancer", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Epilepsy", "Gout", "Hip_Osteoarthritis", "ILD", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "T1D", "T2D")

for(i in 1:length(pheno)){
  FinalScore <- fread(input=paste0("/path/to/PRSscore/files/",pheno[i],"_PRS_chr1.sscore"), data.table=FALSE)
  FinalScore <- as.data.frame(FinalScore$IID)
  colnames(FinalScore) <- "IID"
  print(dim(FinalScore))

  #Iterate over chromosomes
  for(j in 1:22){
    #Read in the dataset
    score <- fread(input=paste0("/path/to/PRSscore/files/",pheno[i],"_PRS_chr",j,".sscore"), data.table=FALSE)
  
    #Drop the FID and NMISS_ALLELE_CT columns as redundant
    score <- score[,-c(1,3)]
    score[[paste0("PRS_Sum_",j)]] <- score$NAMED_ALLELE_DOSAGE_SUM * score$SCORE1_AVG
  
    #Drop allele dosage and average as when next chromosome added in cannot join. Left with the sum and the IID
    score <- score[,-c(2,3)]
    FinalScore <- left_join(FinalScore, score)
    print(dim(FinalScore))
  }

  FinalScore[[paste0(pheno[i],"_PRS")]] <- rowSums(FinalScore[,c(2:23)]) / 22
  print(head(FinalScore))  
  fwrite(FinalScore, paste0("/path/to/PRSscore/files/",pheno[i],"_PRS"), sep="\t")
}

