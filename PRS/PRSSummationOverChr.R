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
  
    #Take IID, ALLELE_CT and SCORE1_AVG
    score <- score[,c("IID","SCORE1_SUM")]
  	colnames(score)[2] <- paste0(pheno[i],"_SUM_",j)
  	
    #Drop average score as when next chromosome added in cannot join. Left join with the sum and the IID
    FinalScore <- left_join(FinalScore, score)
    print(dim(FinalScore))
  }

  FinalScore[["SCORE1_SUM"]] <- rowSums(FinalScore[,c(2:23)]) 
  print(head(FinalScore))  
  fwrite(FinalScore, paste0("/path/to/PRSscore/files/",pheno[i],"_PRS.sscore"), sep="\t")
}
