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
    score <- score[,c("IID","ALLELE_CT","SCORE1_AVG")]
    score[[paste0("PRS_Sum_",j)]] <- score$ALLELE_CT * score$SCORE1_AVG
  	colnames(score)[2] <- paste0("ALLELE_CT_",j)
  	
    #Drop average score as when next chromosome added in cannot join. Left join with the sum and the IID
    score <- score[,-3]
    FinalScore <- left_join(FinalScore, score)
    print(dim(FinalScore))
  }

  FinalScore[[paste0(pheno[i],"_PRS")]] <- rowSums(FinalScore[,c(seq(3,45,by=2))]) / rowSums(FinalScore[,c(seq(2,44,by=2))])
  print(head(FinalScore))  
  fwrite(FinalScore, paste0("/path/to/PRSscore/files/",pheno[i],"_PRS"), sep="\t")
}
