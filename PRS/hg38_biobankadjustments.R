# When converting the adjusted summary statistics to hg38 using Remos dataset, we did not have enough information to know whether A1 was always the reference allele. 
# As such, we have created both possible combinations of the variant_id in the form of CHR_POS_REF_ALT and CHR_POS_ALT_REF.
# This script will identify which one is correct according to your bim file to make sure SNPs are not removed from PRS calculation unecessarily. 

library(data.table)
library(dplyr)

#read bim file
bim <- fread("/path/to/bim/filename.bim", data.table=FALSE)

#give the file column names: I have assumed standard bim format.
colnames(bim) <- c("chr","variant_id","cM","pos","a1","a2") 

#Note: the below code assumes you are using the same filenames as I originally saved
phenotypes <- c("AllCancers", "Appendicitis", "Colorectal_Cancer", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Epilepsy", "Gout", "Hip_Osteoarthritis", "ILD", "Knee_Osteoarthritis", "Lung_Cancer", "MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "T1D", "T2D")

for(i in phenotypes){ 
#read in adjusted mega PRS summary statistics
score <- fread(input=paste0("/path/to/score/file/",i,"_megaPRS_scores_hg38.txt"), data.table=FALSE)

#Use this for reference later
print(dim(score))

#identify which of the two variant_ids are found in your bim file
first_id <- subset(score, Predictor_v1 %in% bim$variant_id)
first_id <- first_id$Predictor_v1
second_id <- subset(score, Predictor_v2 %in% bim$variant_id)
second_id <- second_id$Predictor_v2

score$Predictor <- case_when(score$Predictor_v1 %in% first_id ~ score$Predictor_v1,
                             score$Predictor_v2 %in% second_id ~ score$Predictor_v2,
                             TRUE ~ NA_character_)

score <- score[,c("Predictor", "A1", "A2", "Centre", "Effect_Best")]

#Check to make sure you haven't lost a whole bunch of SNPs from this.
print(dim(score))

#Save adjusted score file so that it can be read by plink
fwrite(score, paste0("/path/to/score/file/",i,"_megaPRS_scores_hg38.txt"), sep="\t")
}
