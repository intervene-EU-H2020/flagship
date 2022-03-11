# When converting the adjusted summary statistics to hg38 using Remos dataset, we did not have enough information to know whether A1 was always the reference allele. 
# As such, we have created both possible combinations of the variant_id in the form of CHR_POS_REF_ALT and CHR_POS_ALT_REF.
# This script will identify which one is correct according to your bim file to make sure SNPs are not removed from PRS calculation unecessarily. 

library(data.table)
library(dplyr)
require(R.utils)

args <- commandArgs(TRUE)
bim_file<-args[1]
map_file<-args[2]

#read bim file
bim <- fread(bim_file, data.table=FALSE)

#read mapping file
if (map_file) {
   map<-fread(map_file,data.table=FALSE)
}

/home/bwolford/scratch/brooke
bcf/PART_09.bim
1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz
#give the file column names: I have assumed standard bim format.
colnames(bim) <- c("chr","variant_id","cM","pos","a1","a2") 

#Note: the below code assumes you are using the same filenames as I originally saved
phenotypes <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "BMI", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Educational_Attainment", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "IPF", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lifespan", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "Pain", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "smoking", "Stroke", "Subarachnoid_Haemmorhage", "TAA", "T1D", "T2D", "Thyroid_Stimulating_Hormone")

for(i in phenotypes){ 
#read in adjusted mega PRS summary statistics
score <- fread(input=paste0("/home/bwolford/scratch/brooke/PRS/",i,"_megaPRS_scores_hg19.txt.gz"), data.table=FALSE)

#Use this for reference later
print(dim(score))

merged<-left_join(score,map,by=c("Predictor"="rsid"))
merged$Predictor_v1<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a1,"/",merged$a2)
merged$Predictor_v2<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a2,"/",merged$a1)


score <- score[,c("Predictor", "A1", "A2", "Centre", "Effect_Best")]

#Check to make sure you haven't lost a whole bunch of SNPs from this.
print(dim(score))

#Save adjusted score file so that it can be read by plink
fwrite(score, paste0("/path/to/score/file/",i,"_megaPRS_scores_hg19.txt"), sep="\t")
}
