# When converting the adjusted summary statistics to hg38 using Remos dataset, we did not have enough information to know whether A1 was always the reference allele. 
# As such, we have created both possible combinations of the variant_id in the form of CHR_POS_REF_ALT and CHR_POS_ALT_REF.
# This script will identify which one is correct according to your bim file to make sure SNPs are not removed from PRS calculation unecessarily. 

library(data.table)
library(dplyr)
require(R.utils)
library(tidyr)

args <- commandArgs(TRUE)
bim_file<-args[1] #bim_file<-"/mnt/scratch/brooke/bcf/PART_09.bim" #HUNT bim are in hg19
map_file<-args[2] #map_file<-"/mnt/scratch/brooke/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz"
score_file_path<-args[3] #"/mnt/scratch/brooke/PRS/"

#read bim file
bim <- fread(bim_file, data.table=FALSE)

#read mapping file
if (file.exists(map_file)) {
   map<-fread(map_file,data.table=FALSE)
}

#give the file column names: I have assumed standard bim format.
colnames(bim) <- c("chrom","variant_id","cM","pos","a1","a2") 

#Note: the below code assumes you are using the same filenames as I originally saved
phenotypes <- c("Alcohol_Use_Disorder", "Alzheimers_Disease", "Asthma", "Atrial_Fibrillation", "BMI", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Educational_Attainment", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "IPF", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lifespan", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "Pain", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "smoking", "Stroke", "Subarachnoid_Haemmorhage", "TAA", "T1D", "T2D", "Thyroid_Stimulating_Hormone")

for(i in phenotypes){ 
  #read in adjusted mega PRS summary statistics
  file<-paste0(score_file_path,i,"_megaPRS_scores_hg19.txt.gz")
  print(file)
  score <- fread(input=file, data.table=FALSE)

  #Use this for reference later
  print(dim(score))

  #hg19 to hg38, rsID to chr:pos_A1/A2
  merged<-left_join(score,map,by=c("Predictor"="rsid"))
  merged$Predictor_v1<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a1,"/",merged$a2)
  merged$Predictor_v2<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a2,"/",merged$a1)

  v1<-left_join(merged,bim,by=c("Predictor_v1"="variant_id")) %>% drop_na %>% mutate(varid=Predictor_v1)
  v2<-left_join(merged,bim,by=c("Predictor_v2"="variant_id")) %>% drop_na %>% mutate(varid=Predictor_v2)
  all<-rbind(v1,v2)
  all_sorted<-arrange(all,chr,pos_hg19)
  
  score2 <- all_sorted[,c("varid", "A1", "A2", "Centre", "Effect_Best")]
 
  #Check to make sure you haven't lost a whole bunch of SNPs from this.
  print(dim(score2))

  #Save adjusted score file so that it can be read by plink
  fwrite(score2, paste0(score_file_path,i,"_megaPRS_scores_hg38_varid.txt"), sep="\t")
}
