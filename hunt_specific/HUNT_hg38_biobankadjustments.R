# When converting the adjusted summary statistics to hg38 using Remos dataset, we did not have enough information to know whether A1 was always the reference allele. 
# As such, we have created both possible combinations of the variant_id in the form of CHR_POS_REF_ALT and CHR_POS_ALT_REF.
# This script will identify which one is correct according to your bim file to make sure SNPs are not removed from PRS calculation unecessarily. 

library(data.table)
library(dplyr)
require(R.utils)
library(tidyr)

args <- commandArgs(TRUE)
print(args)
bim_file<-args[1] #bim_file<-"/mnt/scratch/brooke/bcf/all.log.bim" #HUNT bim are in hg19
map_file<-args[2] #map_file<-"/mnt/scratch/brooke/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz"
score_file_path<-args[3] #score_file_path<-"/mnt/scratch/brooke/PRS_v2/"
snplist_file_path<-args[4] #snplist_file_path<-"/mnt/scratch/brooke/flagship/hunt_specific/"
rsid_TF<-as.logical(args[5]) #TRUE if we want to use rsIDs, FALSE if we want to use chr:pos:A2/A1 or A1/A2
pheno<-args[6] #string if you don't want to use the list 

#read bim file
bim <- fread(bim_file, data.table=FALSE)
#give the file column names: I have assumed standard bim format.
colnames(bim) <- c("chrom","variant_id","cM","pos","a1","a2") 

#read mapping file
if (file.exists(map_file)) {
   map<-fread(map_file,data.table=FALSE)
}

if (pheno!=""){ #provide custom phenotype string
   print(pheno)
   phenotypes<-c(pheno)
} else {
#Note: the below code assumes you are using the same filenames as I originally saved
   phenotypes <- c("AAA","Alcohol_Use_Disorder", "AllCancers", "Alzheimers_Disease", "Appendicitis", "Asthma", "Atrial_Fibrillation", "BMI", "Breast_Cancer", "CHD", "Chronic_Kidney_Disease", "Colorectal_Cancer","Educational_Attainment", "Epilepsy", "Focal_Epilepsy", "Generalised_Epilepsy", "Gout", "Heart_Failure", "Hip_Osteoarthritis", "IPF", "ILD", "Inflammatory_Bowel_Disease", "Knee_Osteoarthritis", "Lifespan", "Lung_Cancer", "MDD", "Melanoma", "Osteoporosis", "Pain", "POAG", "Prostate_Cancer", "Rheumatoid_Arthritis", "Sleep_Apnoea", "smoking", "Stroke", "Subarachnoid_Haemmorhage", "TAA", "T1D", "T2D", "Thyroid_Stimulating_Hormone")
}


for(i in phenotypes){ 
  #read in adjusted mega PRS summary statistics
  file<-paste0(score_file_path,i,"_megaPRS_scores_hg19.txt.gz")
  print(file)
  score <- fread(input=file, data.table=FALSE)

  #Use this for reference later
  print(dim(score))

  #Predictor       A1      A2      Centre  Effect_Best and Predictor is rsID
  if (rsid_TF==TRUE){
    #hg19 to hg38, rsID to chr:pos_A1/A2
    merged<-left_join(score,map,by=c("Predictor"="rsid")) #maybe check to see if matches rsid from another build?
    merged$Predictor_v1<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a1,"/",merged$a2)
    merged$Predictor_v2<-paste0(merged$chr,":",merged$pos_hg19,"_",merged$a2,"/",merged$a1)
  } else {
    #   rsid A1 A2   Centre Effect_Best Predictor chr position and Predictor is chr:pos in hg19
    merged<-left_join(score,map,by=c("chr"="chr","position"="pos_hg19"))
    merged$Predictor_v1<-paste0(merged$chr,":",merged$position,"_",merged$a1,"/",merged$a2)
    merged$Predictor_v2<-paste0(merged$chr,":",merged$position,"_",merged$a2,"/",merged$a1)
  }

  v1<-left_join(merged,bim,by=c("Predictor_v1"="variant_id")) %>% drop_na %>% mutate(varid=Predictor_v1)
  v2<-left_join(merged,bim,by=c("Predictor_v2"="variant_id")) %>% drop_na %>% mutate(varid=Predictor_v2)
  all<-rbind(v1,v2)
  if (rsid_TF==TRUE){
    all_sorted<-arrange(all,chr,pos_hg19)
    score2 <- all_sorted[,c("varid", "A1", "A2", "Centre", "Effect_Best")]
  } else{
    all_sorted<-arrange(all,chr,position) %>%select(-rsid.y) %>% rename(rsid=rsid.x)
    score2 <- all_sorted[,c("rsid", "A1", "A2", "Centre", "Effect_Best","varid","chr","position")]
  }
  
 
  #Check to make sure you haven't lost a whole bunch of SNPs from this.
  print(dim(score2))

  #Save adjusted score file so that it can be read by plink
  fwrite(score2, paste0(score_file_path,i,"_megaPRS_scores_hg19_varid.txt"), sep="\t")
}


#### also create snp list file (addapted from create_snplist.R)

#only want to do once, not for every phenotype since htis is mapping and bim files and not involving score files 
if (!file.exists(paste0(snplist_file_path,"snplist_hg19_rsid")) | !file.exists(paste0(snplist_file_path,"snplist_hg19_varid"))){
  mapping <- fread(map_file,data.table=FALSE)
  mapping$Predictor_v1 <-paste0(mapping$chr,":",mapping$pos_hg19,"_",mapping$a1,"/",mapping$a2)
  mapping$Predictor_v2 <-paste0(mapping$chr,":",mapping$pos_hg19,"_",mapping$a2,"/",mapping$a1)

  snplist <- mapping[,c("Predictor_v1", "Predictor_v2","rsid")]
  print(nrow(snplist))
  #identify which of the two variant_ids are found in your bim file
  first_id <- subset(snplist, Predictor_v1 %in% bim$variant_id)
  first_id <- first_id$Predictor_v1
  second_id <- subset(snplist, Predictor_v2 %in% bim$variant_id)
  second_id <- second_id$Predictor_v2

  snplist$Predictor <- case_when(snplist$Predictor_v1 %in% first_id ~ snplist$Predictor_v1,
                               snplist$Predictor_v2 %in% second_id ~ snplist$Predictor_v2,
                               TRUE ~ NA_character_)
  #Save adjusted snplist file so that it can be read by plink
  print(table(is.na(snplist$Predictor))) #how many match?
  snplist_varid <- snplist[!is.na(snplist$Predictor),]$Predictor
  write.table(snplist_varid, paste0(snplist_file_path,"snplist_hg19_varid"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  #write out rsIDs for comparison with other biobanks SNP coverage
  snplist_rsid <- snplist[!is.na(snplist$Predictor),]$rsid
  write.table(snplist_rsid, paste0(snplist_file_path,"snplist_hg19_rsid"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
