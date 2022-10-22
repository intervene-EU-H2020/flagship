library(data.table)
library(dplyr)
library(metafor)
library(grid)
library(ggplot2)
library(googlesheets4)
library(googledrive)
library(stringr)
library(scales)
library(forcats)
library(RColorBrewer)
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/MetaAnalysis/"
gd_dir<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/"


##### Read in data from Google Drive

#identify folder
folder_id = drive_get(as_id("1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ"))

#find files in folder
files = drive_ls(folder_id)

#loop dirs and download files inside them
for (i in seq_along(files$name)) {
  #list files
  i_dir = drive_ls(files[i, ])
  
  #mkdir
  try({dir.create(paste0(gd_dir,files$name[i]))})
  
  #download files
  for (file_i in seq_along(i_dir$name)) {
    #fails if already exists
    try({
      drive_download(
        as_id(i_dir$id[file_i]),
        path = paste0(gd_dir,files$name[i], "/",i_dir$name[file_i])
      )
    })
  }
}

#estb has weird file structure so do it separately
dir<-drive_ls(as_id("1CySw57_ICg0e-DW3M4wiuRTccmp8QYAn"))
dir.create(paste0(gd_dir,"EstBB_HazardRatios/"))
for (idx in seq_along(dir$name)) {
  #fails if already exists
  try({
    drive_download(
      as_id(dir$id[idx]),
      path = paste0(gd_dir,"EstBB_HazardRatios/",dir$name[idx])
    )
  })
}


#list.files(path=gd_dir,pattern="baseline_summary_stats",recursive=TRUE)
finngen<-fread(paste0(gd_dir,"FinnGen_HazardRatios/biobank_baseline_summary_stats_FinnGen.csv"))
ukb_sas<-fread(paste0(gd_dir,"UKB_HazardRatios/biobank_baseline_summary_stats_UKBSAS.csv"))
ukb_eas<-fread(paste0(gd_dir,"UKB_HazardRatios/biobank_baseline_summary_stats_UKBEAS.csv"))
ukb_eur<-fread(paste0(gd_dir,"UKB_HazardRatios/biobank_baseline_summary_stats_UKBEUR.csv"))
mgb<-fread(paste0(gd_dir,"MGB_HazardRatios/biobank_baseline_summary_stats_removedups_093022.csv"))
bbj<-fread(paste0(gd_dir,"Biobank_Japan_HazardRatios/biobank_baseline_summary_stats_BBJ.csv"))
hunt<-fread(paste0(gd_dir,"HUNT_HazardRatios/HUNT_baseline_summary_stats.csv"))
gs<-fread(paste0(gd_dir,"GenerationScotland_HazardRatios/GS_baseline_summary_stats.csv"))
estbb<-fread(paste0(gd_dir,"EstBB_HazardRatios/ESTBB_baseline_summary_stats.csv"))
ge<-fread(paste0(gd_dir, "GenomicsEngland_HazardRatios/GE_baseline_summary_stats_EUR.csv"))
gnh<-fread(paste0(gd_dir,"GNH_HazardRatios/biobank_baseline_summary_stats_GNH.csv"))

#ukb
ukb_sas$Biobank <- "UK Biobank"
ukb_eas$Biobank<-"UK Biobank"
ukb_eur$Biobank<-"UK Biobank"
ukb_sas$Ancestry<-"SAS"
ukb_eas$Ancestry<-"EAS"
ukb_eur$Ancestry<-"EUR"

#HUNT 
drophunt <-c("T1D","C3_CANCER") #these are weird on first pass 
hunt<-subset(hunt,!(trait %in% drophunt))
hunt$Biobank<-"HUNT"
hunt$Ancestry<-"EUR"

#FinnGen
finngen$Biobank <- "FinnGen"
finngen$Ancestry <- "EUR"

#Genes&Health
gnh$Biobank <- "Genes & Health"
gnh$Ancestry <- "SAS"

#Biobank Japan
dropbbj <- c("I9_AF","C3_BRONCHUS_LUNG","RHEUMA_SEROPOS_OTH","ILD","GOUT")
bbj <- subset(bbj, !(trait %in% dropbbj))
bbj$Biobank <- "Biobank Japan"
bbj$Ancestry <- "EAS"

#Estonian Biobank
estbb$Biobank <- "Estonian Biobank"
estbb$Ancestry <- "EUR"

#Generation Scotland
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(trait %in% dropgs))
gs$Biobank <- "Generation Scotland"
gs$Ancestry <- "EUR"

#Mass General Brigham (probably need to redo for african ancestry)
dropmgb <- c("I9_AF")
mgb <- subset(mgb, !(trait %in% dropmgb))
mgb$Biobank <- "Mass General Brigham"

#Genomics England
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
ge<-ge %>% mutate(trait=ifelse(trait=="Knee_ARTHROSIS","KNEE_ARTHROSIS",trait)) #fix phenotype typo

all <- rbind(gnh, finngen,fill=TRUE) %>%
  rbind(., estbb,fill=TRUE) %>%
  rbind(., gs,fill=TRUE) %>%
  rbind(., mgb,fill=TRUE) %>% 
  rbind(., ukb_eur,fill=TRUE) %>% 
  rbind(., ukb_eas,fill=TRUE) %>% 
  rbind(., ukb_sas,fill=TRUE) %>% 
  rbind(., bbj,fill=TRUE) %>% 
  rbind(., ge,fill=TRUE) %>% 
  rbind(.,hunt,fill=TRUE)

gbd_phenos <- c("Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer", "Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
hr_phenos <- c("ILD", "C3_BRONCHUS_LUNG","C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
labels<-c("ILD","Lung Cancer","All Cancers","Appendicitis", "Asthma", "Atrial Fibrillation", "Breast Cancer","Coronary Heart Disease","Colorectal Cancer","Epilepsy","Gout","Hip Osteoarthritis","Knee Osteoarthritis","Major Depression","Skin Melanoma","Prostate Cancer","Rheumatoid Arthritis","Type 1 Diabetes","Type 2 Diabetes")
label_df<-data.frame(gbd=gbd_phenos,hr=hr_phenos,pretty=labels)
all<-all %>% left_join(label_df,by=c("trait"="hr")) %>% mutate(biobank_anc=paste0(Biobank,"(",Ancestry,")"))
#pretty labels don't have all the traits we did summary statistics for
  
subset<-all %>% filter(trait!="C3_BREAST" & trait!="C3_PROSTATE" & !is.na(pretty) & Ancestry %in% c("EUR","SAS","EAS","AFR"))
pdf(file=paste0(output_dir,"correlations.pdf"),height=10,width=12)
ggplot(subset,aes(y=pretty,x=biobank_anc,fill=sex_corr)) + geom_tile() + theme_bw() +
  scale_fill_gradient2(low="#3C6BAF",mid="white",high="#864684",name="Correlation",midpoint=0,limits=c(-0.25,0.25)) +
  guides(fill = guide_legend(title = "Correlation\nwith male\nsex")) +
  theme(title = element_text(size = 22),
        panel.grid = element_blank(),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24,angle=45,hjust=1),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 28))
dev.off()
