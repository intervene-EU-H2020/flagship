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
bbj<-fread(paste0(gd_dir,"Biobank Japan_HazardRatios/biobank_baseline_summary_stats_BBJ.csv"))
hunt<-fread(paste0(gd_dir,"HUNT_HazardRatios/HUNT_baseline_summary_stats.csv"))
gs<-fread(paste0(gd_dir,"GenerationScotland_HazardRatios/GS_baseline_summary_stats.csv"))
estbb<-fread(paste0(gd_dir,"EstBB_HazardRatios/ESTBB_baseline_summary_stats.csv"))
ge<-fread(paste0(gd_dir, "GenomicsEngland_HazardRatios/GE_baseline_summary_stats_EUR.csv"))
gnh<-fread(paste0(gd_dir,"Genes&Health_HazardRatios/biobank_baseline_summary_stats_GNH.csv"))

#ukb
ukb_sas$Biobank <- "UK Biobank"
ukb_eas$Biobank<-"UK Biobank"
ukb_eur$Biobank<-"UK Biobank"

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
mgb$Ancestry <- "EUR"

#Genomics England
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
ge<-ge %>% mutate(trait=case_when(trait=="Knee_ARTHROSIS"~"KNEE_ARTHROSIS")) #fix phenotype typo
h
            