#Libraries
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(googledrive)
library(purrr)
library(forcats)
library(stringr)

brewer.pal(n = 3, name = 'Dark2')

setwd("/mnt/work/workbench/bwolford/multi_ancestry")
#setwd("/net/snowwhite/home/bwolford/2021_analysis/intervene")

#read in files from google drive 
#issue with authorization but I think we can proceed without
#https://cran.r-project.org/web/packages/gargle/vignettes/auth-from-web.html

folder_url<-"https://drive.google.com/drive/u/0/folders/1W3tqDX6K51TfKZhHabi4Gn5xVcjLQZsB"
folder<-drive_get(as_id(folder_url))
files <- drive_ls(folder$id)

for (i in seq_along(files$name)) {
  #list files
  i_dir = drive_ls(files[i, ])
  
  #mkdir
  dir.create(files$name[i])
  
  #download files
  for (file_i in seq_along(i_dir$name)) {
    #fails if already exists
    try({
      drive_download(
        as_id(i_dir$id[file_i]),
        path = str_c(files$name[i], "/", i_dir$name[file_i])
      )
    })
  }
}

#gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
#hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
#biobanks <-c("BBJ","GNH")
#countries<-c("Japan","England")

code<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/scripts"

#Rscript Sex_LifetimeRiskEstimation_Bootstrap_script.R --country="Japan" --csv="/mnt/work/workbench/bwolford/multi_ancestry/BiobankJapan/HR_MaleSample_BBJ.csv" --sex=Male --biobank="BBJ" --hr_pheno="T2D" --gbd_pheno="Diabetes mellitus type 2" --bootstrap=200 &
#Rscript Sex_LifetimeRiskEstimation_Bootstrap_script.R --country="Japan" --csv="/mnt/work/workbench/bwolford/multi_ancestry/BiobankJapan/HR_FemaleSample_BBJ.csv" --sex=Female --biobank="BBJ" --hr_pheno="T2D" --gbd_pheno="Diabetes mellitus type 2" --bootstrap=200 &

#Rscript Sex_LifetimeRiskEstimation_Bootstrap_script.R --country="England" --csv="/mnt/work/workbench/bwolford/multi_ancestry/Genes_and_Health/HR_FemaleSample_GNH.csv" --sex=Female --biobank="GNH" --hr_pheno="I9_CHD" --gbd_pheno="Ischemic heart disease" --bootstrap=200
#Rscript Sex_LifetimeRiskEstimation_Bootstrap_script.R --country="England" --csv="/mnt/work/workbench/bwolford/multi_ancestry/Genes_and_Health/HR_MaleSample_GNH.csv" --sex=Male --biobank="GNH" --hr_pheno="I9_CHD" --gbd_pheno="Ischemic heart disease" --bootstrap=200 
  