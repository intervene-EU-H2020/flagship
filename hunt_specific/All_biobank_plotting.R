#Libraries
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(googledrive)
library(purrr)

brewer.pal(n = 3, name = 'Dark2')

setwd("/mnt/work/workbench/bwolford/intervene/results")

file<-drive_find("INTERVENE flagship/Manuscript/Figures")
drive_download(file$name,recursive=TRUE)
#will have to authorize your google account, check option to see accounts, wil copy and apste an authorization code 

ukbb_full<-drive_ls(as_id("1ZwFUsYnU8VGcUcEMJ0jG3iq99UulLtXu"),type="csv")
estb_full<-drive_ls(as_id("1dClpNyKv788ouoFoIEq_vZXa04KOg9ai"),type="csv")
finn_full<-drive_ls(as_id("1GmG4GBPgnpq340jAMxlYSxWJw0P5l5gm"),type="csv")
hunt_full<-drive_ls(as_id("1ICKQLXN2JTyqUOhzwVO6iY8WMiLqGyjW"),type="csv")
walk(c(estb_full$id,finn_full$id,hunt_full$id,ukbb_full$id), ~ drive_download(as_id(.x)))



hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")


biobank<-c("EstBB","HUNT","Finngen","UKBiobank")
for (b in 1:length(biobank)){
  for (p in 1:length(hr_phenos){
    fn<-paste(sep="_",hr_phenos[p],"LifetimeRisk","BootstrappedConfidenceIntervals",biobank[b])
    if (file.exts(fn)){
      df<-fread(rn)
      if (b==1 & p==1){
        all<-df
      } else {
        all<-rbind(all,df)
      }
      
    }}}

