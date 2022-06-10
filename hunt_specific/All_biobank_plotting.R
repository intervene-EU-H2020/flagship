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
#setwd("/net/snowwhite/home/bwolford/2021_analysis/intervene")

#download files from google drive into local folder
#file<-drive_find("INTERVENE flagship/Manuscript/Figures")
#drive_download(file$name,recursive=TRUE)
#will have to authorize your google account, check option to see accounts, wil copy and apste an authorization code 

#folder_url<-"https://drive.google.com/drive/u/1/folders/1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ"
#folder<-drive_get(as_id(folder_url))

ukbb_full<-drive_ls(as_id("1ZwFUsYnU8VGcUcEMJ0jG3iq99UulLtXu"),type="csv")
estb_full<-drive_ls(as_id("1dClpNyKv788ouoFoIEq_vZXa04KOg9ai"),type="csv")
finn_full<-drive_ls(as_id("1GmG4GBPgnpq340jAMxlYSxWJw0P5l5gm"),type="csv")
hunt_full<-drive_ls(as_id("1ICKQLXN2JTyqUOhzwVO6iY8WMiLqGyjW"),type="csv")
walk(c(estb_full$id,finn_full$id,hunt_full$id,ukbb_full$id), ~ drive_download(as_id(.x),overwrite=TRUE))

estb_age<-drive_ls(as_id("1ePW6dcW79LGnZ9I5oTzNNK-xB5SJ-hfn"),type="csv")
ukbb_age<-drive_ls(as_id("1LieYp_dhXihEoU3QMDljiAehuVdr3f9o"),type="csv")
hunt_age<-drive_ls(as_id("1wziRjQMGAvfza_5nPeO2uDAQTATKl4CP"),type="csv")
finn_age<-drive_ls(as_id("1OgFSitzUxBeZTDba4z8kcjayCSRBY7Lw"),type="csv")
walk(c(estb_age$id,finn_age$id,hunt_age$id,ukbb_age$id), ~ drive_download(as_id(.x),overwrite=TRUE))

finn_male<-drive_ls(as_id("1-rjyJuUYJIz61w_FtmyPPxBzrx95hKcc"),type="csv")


#hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
phenocols <- c("J10_ASTHMA", "C3_CANCER", "K11_APPENDACUT", "I9_AF", "C3_BREAST", "I9_CHD","C3_COLORECTAL",  "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "I9_SAH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
#prscols <- c("Asthma","AllCancers","Appendicitis", "Atrial_Fibrillation", "Breast_Cancer", "CHD","Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer","Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D","T2D", "ILD", "Lung_Cancer")
#gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")

string<-"LifetimeRisk_AgeStratification_BootstrappedConfidenceIntervals"
age_strat<-all

string<-"LifetimeRisk_BootstrappedConfidenceIntervals"
full<-all

biobank<-c("EstBB","HUNT","Finngen","UKBiobank")
for (b in 1:length(biobank)){
  for (p in 1:length(phenocols)){
    fn<-paste0(phenocols[p],"_",string,"_",biobank[b],".csv")
    if (file.exists(fn)){
      df<-fread(fn)
      df$pheno<-phenocols[p]
      df$biobank<-biobank[b]
      if (b==1 & p==1){
        all<-df
      } else {
        all<-rbind(all,df)
      }
    }}}

age_strat$analysis<-"age_strat"
full$analysis<-"full"
df<-as.data.frame(rbind(age_strat,full))
df$label<-as.factor(paste(sep="_",df$analysis,df$Group))

#df2<-df %>% subset(Group=="Group6" & analysis=="full" & pheno=="I9_CHD")
##### PROSTATE CANCER PLOTS ######
max<-as.numeric(df %>% subset(pheno=="C3_PROSTATE") %>% summarize(max(LifetimeRisk)))
colors<-c("light green","dark green","plum2","orchid4")

#median risk, full sample
df2<-df %>% subset(Group=="Group6" & analysis=="full" & pheno=="C3_PROSTATE")
pdf(file="prostate_cancer_full_median_risk.pdf",width=10,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  ylim(0,max) +
  xlab("Age Range") + 
  ylab("Cumulative Risk (%)") + 
  theme_bw() +   theme(legend.position="bottom",title = element_text(size = 22),
                       legend.text = element_text(size = 16),
                       legend.title = element_text(size = 18),
                       axis.title.x = element_text(size = 18),
                       axis.text.x = element_text(size = 12, angle=45, hjust=1),
                       axis.title.y = element_text(size = 18),
                       axis.text.y = element_text(size = 16))
dev.off()

#median risk, age stratified HR
df2<-df %>% subset(Group=="Group6" & pheno=="C3_PROSTATE") %>% mutate(label= fct_relevel(label,"full_Group6"))
pdf(file="prostate_cancer_age_full_median_risk.pdf",width=10,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Cumulative Risk (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      legend.text = element_text(size = 16),
                      legend.title = element_text(size = 18),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

#add in top risk, full and median
df2<-df %>% subset((Group=="Group6"  | Group =="Group10") & pheno=="C3_PROSTATE") %>% 
  mutate(label=fct_relevel(label,"full_Group6","age_strat_Group6,full_Group10,age_strat_Group10"))
pdf(file="prostate_cancer_age_full_with_top_risk.pdf",width=10,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Cumulative Risk (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      legend.text = element_text(size = 16),
                      legend.title = element_text(size = 18),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

#Considering confidence intervals
#riskwithintervals <- subset(all, Group=="Group1" | Group=="Group6" | Group=="Group11")
#labs(color='PRS Group', fill='PRS Group') +
# scale_color_hue(labels = c("0-1%", "40-60%", "99-100%")) +
#scale_fill_hue(labels = c("0-1%", "40-60%", "99-100%")) +

pdf(file="test.pdf",height=10,width=10)
ggplot(all, aes(Age, LifetimeRisk, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~pheno+biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age Range") + 
  ylab("Cumulative Risk (%)") + 
  theme_bw() +
  
  theme(title = element_text(size = 22),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, angle=-90, hjust=0),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16))
dev.off()




