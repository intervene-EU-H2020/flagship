#Libraries
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(googledrive)
library(purrr)
library(forcats)

brewer.pal(n = 3, name = 'Dark2')

setwd("/mnt/work/workbench/bwolford/intervene/results")
#setwd("/net/snowwhite/home/bwolford/2021_analysis/intervene")

#download files from google drive into local folder
#file<-drive_find("INTERVENE flagship/Manuscript/Figures")
#drive_download(file$name,recursive=TRUE)
#will have to authorize your google account, check option to see accounts, wil copy and apste an authorization code 

#folder_url<-"https://drive.google.com/drive/u/1/folders/1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ"
#folder<-drive_get(as_id(folder_url))

#drive_auth(reset=TRUE)

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
ukbb_male<-drive_ls(as_id("111421mqtOe0b7n5WnncvheyIY-b7_HQq"),type="csv")
hunt_male<-drive_ls(as_id("1EZc6YQDyH7JLtUMtLLN9geshoAHnsOry"),type="csv")
estb_male<-c()
walk(c(finn_male$id,hunt_male$id,ukbb_male$id), ~ drive_download(as_id(.x),overwrite=TRUE))

finn_female<-drive_ls(as_id("1zrgl12J93csdkxes-3e7aS-vZBcmOFVm"),type="csv")
ukbb_female<-drive_ls(as_id("1H0Po0OFsYBiIwDKTImcLMeFUr02WVFGe"),type="csv")
hunt_female<-drive_ls(as_id("13awAoEV0to4ImOcNcwtQTMWmjszIHWH6"),type="csv")
estb_female<-c()
walk(c(finn_female$id,hunt_female$id,ukbb_female$id), ~ drive_download(as_id(.x),overwrite=TRUE))

#hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
phenocols <- c("J10_ASTHMA", "C3_CANCER", "K11_APPENDACUT", "I9_AF", "C3_BREAST", "I9_CHD","C3_COLORECTAL",  "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "I9_SAH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
#prscols <- c("Asthma","AllCancers","Appendicitis", "Atrial_Fibrillation", "Breast_Cancer", "CHD","Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer","Rheumatoid_Arthritis", "Subarachnoid_Haemmorhage", "T1D","T2D", "ILD", "Lung_Cancer")
#gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")

### function to read the files 
combine_files<-function(string){
  biobank<-c("EstBB","HUNT","FinnGen","UKBiobank")
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
      }
   }
  }
  return(all)
}
string<-"LifetimeRisk_AgeStratification_BootstrappedConfidenceIntervals"
age_strat<-combine_files(string)

string<-"LifetimeRisk_BootstrappedConfidenceIntervals"
full<-combine_files(string)

string<-"Male_LifetimeRisk_BootstrappedConfidenceIntervals"
male<-combine_files(string)

string<-"Female_LifetimeRisk_BootstrappedConfidenceIntervals"
female<-combine_files(string)

age_strat$analysis<-"age_strat"
full$analysis<-"full"
female$analysis<-"female"
male$analysis<-"male"

###### CAD PLOTS ######
df<-as.data.frame(rbind(age_strat,full))
df$label<-as.factor(paste(sep="_",df$analysis,df$Group))
max<-as.numeric(df %>% subset(pheno=="I9_CHD") %>% summarize(max(CIpos)))
colors<-c("light green","dark green","plum2","orchid4")

#median risk, full sample
df2<-df %>% subset(Group=="Group6" & analysis=="full" & pheno=="I9_CHD") %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="CHD_full_median_risk.pdf",width=11,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  ylim(0,max) +
  xlab("Age Range") + 
  ylab("Lifetime Risk of Coronary Heart Disease (%)") + 
  theme_bw() +   theme(legend.position="bottom",title = element_text(size = 22),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text = element_text(size = 20, margin = margin()),
                       legend.text = element_text(size = 16),
                       legend.title = element_blank(),
                       axis.title.x = element_text(size = 18),
                       axis.text.x = element_text(size = 12, angle=45, hjust=1),
                       axis.title.y = element_text(size = 18),
                       axis.text.y = element_text(size = 16))
dev.off()

#median risk, age stratified HR
df2<-df %>% subset(Group=="Group6" & pheno=="I9_CHD") %>% 
  mutate(label=fct_relevel(label,"full_Group6")) %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="CHD_age_full_median_risk.pdf",width=11,height=6)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk of Coronary Heart Disease (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

#add in top risk, full and median
df2<-df %>% subset((Group=="Group6"  | Group =="Group10") & pheno=="I9_CHD") %>% 
  mutate(label=fct_relevel(label,"full_Group6","age_strat_Group6","full_Group10","age_strat_Group10")) %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="CHD_age_full_with_top_risk.pdf",width=11,height=6)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk of Coronary Heart Disease (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

###### T2D PLOTS ######
df<-as.data.frame(rbind(age_strat,full))
df$label<-as.factor(paste(sep="_",df$analysis,df$Group))
max<-as.numeric(df %>% subset(pheno=="T2D") %>% summarize(max(CIpos)))
colors<-c("light green","dark green","plum2","orchid4")

#median risk, full sample
df2<-df %>% subset(Group=="Group6" & analysis=="full" & pheno=="T2D") %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="T2D_full_median_risk.pdf",width=11,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  ylim(0,max) +
  xlab("Age Range") + 
  ylab("Lifetime Risk of Type 2 Diabetes (%)") + 
  theme_bw() +   theme(legend.position="bottom",title = element_text(size = 22),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text = element_text(size = 20, margin = margin()),
                       legend.text = element_text(size = 16),
                       legend.title = element_blank(),
                       axis.title.x = element_text(size = 18),
                       axis.text.x = element_text(size = 12, angle=45, hjust=1),
                       axis.title.y = element_text(size = 18),
                       axis.text.y = element_text(size = 16))
dev.off()

#median risk, age stratified HR
df2<-df %>% subset(Group=="Group6" & pheno=="T2D") %>% 
  mutate(label=fct_relevel(label,"full_Group6")) %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="T2D_age_full_median_risk.pdf",width=11,height=6)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk of Type 2 Diabetes (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

#add in top risk, full and median
df2<-df %>% subset((Group=="Group6"  | Group =="Group10") & pheno=="T2D") %>% 
  mutate(label=fct_relevel(label,"full_Group6","age_strat_Group6","full_Group10","age_strat_Group10")) %>%
  mutate(Age=fct_relevel(Age,"1 to 4","5 to 9"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="T2D_age_full_with_top_risk.pdf",width=11,height=6)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk of Type 2 Diabetes (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()


##### PROSTATE CANCER PLOTS ######
df<-as.data.frame(rbind(age_strat,full))
df$label<-as.factor(paste(sep="_",df$analysis,df$Group))
max<-as.numeric(df %>% subset(pheno=="I9_CHD") %>% summarize(max(CIpos)))
colors<-c("light green","dark green","plum2","orchid4")

#median risk, full sample
df2<-df %>% subset(Group=="Group6" & analysis=="full" & pheno=="C3_PROSTATE")
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="prostate_cancer_full_median_risk.pdf",width=11,height=6)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  ylim(0,max) +
  xlab("Age Range") + 
  ylab("Lifetime Risk (%)") + 
  theme_bw() +   theme(legend.position="bottom",title = element_text(size = 22),
                       strip.background = element_rect(color="black", fill="white"),
                       strip.text = element_text(size = 20, margin = margin()),
                       legend.text = element_text(size = 16),
                       legend.title = element_blank(),
                       axis.title.x = element_text(size = 18),
                       axis.text.x = element_text(size = 12, angle=45, hjust=1),
                       axis.title.y = element_text(size = 18),
                       axis.text.y = element_text(size = 16))
dev.off()

#median risk, age stratified HR
df2<-df %>% subset(Group=="Group6" & pheno=="C3_PROSTATE") %>% mutate(label=fct_relevel(label,"full_Group6"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                age_strat_Group10="Age Stratified, >95%",
                full_Group6="Standard, 40-60%",
                full_Group10="Standard, >95%")
pdf(file="prostate_cancer_age_full_median_risk.pdf",width=11,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

#add in top risk, full and median
df2<-df %>% subset((Group=="Group6"  | Group =="Group10") & pheno=="C3_PROSTATE") %>% 
  mutate(label=fct_relevel(label,"full_Group6","age_strat_Group6","full_Group10","age_strat_Group10"))
df2$label<-recode(df2$label,age_strat_Group6="Age Stratified, 40-60%",
                  age_strat_Group10="Age Stratified, >95%",
                  full_Group6="Standard, 40-60%",
                  full_Group10="Standard, >95%")
pdf(file="prostate_cancer_age_full_with_top_risk.pdf",width=11,height=5)
ggplot(df2, aes(Age, LifetimeRisk, color=label, group=label)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point() + facet_wrap(~biobank) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=label), alpha=0.2) +
  ylim(0,max)+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)+
  xlab("Age Range") + 
  ylab("Lifetime Risk (%)") + 
  theme_bw() +  theme(legend.position="bottom",title = element_text(size = 22),
                      strip.background = element_rect(color="black", fill="white"),
                      strip.text = element_text(size = 20, margin = margin()),
                      legend.text = element_text(size = 16),
                      legend.title = element_blank(),
                      axis.title.x = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle=45, hjust=1),
                      axis.title.y = element_text(size = 18),
                      axis.text.y = element_text(size = 16))
dev.off()

df2[df2$Age=="75 to 79",] %>% select(LifetimeRisk,label)
############################

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




