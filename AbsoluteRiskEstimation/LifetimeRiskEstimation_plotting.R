library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

results_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #set this to your directory of choice
files<-list.files(path=results_dir,pattern=".csv")
my_files<-files[grepl("Bootstrapped",files)&!grepl("Male",files)&!grepl("Female",files)&!grepl("AgeStrat",files)] #remove male and female for now
#"","Age","Group","CIneg","CIpos","LifetimeRisk"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/"

results<-c()
for (f in 1:length(my_files)){
  df<-fread(paste0(results_dir,my_files[f]))
  trait<-unlist(strsplit(my_files[f],"_LifetimeRisk_BootstrappedConfidenceIntervals_"))[1]
  biobank<-unlist(strsplit(unlist(strsplit(my_files[f],"_LifetimeRisk_BootstrappedConfidenceIntervals_"))[2],".csv"))[1]
  df$trait<-trait
  df$biobank<-biobank
  results<-rbind(results,df,fill=TRUE)
}
results <- results[!is.na(results$CIpos),]
results$Age <- factor(results$Age,levels=c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79"))
results$Group <- factor(results$Group, levels=c("Group1","Group2","Group3","Group4","Group5","Group6","Group7","Group8","Group9","Group10","Group11"))

traits<-c("T2D","GOUT","I9_CHD","C3_PROSTATE")

riskwithintervals <- results %>% filter(Group %in% c("Group1","Group6","Group11")) %>% filter(trait %in% traits)
pdf(file=paste0(output_dir,"Risks_Facet.pdf"),height=10,width=10,useDingbats=TRUE)
#colors<-c(brewer.pal(11,"RdYlBu")[11],"dark grey",brewer.pal(11,"RdYlBu")[1])
colors<-c("#824A85","#2C6687","#3C884B")
ggplot(riskwithintervals, aes(Age, LifetimeRisk, fill=Group, color=Group, group=Group)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 13), se = FALSE) +
  geom_point(alpha=0.5) +
  geom_ribbon(aes(ymin=CIneg, ymax=CIpos, fill=Group), alpha=0.2) +
  xlab("Age Range") + 
  ylab("Absolute Cumulative Risk (%)") + 
  theme_bw() + facet_wrap(~trait~biobank,ncol=4) +
  labs(color='PRS Group', fill='PRS Group') +
  scale_color_manual(values=colors,guide = guide_legend(reverse = TRUE),labels = c("0-1%", "40-60%", "99-100%")) +
  scale_fill_manual(values=colors,guide = guide_legend(reverse = TRUE),labels = c("0-1%", "40-60%", "99-100%")) +
  theme(title = element_text(size = 22),
        strip.background =element_rect(fill="white"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, angle=-90, hjust=0),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16))
 dev.off()
 
 
############################# comparisons across strata 
 
 hr_phenos <- c("C3_PROSTATE","C3_BREAST", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
 gbd_phenos <- c("Prostate cancer", "Breast cancer", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
 
 #for now remove total cancer because estbb has weird group1 and group 2 lifetime risk flip
 biobanks<-c("EstBB","HUNT","UKB","FinnGen","BBJ","MGB","MGB")
 
 ancestries<-c("EUR","SAS","EAS","AFR")
 
 rel_differences<-data.frame(NULL)
 differences<-data.frame(NULL)
 for(j in 1:length(hr_phenos)){
   for (b in 1:length(biobanks)){
     for (a in 1:length(ancestries)){
       #file<-paste0("/mnt/work/workbench/bwolford/intervene/results/",hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_",biobanks[b],".csv")
       file<-paste0("/home/bwolford/scratch/brooke/results/",hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_",biobanks[b],"_",ancestries[a],".csv")
       if (file.exists(file)){
         df<-fread(file)
         
         #Average PRS cumulative risk at age 80
         avg<-df %>% filter(Group %in% c("Group6") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
         avg$biobank<-biobanks[b]
         avg$pheno<-gbd_phenos[j]
         avg$comparison<-"40-60%"
         avg$ancestry<-ancestries[a]
         
         top<-df %>% filter(Group %in% c("Group11") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
         top$biobank<-biobanks[b]
         top$pheno<-gbd_phenos[j]
         top$comparison<-">99%"
         top$ancestry<-ancestries[a]
         
         bottom<-df %>% filter(Group %in% c("Group1") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
         bottom$biobank<-biobanks[b]
         bottom$pheno<-gbd_phenos[j]
         bottom$comparison<-"<1%"
         bottom$ancestry<-ancestries[a]
         
         if (nrow(bottom)==1 & nrow(top)==1){
           diff<-data.frame(biobank=biobanks[b], pheno=gbd_phenos[j],comparison="top vs bottom", ancestry=ancestries[a], 
                            CIneg=abs(top$CIneg-bottom$CIneg), CIpos=abs(top$CIpos-bottom$CIpos), LifetimeRisk=abs(top$LifetimeRisk-bottom$LifetimeRisk))
         }
         ####need to handle when Group10 is the top percentile, maybe take max of the Group 10 and Group 11 if both exist?
         
         rel_differences<-rbind(rel_differences,avg,top,bottom)
         differences<-rbind(differences,diff)
       }}}}
 
 rel_differences$comparison<-factor(rel_differences$comparison,levels=c("<1%","40-60%",">99%"))
 rel_differences<-rel_differences %>% filter(str_detect(ancestry,"EUR")|str_detect(ancestry,"EAS"))
 png(file="/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/cumulative_risk_differences.png",units="in",res=300,height=6,width=10)
 ggplot(rel_differences,aes(x=LifetimeRisk,y=pheno,color=biobank)) + geom_point(size=3,alpha=0.5) +
   theme_bw() + geom_errorbarh(aes(xmin=CIneg,xmax=CIpos),height=0.4,alpha=0.5) + facet_wrap(~comparison) +
   labs(x="Cumulative Risk at age 80",y='Trait') +
   theme(title = element_text(size = 22),
         legend.text = element_text(size = 16),
         legend.title = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.text.x = element_text(size = 12),
         axis.title.y = element_text(size = 18),
         axis.text.y = element_text(size = 16)) +
   scale_color_manual(values=carto_pal(n=length(biobanks), name="Safe"))
 dev.off()
 
 differences<-differences %>% filter(str_detect(ancestry,"EUR")|str_detect(ancestry,"EAS"))
 png(file="/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/cumulative_risk_top_vs_bottom.png",units="in",res=300,height=6,width=10)
 ggplot(differences,aes(x=LifetimeRisk,y=pheno,color=biobank)) + geom_point(size=3,alpha=0.5) +
   theme_bw() + geom_errorbarh(aes(xmin=CIneg,xmax=CIpos),height=0.4,alpha=0.5) + 
   labs(x="Difference in absolute cumulative risk between >99% and <1% ",y='Trait') +
   theme(title = element_text(size = 22),
         legend.text = element_text(size = 16),
         legend.title = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.text.x = element_text(size = 12),
         axis.title.y = element_text(size = 18),
         axis.text.y = element_text(size = 16)) +
   geom_vline(xintercept=0,linetype="dashed",color="black",alpha=0.8) +
   scale_color_manual(values=carto_pal(n=length(biobanks), name="Safe"))
 dev.off()
 
 
 
 
  