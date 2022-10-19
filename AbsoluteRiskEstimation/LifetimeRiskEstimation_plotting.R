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
  geom_point() +
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
  