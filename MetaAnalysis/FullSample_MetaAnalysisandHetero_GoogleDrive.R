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

#Read in full sample hazard ratios per standard deviation

run_googledrive<-FALSE

if (run_googledrive==TRUE){
  #identify folder
  folder_id = drive_get(as_id("1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ"))
  
  #find files in folder
  files = drive_ls(folder_id)
  
  #loop dirs and download files inside them
  for (i in seq_along(files$name)) {
    #list files
    i_dir = drive_ls(files[i, ])
    
    #mkdir
    try({dir.create(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/",files$name[i]))})
    
    #download files
    for (file_i in seq_along(i_dir$name)) {
      #fails if already exists
      try({
        drive_download(
          as_id(i_dir$id[file_i]),
          path = paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/",files$name[i], "/",i_dir$name[file_i])
        )
      })
    }
  }
  
  #estb has weird file structure so do it separately
  dir<-drive_ls(as_id("1CySw57_ICg0e-DW3M4wiuRTccmp8QYAn"))
  dir.create("/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/")
  for (idx in seq_along(dir$name)) {
    #fails if already exists
    try({
      drive_download(
        as_id(dir$id[idx]),
        path = paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/",dir$name[idx])
      )
    })
  }
}
### full sample, HR per std dev
file<-"HRperSD"
estbb<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/",file,"_EstBB.csv"))
finngen<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/FinnGen_HazardRatios/",file,"_FinnGen.csv"))
hunt<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/HUNT_HazardRatios/",file,"_HUNT.csv"))
bbj<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/Biobank Japan_HazardRatios/",file,"_BBJ.csv"))
gs<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/GenerationScotland_HazardRatios/",file,"_GS.csv"))
gnh<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/Genes&Health_HazardRatios/",file,"_GNH.csv"))
ge<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/GenomicsEngland_HazardRatios/",file,"_GenomicsEngland.csv"))
mgb_eur<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/MGB_HazardRatios/",file,"_MGBB_EUR.csv"))
mgb_afr<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/MGB_HazardRatios/",file,"_MGBB_AFR.csv"))

#UKB
ukb$Biobank <- "UK Biobank"
ukb<-fread(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/UKB_HazardRatios/PRS_HRsperSD_UKBiobank_AllAncestries.csv")) #all ancestries

#HUNT 
drophunt <-c("T1D","C3_CANCER") #these are weird on first pass 
hunt<-subset(hunt,!(Phenotype %in% drophunt))
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
bbj <- subset(bbj, !(Phenotype %in% dropbbj))
bbj$Biobank <- "Biobank Japan"
bbj$Ancestry <- "EAS"

#Estonian Biobank
estbb$Biobank <- "Estonian Biobank"
estbb$Ancestry <- "EUR"

#Generation Scotland
dropgs <- c("C3_COLORECTAL","I9_AF","GOUT")
gs <- subset(gs, !(Phenotype %in% dropgs))
gs$Biobank <- "Generation Scotland"
gs$Ancestry <- "EUR"

#Mass General Brigham
dropmgb <- c("I9_AF")
mgb_eur <- subset(mgb_eur, !(Phenotype %in% dropmgb))
mgb_afr <-subset(mgb_afr, !(Phenotype %in% dropmgb))
mgb_eur$Biobank <- "Mass General Brigham"
mgb_eur$Ancestry <- "EUR"
mgb_afr$Ancestry <- "AFR"
mgb_afr$Biobank <-  "Mass General Brigham"

#Genomics England
ge$Biobank <- "Genomics England"
ge$Ancestry <- "EUR"
names(ge)<-names(hunt)
ge<-ge %>% mutate(Phenotype=case_when(Phenotype=="Knee_ARTHROSIS"~"KNEE_ARTHROSIS")) #fix phenotype typo

#Combine into one dataset
all <- rbind(ukb, finngen) %>%
          rbind(., gnh) %>% 
            rbind(., estbb) %>%
              rbind(., gs) %>%
                rbind(., mgb_eur) %>% 
                  rbind(., mgb_afr) %>% 
                  rbind(., bbj) %>%
                    rbind(.,ge) %>% rbind(., hunt)

full <- subset(all, Sample=="Full Sample")

metaresults <- c()
for(j in c("EUR","SAS","EAS","AFR")){

  for(i in unique(full$Phenotype)){
    print(i)
    #Test with T2D phenotype and EUR ancestry
    if(j=="EAS"){
      disease <- subset(full, Phenotype==i & Phenotype %in% bbj$Phenotype & Ancestry==j & !(is.na(Beta)))
    } else{
      disease <- subset(full, Phenotype==i & Ancestry==j & !(is.na(Beta)))
    }
    if(dim(disease)[1]==0){
      next
    }

    #Meta analysis should be done at the beta level and stratified by ancestry 
    meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
    
    #make new plots: https://www.metafor-project.org/doku.php/plots
    png(file=paste0(output_dir,i,"_",j,"_Forest.png"), width = 7,
        height    = 7,
        units     = "in",
        res       = 300)
    forest(meta, slab=disease$Biobank)
    dev.off()
  
    #Save the forest plot so can compare with the table
    metaresults <- rbind(metaresults, c(j, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
  }
}


metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Ancestry","Phenotype","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(3:7),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))
fwrite(metaresults, paste0(output_dir,"metaanalysisFE.csv"))

#forest plots v2 
#wide<-full %>% mutate(Ancestry2=Ancestry) %>% pivot_wider(names_from=c(Biobank,Ancestry2),values_from=c(HR,Cipos,Cineg,Controls,Cases))
  
df<-metaresults %>% mutate(Biobank="Fixed-effects meta-analysis") %>% rbind.fill(full) %>% mutate(shape=as.factor(case_when(Biobank=="Fixed-effects meta-analysis"~1, Biobank!="Fixed-effects meta-analysis"~0)))
df$Biobank<-as.factor(df$Biobank)
df$Biobank<-relevel(df$Biobank,"Fixed-effects meta-analysis") #puts meta analysis point at bottom
df$pretty<-gsub("_", " ", df$PRS, fixed=TRUE) #make pretty phenotype label
#custom forest plot
by(df,df$Phenotype, function(x){
  by(x,x$Ancestry, function(y){
    name=paste(sep="_",unique(y$Phenotype),unique(y$Ancestry))
    pheno=unique(na.omit(y$pretty))
    png(file=paste0(output_dir,name,"_Forest.png"),height=3,width=4.5,units="in",res=300)
    print(ggplot(y,aes(x=HR,y=Biobank,shape=shape,color=shape)) + geom_point(aes(size=y$shape)) +theme_bw() +
            geom_vline(xintercept=1,color="red",linetype="dashed") + 
            scale_x_continuous(labels = label_number(accuracy = 0.01)) +
            geom_errorbarh(aes(xmin=Cineg,xmax=Cipos),height=0.25) +
            scale_shape_manual(values=c(19,18)) +  
            scale_color_manual(values=c("#000000","#3274B4")) +
            scale_size_manual(values=c(2,4)) +
            theme(title = element_text(size = 18),
                  legend.position="none",plot.title=element_text(hjust=0.5),
                  legend.text = element_text(size = 14),
                  legend.title = element_blank(),
                  axis.title.x = element_text(size = 18),
                  axis.text.x = element_text(size = 12),
                  axis.title.y = element_text(size = 18),
                  axis.text.y = element_text(size = 12)) +
          labs(title=pheno,x="Hazard Ratio (95% CI)",y="Study"))
    dev.off()
  })
})
#could do facet by ancestry? 

#visually compare the heterogeniety pvalues
ggplot(df[df$Ancestry=="EUR",],aes(x=HetPval,y=Phenotype)) + geom_point() + scale_x_log10()

#compare hunt to meta-analysis
compare<-metaresults %>% left_join(hunt,by=c("Ancestry"="Ancestry","Phenotype"="Phenotype")) %>% filter(Sample=="Full Sample")
ggplot(compare,aes(x=HR.x,y=HR.y,label=PRS)) +geom_point() + geom_label() + geom_abline()
 
######################################################################################################################################
#Combine the meta analysed results with each biobank and then compare - reduce to those with differences and save file only
metaresults$Biobank <- "All"
metadiff <- metaresults[,c("Ancestry","Phenotype","Beta","SE","Pval","Biobank")]
##fulldiff <- full[,names(metadiff)] ###what is this line doing?
#diffs <- rbind(fulldiff, metadiff)
diffs<-full %>% select(names(full)[names(full) %in% names(metadiff)]) %>% rbind(metadiff)

diff <- c()
for(i in c("EUR","SAS","EAS","AFR")){
  for(j in unique(diffs$Phenotype)){
    disease <- subset(diffs, Phenotype==j & Ancestry==i)
    disease$BetaDiff <- ifelse(!(is.na(disease$Beta)), disease$Beta - disease$Beta[disease$Biobank=="All"], NA)
    disease$SEDiff <- ifelse(!(is.na(disease$Beta)), sqrt(disease$SE**2 + disease$SE[disease$Biobank=="All"]**2), NA)
    disease$ZDiff <- ifelse(!(is.na(disease$Beta)), disease$BetaDiff/disease$SEDiff, NA)
    disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
    diff <- rbind(diff, disease)
  }
}

#what phenotypes have significant differences?
#115 independent european tests
unique(diffs[diffs$PvalDiff<0.05/115 & diffs$Ancestry=="EUR",]$Phenotype)

diffs <- as.data.frame(diff)
fwrite(diffs, paste0(output_dir,"FullSampleIndBiobankComparisonwMAEffectSize.csv"))

######################################################################################################################################
# Main analysis
meta <- fread(paste0(output_dir,"metaanalysisFE.csv"), data.table=FALSE)
meta$Ancestry <- factor(meta$Ancestry, levels=c("AFR","SAS","EAS","EUR"))
orderPheno <- subset(meta, Ancestry=="EUR")
meta$Phenotype <- factor(meta$Phenotype, levels=c(meta[order(orderPheno$HR),"Phenotype"]))
colors<-brewer.pal(n=4,name="Dark2")

#Main Plot - EUR, SAS, EAS, and AFR
figure2a <- ggplot(meta) +
              geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
              theme_bw() +
              geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
              ylab("Hazard Ratio (95% CI)") +
              xlab("") +
              geom_hline(yintercept = 1.0) +
              scale_color_manual(values=colors)+
              scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
              theme(title = element_text(size = 18),
                    legend.text = element_text(size = 14),
                    legend.title = element_blank(),
                    axis.title.x = element_text(size = 18),
                    axis.text.x = element_text(size = 14),
                    axis.title.y = element_text(size = 18),
                    axis.text.y = element_text(size = 14)) + ylim(0.5,2.5) +
              coord_flip() 
ggsave(paste0(output_dir,"all_ancestry_HR_meta.png"), height=6 , width=8)

### european only
figure2a <- ggplot(meta[meta$Ancestry=="EUR",]) +
  geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=colors[4])+
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) + ylim(0.5,2.5) + 
  coord_flip() 
ggsave(paste0(output_dir,"EUR_HR_meta.png"), height=6 , width=8)

##all except african
figure2a <- ggplot(meta[meta$Ancestry!="AFR",]) +
  geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=colors[2:4]) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) + ylim(0.5,2.5) +
  coord_flip()
ggsave(paste0(output_dir,"EUR_SAS_EAS_HR_meta.png"), height=6 , width=8)



######################################################################################################################################

ukbaframr <- subset(ukb, (Ancestry=="AFR" | Ancestry=="AMR") & Sample=="Full Sample" & !(is.na(HR)))
afronly <- subset(ukbaframr, Ancestry=="AFR")
ukbaframr$Phenotype <- factor(ukbaframr$Phenotype, afronly[order(afronly$HR),"Phenotype"])

#Supplementary Plot - AFR and AMR
ggplot(ukbaframr) +
  geom_point(aes(Phenotype, HR, group=Ancestry, col=Ancestry), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=Ancestry, col=Ancestry), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#88CCEE", "#CC6677")) +
  scale_x_discrete(labels=c("ILD","Appendicitis","Epilepsy","CHD","All Cancers","Lung Cancer","Type 1 Diabetes","Major Depression","Knee Osteoarthritis","Rheumatoid Arthritis","Prostate Cancer","Hip Osteoarthritis","Colorectal Cancer","Atrial Fibrillation","Type 2 Diabetes","Breast Cancer","Gout","Asthma")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigureAFR&AMRAncestriesHRperSDFullSample.png", height=8 , width=8)

######################################################################################################################################

# Leave one out meta-analysis - Only really relevant for European meta-analysis - Each mark can show the effect on the effect size
metaresults <- c()
european <- subset(full, Ancestry=="EUR")

for(i in unique(european$Biobank)){
  
  loo <- subset(european, Biobank!=i)
  
  #Repeat meta-analysis and save down the results but allude to the biobank

  #European ancestry meta-analysis
    
  for(j in unique(loo$Phenotype)){
    print(j)
    
    #Check to see if the biobank has been tested
    if(dim(subset(european, Phenotype==j & Biobank==i))[1]==0){
      next
    }
    
    #Test with T2D phenotype and EUR ancestry
    disease <- subset(loo, Phenotype==j & !(is.na(Beta)))
    
    if(dim(disease)[1]==0){
      next
    }
    
    #Meta analysis should be done at the beta level and stratified by ancestry 
    meta <- rma(yi=Beta, sei=SE, data=disease, method="FE")
    
    forest(meta, slab=disease$Biobank)
    grid.text(paste0(i," - EUR Ancestry"), .5, .9, gp=gpar(cex=2))
    
    metaresults <- rbind(metaresults, c(j, i, meta$b, meta$se, meta$pval, meta$QE, meta$QEp))
  
  }
  
}

metaresults <- as.data.frame(metaresults)
colnames(metaresults) <- c("Phenotype","LeftOutBiobank","Beta","SE","Pval","QHet","HetPval")
metaresults <- metaresults %>% mutate_at(c(3:7),as.numeric)
metaresults$HR <- exp(metaresults$Beta)
metaresults$Cipos <- exp(metaresults$Beta + (1.96*metaresults$SE))
metaresults$Cineg <- exp(metaresults$Beta - (1.96*metaresults$SE))
  
fwrite(metaresults, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv")

loometa <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv", data.table=FALSE)
allmeta <- fread("/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/metaanalysisFE.csv", data.table=FALSE)
allmeta <- subset(allmeta, Ancestry=="EUR")
allmeta$LeftOutBiobank <- "All"
allmeta <- allmeta[,c(names(loometa))]
meta <- rbind(loometa, allmeta)
meta$LeftOutBiobank <- factor(meta$LeftOutBiobank, levels=c("All","FinnGen","Estonian Biobank","UK Biobank","Mass General Brigham","Generation Scotland"))
meta$Phenotype <- factor(meta$Phenotype, levels=c(allmeta[order(allmeta$HR),"Phenotype"]))

metaupdate <- c()
#Test to see if the individual biobank association is significantly different to that of the meta-analysed association
for(i in unique(meta$Phenotype)){
  disease <- subset(meta, Phenotype==i)
  disease$BetaDiff <- disease$Beta - disease$Beta[disease$LeftOutBiobank=="All"]
  disease$SEDiff <- sqrt(disease$SE**2 + disease$SE[disease$LeftOutBiobank=="All"]**2)
  disease$ZDiff <- disease$BetaDiff/disease$SEDiff
  disease$PvalDiff <- 2*pnorm(abs(disease$ZDiff), lower.tail=FALSE)
  metaupdate <- rbind(metaupdate, disease)
}

meta <- as.data.frame(metaupdate)

#Plot leave one out results
ggplot(meta) +
  geom_point(aes(Phenotype, HR, group=LeftOutBiobank, col=LeftOutBiobank), position=position_dodge(width=0.7)) +
  theme_bw() +
  geom_errorbar(aes(x=Phenotype, ymin = Cineg, ymax = Cipos, group=LeftOutBiobank, col=LeftOutBiobank), size=0.5, width=0.5, position=position_dodge(width=0.7)) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  geom_hline(yintercept = 1.0) +
  scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")) +
  scale_x_discrete(labels=c("Appendicitis","Epilepsy","All Cancers","Lung Cancer","Major Depression","ILD","Knee Osteoarthritis","Skin Melanoma","Hip Osteoarthritis","CHD","Asthma","Colorectal Cancer","Atrial Fibrillation","Breast Cancer","Rheumatoid Arthritis","Gout","Type 2 Diabetes","Prostate Cancer","Type 1 Diabetes")) + 
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  coord_flip()
ggsave("/Users/jermy/Documents/INTERVENE/Write-up/SupplementaryFigure_LOOMetaAnalysisEUR.png", height=8 , width=10)

################################################################################################################################################################################################################################################################################################################################
#Drop the All results and resave the file
meta <- subset(meta, Phenotype!="All")
fwrite(meta, "/Users/jermy/Documents/INTERVENE/Results/MetaAnalysis/LOOmetaanalysisFE.csv")

################################################################################################################################################################################################################################################################################################################################







