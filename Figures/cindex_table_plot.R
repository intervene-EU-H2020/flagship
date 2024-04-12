library(data.table)
library(tidyverse)

setwd("~/Documents/MyStuff/NTNU/Research/INTERVENE/cindex")
fg<-fread("cindex_output_finngen.tab",header=TRUE)
fg$biobank<-"FinnGen"
mg<-fread("cindex_output_MGBB_EUR.tab",header=TRUE)
mg$biobank<-"Mass Gen Biobank"
ge<-fread("cindex_output_GE.tab",header=TRUE)
ge$biobank<-"Genomics England"
es<-fread("cindex_output_ESTBB.tab",header=TRUE)
es$biobank<-"Estonian Biobank"
gs<-fread('cindex_output_GS.tab',header=TRUE)
gs$biobank<-"Generation Scotland"
hnt<-fread("cindex_output_HUNT.tab",header=TRUE)
hnt$biobank<-"HUNT Study"

all<-all<-rbind(fg,mg,ge,es,gs,hnt) %>% drop_na()  %>% filter(pheno!="ILD") %>%
 filter(model=="null"|model=="PRS_percentile_group"|model=="PRS_percentile_group_and_birthyear"|model=="PRS_percentile_group_and_birthyear_and_sex")


all$model[all$model=="null"]<-"Baseline"
all$model[all$model=="PRS_percentile_group"]<-"PRS"
all$model[all$model=="PRS_percentile_group_and_birthyear"]<-"PRS and birth year"
all$model[all$model=="PRS_percentile_group_and_birthyear_and_sex"]<-"PRS, birth year, and sex"
all$pheno[all$pheno=="Knee_ARTHROSIS"]<-"KNEE_ARTHROSIS"


all$ub<-all$cindex+(1.96*all$cindex_SE)
all$lb<-all$cindex-(1.96*all$cindex_SE)

jpeg(file="cindex_facet_phenos.jpg",height=1500,width=1500,res=200)
ggplot(all,aes(x=cindex,y=pheno,color=model)) + geom_point(alpha=0.5) + facet_wrap(~biobank) + 
  theme_bw() + theme(legend.position="bottom") +
  geom_errorbarh(aes(xmin=all$lb,xmax=all$ub)) +
  labs(xlab="Harrell's C-statistic",y="Phenotype")
dev.off()

#phenotypes and the "best" model per each
#hr_phenos <- c("C3_BRONCHUS_LUNG","C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
#$best<-c()

