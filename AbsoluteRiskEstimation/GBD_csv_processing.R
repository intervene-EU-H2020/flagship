#prepare GBD stats

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
gbd_bcpc <-c("Breast cancer","Prostate cancer")
  
path<-"/mnt/work/workbench/bwolford/hunt_flagship/GBD/"
output_dir<-"/mnt/work/workbench/bwolford/hunt_flagship/GBD/"

dat1<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-1.csv"))
dat2<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-2.csv"))
dat3<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-3.csv"))
dat4<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-4.csv"))
dat5<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-5.csv"))
dat6<-fread(paste0(path,"IHME-GBD_2019_DATA-da65912d-6.csv"))
dat_all<-rbind(dat1,dat2,dat3,dat4,dat5,dat6)

#### fill in the gaps
dat<-complete(dat_all,location,age,sex,cause,measure,fill=list(0))

####### Mortality


dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>% filter(str_detect(measure,"Deaths")) %>% filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Mortality.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Deaths")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(!str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Mortality.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Deaths")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_bcpc) %>% 
  filter(str_detect(cause,"Breast cancer") & str_detect(sex,"Female") | str_detect(cause,"Prostate cancer") & str_detect(sex,"Male")) %>%
select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Mortality.csv"))


######### Prevalence

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>% filter(str_detect(measure,"Prevalence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Prevalence.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Prevalence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(!str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Prevalence.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Prevalence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_bcpc) %>% 
  filter(str_detect(cause,"Breast cancer") & str_detect(sex,"Female") | str_detect(cause,"Prostate cancer") & str_detect(sex,"Male")) %>%
  select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Prevalence.csv"))

##### Incidence

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Incidence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Incidence.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Incidence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_phenos) %>% filter(!str_detect(sex,"Both")) %>% select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Incidence.csv"))

dat %>% filter(str_detect(metric,"Rate")|str_detect(metric,"Number")) %>%filter(str_detect(measure,"Incidence")) %>%filter(!is.na(val)) %>%
  filter(cause %in% gbd_bcpc) %>% 
  filter(str_detect(cause,"Breast cancer") & str_detect(sex,"Female") | str_detect(cause,"Prostate cancer") & str_detect(sex,"Male")) %>%
  select(measure,location,sex,age,cause,metric,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Incidence.csv"))
