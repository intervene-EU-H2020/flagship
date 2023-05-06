#To investigate - Prostate Cancer and the Epilepsys and maybe heart failure
library(data.table)
library(dplyr)

phenotypes <- c("smk","scz","sbp","risk","pain","bmi","drk","ea","externalizing","ldl")

for(i in phenotypes){
# Calculate Per-Predictor Heritabilities
system(paste0("/Users/jermy/Software/megaPRS/ldak5.1.mac --sum-hers /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak --tagfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak.tagging --summary /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/",i,"_megaPRS.tsv --matrix /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak.matrix --check-sums NO"))

# Create pseudo summaries using the reference panel
system(paste0("/Users/jermy/Software/megaPRS/ldak5.1.mac --pseudo-summaries /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_sumstats.pseudo --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --summary /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/",i,"_megaPRS.tsv --training-proportion .9 --keep /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepa --allow-ambiguous YES --extract /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/",i,"_megaPRS.tsv"))

#################################################################
# Estimate effect sizes for training and full prediction models.
#################################################################

#Training
system(paste0("/Users/jermy/Software/megaPRS/ldak5.1.mac --mega-prs /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full --model mega --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --cors /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/cors_full --ind-hers /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak.ind.hers --summary /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_sumstats.pseudo.train.summaries --one-sums YES --window-cm 1 --allow-ambiguous YES --extract /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/",i,"_megaPRS.tsv"))

#Test 
system(paste0("/Users/jermy/Software/megaPRS/ldak5.1.mac --calc-scores /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --scorefile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full.effects --summary /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_sumstats.pseudo.test.summaries --power 0 --final-effects /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full.effects --keep /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepc --allow-ambiguous YES --exclude /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/highld/genes.predictors.used")) 

#Work out best model and extract parameters

system(paste0("mv /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full.effects.best /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_megaPRS_scores.txt"))

bestmodel <- fread(input=paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_megaPRS_scores.txt"), data.table=FALSE) 
bestmodel$Effect_Best <- signif(bestmodel$Effect_Best, digits=2)

allmodels <- fread(input=paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full.effects"), data.table=FALSE) 

#Put each result into the same format
allmodels[,c(5:240)] <- signif(allmodels[,c(5:240)], digits=2)

#Find which model has been taken to be the best model
identitycheck <- apply(allmodels[,c(5:240)], 2, function(x) sum(x==bestmodel$Effect_Best))
print(max(identitycheck))
index <- match(max(identitycheck),identitycheck)

#Read in parameter file
parameters <- fread(input=paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_mega_full.parameters"), data.table=FALSE)

#Select best model
bestparams <- subset(parameters, Model==index)
fwrite(bestparams, paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_bestmodel_parameters"), sep="\t")
  
#Move important files to relevant location and remove unecessary files
system(paste0("mv /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_bestmodel_parameters /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/ModelParameters/",i,"_bestmodel_parameters"))
system(paste0("mv /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"_megaPRS_scores.txt /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/hg19_sumstats/",i,"_megaPRS_scores.txt"))

system(paste0("rm /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/",i,"*"))

#Adjust summary statistics to be suitable for hg38
#Convert files to Chr:BP for hg38
scores <- fread(input=paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/hg19_sumstats/",i,"_megaPRS_scores.txt"), data.table=FALSE)

hg38conv <- fread(input=paste0("/Users/jermy/Documents/INTERVENE/Sumstats/1KGPhase3_hm3_hg19_hg38_mapping.tsv"), data.table=FALSE)
hg38conv$Predictor <- paste(hg38conv$chr, hg38conv$pos_hg19, sep=":")
hg38conv$chr <- paste("chr", hg38conv$chr, sep="")
hg38conv$variant_id_hg38_v1 <- paste(hg38conv$chr, hg38conv$pos_hg38, hg38conv$a1, hg38conv$a2, sep="_")
hg38conv$variant_id_hg38_v2 <- paste(hg38conv$chr, hg38conv$pos_hg38, hg38conv$a2, hg38conv$a1, sep="_")

hg38conv <- hg38conv[,c("Predictor","variant_id_hg38_v1","variant_id_hg38_v2","a1","a2")]
colnames(hg38conv) <- c("Predictor","variant_id_hg38_v1","variant_id_hg38_v2","A1","A2")

scores$IUPAC[scores$A1 == 'A' & scores$A2 =='T' | scores$A1 == 'T' & scores$A2 =='A'] <- 'W'
scores$IUPAC[scores$A1 == 'C' & scores$A2 =='G' | scores$A1 == 'G' & scores$A2 =='C'] <- 'S'
scores$IUPAC[scores$A1 == 'A' & scores$A2 =='G' | scores$A1 == 'G' & scores$A2 =='A'] <- 'R'
scores$IUPAC[scores$A1 == 'C' & scores$A2 =='T' | scores$A1 == 'T' & scores$A2 =='C'] <- 'Y'
scores$IUPAC[scores$A1 == 'G' & scores$A2 =='T' | scores$A1 == 'T' & scores$A2 =='G'] <- 'K'
scores$IUPAC[scores$A1 == 'A' & scores$A2 =='C' | scores$A1 == 'C' & scores$A2 =='A'] <- 'M'

scores <- scores[(scores$IUPAC %in% c('R', 'Y', 'K', 'M')),]

hg38conv$IUPAC[hg38conv$A1 == 'A' & hg38conv$A2 =='T' | hg38conv$A1 == 'T' & hg38conv$A2 =='A'] <- 'W'
hg38conv$IUPAC[hg38conv$A1 == 'C' & hg38conv$A2 =='G' | hg38conv$A1 == 'G' & hg38conv$A2 =='C'] <- 'S'
hg38conv$IUPAC[hg38conv$A1 == 'A' & hg38conv$A2 =='G' | hg38conv$A1 == 'G' & hg38conv$A2 =='A'] <- 'R'
hg38conv$IUPAC[hg38conv$A1 == 'C' & hg38conv$A2 =='T' | hg38conv$A1 == 'T' & hg38conv$A2 =='C'] <- 'Y'
hg38conv$IUPAC[hg38conv$A1 == 'G' & hg38conv$A2 =='T' | hg38conv$A1 == 'T' & hg38conv$A2 =='G'] <- 'K'
hg38conv$IUPAC[hg38conv$A1 == 'A' & hg38conv$A2 =='C' | hg38conv$A1 == 'C' & hg38conv$A2 =='A'] <- 'M'

hg38conv <- hg38conv[(hg38conv$IUPAC %in% c('R', 'Y', 'K', 'M')),]
hg38conv <- hg38conv[,c(1,2,3,6)]

scores <- inner_join(scores, hg38conv)

scores <- scores[,c(7,8,2:5)]
colnames(scores)[c(1,2)] <- c("Predictor_v1","Predictor_v2")

fwrite(scores, paste0("/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/hg38_sumstats_INTERVENE/",i,"_megaPRS_scores_hg38.txt"), sep="\t")
}

