library(data.table)

mapping <- fread("/path/to/mapping/file/1KGPhase3_hm3_hg19_hg38_mapping.tsv", data.table=FALSE)

mapping$chr <- paste('chr', mapping$chr, sep="")
mapping$Predictor_v1 <- paste(mapping$chr, mapping$pos_hg38, mapping$a1, mapping$a2, sep="_")
mapping$Predictor_v2 <- paste(mapping$chr, mapping$pos_hg38, mapping$a2, mapping$a1, sep="_")

snplist <- mapping[,c("Predictor_v1", "Predictor_v2")]

#Read in bim file
#read bim file
bim <- fread("/path/to/bim/filename.bim", data.table=FALSE)

#give the file column names: I have assumed standard bim format.
colnames(bim) <- c("chr","variant_id","cM","pos","a1","a2") 

#identify which of the two variant_ids are found in your bim file
first_id <- subset(snplist, Predictor_v1 %in% bim$variant_id)
first_id <- first_id$Predictor_v1
second_id <- subset(snplist, Predictor_v2 %in% bim$variant_id)
second_id <- second_id$Predictor_v2
  
snplist$Predictor <- case_when(snplist$Predictor_v1 %in% first_id ~ snplist$Predictor_v1,
                               snplist$Predictor_v2 %in% second_id ~ snplist$Predictor_v2,
                               TRUE ~ NA_character_)
  
snplist <- snplist$Predictor
  
#Save adjusted snplist file so that it can be read by plink
write.table(snplist, "/path/to/snplist/snplist_file", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

