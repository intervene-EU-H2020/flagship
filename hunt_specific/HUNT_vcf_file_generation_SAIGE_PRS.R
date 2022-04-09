library(data.table)
library(dplyr)

#Read in PRS score file
PRS <- fread("path/to/prs/file", data.table=FALSE)

#Subset columns to the IDs and score only
PRS <- PRS[,c(2,5)]

#Rename columns
colnames(PRS) <- c("ID", "PRS")

#scale PRS
PRS$PRS <- scale(PRS$PRS)

#IGNORE lines 17 to 35 IF NOT GROUPING INTO PERCENTILES

#Assign into percentiles 
q <- quantiles(PRS$PRS, probs=c(0,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99,1))

PRS$PRS_Group <- cut(PRS$PRS, q, include.lowest=TRUE,
                     labels=paste("Group",1:11))

#Keep relevant columns
PRS <- PRS[,c("ID","PRS_Group")]

#For each group - assign a column which has a 1 if the individual belongs to this group. If the individual belongs to the reference group (group 6 - 40-60%) then assign a zero. 
for(i in c(1:5,7:11)){
  PRS[[paste0("Group",i)]] <- case_when(PRS[["PRS_Group"]] == paste0("Group ", i) ~ 1,
                                          PRS[["PRS_Group"]] == paste0("Group 6") ~ 0,
                                          TRUE ~ NA_real_)
}

#Remove grouping
PRS <- PRS[,-2]

#Assign rownames to ID so that transposing the dataframe assigns these as column names. Then remove ID column 
rownames(PRS) <- PRS$PRS
PRS <- PRS[,-1]

#Convert NAs to "." in line with vcf requirements
PRS[is.na(PRS)] <- "."

#Transpose dataframe so that samples are columns (required for vcf)
PRS_transpose <- t(PRS)
colnames(PRS_transpose) <- rownames(PRS)

#FOR PRS GROUPS
#Create fake columns for vcf file 
CHROM <- rep("chr1", 10)
POS <- seq(1,10)
ID <- c("rs1", "rs15", "rs510", "rs1020", "rs2040", "rs6080", "rs8090", "rs9095", "rs9599", "rs99")
REF <- rep("A", 10)
ALT <- rep("G", 10)
QUAL <- rep(".", 10)
FILTER <- rep(".", 10)
INFO <- rep(".", 10)
FORMAT <- rep("DS", 10)
vcf <- data.frame(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(vcf)[1] <- #CHROM
  
#FOR CONTINUOUS PRS
#Create fake columns for vcf file 
CHROM <- "chr1"
POS <- 1
ID <- "rs1"
REF <- "A"
ALT <- "G"
QUAL <- "."
FILTER <- "."
INFO <- "."
FORMAT <- "DS"
vcf <- data.frame(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(vcf)[1] <- #CHROM
  
#Join samples to required VCF columns
vcf <- cbind(vcf, PRS_transpose)
rownames(vcf) <- NULL

fwrite(vcf, "specify/file/path", sep="\t")

#Write meta data to be appended onto vcf file above
meta <- '##fileformat=VCFv4.2\n##fileDate=20171104\n##source=PLINKv1.90\n##contig=<ID=1,length=100001>\n##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">'

#Append to current vcf
con <- file('file/path/to/vcf','wt')

cat(paste0(meta, sep='\n'), file = con)
write.table(vcf,
            con,
            append = TRUE,
            sep = "\t",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

close(con)

#bgzip vcf file
bgzip -c your_vcf_file > your_vcf_file.gz

#create index file using gatk
gatk IndexFeatureFile -I your_vcf_file.gz
