#Formatting information
## "" are not to be included - merely used to designate where changes are required.
## plinkFile for FinnGen is a single file. May need adapting if split into different chromosomes. 
## phenoFile is a gzipped file containing phenotype and covariate data in a tab delimited format. phenoFile should include all IDs from sampleFile even if they are all NAs. 
## sampleFile is a txt file with a single column of IDs - one line per ID. 
## GMMATmodelFile is the output from step1 with the rda prefix.
## varianceRatioFile is similar to above.
## vcfFile created for PRS using previously sent script. Same goes for vcfFileIndex

step1_fitNULLGLMM.R \
      --plinkFile="path/to/plink/file" \
      --phenoFile="path/to/pheno/file" \
      --phenoCol="name of column in phenoFile" \
      --covarColList="name of columns for first ten PCs, i.e. PC1,PC2,...,PC10" \
      --sampleIDColinphenoFile="name of ID column" \
      --traitType=binary \
      --outputPrefix="whatever you want your result to be outputted as" \
      --nThreads="we chose 32" \
      --LOCO=FALSE \
      --traceCVcutoff 0.0025 \
      --ratioCVcutoff 0.001 \
      --minCovariateCount 10 
      
step2_SPAtests.R \
      --minMAC=0 \
      --sampleFile="/path/to/samplefile" \
      --GMMATmodelFile="path/to/modelfile/filename.rda" \
      --varianceRatioFile="path/to/varianceRatioFile/filename.varianceRatio.txt" \
      --numLinesOutput=1000 \
      --IsOutputAFinCaseCtrl=TRUE \
      --IsDropMissingDosages=TRUE \
      --LOCO=FALSE \
      --analysisType=additive \
      --vcfFile="path/to/vcffile" \
      --vcfFileIndex="path/to/vcffileindex" \
      --vcfField=DS \
      --chrom=chr1 \
      --SAIGEOutputFile="whatever you want your result to be outputted as" \
