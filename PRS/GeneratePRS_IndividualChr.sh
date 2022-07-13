#!/bin/bash

#command read-freq used to stop plink calculating allele frequencies every time the PRS is calculated.
#if you do not have an allele frequency file, i suggest creating one and then referring to this as it will save you a bunch of time. 

pheno=(AllCancers Appendicitis Colorectal_Cancer Asthma Atrial_Fibrillation Breast_Cancer CHD Epilepsy Gout Hip_Osteoarthritis ILD Knee_Osteoarthritis Lung_Cancer MDD Melanoma Prostate_Cancer Rheumatoid_Arthritis T1D T2D)

output=/path/to/output
score_directory=/path/to/score/file
frequency_directory=/path/to/allele/frequency
snplist_directory=/path/to/snplist #Optional
genotype_directory=/path/to/genotype

for i in ${!pheno[@]}; do

pheno_i=${pheno[i]}

for j in {1..22}; do
path/to/plink2 \
--bfile ${genotype_directory}/genotype_plink_files \ #Include reference to chromosome number - i.e. ukb_imp_chr${j}_v3
--extract ${snplist_directory}/snplist_hg19 \ #Optional
--read-freq ${frequency_directory}/frequency_file \
--out ${output}/${pheno_i}_PRS \
--score ${score_directory}/${pheno_i}_megaPRS_scores_hg[19/38]_chr${j}.txt.gz 1 2 5 header list-variants

done

done
