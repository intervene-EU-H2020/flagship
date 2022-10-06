#!/bin/bash

#command read-freq used to stop plink calculating allele frequencies every time the PRS is calculated.
#if you do not have an allele frequency file, i suggest creating one and then referring to this as it will save you a bunch of time. 

#pheno=(AllCancers Appendicitis Colorectal_Cancer Asthma Atrial_Fibrillation Breast_Cancer CHD Epilepsy Gout Hip_Osteoarthritis ILD Knee_Osteoarthritis Lung_Cancer MDD Melanoma Prostate_Cancer Rheumatoid_Arthritis T1D T2D)

pheno=$1
output=/home/bwolford/bluebox/results
score_directory=/home/bwolford/bluebox/data/
frequency_directory=/home/bwolford/bluebox/data
snplist_directory=/home/bwolford/bluebox/hunt_specific
genotype_directory=/home/bwolford/bluebox/data
plink_path=/home/bwolford/miniconda/bin/plink2

${plink_path} \
--bfile ${genotype_directory}/all.log \
--extract ${snplist_directory}/snplist_hg19_varid \
--read-freq ${frequency_directory}/all.frq \
--score ${score_directory}/${pheno_i}_megaPRS_scores_hg[19/38].txt.gz 1 2 5 header cols=+scoresums list-variants \
--out ${output}/${pheno}_PRS



