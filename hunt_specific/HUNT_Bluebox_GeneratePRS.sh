#!/bin/bash

#command read-freq used to stop plink calculating allele frequencies every time the PRS is calculated.
#if you do not have an allele frequency file, i suggest creating one and then referring to this as it will save you a bunch of time. 

#pheno=(Alcohol_Use_Disorder Alzheimers_Disease Asthma Atrial_Fibrillation BMI Breast_Cancer CHD Chronic_Kidney_Disease Educational_Attainment Epilepsy Focal_Epilepsy Generalised_Epilepsy Gout Heart_Failure Hip_Osteoarthritis IPF ILD Inflammatory_Bowel_Disease Knee_Osteoarthritis Lifespan Lung_Cancer MDD Melanoma Osteoporosis Pain POAG Prostate_Cancer Rheumatoid_Arthritis Sleep_Apnoea smoking Stroke Subarachnoid_Haemmorhage TAA T1D T2D Thyroid_Stimulating_Hormone)
pheno=$1
output=/home/bwolford/bluebox/results
score_directory=/home/bwolford/bluebox/data/PRS_v2
frequency_directory=/home/bwolford/bluebox/data
snplist_directory=/home/bwolford/bluebox/hunt_specific
genotype_directory=/home/bwolford/bluebox/data
plink_path=/home/bwolford/miniconda/bin/plink2

${plink_path} \
--bfile ${genotype_directory}/all.log \
--extract ${snplist_directory}/snplist_hg19_varid \
--read-freq ${frequency_directory}/all.frq \
--score ${score_directory}/${pheno}_megaPRS_scores_hg19_varid.txt.gz 6 2 5 header list-variants \
--out ${output}/${pheno}_PRS
done


#--score my.scores 3 2 1 header 
#reads variant IDs from column 3, allele codes from column 2, and scores from column 1 
