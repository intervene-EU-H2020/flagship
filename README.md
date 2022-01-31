# flagship
Code for the flagship project (WP2)

## ESHG Abstract Submission - PRS Association Analysis Across Biobanks - Very much a test run...

### Step 0: Requirements and Set-up

* These scripts assume you have plink and R (libraries: data.table, tidyverse/dplyr) installed on your biobanks system. 

* In each step, I have highlighted what parts of the scripts need to be adjusted.

* There are likely to be errors in the scripts. As and when they pop up, please let me know on slack and I can amend.

### Step 1: Download Adjusted Summary Statistics 

* All pre-adjusted summary statistics can be found [here](https://figshare.com/account/home#/projects/131369).

* Download the pre-adjusted summary statistics created using MegaPRS that correspond to the build of the genome for your biobank.

* hg19 files found [here](https://figshare.com/account/projects/131369/articles/19093304) / hg38 files found [here](https://figshare.com/account/projects/131369/articles/19093313). Hg19 contains rsids whereas hg38 contains variant IDs in the format CHR_POS_REF_ALT. 

### Step 2: Merge variant IDs to match those within the .bim file. 

#### Note - only to be performed if the genome build is hg38 and the variant ID structure is CHR_POS_REF_ALT.

* Run script hg38_biobankadjustments.R

* This script is required for all biobanks which have a genome build of hg38/where the variant ID structure within the bim file is in the form CHR_POS_REF_ALT. 

* Due to the nature of the conversion to hg38 within the pre-adjusted summary statistics, we cannot be sure what the REF or ALT allele is within your bim file as we have relied on those within the summary statistics themselves. As such, both combinations are currently saved within the sum stats. 

* The script reads in your bim file and identifies which variant ID is contained within your bim file and the summary statistics and creates a final snp list for the adjusted summary statistics. 

* For this script to work you will have to:
    1. Line 9 - Adjust the path within the 'bim' variable to read in your biobanks bim file.
    2. Line 19 - Adjust the path to the location of the adjusted summary statistics downloaded during step 1. Note: I have assumed you have mainted the same filename structure as when you downloaded the files. 
    3. Line 40 - Adjust the path to the location where you wish to save the amended summary statistics. Recommended to just overwrite the original summary statistics so keep the same path as specified in line 19. You may want to test the code before saving initially just to make sure it is behaving as expected. 

### Step 3: Compute PRS

* **Note: if you use a job scheduler which allows multiple jobs to be submitted, please adjust the script below and create a single job for each phenotype. This script assumes an interactive job and loops over phenotypes which will be much much slower than running in parallel.**

* Run script GeneratePRS.sh

* For this script to work you will have to:
    1. Line 6 - Specify the path where PRS are to be saved.
    2. Line 7 - Specify the path of the downloaded summary statistics from Step 1.
    3. Line 8 - Specify the path to the directory where the allele frequencys are stored.
    4. Line 9 - Specify the path to the directory where the genotypes are stored (it is likely line 8 and line 9 will have the same path). 
    5. Line 16 - Specify the path to plink. If you do not have plink2 installed, also change to the version used by your biobank.
    6. Line 17 - Change 'genotype_plink_files' to the name of your genotype files
    7. Line 18 - Change 'frequency_file' to the name of your file containing allele frequencies.
    8. Line 19 - Select 19 or 38 depending on the build of your biobanks genome. Also remove square brackets surrounding the number.
    
* *Note: If you do not have a file containing allele frequencies, we recommend producing one before computing PRS as otherwise plink will redo this step for every phenotype.*

* *Note: If the number of SNPs within your genotype file is very large and the run is taking a long time due to plink having to read in all the SNPs, consider using an argument that only uses Remos SNP list (hapmap + well-imputed 1000G). Remos SNP list can be found [here](https://github.com/intervene-EU-H2020/prspipe/tree/main/resources/1kg).*

### Step 4: Calculate associations between PRS and Phenotype - logistic regression

* **Note: this script has four assumptions.** 
    1) You have a phenotype file with case control assignments for each phenotype.
    2) You have been able to allocate genetic ancestry for the participants within your biobank. This does not have to be harmonized, any approach taken by your biobank will suffice for this analysis.
    3) You have kept the same short hand names for the phenotypes as within [FinnGen](https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519) (column B).
    4) You have kept the same naming structure for the PRS files as when you downloaded them.

* Run script OddsRatioCalculation.R

* For this script to work you will have to:
    1. Line 9 - Specify the path to the phenotype file. *See assumption 1.*
    2. Line 15 - Specify the path to the files containing the individual level PRS for each participant. 
    3. Line 21 - Change 'ENTER_ID' to correspond to the ID column within your biobank, i.e. FINNGENID. *Note: there may be redundant info taken by choosing the two columns. If this is the case, just select one of the IDs.*
    4. Line 28 - Specify the path to file which contains the genetic principal components calculated within participants of european ancestries.
    5. Line 31 - Amend if the column numbers selected do not automatically select the ID column and the first ten principal components. If necessary, also amend the column names of the file to be in the format 'PC1, PC2, ... PCN'.
    6. Line 38 - If this does not subset individuals to those of european ancestries, amend to do so with the files you have available to you.
    7. Line 88 - Change ENTER_BIOBANK_NAME to your biobank.

**Send results to bradley.jermy@helsinki.fi :)**

### UKBBPhenotyper.R
v1.0 under development :warning:  
Script by Bradley Jermy can work across biobanks as long as the ICD code individual level data is in long format.  
It extracts information from UKBB_definitions_demo_TEST.csv to extract all the endpoints for the flagship project.  
dependencies: data.table, tidyverse. 

Please check your ICD codes are in the correct format within the individual level data before applying the code. 

To make sure the script can separate between ICD 10 and ICD 9 codes, it places a '10x' or a '9x' at the start of the code according to whether the code is ICD9 or 10. 

The regex pattern then searches for strings starting with 10x or 9x before identifying the codes themselves. 

### biobank_summary_plot.R
Script by Brooke Wolford to analyze endpoint metrics collected from biobanks. 
Reads in data from google drive excel file. 
