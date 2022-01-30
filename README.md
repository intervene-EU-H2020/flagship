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
#### Note - only to be performed if the genome build is hg38 and the variant ID structure is CHR_POS_REF_ALT.####
* This script is required for all biobanks which have a genome build of hg38/where the variant ID structure within the bim file is in the form CHR_POS_REF_ALT. 
* Due to the nature of the conversion to hg38 within the pre-adjusted summary statistics, we cannot be sure what the REF or ALT allele is within your bim file as we have relied on those within the summary statistics themselves. As such, both combinations are currently saved within the sum stats. 
* The script reads in your bim file and identifies which variant ID is contained within your bim file and the summary statistics and creates a final snp list for the adjusted summary statistics. 
* Run script hg38_biobankadjustments.R
* For this script to work you will have to:
    1. Line 9 - Adjust the path within the 'bim' variable to read in your biobanks bim file.
    2. Line 19 - Adjust the path to the location of the adjusted summary statistics downloaded during step 1. Note: I have assumed you have mainted the same filename structure as when you downloaded the files. 
    3. Line 40 - Adjust the path to the location where you wish to save the amended summary statistics. Recommended to just overwrite the original summary statistics so keep the same path as specified in line 19. You may want to test the code before saving initially just to make sure it is behaving as expected. 

### Step 3: Compute PRS
* **note: if you have a system similar to a high-performance cluster which can submit multiple jobs in parallel, please adjust the below script and create multiple jobs for each phenotype. This script assumes an interactive job is in play and simply loops over phenotypes which will be much slower than the method highlighted above**
* Run script GeneratePRS.sh
* For this script to work you will have to:
    1. Line 6 - Specify the path where PRS are to be saved.
    2. Line 7 - Specify the path ...
    3. Line 8 - Specify the path ...
    4. Line 9 - Specify the path ... 

### Step 4: Calculate associations between PRS and Phenotype - logistic regression

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
