# flagship
Code for the flagship project (WP2)

### UKBBPhenotyper.R
v1.0 under development :warning:  
Script by Bradley Jermy can work across biobanks as long as the ICD code individual level data is in long format.  
It extracts information from UKBB_definitions_demo_TEST.csv to extract all the endpoints for the flagship project.  
dependencies: data.table, tidyverse. 

Please check your ICD codes are in the correct format within the individual level data before applying the code. 

To make sure the script can separate between ICD 10 and ICD 9 codes, it places a '10x' or a '9x' at the start of the code according to whether the code is ICD9 or 10. 

The regex pattern then searches for strings starting with 10x or 9x before identifying the codes themselves. 

