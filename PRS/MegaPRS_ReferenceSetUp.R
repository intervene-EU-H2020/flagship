############################################################################################################
#Per predictor heritabilities 
############################################################################################################

# Calculate LDAK weights
/Users/jermy/Software/megaPRS/ldak5.1.mac --cut-weights /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/sections --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref
/Users/jermy/Software/megaPRS/ldak5.1.mac --calc-weights-all /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/sections --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref
mkdir /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/bld
mv /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/sections/weights.short /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld/bld65

# Calculate taggings
/Users/jermy/Software/megaPRS/ldak5.1.mac --calc-tagging /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld/bld --window-cm 1 --save-matrix YES

# Calculate Per-Predictor Heritabilities - First stage to use T2D summary statistics
/Users/jermy/Software/megaPRS/ldak5.1.mac --sum-hers /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak --tagfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak.tagging --summary /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/T2D_megaPRS.tsv --matrix /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/bld.ldak.matrix --check-sums NO

# Identify SNPs in high LD regions
/Users/jermy/Software/megaPRS/ldak5.1.mac --cut-genes /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/highld --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --genefile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/highld.txt

#################################################################
#Split reference into three
#################################################################

awk < /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref.fam '(NR%3==1){print $0 > "/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepa"}(NR%3==2){print $0 > "/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepb"}(NR%3==0){print $0 > "/Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepc"}'

# Create pseudo summaries using the reference panel
/Users/jermy/Software/megaPRS/ldak5.1.mac --pseudo-summaries /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/T2D_sumstats.pseudo --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --summary /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/T2D_megaPRS.tsv --training-proportion .9 --keep /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepa --allow-ambiguous YES --extract /Users/jermy/Documents/INTERVENE/Sumstats/MegaPRS/T2D_megaPRS.tsv 

#Calculate predictor-predictor correlations
/Users/jermy/Software/megaPRS/ldak5.1.mac --calc-cors /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/cors_full --bfile /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/ref --window-cm 3 --keep /Users/jermy/Documents/INTERVENE/Results/FinnGen/MegaPRS/WorkingFiles/keepb 
