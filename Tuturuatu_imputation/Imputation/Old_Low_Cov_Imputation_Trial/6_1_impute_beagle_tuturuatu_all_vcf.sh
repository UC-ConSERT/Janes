#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Imputing with Beagle 5.4 (beagle.22Jul22.46e.jar)
# Imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals, and some moderate coverage wild indv.
# Imputation will be on the variant call vcf that has been filtered for 
#   0x,4x,5x coverage, no site missingness, bi-allelic sites only, max depth 50, minGQ 10, minQ 20, SP<60.

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
mkdir -p ${impdir}beagle_imputations ${impdir}stats
mkdir -p ${impdir}beagle_imputations/impute_trials
impoutdir=${impdir}beagle_imputations/impute_trials/


# Imputation
    for file in ${finaldir}*_study.vcf.gz
    do
        base=$(basename ${file} _study.vcf.gz)
        echo "Imputing ${base}"
        refvcf=${finaldir}${base}_ref_phased.vcf.gz

        for test_ne in {50,100,500}
        do
            java -jar ${beaglejar} gt=${file} impute=true gp=true ne=${test_ne} em=false nthreads=16 \
                ref=${refvcf} out=${impoutdir}${base}_${test_ne}ne_beagle_imp
            echo ""; echo "Indexing ${base}"
            bcftools index -f --threads 16 ${impoutdir}${base}_${test_ne}ne_beagle_imp.vcf.gz
        done
    done

echo ""
echo "Imputation script is complete."