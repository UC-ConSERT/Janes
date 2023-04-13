#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
statsdir=${impdir}stats


# Stats
    # To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/tuturuatu_beagle_imp.vcf.gz

    for file in ${finaldir}*_study.vcf.gz
    do
        for test_ne in {50,100,500}
        do
            base=$(basename ${file} _study.vcf.gz)
            echo ""; echo "Calculating stats for ${base}_${test_ne}ne_beagle_imp.vcf.gz"

            # Concordance
            vcf-compare ${file} ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_concordance.txt

            # Allelic/Dosage r^2
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' \
                ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_r2.txt

        done
    done

echo ""
echo "Imputation stats script is complete."