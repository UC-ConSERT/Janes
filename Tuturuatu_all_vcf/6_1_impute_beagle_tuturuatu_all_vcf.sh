#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Imputing with Beagle 5.2
# Imputing tuturuatu_all bcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call bcf that has been filtered for 
#   4x coverage, 0.2 site missingness, 10 minGQ and no LD.

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
mkdir -p ${impdir}beagle_imputations ${impdir}stats

# Imputation
    for file in ${finaldir}tuturuatu_{1..3}_study.vcf.gz
    do
        echo "Imputing ${file}"
        base=$(basename ${file} _study.vcf.gz)
        refvcf=${finaldir}${base}_ref.vcf.gz
        beagle gt=${file} ref=${refvcf} impute=true gp=true out=${impdir}beagle_imputations_ne100/${base}_beagle_imp ne=100
        bcftools index -f --threads 16 ${impdir}beagle_imputations/${base}_beagle_imp.vcf.gz
    done

# Stats
    # To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/tuturuatu_beagle_imp.vcf.gz

    for file in ${impdir}beagle_imputations/tuturuatu_{1..3}_beagle_imp.vcf.gz
    do
        echo "Calculating stats for ${file}"
        base=$(basename ${file} _beagle_imp.vcf.gz)
        # Concordance?
        vcf-compare ${finaldir}${base}_study.vcf.gz ${file} > ${impdir}stats/${base}_beagle_imp_concordance.txt

        # r2?
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' \
            ${file} > ${impdir}stats/${base}_beagle_imp_r2.txt
    done

echo ""
echo "Imputation script is complete."