#!/bin/bash -e
set -e

# 24 April 2023

# Olivia Janes
# Imputation stats
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    sppdir=~/data/tuturuatu_all_vcf/
    ## Edit to be run specific


    # Setting variables.
    beagledir=${sppdir}impute/beagle_imputations/
    mkdir -p ${beagledir}filtered/MAF_filtered_only/
    finaldir=${beagledir}filtered/MAF_filtered_only/


#Filter imputed files
    for vcf in ${beagledir}merged/*_imp_merged.vcf.gz
    do
        base=$(basename ${vcf} _imp_merged.vcf.gz)
        echo "Filtering SNPs for ${base}...." 
        vcftools --gzvcf ${vcf} \
            --out ${finaldir}${base}_0.05MAF.vcf \
            --maf 0.05 \
            --recode \
            --recode-INFO-all &
    done
    wait


# Removing '.recode.vcf' from filter file names.
    echo "Renaming filtered files to remove '.recode.vcf'"
    for vcf in ${finaldir}*recode.vcf
    do
        echo ""
        echo "Renaming ${vcf}"
        name=$(basename ${vcf} .recode.vcf)
        mv -i ${vcf} ${finaldir}${name}
    done 

# Convert filter files to vcf.gz format and index
    for vcf in ${finaldir}*vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${finaldir}${base}.vcf.gz --threads 16
        echo "Indexing ${vcf}.gz"
        bcftools index ${finaldir}${base}.vcf.gz --threads 16
        echo ""
    done

echo ""
echo "Filtering imputed vcfs is complete."