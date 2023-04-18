#!/bin/bash -e
set -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: samtools

sppdir=~/data/tuturuatu_all_vcf/impute/truth/
## Edit to be run specific


# Setting variables.
beagledir=${sppdir}impute/beagle_imputations/
mkdir -p ${beagledir}filtered/
finaldir=${beagledir}filtered/


#Filter imputed files
    for vcf in ${beagledir}merged/*_imp_merged.vcf.gz
    do
        base=$(basename ${vcf} _imp_merged.vcf.gz)
        echo "Filtering SNPs for ${base}...." 
        vcftools --gzvcf ${vcf} \
            --out ${finaldir}${base}_filtered.vcf \
            --max-missing 0.9 \
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