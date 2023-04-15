#!/bin/bash -e
set -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: samtools

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific


# Setting variables.
impdir=${sppdir}impute/
mkdir -p ${impdir}beagle_imputations/filtered
finaldir=${impdir}beagle_imputations/filtered/


#Filter imputed files
for vcf in ${impdir}beagle_imputations/*_imp.vcf.gz
do
    base=$(basename ${vcf} _beagle_imp.vcf.gz)
    echo "Filtering SNPs for ${base}...." 
    vcftools --gzvcf ${vcf} \
        --out ${finaldir}${base}_filtered.vcf \
        --max-missing 0.9 \
        --maf 0.05 \
        --recode \
        --recode-INFO-all &
    done

# Removing '.recode.vcf' from filter file names.
    echo "Renaming noLD files to remove '.recode.vcf'"
    for vcf in ${finaldir}*recode.vcf
    do
        echo ""
        echo "Renaming ${vcf}"
        base=$(basename ${vcf} .recode.vcf)
        mv -i ${vcf} ${finaldir}${base}
        echo "Indexing ${base.vcf}"
        bcftools index ${finaldir}${base}.vcf --threads 16
    done 

# Convert filter files to vcf.gz format and index
for vcf in ${finaldir}*vcf
do
	base=$(basename ${vcf} .vcf)
    echo "Converting ${vcf} to vcf.gz format"
	bcftools view ${vcf} -O z -o ${finaldir}${base}.vcf.gz --threads 16
    echo "Indexing ${vcf}"
    bcftools index ${finaldir}${base}.vcf.gz --threads 16
	echo ""
done

echo ""
echo "Filtering imputed vcfs is complete."