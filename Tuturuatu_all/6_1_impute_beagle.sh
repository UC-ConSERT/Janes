#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Imputing with Beagle 5.2 - a haplotype likelihood based imputation
# Imputing tuturuatu_all bcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call bcf that has been filtered for 
#   4x coverage, 0.2 site missingness, 10 minGQ and no LD.

#Environment: impute

sppdir=~/data/tuturuatu_all/
## Edit to be run specific

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
subsetdir=${impdir}bcf_subsets/
mkdir -p ${impdir}imputations ${impdir}stats

# Imputation
beagle gt=${subsetdir}tuturuatu_study_phased.vcf.gz ref=${subsetdir}tuturuatu_ref_phased.vcf.gz impute=true gp=true out=${impdir}imputations/tuturuatu_beagle
bcftools index -f ${impdir}imputations/tuturuatu_beagle.vcf.gz

# To have a look at the imputation
zless -S ${impdir}imputations/tuturuatu_beagle.vcf.gz

# Concordance?
vcf-compare ${subsetdir}tuturuatu_study.vcf.gz ${impdir}imputations/tuturuatu_beagle.vcf.gz > ${impdir}stats/concordance_beagle5.txt

# r2?
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' ${impdir}imputations/tuturuatu_beagle.vcf.gz > ${impdir}stats/r2_beagle5.txt