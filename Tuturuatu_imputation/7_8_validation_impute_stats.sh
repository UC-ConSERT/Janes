#!/bin/bash -e

# 14 April 2023
# Olivia Janes
# Imputation with Beagle 5.4 (beagle.22Jul22.46e.jar)
# Running stats on the validation imputation file to compare it to known SNPs for the validation (high coverage) individuals

#Environment: impute
run=tuturuatu_all_vcf
sppdir=~/data/${run}/impute/validation/
## Edit to be run specific

#Setting variables
vcfdir=~/data/${run}/bcf/filter_strand_bias/
impdir=${sppdir}impute/
finaldir=${impdir}beagle_imputations/
mkdir -p ${impdir}stats

#Concordance
#Compare between imputed and confirmed SNPs
for file in ${finaldir}*vcf.gz
do
    base=$(basename ${file} _beagle_imp.vcf.gz)
    vcf-compare ${vcfdir}Tuturuatu_VariantCalls_5x_coverage_0.1site_missing_MinGQ10_0.6SP.vcf.gz \
    ${file} > ${impdir}stats/${base}_concordance.txt

#Investigate SNPs
    #Extract information on the haplotypes at each TLR SNP for each individual
        bcftools query -R ~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' \
            ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_4x_0.6SP.vcf.gz >> ~/data/tuturuatu_all_vcf/scripts/tlr_haps_4x.txt
    #Extract the headers to add to the above
        bcftools view -h ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_0x_0.6SP.vcf.gz \
            | tail -n 1 >> ~/data/tuturuatu_all_vcf/scripts/tlr_haps_4x_header.txt
    #Download these and extract into a spreadsheet to analyse
