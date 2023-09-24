#!/bin/bash -e

#19 Sep 2022
#Olivia Janes, from Molly Magid
#Filtering the "chosen" filtered bcfs, to remove any sites where strand bias adjusted phred score is <60

sppdir=~/data/tuturuatu/
filterdir=${sppdir}bcf/filter_trial/noLD/

mkdir -p ${sppdir}bcf_final/

<<"COMMENTS"

#Filtering chosen filtered file to export only sites where strand bias adjusted phred score is <60
bcftools query -i 'FORMAT/SP<60' -f %CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\n \
    ${filterdir}Tuturuatu_VariantCalls_4x_coverage_0.1site_missing_MinGQ10.bcf.recode.bcf \
    -o ${filterdir}Tuturuatu_VariantCalls_4x_coverage_0.1site_missing_MinGQ10_sp60.bcf

bcftools query -i 'FORMAT/SP<60' -f %CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\n \
    ${filterdir}Tuturuatu_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode.bcf \
    -o ${filterdir}Tuturuatu_VariantCalls_5x_coverage_0.1site_missing_MinGQ10_sp60.bcf

#[%GT:PL:DP:SP:ADF:ADR:AD:GQ]
#Error: no such tag defined in the VCF header: INFO/FORMATn

COMMENTS
#How is this different to above? Not sure if molly used both or just below...? 
    #I think Molly used just below according to her MSc Manuscript.
    #This sets individual sites with SP>60 to...?missing?
#filter bcf file for depth and strand bias at individual sites
## Final filtered vcf must be selected
bcftools +setGT ${finaldir}#chosen_filtered_vcf_to_be_added_here -- -t q -n .##[perhaps this needs to be ./.] -i 'FORMAT/SP>60' > ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.vcf
echo "Filtering for strand bias is complete."


<<"COMMENTS2"
#Pick one: 

echo "GZipping and indexing vcf file."
bgzip ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.vcf
tabix -f -p vcf ${finaldir}Tuturuatu_VariantCalls_final_variants.vcf.gz
echo "GZipping and indexing vcf file is complete."

echo "GZipping and indexing bcf file."
bgzip ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.bcf
bcftools index ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.bcf.gz
echo "GZipping and indexing bcf file is complete."

COMMENTS2