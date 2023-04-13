#!/bin/bash -e 

#13 April 2023
#Olivia Janes
#Calculating Ne (effective population size) using SNeP.
#SNeP should be installed already.

## NOTE: Do not use SNeP. SNeP not working (it started investigating the first chr, and never stopped). Try a different program to estimate Ne.

sppdir=~/data/tuturuatu_all_vcf/
bcfdir=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)
vcffile=${bcfdir}filter_strand_bias/Tuturuatu_VariantCalls_5x_coverage_0.1site_missing_MinGQ10_0.6SP.vcf.gz
    #filtered variant calls file
snep=~/data/SNeP_111
    #location of SNeP installation


mkdir -p ${bcfdir}stats/ne/
nedir=${bcfdir}stats/ne/


# Converting the vcf.gz file into a .ped and .map for use in SNeP.
    echo "Converting ${vcffile} to .ped and .map using PLINK"
    name=$(basename ${vcffile} .vcf.gz)
    plink --vcf ${vcffile} --recode --double-id --allow-extra-chr --out ${nedir}${name}

# Calculating Ne using SNeP
    echo ""; echo "Calculating Ne using SNeP"
    ${snep} -ped ${nedir}${name}.ped -map ${nedir}${name}.map -out ${nedir}Tuturuatu_SNeP -threads 16

echo "Script has finished. Find output files as ${nedir}Tuturuatu_SNeP...


