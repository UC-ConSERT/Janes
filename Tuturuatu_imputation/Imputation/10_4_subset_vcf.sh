#!/bin/bash -e
set -e
# 29 April 2023
# Olivia Janes
# Subsetting the final vcfs into wild (both with and without 2019 Rangatira samples) and captive populations.

#Environment: bcftools

sppdir=~/data/tuturuatu_all_vcf/
    ## Edit to be run specific
wild_list_2021="CR01|CR02|CR03|CR04|CR05|CR06|CR07|CR08|CR10|CR11|CR12|CR13|CR14|CR15|CR16|CR18|CR19|CR20"
wild_list_2019="A|B|C0|C1|D|E|F|G|H"

# Defining directories.
impdir=${sppdir}impute/beagle_imputations/filtered/miss_filter_trial/
mkdir -p ${impdir}final_subsetted/
subsetdir=${impdir}final_subsetted/


# Setting population subsets
    echo ""; echo "Setting population subsets"
    for file in ${impdir}*_0.5miss_low_cov_final.vcf.gz
    do
        # Extracting individual IDs
        bcftools query -l ${file} > ${subsetdir}indv.ID
        
        # Subsetting wild populations
        grep -E "${wild_list_2021}" ${subsetdir}indv.ID > ${subsetdir}wild_2021.ID
            grep -E "${wild_list_2021}" ${subsetdir}indv.ID > ${subsetdir}wild_all.ID
        grep -E "${wild_list_2019}" ${subsetdir}indv.ID > ${subsetdir}wild_2019.ID
            grep -E "${wild_list_2019}" ${subsetdir}indv.ID >> ${subsetdir}wild_all.ID

        # Subsetting captive population
        grep -v -f ${subsetdir}wild_all.ID ${subsetdir}indv.ID > ${subsetdir}captive.ID
            # This extracts indvs not present in the wild populations.
        break
            # This stops the loop after one iteration.
    done


# Subsetting the final missing trial vcfs into populations
    for file in ${impdir}*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Subsetting ${base} vcf file into popls."
        # Captive population
        bcftools view -O z -o ${subsetdir}${base}_captive.vcf.gz -S ${subsetdir}captive.ID ${file}
        # Index the captive contig vcf
        bcftools index ${subsetdir}${base}_captive.vcf.gz -f --threads 16

        # Wild 2021 population
        bcftools view -O z -o ${subsetdir}${base}_wild_2021.vcf.gz -S ${subsetdir}wild_2021.ID ${file}
        # Index the study TLR contig vcf
        bcftools index ${subsetdir}${base}_wild_2021.vcf.gz -f --threads 16

        # Wild 2019 population
        bcftools view -O z -o ${subsetdir}${base}_wild_2019.vcf.gz -S ${subsetdir}wild_2019.ID ${file}
        # Index the study TLR contig vcf
        bcftools index ${subsetdir}${base}_wild_2019.vcf.gz -f --threads 16

        #Wild all population
        bcftools view -O z -o ${subsetdir}${base}_wild_all.vcf.gz -S ${subsetdir}wild_all.ID ${file}
        # Index the study TLR contig vcf
        bcftools index ${subsetdir}${base}_wild_all.vcf.gz -f --threads 16        
    done


echo ""
echo "Script has finished subsetting final filtered variant call files."
echo "Your final wild and captive subsetted files ${subsetdir}"
echo "Now ready for haplotyping!!"