#!/bin/bash -e
set -e

# 13 April 2023
# Olivia Janes
# Merge the vcfs for the imputed individuals back with the high coverage individuals.

#Environment: samtools

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific


# Setting variables.
impdir=${sppdir}impute/
mkdir -p ${impdir}beagle_imputations/merged
impoutdir=${impdir}beagle_imputations/impute_trials/
mergedir=${impdir}beagle_imputations/merged/
mkdir -p ${impdir}beagle_imputations/merged/intermediate_merged


# Merge ref and study individuals back together
    # Ref: ~/data/tuturuatu_all_vcf/impute/vcf_finals/Tuturuatu_VariantCalls_[0,4,5]x_[1,2,3]_ref_phased.vcf.gz
    # Study: ~/data/tuturuatu_all_vcf/impute/beagle_imputations/impute_trials/Tuturuatu_VariantCalls_[0,4,5]x_[1,2,3]_[50,100,500]ne_beagle_imp.vcf.gz
    for test_ne in {50,100,500}
    do
        echo ""
        for vcf in ${impoutdir}*_${test_ne}ne_beagle_imp.vcf.gz
        do
            base=$(basename ${vcf} _${test_ne}ne_beagle_imp.vcf.gz)
            study=$(basename ${vcf} _beagle_imp.vcf.gz)
            ref=${impdir}vcf_finals/${base}_ref_phased.vcf.gz

            #Creating file list
            echo "${vcf}" > ${mergedir}intermediate_merged/merge_list.txt
            echo "${ref}" >> ${mergedir}intermediate_merged/merge_list.txt

            echo ""; echo "Merging individuals in ${base}_${test_ne}ne_beagle_imp.vcf.gz and ${ref}"
            bcftools merge -O z --threads 16 --file-list ${mergedir}intermediate_merged/merge_list.txt -o ${mergedir}intermediate_merged/${study}_imp_merged.vcf.gz
            echo "Indexing ${study}_imp_merged.vcf.gz"
            bcftools index -f --threads 16 ${mergedir}intermediate_merged/${study}_imp_merged.vcf.gz
        done
    done


# Merge the 3 seperated TLR contig files back together
    for test_ne in {50,100,500}
    do
        for subset_1 in ${mergedir}intermediate_merged/*1_${test_ne}ne_imp_merged.vcf.gz
        do
            base=$(basename ${subset_1} _1_${test_ne}ne_imp_merged.vcf.gz)
            subset_2=${mergedir}intermediate_merged/${base}_2_${test_ne}ne_imp_merged.vcf.gz
            subset_3=${mergedir}intermediate_merged/${base}_3_${test_ne}ne_imp_merged.vcf.gz

            #Create a file list of files to merge
            echo ""; echo "Merging subsetted tlr contigs in:"
            echo "${subset_1}"; echo "${subset_1}" > ${mergedir}merge_list.txt
            echo "${subset_2}"; echo "${subset_2}" >> ${mergedir}merge_list.txt
            echo "${subset_3}"; echo "${subset_3}" >> ${mergedir}merge_list.txt

            #Merge and index
            bcftools concat -O z --threads 16 -f ${mergedir}merge_list.txt -o ${mergedir}${base}_${test_ne}ne_imp_merged.vcf.gz
            echo "Indexing ${base}_${test_ne}ne_imp_merged.vcf.gz"
            bcftools index -f --threads 16 ${mergedir}${base}_${test_ne}ne_imp_merged.vcf.gz
        done
    done

echo ""; echo "Merging files with just a click,"
echo "No more scattered data, it's slick."
echo "..."; echo "Script is complete, if you're not poetically inclined."