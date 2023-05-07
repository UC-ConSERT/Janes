#!/bin/bash -e

# 24 April 2023
# Olivia Janes
# Preparing for imputation:
# Filtering the preimpute high cov reference panel (5x, 100Ne) for missingness <0.1, 0.2 and 0.5.
#    Then, removing these sites in the imputed study and ref populations (5x, 100Ne).

#Environment: samtools

## Edit to be run specific:
bcfdir=~/data/tuturuatu_all_vcf/bcf/
sppdir=~/data/tuturuatu_all_vcf/
missdir=${sppdir}impute/missingness_trial/
evaldir=${sppdir}impute/
    mkdir -p ${evaldir}beagle_imputations/filtered/miss_filter_trial/
    filterdirs=beagle_imputations/filtered/miss_filter_trial/
        #This is just a suffix to be used in conjunction with any run level directory
run=low_cov

#Defining folders
subsetdir=${sppdir}impute/vcf_subsets/
    # Directory that holds the preimputation filtered variant calls for the reference panel individuals (phased)
mkdir -p ${missdir}vcf_ref_merged/ ${missdir}removed_sites/miss_trial
mergedir=${missdir}vcf_ref_merged/
removedir=${missdir}removed_sites/miss_trial



#Merge the Pre-imputation, TLR contig-seperated, NOT phased, reference vcf back into one file for all TLR contigs
    # I will only be needing the merged 5x vcf.

    for subset_1 in ${subsetdir}*5x_1_ref.vcf.gz
    do
        base=$(basename ${subset_1} _1_ref.vcf.gz)
        subset_2=${subsetdir}${base}_2_ref.vcf.gz
        subset_3=${subsetdir}${base}_3_ref.vcf.gz

        #Create a file list of files to merge
        echo ""; echo "Merging subsetted tlr contigs in:"
        echo "${subset_1}"; echo "${subset_1}" > ${mergedir}merge_list.txt
        echo "${subset_2}"; echo "${subset_2}" >> ${mergedir}merge_list.txt
        echo "${subset_3}"; echo "${subset_3}" >> ${mergedir}merge_list.txt

        #Merge and index
        bcftools concat -O z --threads 16 -f ${mergedir}merge_list.txt -o ${mergedir}${base}_ref_merged.vcf.gz
        echo "Indexing ${base}_study_merged.vcf.gz"
        bcftools index -f --threads 16 ${mergedir}${base}_ref_merged.vcf.gz
    done


# Filter reference panel for missingness, and output removed sites into a text file. 
    for i in {0.9,0.8,0.5,0}
    do
        for file in ${mergedir}*5x_ref_merged.vcf.gz
        do
            echo ""; echo "Printing sites with present in less than ${i} of indv."
            base=$(basename ${file} _ref_merged.vcf.gz)
            vcftools --gzvcf ${file} \
                --max-missing ${i} \
                --out ${removedir}${base}_ref_${i}miss \
                --removed-sites
            echo "Finished printing missing sites, output found at: ${removedir}"
        done
    done


# Remove the missing sites from the imputed, merged final vcfs.
    for i in {0.9,0.8,0.5,0}
    do
        for file in ${evaldir}beagle_imputations/filtered/MAF_filtered_only/*5x_100ne_0.05MAF.vcf.gz
        do
            base=$(basename ${file} _5x_100ne_0.05MAF.vcf.gz)
            vcftools --gzvcf ${file} \
                --exclude-positions ${removedir}Tuturuatu_VariantCalls_5x_ref_${i}miss.removed.sites \
                --recode \
                --recode-INFO-all \
                --out ${evaldir}${filterdirs}${base}_${i}miss_${run}_final.vcf

            echo "Renaming filter file to remove '.recode.vcf'"
            mv -i ${evaldir}${filterdirs}${base}_${i}miss_${run}_final.vcf.recode.vcf ${evaldir}${filterdirs}${base}_${i}miss_${run}_final.vcf
        done
    done


# Convert filter files to vcf.gz format and index
    for vcf in ${evaldir}${filterdirs}*vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${evaldir}${filterdirs}${base}.vcf.gz --threads 16

        echo "Indexing ${vcf}.gz"
        bcftools index ${evaldir}${filterdirs}${base}.vcf.gz --threads 16

        if [ -e "${evaldir}${filterdirs}${base}.vcf.gz" ]; then
            echo "Removing ${evaldir}${filterdirs}${base}.vcf"
            rm "${evaldir}${filterdirs}${base}.vcf"
        fi
        echo ""
    done


echo ""
echo "Script has finished preparing filtered variant call files in ${truthdir}${filterdirs}"
