#!/bin/bash -e

# 13 April 2023

# Olivia Janes
# Preparing for imputation:
# Filtering high cov reference panel indv for missingness <0.2, and then removing these sites in the study populations:
#   This will be done for the validation and truth runs first, to finalise the vcfs and prepare them for comparisons.
# Preparing for imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters.
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    ## Edit to be run specific:
    bcfdir=~/data/tuturuatu_all_vcf/bcf/
    sppdir=~/data/tuturuatu_all_vcf/
    missdir=${sppdir}impute/missingness_trial/
    truthdir=${sppdir}impute/truth/
    valdir=${sppdir}impute/validation/
        filterdirs=impute/beagle_imputations/filtered/
            #This is just a suffix to be used in conjunction with any run level directory
    run1=truth
    run2=validation

    #Defining folders
    subsetdir=${sppdir}impute/vcf_subsets/
        # Directory that holds the preimputation filtered variant calls for the reference panel individuals (phased)
    mkdir -p ${missdir}vcf_ref_merged/ ${missdir}removed_sites/
    mergedir=${missdir}vcf_ref_merged/
    removedir=${missdir}removed_sites/


#Merge the Pre-imputation, TLR contig-seperated, NOT phased, reference vcf back into one file for all TLR contigs
    # For the comparisons (validation vs truth), I will only be needing the merged 5x vcf.
    for dp in {0,4,5}
    do
        for subset_1 in ${subsetdir}*${dp}x_1_ref.vcf.gz
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
    done


# Filter reference panel for missingness, and output removed sites into a text file. 
    for file in ${mergedir}*5x_ref_merged.vcf.gz
    do
        echo ""; echo "Printing sites with missingness over 10%"
        base=$(basename ${file} _ref_merged.vcf.gz)
        vcftools --gzvcf ${file} \
            --max-missing 0.9 \
            --out ${removedir}${base}_ref_0.1miss \
            --removed-sites
        echo "Finished printing missing sites, output found at: ${removedir}"
    done


# Remove the missing sites from the imputed, merged final vcfs.
    #This will be done for both truth and validation trials, to allow a comparison of the same high quality sites.
    #Validation trial
    for file in ${valdir}${filterdirs}MAF_filtered_only/*_0.05MAF.vcf.gz
    do
        base=$(basename ${file} _0.05MAF.vcf.gz)
        vcftools --gzvcf ${file} \
            --exclude-positions ${removedir}Tuturuatu_VariantCalls_5x_ref_0.1miss.removed.sites \
            --recode \
            --recode-INFO-all \
            --out ${valdir}${filterdirs}${base}_${run2}_final.vcf

        echo "Renaming filter file to remove '.recode.vcf'"
        mv -i ${valdir}${filterdirs}${base}_${run2}_final.vcf.recode.vcf ${valdir}${filterdirs}${base}_${run2}_final.vcf
    done

    #Truth trial
    for file in ${truthdir}${filterdirs}MAF_filtered_only/*_0.05MAF.vcf.gz
    do
        base=$(basename ${file} _0.05MAF.vcf.gz)
        vcftools --gzvcf ${file} \
            --exclude-positions ${removedir}Tuturuatu_VariantCalls_5x_ref_0.1miss.removed.sites \
            --recode \
            --recode-INFO-all \
            --out ${truthdir}${filterdirs}${base}_${run1}_final.vcf

        echo "Renaming filter file to remove '.recode.vcf'"
        mv -i ${truthdir}${filterdirs}${base}_${run1}_final.vcf.recode.vcf ${truthdir}${filterdirs}${base}_${run1}_final.vcf
    done


# Convert filter files to vcf.gz format and index
    #Validation trial
    for vcf in ${valdir}${filterdirs}*vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${valdir}${filterdirs}${base}.vcf.gz --threads 16

        echo "Indexing ${vcf}.gz"
        bcftools index ${valdir}${filterdirs}${base}.vcf.gz --threads 16

        if [ -e "${valdir}${filterdirs}${base}.vcf.gz" ]; then
            echo "Removing ${valdir}${filterdirs}${base}.vcf"
            rm "${valdir}${filterdirs}${base}.vcf"
        fi
        echo ""
    done

    #Truth trial
    for vcf in ${truthdir}${filterdirs}*vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${truthdir}${filterdirs}${base}.vcf.gz --threads 16

        echo "Indexing ${vcf}.gz"
        bcftools index ${truthdir}${filterdirs}${base}.vcf.gz --threads 16

        if [ -e "${truthdir}${filterdirs}${base}.vcf.gz" ]; then
            echo "Removing ${truthdir}${filterdirs}${base}.vcf"
            rm "${truthdir}${filterdirs}${base}.vcf"
        fi
        echo ""
    done

echo ""
echo "Script has finished preparing filtered variant call files in ${truthdir}${filterdirs}"
