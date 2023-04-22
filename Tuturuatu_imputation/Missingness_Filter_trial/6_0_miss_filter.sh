#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Preparing for imputation:
# Filtering high cov reference panel indv for missingness <0.2, and then removing these sites in the study populations:
#   This will be done for the validation and truth runs first, to finalise the vcfs and prepare them for comparisons.
# Preparing for imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters.

#Environment: samtools

## Edit to be run specific:
bcfdir=~/data/tuturuatu_all_vcf/bcf/
sppdir=~/data/tuturuatu_all_vcf/
missdir=${sppdir}impute/missingness_trial/
truthdir=${sppdir}impute/truth/
valdir=${sppdir}impute/validation/
    filterdirs=impute/beagle_imputations/filtered/
        #This is just a suffix to be used in conjunction with any run level directory

beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

#Defining folders
subsetdir=${sppdir}impute/vcf_subsets/
    # Directory that holds the preimputation filtered variant calls for the reference panel individuals (phased)
mkdir -p ${missdir}vcf_ref_merged/ ${missdir}removed_sites/
mergedir=${missdir}vcf_ref_merged/
removedir=${missdir}removed_sites/


<<"COMMENTS"
# Making a directory to hold the imputation work.
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}vcf_finals/
subsetdir=${impdir}vcf_subsets/
finaldir=${impdir}vcf_finals/
COMMENTS


<<"DONE"
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
            --out ${valdir}${filterdirs}${base}_validation_final.vcf

        echo "Renaming filter file to remove '.recode.vcf'"
        mv -i ${valdir}${filterdirs}${base}_validation_final.vcf.recode.vcf ${valdir}${filterdirs}${base}_validation_final.vcf
    done
DONE
    #Truth trial
    for file in ${truthdir}${filterdirs}MAF_filtered_only/*_0.05MAF.vcf.gz
    do
        base=$(basename ${file} _0.05MAF.vcf.gz)
        vcftools --gzvcf ${file} \
            --exclude-positions ${removedir}Tuturuatu_VariantCalls_5x_ref_0.1miss.removed.sites \
            --recode \
            --recode-INFO-all \
            --out ${truthdir}${filterdirs}${base}_truth_final.vcf

        echo "Renaming filter file to remove '.recode.vcf'"
        mv -i ${truthdir}${filterdirs}${base}_truth_final.vcf.recode.vcf ${truthdir}${filterdirs}${base}_truth_final.vcf
    done



echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
echo "Now ready for imputing!!"