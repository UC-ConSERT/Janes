#!/bin/bash -e
set -e


# 13 April 2023
# Olivia Janes
# Preparing for imputation:
# Filtering high cov reference panel indv for missingness <0.2, and then removing these sites in the study populations:
#   This will be done for the validation and truth runs first, to finalise the vcfs and prepare them for comparisons.
# Preparing for imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters.

#Environment: samtools

## Edit to be run specific:
sppdir=~/data/tuturuatu_all_vcf/
missdir=${sppdir}impute/missingness_trial/
run=downsampling

for ds in {0.05,0.1,0.2,0.3}
do

    echo "Beginning script for ${ds}"



    evaldir=${sppdir}impute/downsampling_trial/${ds}_downsample/
        filterdirs=impute/beagle_imputations/filtered/
            #This is just a suffix to be used in conjunction with any run level directory


    #Defining folders
    removedir=${missdir}removed_sites/


    # Remove the missing sites from the imputed, merged final vcfs.
        for file in ${evaldir}${filterdirs}MAF_filtered_only/*_0.05MAF.vcf.gz
        do
            base=$(basename ${file} _0.05MAF.vcf.gz)
            vcftools --gzvcf ${file} \
                --exclude-positions ${removedir}Tuturuatu_VariantCalls_5x_ref_0.1miss.removed.sites \
                --recode \
                --recode-INFO-all \
                --out ${evaldir}${filterdirs}${base}_${run}_final.vcf

            echo "Renaming filter file to remove '.recode.vcf'"
            mv -i ${evaldir}${filterdirs}${base}_${run}_final.vcf.recode.vcf ${evaldir}${filterdirs}${base}_${run}_final.vcf
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

done

echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
