#!/bin/bash -e

# 13 April 2023

# Olivia Janes
# Imputation stats: comparing all of the trials (depth and effective popl size) within the truth imputation to the
#   chosen truth imputation (5x and 100ne) to ensure that they are all mostly the same. This ensures the chosen truth
#   imputation is the most accurate representation of the truth, regardless of depth or Ne.
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    ## Edit to be run specific
    sppdir=~/data/tuturuatu_all_vcf/impute/truth/
    impstats=~/data/tuturuatu_all_vcf/impute/
        #Location of imputation stats folder
    tlr_regions=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
        #Define location of TLR regions bed file
    run=truth

    # Defining directories
    impdir=${sppdir}impute/
    finaldir=${impdir}vcf_finals/
    mkdir -p ${impdir}vcf_finals/vcf_merged/
    mergedir=${impdir}vcf_finals/vcf_merged/

    mkdir -p ${impstats}stats/${run}_stats/
    statsdir=${impstats}stats/${run}_stats/


# To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/filtered/[imputation.vcf.gz]

#Merge the Pre-imputation, TLR contig-separated, study vcf back into one file for all TLR contigs
    for dp in {0,4,5}
    do
        for subset_1 in ${impdir}vcf_finals/*${dp}x_1_study.vcf.gz
        do
            base=$(basename ${subset_1} _1_study.vcf.gz)
            subset_2=${impdir}vcf_finals/${base}_2_study.vcf.gz
            subset_3=${impdir}vcf_finals/${base}_3_study.vcf.gz

            #Create a file list of files to merge
            echo ""; echo "Merging subsetted tlr contigs in:"
            echo "${subset_1}"; echo "${subset_1}" > ${mergedir}merge_list.txt
            echo "${subset_2}"; echo "${subset_2}" >> ${mergedir}merge_list.txt
            echo "${subset_3}"; echo "${subset_3}" >> ${mergedir}merge_list.txt

            #Merge and index
            bcftools concat -O z --threads 16 -f ${mergedir}merge_list.txt -o ${mergedir}${base}_study_${run}_merged.vcf.gz
            echo "Indexing ${base}_study_${run}_merged.vcf.gz"
            bcftools index -f --threads 16 ${mergedir}${base}_study_${run}_merged.vcf.gz
        done
    done



#Investigating TLR Genotypes: Extract the TLR Genotypes out of the pre-imputed files and imputed at Ne=100 files
    for dp in {0,4,5}
    do
        echo ""; echo "Extracting TLR Genotypes for ${dp}x files, pre and post impute"

        #All individual's Genotypes:
            file0=~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/*VariantCalls_${dp}x_0.6SP.vcf.gz
            #Extract information on the Genotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file0} > ${statsdir}tlr_genotypes_preimpute_${dp}x_all_indv_${run}.txt
            #Extract the headers to add to the above
                bcftools view -h ${file0} | tail -n 1 > ${statsdir}tlr_genotypes_preimpute_header_all_indv_${run}.txt
            #Download these and extract into a spreadsheet to analyse

        #Preimpute Genotypes
            file=${mergedir}*VariantCalls_${dp}x_study_${run}_merged.vcf.gz
            #Extract information on the genotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file} > ${statsdir}tlr_genotypes_preimpute_${dp}x_${run}.txt
            #Extract the headers to add to the above
                bcftools view -h ${file} | tail -n 1 > ${statsdir}tlr_genotypes_preimpute_header_${run}.txt
            #Download these and extract into a spreadsheet to analyse

        #Imputed Genotypes (for 100ne only)
            file2=${impdir}beagle_imputations/filtered/*VariantCalls_${dp}x_100ne_filtered.vcf.gz
            #Extract information on the Genotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file2} > ${statsdir}tlr_genotypes_impute_${dp}x_${run}.txt
            #Extract the headers to add to the above
                bcftools view -h ${file2} | tail -n 1 > ${statsdir}tlr_genotypes_impute_header_${run}.txt
            #Download these and extract into a spreadsheet to analyse

    done
    echo "TLR Genotype files can be found at: ${statsdir}tlr_genotypes..."



echo ""; echo "To download all of the stats, navigate to the right directory on your desktop: ~/Documents/Tuturuatu_resources/tuturuatu_all_vcf/impute/impute_stats/"
echo "Enter code (edited for the right run):"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/impute/stats/${run}_stats/* ./"

echo ""
echo "Imputation stats script is complete."