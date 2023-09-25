#!/bin/bash -e
set -e

# 16 April 2023

# Olivia Janes
# Post-imputation filtering stats for Missingness Trial
# Just a quick wee script to look at snp counts pre- and post- filtering (postimputation)
# From: tuturuatu_imputation

## Environment: samtools or bcftools

# Setting up
    sppdir=~/data/tuturuatu_all_vcf/
    ## Edit to be run specific
    run=low_cov

    # Setting variables.
    bcfdir=${sppdir}bcf/
        #Define location of tlr_regions.bed file in script. Should be in bcf/
    beagledir=${sppdir}impute/beagle_imputations/
    mergedir=${beagledir}merged/
        #Directory of pre-filtered, merged imputation vcfs
    filterdir=${beagledir}filtered/miss_filter_trial/
        #Directory of filtered, merged imputation vcfs
    mkdir -p ${sppdir}impute/stats/
    mkdir -p ${sppdir}impute/stats/beagle_imp_stats/miss_filter_trial/
    statsdir=${sppdir}impute/stats/beagle_imp_stats/miss_filter_trial/


# (5_2) SNP Counts: Compiling all SNP counts from filtering, to compare between filtering methods.
echo "(5_2) SNP Counts beginning. Please fasten your seatbelts."
    
    snpfile=${statsdir}SNP_counts_${run}_miss.txt
    tlrsnpfile=${statsdir}TLR_SNP_counts_${run}_miss.txt

    ## Total SNP counts ##
    echo "Counting TOTAL SNPs in TLR contigs"
    echo "TOTAL SNPs in TLR contigs,Number" >> ${snpfile}    

    # Prefiltered SNP counts
    #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script.
        echo "Prefiltered SNP counts for 5x and 100ne"

        for file in ${mergedir}Tuturuatu_VariantCalls_5x_100ne_imp_merged.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${snpfile}
        done
        wait

    ## TLR SNP counts ##
    echo "Counting TLR SNPs"
    echo "TLR SNPs,Number" >> ${tlrsnpfile}

    #Prefilter TLR SNP counts
        echo "Prefilter TLR SNP counts for 5x and 100ne"

        for file in ${mergedir}Tuturuatu_VariantCalls_5x_100ne_imp_merged.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${tlrsnpfile}
        done
        wait

    for i in {0.9,0.8,0.5,0}
    do

        # Filtered SNP counts
            echo "Filtered SNP counts for 5x, 100ne and 1-${i} missingness"

            for file in ${filterdir}Tuturuatu_VariantCalls_${i}miss_low_cov_final.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${snpfile}
            done
            wait

        # Filtered TLR SNP counts
            echo "Filtered TLR SNP counts for 5x, 100ne and 1-${i} missingness"

            for file in ${filterdir}Tuturuatu_VariantCalls_${i}miss_low_cov_final.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${tlrsnpfile}
            done
    done


echo "Script is complete."