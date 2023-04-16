#!/bin/bash -e
set -e

# 16 April 2023
# Olivia Janes
# Post-imputation filtering stats
# Just a quick wee script to look at snp counts pre- and post- filtering (postimputation)


sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific


# Setting variables.
bcfdir=${sppdir}bcf/
    #Define location of tlr_regions.bed file in script. Should be in bcf/
beagledir=${sppdir}impute/beagle_imputations/
mergedir=${beagledir}merged/
    #Directory of pre-filtered, merged imputation vcfs
filterdir=${beagledir}filtered/
    #Directory of filtered, merged imputation vcfs
mkdir -p ${beagledir}filtered/filter_stats
statsdir=${beagledir}filtered/filter_stats/

# (5_2) SNP Counts: Compiling all SNP counts from filtering, to compare between filtering methods.
echo "(5_2) SNP Counts beginning. Please fasten your seatbelts."

for dp in {0,4,5}
do
    for test_ne in {50,100,500}
    do

        snpfile=${statsdir}SNP_counts_${dp}x.txt
        tlrsnpfile=${statsdir}TLR_SNP_counts_${dp}x.txt

            ## Total SNP counts ##
        echo "Counting TOTAL SNPs in TLR contigs"
        echo "TOTAL SNPs in TLR contigs,Number" >> ${snpfile}

        # Prefiltered SNP counts
        #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script.
            echo "Prefiltered SNP counts for ${dp}x and ${test_ne}ne"

            for file in ${mergedir}Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_imp_merged.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${snpfile}
            done
            wait

        # Filtered SNP counts
            echo "Filtered SNP counts for ${dp}x and ${test_ne}ne"

            for file in ${filterdir}Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_filtered.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${snpfile}
            done
            wait


    ## TLR SNP counts ##
        echo "Counting TLR SNPs"
        echo "TLR SNPs,Number" >> ${tlrsnpfile}

        #Prefilter SNP counts
            echo "Prefilter TLR SNP counts for ${dp}x and ${test_ne}ne"
    
            for file in ${mergedir}Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_imp_merged.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${tlrsnpfile}
            done
            wait

        # Filtered SNP counts
            echo "Filtered SNP counts for ${dp}x and ${test_ne}ne"

            for file in ${filterdir}Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_filtered.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
                echo "${base},${snp}" >> ${tlrsnpfile}
            done
    done
done

echo "Script is complete."