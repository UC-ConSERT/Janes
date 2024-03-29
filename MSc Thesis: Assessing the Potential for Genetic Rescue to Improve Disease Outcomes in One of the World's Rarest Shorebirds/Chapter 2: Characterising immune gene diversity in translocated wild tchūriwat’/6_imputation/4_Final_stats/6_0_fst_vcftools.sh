#!/bin/bash -e 
set -e

# 22 May 2023

# Olivia Janes
# Pairwise Fst with vcftools
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    ##  Needs to be edited to be run specific   ##
        sppdir=~/data/tuturuatu_all_vcf/
        impvcf=${sppdir}final_outputs/final_vcfs/Tuturuatu_VariantCalls_0.9miss_low_cov_final.vcf.gz
        wgsvcf=${sppdir}final_outputs/final_vcfs/whole_genome/Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz
        ogvcf=${sppdir}bcf/filter_trial/noLD/Tuturuatu_VariantCalls_5x_coverage_0.2site_missing_MinGQ10.vcf.gz

    #Setting directories
        mkdir -p ${sppdir}final_outputs/final_stats/fst/
        statsdir=${sppdir}final_outputs/final_stats/fst/


#Setting populations
    # Extracting individual IDs
        bcftools query -l ${ogvcf} > ${statsdir}indv.ID

        # Subsetting captive population
        grep -E "CT|I164" ${statsdir}indv.ID > ${statsdir}captive.ID
            # This extracts CT... or I164...

        # Subsetting wild population
        grep -v -f ${statsdir}captive.ID ${statsdir}indv.ID > ${statsdir}wild_all.ID
            # This extracts indvs not present in the captive population.

        # Subsetting wild_2021
        grep -E "CR" ${statsdir}indv.ID > ${statsdir}wild_2021.ID
            # This extracts indvs CR01-20.
        
        # Subsetting wild_2019
        grep -E "A|B|C0|C1|D|E|F|G|H" ${statsdir}indv.ID > ${statsdir}wild_2019.ID

# Running vcftools pairwise fst
    #WGS
        #Captive vs Wild_all
            for kb in {10,500}
            do
                echo ""; echo "Running fst calculations on wgs, ${kb}kb"
                vcftools --gzvcf ${ogvcf} \
                    --weir-fst-pop ${statsdir}captive.ID \
                    --weir-fst-pop ${statsdir}wild_all.ID \
                    --fst-window-size ${kb}000 \
                    --out ${statsdir}wgs_captive_wild_all_fst_${kb}kb
                echo ""; echo "Finished wgs, ${kb}kb"
            done

        #Captive vs Wild_2019 vs Wild_2021
            for kb in {10,500}
            do
                echo ""; echo "Running fst calculations on wgs, ${kb}kb, splitting wild popls"
                vcftools --gzvcf ${ogvcf} \
                    --weir-fst-pop ${statsdir}captive.ID \
                    --weir-fst-pop ${statsdir}wild_2019.ID \
                    --weir-fst-pop ${statsdir}wild_2021.ID \
                    --fst-window-size ${kb}000 \
                    --out ${statsdir}wgs_captive_wild_years_fst_${kb}kb
                echo ""; echo "Finished wgs, ${kb}kb, splitting wild popls"
            done

    #TLRs
        #Captive vs Wild_2019 vs Wild_2021
            echo ""; echo "Running fst calculations on tlr, splitting wild popls"
            vcftools --gzvcf ${impvcf} \
                --weir-fst-pop ${statsdir}captive.ID \
                --weir-fst-pop ${statsdir}wild_2019.ID \
                --weir-fst-pop ${statsdir}wild_2021.ID \
                --out ${statsdir}tlr_captive_wild_years_fst
            echo ""; echo "Finished tlr, splitting wild popls"

    #WGS trying subpopl by subpopl
        #Captive vs Wild
            for wild in {wild_2019.ID,wild_2021.ID}
            do
                echo ""; echo "Running fst calculations on wgs, 10kb, captive vs ${wild}"
                vcftools --gzvcf ${ogvcf} \
                    --weir-fst-pop ${statsdir}captive.ID \
                    --weir-fst-pop ${statsdir}${wild} \
                    --fst-window-size 10000 \
                    --out ${statsdir}wgs_captive_${wild}_fst_10kb
                echo ""; echo "Finished wgs, 10kb"
            done

        #Wild 2019 vs Wild 2021
            echo ""; echo "Running fst calculations on wgs, 10kb, wild_2019 vs wild_2021"
            vcftools --gzvcf ${ogvcf} \
                --weir-fst-pop ${statsdir}wild_2019.ID \
                --weir-fst-pop ${statsdir}wild_2019.ID\
                --fst-window-size 10000 \
                --out ${statsdir}wgs_wild_2019_2021_fst_10kb
            echo ""; echo "Finished wgs, 10kb"


    #TLR trying subpopl by subpopl
        #Captive vs Wild
            for wild in {wild_2019.ID,wild_2021.ID}
            do
                echo ""; echo "Running fst calculations on tlr, 10kb, captive vs ${wild}"
                vcftools --gzvcf ${impvcf} \
                    --weir-fst-pop ${statsdir}captive.ID \
                    --weir-fst-pop ${statsdir}${wild} \
                    --fst-window-size 10000 \
                    --out ${statsdir}tlr_captive_${wild}_fst_10kb
                echo ""; echo "Finished tlr, 10kb"
            done

        #Wild 2019 vs Wild 2021
            echo ""; echo "Running fst calculations on tlr, 10kb, wild_2019 vs wild_2021"
            vcftools --gzvcf ${impvcf} \
                --weir-fst-pop ${statsdir}wild_2019.ID \
                --weir-fst-pop ${statsdir}wild_2019.ID\
                --fst-window-size 10000 \
                --out ${statsdir}tlr_wild_2019_2021_fst_10kb
            echo ""; echo "Finished tlr, 10kb"


echo ""; echo "To download all of the stats, navigate to the right directory on your desktop: ~/Documents/Tuturuatu_resources/tuturuatu_all_vcf/final_outputs/final_stats/"
echo "Enter code (edited for the right run):"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/fst ./"

echo ""
echo "Whole genome genotype script is complete."       
