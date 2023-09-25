#!/bin/bash -e
set -e

# 25 April 2023

# Olivia Janes
# Formatting the output meg files to enable analysis in DnaSP
# CAUTION: only run this once on each file otherwise it might edit out important rows
# From: tuturuatu_imputation

## Environment: N/A

# Setting up:
    #Edit to be run specific:
    sppdir=~/data/tuturuatu_all_vcf/
    run=low_cov

    #Setting directories
    megdir=${sppdir}final_outputs/subset_popls/meg_subsets/
    tlrlist=${sppdir}final_outputs/tlr_list.txt

#Editing meg files to suit DnaSP
    while read tlr_pos tlr_name
        #Cycle through TLRs listed in ${tlrlist}
    do
        
        for i in {0.9,0.8,0.5,0}
        do
            for popl in {captive,wild_2021,wild_2019,wild_all}
            do

            echo ""; echo "Editing tlr_haplotypes_${i}miss_${tlr_name}_${popl}.meg"
            
            #Remove TLR chr/pos lines (beginning with ">")
            grep -v "^>" -i ${megdir}tlr_haplotypes_${i}miss_${tlr_name}_${popl}.meg > tmp && mv tmp ${megdir}tlr_haplotypes_${i}miss_${tlr_name}_${popl}.meg
            wait

            #Add title line as second line
            sed -i "2s/^/TITLE:${tlr_name}\n/" ${megdir}tlr_haplotypes_${i}miss_${tlr_name}_${popl}.meg

            done
        done
    done < ${tlrlist}

echo ""; echo "Script is complete."
