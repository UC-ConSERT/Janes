#!/bin/bash -e
set -e

#25 April 2023
#Olivia Janes
#Formatting the output meg files to enable analysis in DnaSP
#CAUTION: only run this once on each file otherwise it might edit out important rows

#Environment: N/A

#Edit to be run specific:
sppdir=~/data/tuturuatu_all_vcf/
run=low_cov

#Setting directories
outdir=${sppdir}final_outputs/
tlrlist=${sppdir}final_outputs/tlr_list.txt

#Editing meg files to suit DnaSP
while read tlr_pos tlr_name
    #Cycle through TLRs listed in ${tlrlist}
do
    
    for i in {0.9,0.8,0.5,0}
    do
        echo ""; echo "Editing tlr_haplotypes_${i}miss_${tlr_name}.meg"
        
        #Remove TLR chr/pos lines (beginning with ">")
        grep -v "^>" -i ${outdir}meg/tlr_haplotypes_${i}miss_${tlr_name}.meg > tmp && mv tmp ${outdir}meg/tlr_haplotypes_${i}miss_${tlr_name}.meg
        wait

        #Add title line as second line
        sed -i "2s/^/TITLE:${tlr_name}\n/" ${outdir}meg/tlr_haplotypes_${i}miss_${tlr_name}.meg

    done
done < ${tlrlist}

echo ""; echo "Script is complete."
