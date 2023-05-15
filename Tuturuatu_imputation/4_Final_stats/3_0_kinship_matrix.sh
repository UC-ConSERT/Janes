#!/bin/bash -e 
set -e

#15 May 2023

# Olivia Janes
# Calculating a kinship matrix based on TLR imputed data.

#Environment: plink

##  Needs to be edited to be run specific   ##
    sppdir=~/data/tuturuatu_all_vcf/
    run=low_cov
    tlr_bed=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
        #Define the tlr regions bed location


#Setting directories
    mkdir -p ${sppdir}final_outputs/final_stats/kinship
    statsdir=${sppdir}final_outputs/final_stats/kinship/
    workdir=${sppdir}final_outputs/final_vcfs/


#Calculate kinship matrixes
    for file in ${workdir}*vcf.gz
    do
        base=$(basename ${file} _${run}_final.vcf.gz)
        echo "Making kinship matrix for ${base}"
        plink2 --vcf ${file} --make-king --allow-extra-chr --out ${statsdir}${base}_kinship > ${statsdir}${base}_kinship.log
        echo "Making kinship table for ${base}"
        plink2 --vcf ${file} --make-king-table --allow-extra-chr --out ${statsdir}${base}_kinship >> ${statsdir}${base}_kinship.log
    done



#Downloading stats to home computer. Navigate to final_outputs folder.
#mkdir -p final_stats/kinship/
#rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/kinship/* final_stats/kinship/

echo "Script has finished"