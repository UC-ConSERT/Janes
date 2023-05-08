#!/bin/bash -e 
set -e

#07 May 2023

# Olivia Janes
# Extracting Allele Frequencies for each individual over the various populations.

#Environment: samtools

##  Needs to be edited to be run specific   ##
sppdir=~/data/tuturuatu_all_vcf/
run=low_cov
tlr_bed=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
    #Define the tlr regions bed location
finaldir=${sppdir}impute/beagle_imputations/filtered/miss_filter_trial/
    #Directory holding the final (missingness trial) vcfs.
    #Will run the stats on all to compare, see if there is a difference in them.

#Setting directories
    mkdir -p ${sppdir}final_outputs/final_stats/subsetted/ ${sppdir}final_outputs/final_vcfs/final_subsetted/
    statsdir=${sppdir}final_outputs/final_stats/
    workdir=${sppdir}final_outputs/final_vcfs/


#Extract allele frequencies for final imputed (TLR) vcfs, all individuals.
    #No need to compare these to each other, they are all the same and only calculated just in case there's a need.
    echo "Extracting Allele Frequencies, all individuals"
    for file in ${finaldir}*_${run}_final.vcf.gz
    do
        base=$(basename ${file} _${run}_final.vcf.gz)
        echo ""; echo "Calculating AF for ${base}"
        bcftools +fill-tags ${file} \
            -O z -o ${workdir}${base}_${run}_final.vcf.gz \
            -- -t AN,AC,AF 

        echo "Indexing ${base}"
        bcftools index ${workdir}${base}_${run}_final.vcf.gz --threads 16

        echo "Extracting AF for ${base}..."
        bcftools query -f '%CHROM %POS %AC %AN %AF\n' \
            -R ${tlr_bed} -H \
            ${workdir}${base}_${run}_final.vcf.gz \
            -o ${statsdir}${base}_AF.txt
    done

#Extract allele frequencies for subsetted populations
    #These populations should have been subsetted in 10_4_subset_vcf.sh
    #This will be done on the final chosen filter and impute parameters:
        #5x, 100ne, 0.1 max missingness
    echo ""; echo "Extracting Allele Frequencies, subsetted populations"
    for file in ${finaldir}final_subsetted/*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Calculating AF for ${base}"
        bcftools +fill-tags ${file} \
            -O z -o ${workdir}final_subsetted/${base}.vcf.gz \
            -- -t AN,AC,AF 

        echo "Indexing ${base}"
        bcftools index ${workdir}final_subsetted/${base}.vcf.gz --threads 16

        echo "Extracting AF for ${base}..."
        bcftools query -f '%CHROM %POS %AC %AN %AF\n' \
            -R ${tlr_bed} -H \
            ${workdir}final_subsetted/${base}.vcf.gz \
            -o ${statsdir}subsetted/${base}_AF.txt
    done

#Downloading stats to home computer. Navigate to final_outputs folder.
#mkdir -p final_stats/subsetted/
#rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/* final_stats/

echo "Script has finished"