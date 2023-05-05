#!/bin/bash -e 
set -e

#04 May 2023
#Olivia Janes
#Calculating Runs of Homozygosity (ROH) using filtered vcf file.

#Environment: plink

sppdir=~/data/tuturuatu_roh/
    ## Must be edited to be run specific

#Defining directories
    filterdir=${sppdir}bcf/filter_trial/
        #Directory containing the filtered vcf(s)
    mkdir -p ${sppdir}roh/
    rohdir=${sppdir}roh/

#Calculating ROH with Plink v1.90. Using sliding window sizes of 300kb and 1000kb.
    for dp in {5,6,8}
    do
        for miss in {0.1,0.2}
        do
            for file in ${filterdir}*${dp}x_coverage_${miss}site_missing_0.6SP.vcf.gz
            do
                base=$(basename ${file} .vcf.gz)
                echo ""; echo "Calculating ROH for ${base}.vcf.gz"

                echo "300kb sliding window size"
                plink --vcf ${file} --homozyg \
                    --out ${rohdir}tuturuatu_roh_${dp}x_${miss}miss_300kb_window \
                    --homozyg-kb 300 \
                    --homozyg-snp 50 \
                    --homozyg-window-snp 50 \
                    --homozyg-density 50 \
                    --homozyg-gap 1000 \
                    --homozyg-window-het 3 \
                    --homozyg-window-missing 5 \
                    --allow-extra-chr \
                    --double-id
                
                
                echo "1000kb sliding window size"
                plink --vcf ${file} --homozyg \
                    --out ${rohdir}tuturuatu_roh_${dp}x_${miss}miss_1000kb_window \
                    --homozyg-kb 1000 \
                    --homozyg-snp 50 \
                    --homozyg-window-snp 50 \
                    --homozyg-density 50 \
                    --homozyg-gap 1000 \
                    --homozyg-window-het 3 \
                    --homozyg-window-missing 5 \
                    --allow-extra-chr \
                    --double-id
            done
        done
    done

echo "Runs of homozygosity script has finished running, homozygositily."