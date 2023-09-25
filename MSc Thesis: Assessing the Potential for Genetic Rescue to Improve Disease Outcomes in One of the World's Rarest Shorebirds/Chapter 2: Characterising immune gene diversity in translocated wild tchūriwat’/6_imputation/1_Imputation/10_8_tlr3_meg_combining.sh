#!/bin/bash -e
set -e

# 02 May 2023

# Olivia Janes
# Merging TLR3_1, TLR3_2 and TLR3_3 meg files together for downstream analyses.

## DOESN'T WORK - it doesn't join the files for all 3 tlr3s together nicely - leaves a line break ##

# From: tuturuatu_imputation

## Environment: N/A

# Setting up:
    sppdir=~/data/tuturuatu_all_vcf/

    mkdir -p ${sppdir}final_outputs/meg/meg_test/final_tlr3/ ${sppdir}final_outputs/meg/tlr3_combined/
    megdir=${sppdir}final_outputs/meg/meg_test/
    tlr3dir=${sppdir}final_outputs/meg/tlr3_combined/


for miss in {0,0.5,0.8,0.9}
do
    for i in {1..84}
    do
        echo "Processing individual: ${i}-1"
        echo "#${i}_1" > ${miss}miss_TLR3_indv_${i}_1.meg

        for file in ${megdir}*_${miss}miss_TLR3*
        do
            echo "Processing file: ${file}"
            awk -v start=$(awk "/^#${i}-1/{print NR+1; exit}" ${file}) "/^#${i}-2/{exit} NR>=start" ${file} >> ${miss}miss_TLR3_indv_${i}_1.meg

        done
    done
done

for miss in {0,0.5,0.8,0.9}
do
    for i in {1..84}
    do
        echo "Processing individual: ${i}-2"
        echo "#${i}_2" > ${miss}miss_TLR3_indv_${i}_2.meg

        for file in ${megdir}*_${miss}miss_TLR3*
        do
            echo "Processing file: ${file}"

            c=$(expr ${i} + 1)
            awk -v start=$(awk "/^#${i}-2/{print NR+1; exit}" ${file}) "/^#${c}-1/{exit} NR>=start" ${file} >> ${miss}miss_TLR3_indv_${i}_2.meg

        done
    done
done


for miss in {0,0.5,0.8,0.9}
do
    echo "Processing missingness ${miss}"
    for file in ${megdir}${miss}miss_TLR3_indv*
    do
        echo "Processing file: ${file}"
        cat ${file} >> ${megdir}final_tlr3/tlr_haplotypes_${miss}miss_TLR3.meg
    done
done

echo "Moving final files"
mv ${megdir}final_tlr3/* ${tlr3dir}

echo ""; echo "Script complete"
