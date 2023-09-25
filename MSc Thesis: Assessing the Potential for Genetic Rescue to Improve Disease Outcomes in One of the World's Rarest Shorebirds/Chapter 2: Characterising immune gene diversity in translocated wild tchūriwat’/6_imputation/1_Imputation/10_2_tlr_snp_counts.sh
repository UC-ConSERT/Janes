#!/bin/bash -e
set -e

# 25 April 2023

# Olivia Janes
# SNP counts for within each TLR, for each missingness filtered final vcf
# From: tuturuatu_imputation

## Environment: N/A

# Setting up:
    sppdir=~/data/tuturuatu_all_vcf/
    finaldir=${sppdir}impute/beagle_imputations/filtered/miss_filter_trial/
        #Directory holding the final phased, imputed and filtered vcf
    run=low_cov

    #Setting up directories for use
    mkdir -p ${sppdir}final_outputs/tlr_snp_counts/
    outdir=${sppdir}final_outputs/tlr_snp_counts/


#Creating the TLR file list
    tlrlist=${sppdir}final_outputs/tlr_list.txt
    echo "jcf7180002669510:1-2191 TLR1A" > ${tlrlist}
    echo "jcf7180002696225:79300-80387 TLR1B" >> ${tlrlist}
    echo "jcf7180002687310:7857-10285 TLR2A" >> ${tlrlist}
    echo "jcf7180002686685:232999-235155 TLR2B" >> ${tlrlist}
    echo "jcf7180002693511:279461-279884 TLR3_1" >> ${tlrlist}
    echo "jcf7180002693511:284293-286123 TLR3_2" >> ${tlrlist}
    echo "jcf7180002693511:286306-286420 TLR3_3" >> ${tlrlist}
    echo "jcf7180002688524:101290-103549 TLR4" >> ${tlrlist}
    echo "jcf7180002688267:101472-103857 TLR5" >> ${tlrlist}
    echo "jcf7180002696332:371554-374695 TLR7" >> ${tlrlist}
    echo "jcf7180002694481:2000-3649 TLR21" >> ${tlrlist}
    echo "jcf7180002693589:2103-0 TLR15_part" >> ${tlrlist}
        #Please note: TLR15_part is reversed


#Counting SNPs
    while read tlr_pos tlr_name
        #Cycle through TLRs listed in ${tlrlist}
    do
        echo ""; echo "Working on ${tlr_name}"
        
        
        for i in {0.9,0.8,0.5,0}
        do
            echo ""; echo "Outputting counts for ${i} missingness"
            
            tlrsnpfile=${outdir}${i}miss_SNP_counts.txt
            file=${finaldir}Tuturuatu_VariantCalls_${i}miss_low_cov_final.vcf.gz

            ## TLR SNP counts ##
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -r ${tlr_pos} -f '%POS\n' ${file} | wc -l)
            echo "${tlr_name},${snp}" >> ${tlrsnpfile}
        
        done
    done < ${tlrlist}

echo ""; echo "Script is complete."