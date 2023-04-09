#!/bin/bash -e 

#31 Mar 2023
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu_all/
bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
sbiasdir=${sppdir}bcf/filter_strand_bias/
#Define location of tlr_regions.bed file in script. Should be in bcf/

# Prefiltered SNP counts
#For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script.
    cd ${bcfdir}
    for file in ${bcfdir}*VariantCalls_concat.bcf
    do
        base=$(basename ${file} .bcf)
    #    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
        wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    done

# LD filtered SNP counts
    cd ${filterdir}LD_filter/

    for file in ${filterdir}LD_filter/*.bcf
    do
        base=$(basename ${file} .bcf)
    #    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
        wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    done

# no LD filtered SNP counts
    cd ${filterdir}noLD/

    for file in ${filterdir}noLD/*.bcf
    do
        base=$(basename ${file} .bcf)
    #    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
        wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    done

# Imputation filtered SNP counts
    cd ${filterdir}impute/

    for file in ${filterdir}impute/*.bcf
    do
        base=$(basename ${file} .bcf)
    #    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
        wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    done

# Strand bias filtered SNP counts
    cd ${sbiasdir}

    for file in ${sbiasdir}*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
    #    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
        wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    done

echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/SNP_counts.txt"