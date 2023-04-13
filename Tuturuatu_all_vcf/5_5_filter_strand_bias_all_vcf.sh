#!/bin/bash -e

#07 March 2023
#Olivia Janes
#Filtering the filtered bcfs, to remove any sites where strand bias adjusted phred score is <60
#All filtered files must be in vcf.gz format.

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific

filterdir=${sppdir}bcf/filter_trial/

# Making a directory to hold the strand bias filtered bcfs.
mkdir -p ${sppdir}bcf/filter_strand_bias/
vcfdir=${sppdir}bcf/filter_trial/
sbiasdir=${sppdir}bcf/filter_strand_bias/


# Filter vcf files for strand bias.
    #This sets individual sites with SP <60 to "."
    #Filtering noLD files
    for file in ${vcfdir}noLD/*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Filtering ${base} for SP <60"
        bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.vcf
    done

    #Filtering LD files
    for file in ${vcfdir}LD_filter/*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Filtering ${base} for SP <60"
        bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.vcf
    done

# Zipping and indexing strand bias filtered vcfs.
    for file in ${sbiasdir}*.vcf
    do
        echo "Zipping and indexing ${file}"
        base=$(basename ${file} .vcf)
        bcftools view ${file} -O z -o ${sbiasdir}${base}.vcf.gz --threads 16
        bcftools index ${sbiasdir}${base}.vcf.gz --threads 16
    done


echo "Filtering for strand bias is complete."
echo "Strand bias filtered files can be found at ${sbiasdir}...vcf.gz"

