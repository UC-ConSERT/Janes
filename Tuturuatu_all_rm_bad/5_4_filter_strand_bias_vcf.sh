#!/bin/bash -e

#19 Mar 2023
#Olivia Janes, from Molly Magid
#Filtering the filtered bcfs, to remove any sites where strand bias adjusted phred score is <60

sppdir=~/data/tuturuatu_all_rm_bad/
## Edit to be run specific

filterdir=${sppdir}bcf/filter_trial/

# Making a directory to hold the strand bias filtered bcfs.
mkdir -p ${sppdir}bcf/filter_strand_bias/vcf/
sbiasdir=${sppdir}bcf/filter_strand_bias/vcf/


# Filter bcf files for strand bias.
    #This sets individual sites with SP <60 to "."

for file in ${filterdir}noLD/vcf/*.vcf.gz
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .vcf.gz)
    bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.vcf
done


echo "Filtering for strand bias is complete."

