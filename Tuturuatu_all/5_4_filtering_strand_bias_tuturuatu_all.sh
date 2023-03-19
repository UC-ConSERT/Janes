#!/bin/bash -e

#19 Mar 2023
#Olivia Janes, from Molly Magid
#Filtering the filtered bcfs, to remove any sites where strand bias adjusted phred score is <60

sppdir=~/data/tuturuatu_all/
## Edit to be run specific

filterdir=${sppdir}bcf/filter_trial/

# Making a directory to hold the strand bias filtered bcfs.
mkdir -p ${sppdir}bcf/filter_strand_bias/
sbiasdir=${sppdir}bcf/filter_strand_bias/


# Filter bcf files for strand bias.
    #This sets individual sites with SP <60 to "."

for file in ${filterdir}noLD/*.bcf
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .bcf)
    bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.bcf
done

for file in ${filterdir}LD_filter/*.bcf
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .bcf)
    bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.bcf
done

echo "Filtering for strand bias is complete."

