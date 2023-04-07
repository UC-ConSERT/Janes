#!/bin/bash -e

#19 Mar 2023
#Olivia Janes, from Molly Magid
#Filtering the "chosen" filtered bcfs, to remove any sites where strand bias adjusted phred score is <60.
## This scripts uses bcftools query as opposed to bcftools +setGT.
##  This is to see if it works, and to see if it might fix the problem of +setGT setting all filtered sites to ".", 
##  as I'm not sure if this is problematic down the line. Count the outputs of this script as a backup to +setGT.

sppdir=~/data/tuturuatu_all/
## Edit to be run specific

filterdir=${sppdir}bcf/filter_trial/

# Making a directory to hold the strand bias filtered bcfs.
mkdir -p ${sppdir}bcf/filter_strand_bias_bcftools_query/
sbiasdir=${sppdir}bcf/filter_strand_bias_bcftools_query/
#####   If this is chosen script, edit these to remove the bcftools query part of the folder name


#Filtering chosen filtered file to export only sites where strand bias adjusted phred score is <60

for file in ${filterdir}noLD/*.bcf
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .bcf)
    bcftools query -i 'FORMAT/SP<60' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\n' \
        ${file} -o ${sbiasdir}${base}_0.6SP.bcf
done

for file in ${filterdir}LD_filter/*.bcf
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .bcf)
    bcftools query -i 'FORMAT/SP<60' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\n' \
        ${file} -o ${sbiasdir}${base}_0.6SP.bcf
done

#[%GT:PL:DP:SP:ADF:ADR:AD:GQ]
#Error: no such tag defined in the VCF header: INFO/FORMATn

<<"COMMENTS"
echo "GZipping and indexing bcf file."
bgzip ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.bcf
bcftools index ${sppdir}bcf_final/Tuturuatu_VariantCalls_final_variants.bcf.gz
echo "GZipping and indexing bcf file is complete."

COMMENTS

echo "Filtering for strand bias is complete."

