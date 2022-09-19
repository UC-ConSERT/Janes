#!/bin/bash -e 

#12 Sep 2022
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu/
filterdir=${sppdir}bcf/filter_trial/


>>"COMMENTS"
#This is for a vcf file: not indexed or compressed. However, may be working with a bcf file (already compressed by auto, indexed in 5_1....sh).
for file in ${bcfdir}*.vcf 
do 
    base=$(basename $file .vcf) 
    bgzip $file 
    bcftools index ${base}.vcf.gz 
    echo $base > TLR_SNP_counts.txt 
    bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l > TLR_SNP_counts.txt 
done 

COMMENTS

#For bcf file that was indexed in previous script (5_1_prefiltering_stats)
for file in ${bcfdir}*VariantCalls_concat.bcf 
do 
    base=$(basename $file _concat.bcf)  
    echo $base > TLR_SNP_counts.txt 
    bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l > TLR_SNP_counts.txt 
done

for file in *.vcf 
do 
base=$(basename $file .vcf) 
bgzip $file 
bcftools index ${base}.vcf.gz 
echo $base > TLR_SNP_counts.txt 
bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l > TLR_SNP_counts.txt 
done 