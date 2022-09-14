#!/bin/sh

#12 Sep 2022
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu/
work=${sppdir}bcf/


for file in *.vcf 
do 
    base=$(basename $file .vcf) 
    bgzip $file 
    bcftools index ${base}.vcf.gz 
    echo $base > TLR_SNP_counts.txt 
    bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l > TLR_SNP_counts.txt 
done 