#!/bin/bash -e 

#12 Sep 2022
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu/
filterdir=${sppdir}bcf/filter_trial/
#Define location of tlr_regions.bed file in script


>>"COMMENTS"
#This is for a vcf file: not indexed or compressed. However, may be working with a bcf file (already compressed by auto, indexed in 5_1....sh).
#Change to >>?
for file in ${bcfdir}*.vcf 
do 
    base=$(basename $file .vcf) 
    bgzip $file 
    bcftools index ${base}.vcf.gz 
    echo $base > TLR_SNP_counts.txt 
    bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l > TLR_SNP_counts.txt 
done 

COMMENTS

echo "##### Before filtering SNP count" >> TLR_SNP_counts.txt

#For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
for file in ${bcfdir}*VariantCalls_concat.bcf 
do 
    base=$(basename $file .bcf)  
    echo ${base} >> TLR_SNP_counts.txt 
    bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> TLR_SNP_counts.txt 
done

echo "##### LD Filtered SNP counts" >> TLR_SNP_counts.txt

for file in ${filterdir}LD_filter/*.bcf 
do 
base=$(basename $file .bcf) 
echo ${base} >> TLR_SNP_counts.txt 
bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> TLR_SNP_counts.txt 
done 

echo "##### no LD Filtered SNP counts" >> TLR_SNP_counts.txt

for file in ${filterdir}noLD/*.bcf 
do 
base=$(basename $file .bcf) 
echo ${base} >> TLR_SNP_counts.txt 
bcftools query -R tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> TLR_SNP_counts.txt 
done 

echo "SNP counting is complete. Yay!"