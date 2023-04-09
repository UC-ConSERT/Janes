#!/bin/bash -e 

#12 Sep 2022
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu_all_vcf/

bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
#Define location of tlr_regions.bed file in script. Should be in bcf/

echo "Beginning SNP counting."

echo "##### Before filtering SNP count" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${bcfdir}

#For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
for file in ${bcfdir}*VariantCalls_concat.vcf.gz
do
    base=$(basename ${file} .vcf.gz)
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done

# Filtered files should be indexed in 5_0_filtering.sh
echo "##### LD Filtered SNP counts" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${filterdir}LD_filter/

for file in ${filterdir}LD_filter/*.vcf.gz
do
    base=$(basename ${file} .vcf.gz)
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done


echo "##### no LD Filtered SNP counts" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${filterdir}noLD/

for file in ${filterdir}noLD/*.vcf.gz
do
    base=$(basename ${file} .vcf.gz)
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done

echo "##### Imputation Filtered SNP counts" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${filterdir}impute/

for file in ${filterdir}impute/*.vcf.gz
do
    base=$(basename ${file} .vcf.gz)
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.vcf.gz | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done


echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/TLR_SNP_counts.txt"