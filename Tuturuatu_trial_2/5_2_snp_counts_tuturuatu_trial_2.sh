#!/bin/bash -e 

#12 Sep 2022
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu_trial_2/
bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
#Define location of tlr_regions.bed file in script. Should be in bcf/


echo "##### Before filtering SNP count" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${bcfdir}

#For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
for file in ${bcfdir}*VariantCalls_concat.bcf
do
    base=$(basename ${file} .bcf)
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done


echo "##### LD Filtered SNP counts" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${filterdir}LD_filter/

for file in ${filterdir}LD_filter/*.bcf
do
    base=$(basename ${file} .bcf)
    bcftools index ${base}.bcf
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done


echo "##### no LD Filtered SNP counts" >> ${bcfdir}stats/TLR_SNP_counts.txt

cd ${filterdir}noLD/

for file in ${filterdir}noLD/*.bcf
do
    base=$(basename ${file} .bcf)
    bcftools index ${base}.bcf
    echo ${base} >> ${bcfdir}stats/TLR_SNP_counts.txt
    bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${bcfdir}stats/TLR_SNP_counts.txt
done

echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/TLR_SNP_counts.txt"