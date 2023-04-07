#!/bin/bash -e 

#31 Mar 2023
#Molly Magid adapted by Olivia Janes
#Compiling all SNP counts from filtering, to compare between filtering methods.

sppdir=~/data/tuturuatu_all_rm_bad/
bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
sbiasdir=${sppdir}bcf/filter_strand_bias/
#Define location of tlr_regions.bed file in script. Should be in bcf/


echo "##### Before filtering SNP count" >> ${bcfdir}stats/SNP_counts.txt

cd ${bcfdir}

# Prefiltered SNP counts
#For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
for file in ${bcfdir}*VariantCalls_concat.bcf
do
    base=$(basename ${file} .bcf)
    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
    wait
    wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    echo ""
done

# LD filtered SNP counts
echo "##### LD Filtered SNP counts" >> ${bcfdir}stats/SNP_counts.txt

cd ${filterdir}LD_filter/

for file in ${filterdir}LD_filter/*.bcf
do
    base=$(basename ${file} .bcf)
    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
    wait
    wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    echo ""
done

# no LD filtered SNP counts
echo "##### no LD Filtered SNP counts" >> ${bcfdir}stats/SNP_counts.txt

cd ${filterdir}noLD/

for file in ${filterdir}noLD/*.bcf
do
    base=$(basename ${file} .bcf)
    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
    wait
    wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    echo ""
done

# Imputation filtered SNP counts
echo "##### Imputation Filtered SNP counts" >> ${bcfdir}stats/SNP_counts.txt

cd ${filterdir}impute/

for file in ${filterdir}impute/*.bcf
do
    base=$(basename ${file} .bcf)
    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
    wait
    wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    echo ""
done

# Strand bias filtered SNP counts
echo "##### Strand Bias Filtered SNP counts" >> ${bcfdir}stats/SNP_counts.txt

cd ${sbiasdir}

for file in ${sbiasdir}*.vcf.gz
do
    base=$(basename ${file} .vcf.gz)
    echo ${base} >> ${bcfdir}stats/SNP_counts.txt
    wait
    wc -l ${file} >> ${bcfdir}stats/SNP_counts.txt
    echo ""
done

echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/SNP_counts.txt"