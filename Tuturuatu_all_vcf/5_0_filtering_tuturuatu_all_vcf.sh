#!/bin/bash -e 


#09 April 2023
#Molly Magid adapted by Olivia Janes
#Filtering trials for Tuturuatu variant calls, from [Molly's github](https://github.com/UC-ConSERT/Magid_et_al/blob/main/3_filtering.sh)

sppdir=~/data/tuturuatu_all_vcf/
    ## Must be edited to be run specific

bcfdir=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

mkdir -p ${bcfdir}filter_trial/ 
mkdir -p ${bcfdir}filter_trial/noLD/ ${bcfdir}filter_trial/LD_filter/ ${bcfdir}/stats/

vcf_out=${bcfdir}filter_trial/
noLD=${vcf_out}noLD/
LD=${vcf_out}LD_filter/


# First, the variant calls bcf must be transformed into a gzipped vcf and indexed
base=$(basename ${bcfdir}*concat.bcf .bcf)
echo "Converting Variant Calls bcf to vcf.gz format"
    bcftools view ${bcfdir}*concat.bcf -O z -o ${bcfdir}${base}.vcf.gz --threads 16
    bcftools index ${bcfdir}${base}.vcf.gz --threads 16
echo ""

#for loop to filter file with different values for parameters including
#missingness, depth, and GQ
for vcf in ${bcfdir}*_concat.vcf.gz
do
    base=$(basename ${vcf} _concat.vcf.gz)
    for i in {4..5} #filtering files for 4x and 5x depth
    do
        echo "Filtering SNPs for ${base}...." 
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_noMinGQ.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all &
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_noMinGQ.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all &
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ10.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all &
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ10.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all &
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ20.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all &
        vcftools --gzvcf ${vcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ20.vcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode \
            --recode-INFO-all
    done;
done


# Removing '.recode.vcf' from noLD file names.
echo "Renaming noLD files to remove '.recode.vcf'"
for vcf in ${noLD}*recode.vcf
do
    echo ""
        echo "Renaming ${vcf}"
    name=$(basename ${vcf} .recode.vcf)
    mv -i ${vcf} ${noLD}${name}
done 


echo "Filtering for Linkage parameters..."
#for loop to filter previous filtered files for linkage
for vcf in ${noLD}*.vcf
do
  base=$(basename ${vcf} .vcf)
   echo "Running light LD pruning at 0.8 for ${base}...."
   bcftools +prune \
        -m 0.8 \
        -w 1000 \
        -O v \
        -o ${LD}${base}_0.8LD.vcf.gz \
        ${vcf} &
    echo "Running moderate LD pruning at 0.6 for ${base}...."
    bcftools +prune \
        -m 0.6 \
        -w 1000 \
        -O v \
        -o ${LD}${base}_0.6LD.vcf.gz \
        ${vcf} &
    echo "Running strong LD pruning at 0.4 for ${base}...."
    bcftools +prune \
        -m 0.4 \
        -w 1000 \
        -O v \
        -o ${LD}${base}_0.4LD.vcf.gz \
        ${vcf}
done

# Convert no LD filter files to vcf.gz format and index
for vcf in ${noLD}*vcf
do
	base=$(basename ${vcf} .vcf)
    echo "Converting ${vcf} to vcf.gz format"
	bcftools view ${vcf} -O z -o ${noLD}${base}.vcf.gz --threads 16
    echo "Indexing ${vcf}"
    bcftools index ${noLD}${base}.vcf.gz --threads 16
	echo ""
done

# Index LD filter files (already gzipped)
for vcf in ${LD}*vcf.gz
do
    echo "Indexing ${vcf}"
    bcftools index ${vcf} --threads 16
	echo ""
done

echo "Script has finished running. Now whether it worked or not is another question..."