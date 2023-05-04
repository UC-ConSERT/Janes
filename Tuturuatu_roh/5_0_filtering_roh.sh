#!/bin/bash -e 
set -e

#04 May 2023
#Molly Magid adapted by Olivia Janes
#Filtering trials for Tuturuatu variant calls, from [Molly's github](https://github.com/UC-ConSERT/Magid_et_al/blob/main/3_filtering.sh)
#Altered for filtering for ROH

sppdir=~/data/tuturuatu_roh/
    ## Must be edited to be run specific

bcfdir=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

mkdir -p ${bcfdir}filter_trial/ 
mkdir -p ${bcfdir}filter_trial/intermediate_filters/ ${bcfdir}/stats/

vcf_out=${bcfdir}filter_trial/

<<"COMMENTS"
# First, the variant calls bcf must be transformed into a gzipped vcf and indexed
base=$(basename ${bcfdir}*concat.bcf .bcf)
echo "Converting Variant Calls bcf to vcf.gz format"
    bcftools view ${bcfdir}*concat.bcf -O z -o ${bcfdir}${base}.vcf.gz --threads 16
    bcftools index ${bcfdir}${base}.vcf.gz --threads 16
echo ""
COMMENTS


#Filter file with different values of missingness and depth
    for vcf in ${bcfdir}*_concat.vcf.gz
    do
        base=$(basename ${vcf} _concat.vcf.gz)
        for i in {5,6,8} #filtering files for 5x, 6x and 8x depth
        do
            echo "Filtering SNPs for ${base}...." 
            vcftools --gzvcf ${vcf} \
                --out ${vcf_out}intermediate_filters/${base}_${i}x_coverage_0.1site_missing.vcf \
                --minDP ${i} \
                --maxDP 50  \
                --max-missing 0.9 \
                --maf 0.05 \
                --minQ 20 \
                --minGQ 10 \
                --max-alleles 2 \
                --remove-indels \
                --remove-filtered-all \
                --recode \
                --recode-INFO-all &
            vcftools --gzvcf ${vcf} \
                --out ${vcf_out}intermediate_filters/${base}_${i}x_coverage_0.2site_missing.vcf \
                --minDP ${i} \
                --maxDP 200  \
                --max-missing 0.8 \
                --maf 0.05 \
                --minQ 20 \
                --minGQ 10 \
                --max-alleles 2 \
                --remove-indels \
                --remove-filtered-all \
                --recode \
                --recode-INFO-all &
        done;
    done
    wait

# Removing '.recode.vcf' from filtered file names.
    echo "Renaming filter files to remove '.recode.vcf'"
    for vcf in ${vcf_out}intermediate_filters/*.vcf.recode*
    do
        echo ""
        echo "Renaming ${vcf}"
        name=$(echo ${vcf} | sed 's/.vcf.recode//g')
        #rename "s/${vcf}/${name}/g" ${vcf}
        mv -i ${vcf} ${name}

    done
    wait


# Filter vcf files for strand bias.
    #This sets individual sites with SP <60 to "."
    for file in ${vcf_out}intermediate_filters/*.vcf
    do
        base=$(basename ${file} .vcf)
        echo "Filtering ${base} for SP <60"
        bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${impdir}intermediate_filters/${base}_0.6SP.vcf
    done
    wait


# Convert imputation filter files to vcf.gz format and index
    for vcf in ${vcf_out}intermediate_filters/*SP.vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${vcf_out}${base}.vcf.gz --threads 16
        wait
        echo "Indexing ${vcf}"
        bcftools index ${vcf_out}${base}.vcf.gz --threads 16
        echo ""
    done


echo "Script has finished running. Now whether it worked or not is another question..."