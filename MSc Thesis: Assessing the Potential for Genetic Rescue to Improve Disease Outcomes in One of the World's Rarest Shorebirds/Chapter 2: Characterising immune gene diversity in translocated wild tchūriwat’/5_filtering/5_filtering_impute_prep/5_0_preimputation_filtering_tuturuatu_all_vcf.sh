#!/bin/bash -e 

# 31 Aug 2022

# Olivia Janes
# Filtering trials for Tuturuatu variant calls for imputation
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    sppdir=~/data/tuturuatu_all_vcf/
    work=${sppdir}bcf/
        #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

    mkdir -p ${work}filter_trial/ 
    mkdir -p ${work}/stats/ ${work}filter_trial/impute/
    mkdir -p ${work}/stats/impute_stats/ ${work}filter_trial/impute/intermediate_filters/

    impdir=${work}filter_trial/impute/

# First, the variant calls bcf must be transformed into a gzipped vcf and indexed, if not already done
    base=$(basename ${bcfdir}*concat.bcf .bcf)
    echo "Converting Variant Calls bcf to vcf.gz format"
        bcftools view ${bcfdir}*concat.bcf -O z -o ${bcfdir}${base}.vcf.gz --threads 16
        bcftools index ${bcfdir}${base}.vcf.gz --threads 16
    echo ""


# Filtering for imputation, with various filter trials including:
    for vcf in ${work}*_concat.vcf.gz
    do
        base=$(basename ${vcf} _concat.vcf.gz)
        for dp in {0,4,5}
        do
            echo "Filtering SNPs for ${base}, depth=${dp}x, no missingness filter...."
            # no missingness filter, 10 minGQ, bi-allelic only, no MAF filter
            vcftools --gzvcf ${vcf} \
                --out ${impdir}intermediate_filters/${base}_${dp}x.vcf \
                --minDP ${dp} \
                --max-alleles 2 \
                --maxDP 50  \
                --minGQ 10 \
                --minQ 20 \
                --remove-indels \
                --remove-filtered-all \
                --recode \
                --recode-INFO-all &
        done
    done
    wait

# Removing '.recode.vcf' from imputation file names.
    echo "Renaming impute filter files to remove '.recode.vcf'"
    for vcf in ${impdir}intermediate_filters/*.vcf.recode*
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
    for file in ${impdir}intermediate_filters/*.vcf
    do
        base=$(basename ${file} .vcf)
        echo "Filtering ${base} for SP <60"
        bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${impdir}intermediate_filters/${base}_0.6SP.vcf
    done
    wait


# Convert imputation filter files to vcf.gz format and index
    for vcf in ${impdir}intermediate_filters/*SP.vcf
    do
        base=$(basename ${vcf} .vcf)
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${impdir}${base}.vcf.gz --threads 16
        wait
        echo "Indexing ${vcf}"
        bcftools index ${impdir}${base}.vcf.gz --threads 16
        echo ""
    done


wait
echo ""; echo "Script has finished running. Now whether it worked or not is another question..."
echo "Find final filtered files ready for imputation prep at ${impdir}"