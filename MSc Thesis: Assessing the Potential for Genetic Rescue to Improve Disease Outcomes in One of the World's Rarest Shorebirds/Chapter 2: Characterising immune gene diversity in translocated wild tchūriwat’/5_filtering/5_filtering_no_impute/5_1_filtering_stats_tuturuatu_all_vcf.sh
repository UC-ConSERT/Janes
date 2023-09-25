#!/bin/bash -e 

# 12 Sep 2022

# Adapted by Olivia Janes from https://speciationgenomics.github.io/filtering_vcfs/
# Pre filtering statistics calculated to guide vcf filtering and compare to post-filtering statistics.
# From: tuturuatu_all_vcf

## Environment: samtools

# Setting up
    sppdir=~/data/tuturuatu_all_vcf/
    bcfdir=${sppdir}bcf/
        #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)
    noLD=${bcfdir}filter_trial/noLD/
    LD=${bcfdir}filter_trial/LD_filter/

    mkdir -p ${bcfdir}stats/
    mkdir -p ${bcfdir}stats/stats_raw_files/

# Calculating stats for prefilter variant calls vcf
    varcallbcf=${bcfdir}Tuturuatu_VariantCalls_concat.vcf.gz
        #output concatenated file from variant calling
        ##### Must be edited to be sample specific #####

    outfile=${bcfdir}stats/stats_raw_files/Tuturuatu_VariantCalls_prefilter
        #defining the output file prefix

    #Calculate stats
    echo "Calculating pre-filter stats"

    #Calculate allele frequency for each variant
    vcftools --gzvcf ${varcallbcf} --freq2 --out ${outfile} --max-alleles 2
        ##### Not incl in Mollys #####

    #Calculate site quality
    vcftools --gzvcf ${varcallbcf} --site-quality --out ${outfile}
        ##### Not incl in Mollys #####

    #Calculate mean depth per individual
    vcftools --gzvcf ${varcallbcf} --depth --out ${outfile}

    #Calculate site depth
    vcftools --gzvcf ${varcallbcf} --site-depth --out ${outfile}
        
    #Calculate proportion of missing data per individual
    vcftools --gzvcf ${varcallbcf} --missing-indv --out ${outfile}

    #Calculate proportion of missing data per site
    vcftools --gzvcf ${varcallbcf} --missing-site --out ${outfile}
        
    #Calculate heterozygosity
    vcftools --gzvcf ${varcallbcf} --het --out ${outfile}


# Calculating statistics for no linkage filtered files
    for file in ${noLD}*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Calculating depth for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --site-depth &
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --missing-site &
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --het
    done

# Calculating statistics for linkage filtered files
    for file in ${LD}*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Calculating depth for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --site-depth &
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --missing-site &
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${bcfdir}stats/stats_raw_files/${base} \
            --het
    done

wait
#To view these stats in R, see link above.
echo "Prefilter stats has finished."
