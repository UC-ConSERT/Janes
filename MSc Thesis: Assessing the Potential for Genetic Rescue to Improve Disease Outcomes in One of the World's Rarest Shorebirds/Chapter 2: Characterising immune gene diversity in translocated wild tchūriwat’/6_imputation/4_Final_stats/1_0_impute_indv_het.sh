#!/bin/bash -e 
set -e

# 06 May 2023

# Olivia Janes
# Running F stats for the imputed TLR regions to compare back to stats created for TLR regions, not imputed.
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    ##  Needs to be edited to be run specific   ##
        sppdir=~/data/tuturuatu_all_vcf/
        run=low_cov
        vcf=Tuturuatu_VariantCalls_4x_coverage_0.2site_missing_MinGQ10_0.6SP.vcf.gz
            #Chosen non-imputed vcf to compare imputed het to.

    #Setting directories
        mkdir -p ${sppdir}impute/stats/beagle_imp_stats/
        statsdir=${sppdir}impute/stats/beagle_imp_stats/
        mkdir -p ${sppdir}bcf/filter_strand_bias/tlr_subsetted/
        subsetdir=${sppdir}bcf/filter_strand_bias/tlr_subsetted/

#Calculate .het stats for imputed TLR contigs
    echo "Calculating filter stats"
    for file in ${sppdir}impute/beagle_imputations/filtered/*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}${base} \
            --het
        wait
        rm ${statsdir}${base}.log
    done

#Calculate .het stats for non-imputed TLR contigs
    #Extracting TLR contigs
        for vcf in ${sppdir}bcf/filter_strand_bias/${vcf}
        do
            base=$(basename ${vcf} .vcf.gz)
            # Extract the TLR contigs from the filtered vcf files
            # This is done in 3 files as bcftools only seems to be able to handle 4 "windows"/contigs at a time
            echo "Extracting TLR contigs"
            bcftools view --threads 16 ${vcf} -r jcf7180002669510,jcf7180002696225,jcf7180002687310,jcf7180002686685 \
                -O z -o ${subsetdir}${base}_1_tlr_contigs.vcf.gz
            bcftools view --threads 16 ${vcf} -r jcf7180002693511,jcf7180002688524,jcf7180002688267,jcf7180002696332 \
                -O z -o ${subsetdir}${base}_2_tlr_contigs.vcf.gz
            bcftools view --threads 16 ${vcf} -r jcf7180002694481,jcf7180002693589 \
                -O z -o ${subsetdir}${base}_3_tlr_contigs.vcf.gz
        done

    #Bringing the TLR contigs back into one file
        for subset_1 in ${subsetdir}*1_tlr_contigs.vcf.gz
        do
            base=$(basename ${subset_1} _1_tlr_contigs.vcf.gz)
            subset_2=${subsetdir}${base}_2_tlr_contigs.vcf.gz
            subset_3=${subsetdir}${base}_3_tlr_contigs.vcf.gz

            #Create a file list of files to merge
            echo ""; echo "Merging subsetted tlr contigs in:"
            echo "${subset_1}"; echo "${subset_1}" > ${subsetdir}merge_list.txt
            echo "${subset_2}"; echo "${subset_2}" >> ${subsetdir}merge_list.txt
            echo "${subset_3}"; echo "${subset_3}" >> ${subsetdir}merge_list.txt

            #Merge and index
            bcftools concat -O z --threads 16 -f ${subsetdir}merge_list.txt -o ${subsetdir}${base}_tlr_contigs_merged.vcf.gz
            echo "Indexing ${base}_tlr_contigs_merged.vcf.gz"
            bcftools index -f --threads 16 ${subsetdir}${base}_tlr_contigs_merged.vcf.gz
        done

    #Calculating het
        for file in ${subsetdir}*_tlr_contigs_merged.vcf.gz
        do
            base=$(basename ${file} _merged.vcf.gz)
            echo "Calculating individual heterozygosity for ${base}..."
            vcftools --gzvcf ${file} \
                --out ${statsdir}${base} \
                --het
            wait
            rm ${statsdir}${base}.log
        done       


echo "Script has finished"