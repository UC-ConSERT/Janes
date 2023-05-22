#!/bin/bash -e 
set -e

#09 May 2023

# Olivia Janes
# Genotypes for WGS.

#Environment: samtools  

##  Needs to be edited to be run specific   ##
    sppdir=~/data/tuturuatu_all_vcf/
    run=low_cov
    impvcf=Tuturuatu_VariantCalls_0.9miss_low_cov_final.vcf.gz
    wgsvcf=Tuturuatu_VariantCalls_5x_coverage_0.2site_missing_MinGQ10.vcf.gz
    sampleorder=~/data/tuturuatu_all/nodup_bam/Tuturuatu_bam_list.txt
    tlr_regions=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
        #Define location of TLR regions bed file

#Setting directories
    mkdir -p ${sppdir}final_outputs/final_stats/wgs_genotypes/ ${sppdir}final_outputs/final_vcfs/whole_genome/
    statsdir=${sppdir}final_outputs/final_stats/wgs_genotypes/
    impdir=${sppdir}final_outputs/final_vcfs/
    wgsdir=${sppdir}bcf/filter_trial/noLD/
    wgsout=${sppdir}final_outputs/final_vcfs/whole_genome/

#Reorder the final vcf
    echo "Reordering the final imputed vcf so that samples are in the same order as the non-imp vcf"
    bcftools view -S ${sampleorder} ${impdir}${impvcf} -O z --threads 16 > ${impdir}temp.vcf.gz
    bcftools sort ${impdir}temp.vcf.gz -O z > ${impdir}Tuturuatu_VariantCalls_0.9miss_low_cov_final_reordered.vcf.gz
    bcftools index --threads 16 ${impdir}Tuturuatu_VariantCalls_0.9miss_low_cov_final_reordered.vcf.gz
    rm ${impdir}temp.vcf.gz

#Remove the TLR contigs from the WGS vcf, to be replaced with the imputed TLR contigs.
    echo ""; echo "Converting TLR bed file into chromosome file"
    cut ${tlr_regions} -f1 > ${sppdir}bcf/chroms.txt
    echo ""; echo "Removing TLR contigs from WGS vcf."
    wgsbase=$(basename ${wgsdir}${wgsvcf} .vcf.gz)
    bcftools view -T ^${sppdir}bcf/chroms.txt ${wgsdir}${wgsvcf} -O z -o ${wgsout}${wgsbase}_no_tlrs.vcf.gz --threads 16
    vcftools --gzvcf ${wgsdir}${wgsvcf} \
        --out ${wgsout}${wgsbase}_no_tlrs \
        --exclude-bed ${tlr_regions} \
        --recode \
        --recode-INFO-all

#Renaming and compressing
    mv ${wgsout}${wgsbase}_no_tlrs.recode.vcf ${wgsout}${wgsbase}_no_tlrs.vcf
    
    echo "Converting merged vcf to vcf.gz format"
    bcftools view ${wgsout}${wgsbase}_no_tlrs.vcf -O z -o ${wgsout}${wgsbase}_no_tlrs.vcf.gz --threads 16
    echo "Indexing ${wgsbase}_no_tlrs.vcf.gz"
    bcftools index ${wgsout}${wgsbase}_no_tlrs.vcf.gz --threads 16
    echo ""

#Concatenating the files together
    echo ""; echo "Concatenating wgs vcf and imputed TLR vcf together."

    echo "${wgsout}${wgsbase}_no_tlrs.vcf.gz"; echo "${wgsout}${wgsbase}_no_tlrs.vcf.gz" > ${wgsout}merge_list.txt
    echo "${impdir}Tuturuatu_VariantCalls_0.9miss_low_cov_final_reordered.vcf.gz"
    echo "${impdir}Tuturuatu_VariantCalls_0.9miss_low_cov_final_reordered.vcf.gz" >> ${wgsout}merge_list.txt

    echo "Concatenating now"
    bcftools concat -O z --threads 16 -f ${wgsout}merge_list.txt -o ${wgsout}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz
    echo "Indexing ${wgsout}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz"
    bcftools index -f --threads 16 ${wgsout}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz

#Investigating WGS Genotypes: Extract the Genotypes out of the merged file.
#This is not necessary if doing a plink PCA, only needed for an adegenet PCA in R.
    echo ""; echo "Extracting whole genome genotypes"

    #Extract information on the Genotypes at each SNP for each individual
        bcftools query --format '%CHROM\t%POS[\t%TGT]\n' ${wgsout}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz > ${statsdir}wgs_tlrimp_merged_genotypes.txt

    #Extract the headers to add to the above
        bcftools view -h ${wgsout}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz | tail -n 1 > ${statsdir}wgs_tlrimp_merged_genotypes_header.txt
    #Download these and extract into a spreadsheet to analyse

    done
    echo "TLR Genotype files can be found at: ${statsdir}wgs_tlrimp_merged_genotypes..."



echo ""; echo "To download all of the stats, navigate to the right directory on your desktop: ~/Documents/Tuturuatu_resources/tuturuatu_all_vcf/final_outputs/final_stats/"
echo "Enter code (edited for the right run):"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/wgs_genotypes ./"

echo ""
echo "Whole genome genotype script is complete."       
