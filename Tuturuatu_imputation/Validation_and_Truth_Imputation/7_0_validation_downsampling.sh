#!/bin/bash -e 


#13 April 2023
#Olivia Janes
#Downsampling a selection of individuals to ~4x coverage to validate the imputation. 


sppdir=~/data/tuturuatu_all/
sppdir2=~/data/tuturuatu_all_vcf/

mkdir -p ${sppdir2}impute/
mkdir -p ${sppdir2}impute/validation/
mkdir -p ${sppdir2}impute/validation/validation_bams/

ref=${sppdir}ref_genome/Maui_merged_assembly.fa
         #reference genome for alignment
         ##### Must be edited to be sample specific #####
nodupbamdir=${sppdir}nodup_bam/
        #directory that holds the merged bam files that have been sorted, fixed and had duplicates removed.
valbamdir=${sppdir2}impute/validation/validation_bams/
        #a directory to hold the chunked bam files

<<"COMMENTS_Deciding_on_individuals"
#Decide on individuals to downsample by selecting a mix of captive and wild, with rare SNPs, or with many ALT SNPs:
    #Extract information on the Genotypes at each TLR SNP for each individual
        bcftools query -R ~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' \
            ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_4x_0.6SP.vcf.gz >> ~/data/tuturuatu_all_vcf/scripts/tlr_genotypes_4x.txt
    #Extract the headers to add to the above
        bcftools view -h ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_4x_0.6SP.vcf.gz \
            | tail -n 1 >> ~/data/tuturuatu_all_vcf/scripts/tlr_genotypes_4x_header.txt
    #Download these and extract into a spreadsheet to analyse
COMMENTS_Deciding_on_individuals

# Define the list of files to downsample
file_list="A09_nodup.bam A11_nodup.bam B10_nodup.bam CR20_nodup.bam CT07_nodup.bam CT11_nodup.bam E10_nodup.bam \
    F09_nodup.bam I16468_nodup.bam I16476_nodup.bam"

#Downsample with picard or samtools
    for file in ${file_list}
    do
        echo "Downsampling ${file}"
        name=$(basename ${file} _nodup.bam)
        #picard DownsampleSam -I ${nodupbamdir}${file} -O ${valbamdir}${name}_downsampled.bam \
        #    -P 0.5 --CREATE_INDEX true
        samtools view -@ 16 -b -s 0.4 ${nodupbamdir}${file} -o ${valbamdir}${name}_downsampled.bam
        echo ""; echo ""
    done

# Indexing the downsampled bams
    for file in ${valbamdir}*.bam
    do
        echo "Indexing ${file}"
        samtools index -@ 16 -b ${file}
    done

echo "Script has finished. Find downsampled bams at: ${valbamdir}"
