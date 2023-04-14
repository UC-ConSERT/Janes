#!/bin/bash -e 
set -e

#14 April 2023

#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu variant calling and preparing files for filtering
#Used downsampled files for imputation validation

##### Must be edited to be sample specific #####
sppdir=~/data/tuturuatu_all_vcf/
nodupbamdir=~/data/tuturuatu_all/nodup_bam/
        #Location of non-downsampled bam files to variant call alongside the downsampled bams
ref=~/data/tuturuatu_all/ref_genome/Maui_merged_assembly.fa
         #reference genome for alignment
        
#Preparing the directories
        valdir=${sppdir}impute/validation/
        mkdir -p ${valdir}chunks_tlr/ ${valdir}bcf_tlr/ ${nodupbamdir}not_var_calling ${valdir}validation_bams_tlr/
        valbamdir=${valdir}validation_bams/
                #directory that holds the downsampled bams.
        valbamtlrdir=${valdir}validation_bams_tlr/
                #directory that holds the downsampled bams, only tlr contig regions.
        scriptdir=~/data/general_scripts/
        chunksdir=${valdir}chunks_tlr/
                #a directory to hold the chunked bam files
        bcf_file=${valdir}bcf_tlr/
                #bcf file output
        species="Tuturuatu"

<<"COMMENTS"
#Remove the downsample original bams from the nodup folder to ensure they are not variant called
        # Define the list of files to downsample
        file_list="A09_nodup.bam A11_nodup.bam B10_nodup.bam CR20_nodup.bam CT07_nodup.bam CT11_nodup.bam E10_nodup.bam \
                F09_nodup.bam I16468_nodup.bam I16476_nodup.bam"
        echo "Moving files not to be variant called"
        for file in ${file_list}
        do
                mv ${nodupbamdir}${file}* ${nodupbamdir}not_var_calling/ || true
                        #|| true : this means that if this doesn't find a file to move, it won't stop the whole script
        done

#Subset the bam files so that they only contain the TLR contigs
        #Extracting TLR region for downsampled bam
        for file in ${valbamdir}*.bam
        do
                echo "Extracting TLR region for ${file}"
                base=$(basename ${file} .bam)
                samtools view -@ 16 -b -h ${file} -o ${valbamtlrdir}${base}_tlrs.bam -M -L ${valbamtlrdir}tlr_contigs.bed
        done
        
        #Extracting TLR region for all other samples
        for file in ${nodupbamdir}*.bam
        do
                echo "Extracting TLR region for ${file}"
                base=$(basename ${file} nodup.bam)
                samtools view -@ 16 -b -h ${file} -o ${valbamtlrdir}${base}_tlrs.bam -M -L ${valbamtlrdir}tlr_contigs.bed
        done
COMMENTS


#chunk bam files for mpileup
        echo ""; echo "Chunking files for mpileup"
        ls ${valbamtlrdir}*.bam > ${valbamtlrdir}${species}_bam_list.txt
        perl ${scriptdir}split_bamfiles_tasks.pl \
                -b ${valbamtlrdir}${species}_bam_list.txt \
                -g $ref -n 16 -o ${chunksdir} | parallel -j 16 {}

#run mpileup on chunks of bam files
        echo "Running mpileup on chunks of bam files"
        for ((i=1; i<=16; i++))
        do
                bcftools mpileup \
                        --threads 16 \
                        -f $ref \
                        -a AD,ADF,ADR,DP,SP \
                        -O b -o ${bcf_file}${species}_${i}_raw.bcf \
                        ${chunksdir}${i}/* &
        done
        wait
        echo "mpileup is done running. Beginning variant calling..."


#variant calling on bcf files
        for file in ${bcf_file}*.bcf
        do
        base=$(basename $file .bcf)
        bcftools call --threads 16 $file -mv -O b -f GQ -o ${bcf_file}${base}_tlr_VariantCalls.bcf &    
        done
        wait
        echo "Variant calling is complete. Preparing files for filtering..."


#prepare files for filtering
        for file in ${bcf_file}*Calls.bcf
        do
                base=$(basename $file _VariantCalls.bcf)
                #reheader each chunked bcf so it has the same sample names
                bcftools reheader -s ${valbamtlrdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf  
                wait
                #put bcf files names into a list for concatenation
                ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
        done
        echo "Preparing files complete. Concatenating chunked bcf files"

#concatenate the chunked bcf files
        bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O b -o ${bcf_file}${species}_VariantCalls_concat.bcf --threads 16
        echo "bcf file is ready for filtering!"