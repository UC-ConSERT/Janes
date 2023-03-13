#!/bin/bash -e 

#26 Jan 2023
## For CAPTIVE samples from 2019 (IKMB) and 2021 (LIC) that need to be merged before variant calling with the rest. ##
## TRIAL 2 as for some reason made one file for each indv with 2 samples in them, rather than one file with one merged sample

#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu alignment from bwa_alignment_tara_iti_oj.sh


sppdir=~/data/tuturuatu_all/

mkdir -p ${sppdir}merged_bam_files_trial_2/

datadir=${sppdir}to_merge_bam/
         #directory with trimmed fastq data
mergedbamdir=${sppdir}merged_bam_files_trial_2/
                #directory that holds the aligned, sorted and merged bam files
species="Tuturuatu"

<<"COMMENTS"

#rename files to remove unnecessary text
        ##### Must be edited to be sample specific #####
for sample in ${datadir}*_merged.bam
do 
        echo $sample
        base=$(basename $sample .bam)
        echo $base
        name=$(echo $base | sed 's/_merged//g')
        echo $name
        rename "s/${base}/${name}/g" ${datadir}/${base}* 
done

#Here I manually changed each of the I164... file names to their corresponding CT.. file name plus "d" for duplicate

COMMENTS

#Merging two samples over two lanes of the same individual (L001 & L002).
        ######### Must be edited to be sample specific ######
for file in ${datadir}*m.bam
do
        base=$(basename $file m.bam) 
        echo "Merging file $base"
        samtools merge -@ 32 ${mergedbamdir}${base}_merged.bam \
                ${datadir}${base}.bam \
                ${datadir}${base}m.bam
done
echo "Merging is complete"

#Indexing the merged bam file
for file in ${mergedbamdir}*_merged.bam
do
        base=$(basename $file _merged.bam)
        echo "Indexing merged bam file $base"
        samtools index -@ 16 -b ${mergedbamdir}${base}_merged.bam
done

echo "Indexing merged bam files is complete"

