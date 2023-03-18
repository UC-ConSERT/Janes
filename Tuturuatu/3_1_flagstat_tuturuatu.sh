#!/bin/bash -e 

# 18 March 2023
# Script to check before and after duplicate removal bams.


sppdir=~/data/tuturuatu/

mergedbamdir=${sppdir}merged_bam_files/
    #directory that holds the aligned, sorted and merged bam files


echo "Running flagstat on dup bams"
for bam in ${mergedbamdir}*_merged.bam
do
    samtools flagstat $bam >> ${sppdir}dup_bam_stats/flagstat_dup.txt
done

echo "Running flagstat on nodup bams"
for bam in ${sppdir}nodup_bam/*_nodup.bam
do
    samtools flagstat $bam >> ${sppdir}nodup_bam_stats/flagstat_nodup.txt
done

echo ""
echo "Script complete."