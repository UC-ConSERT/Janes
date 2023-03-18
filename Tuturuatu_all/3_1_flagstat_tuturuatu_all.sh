#!/bin/bash -e 

# 18 March 2023
# Script to check before and after duplicate removal bams.


sppdir=~/data/tuturuatu_all/

mergedbamdir=${sppdir}merged_bam_files/
    #directory that holds the aligned, sorted and merged bam files


echo "Running flagstat on dup bams"
for bam in ${mergedbamdir}*_merged.bam
do
    echo "Running flagstat on ${bam}"
    base=$(basename ${bam} _merged.bam)
    echo "##    ${base} dup bam" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
        echo "" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
    samtools flagstat $bam >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
done

echo "Running flagstat on nodup bams"
for bam in ${sppdir}nodup_bam/*_nodup.bam
do
    echo "Running flagstat on ${bam}"
    base=$(basename ${bam} _nodup.bam)
    echo "" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
        echo "" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
        echo "##    ${base} NOdup bam" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
        echo "" >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
    samtools flagstat $bam >> ${sppdir}nodup_bam_stats/${base}_flagstat.txt
done

echo ""
echo "Script complete."