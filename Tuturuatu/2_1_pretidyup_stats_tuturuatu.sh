#!/bin/bash -e 

#21 July 2022
#Tuturuatu alignment stats script adapted from https://github.com/janawold1/PhD_thesis/blob/main/Chapter_4/3_align_stats.md
#To be used before running tidyup on the merged bam, and compared to the stats created after tidyup.
#Jana Wold, adapted by Olivia Janes

#####   Must be edited to be run specific   #####
sppdir=~/data/tuturuatu/

mergedbamdir=${sppdir}merged_bam_files/
    #directory that holds the aligned, sorted and merged bam files

mkdir -p ${sppdir}dup_bam_stats/
    #Making directories to hold the merged bam stats generated.

#Samtools stats are run for merged/aligned bam files..
for bam in ${mergedbamdir}*_merged.bam
do
    base=$(basename ${bam} _merged.bam)
    echo "Now preparing to create stats for ${base}..."
    samtools stats ${bam} > ${sppdir}dup_bam_stats/${base}_merged_dup.stats
done
echo "Samtools stats has finished. Time to run stats, stat!"

#Running mosdepth and qualimap stats to observe the merged bam files.
for bam in ${mergedbamdir}*_merged.bam
    do
    base=$(basename ${bam} _merged.bam)
    echo "Running Qualimap for ${base}..."
    qualimap bamqc \
        -bam ${bam} \
        -nw 10000 \
        -nt 16 -c \
        -outdir ${sppdir}dup_bam_stats/${base}_merged_dup.graphmap \
        --java-mem-size=8G
    echo "Running calculating stats for ${base}..."
    mosdepth --threads 24 --fast-mode --by 50 ${sppdir}dup_bam_stats/${base} ${bam}
done
echo "Stats have finished."

#Plotting mosdepth outputs.
python ~/data/general_scripts/plot-dist.py ${sppdir}dup_bam_stats/*.global.dist.txt

echo "Script has finished."