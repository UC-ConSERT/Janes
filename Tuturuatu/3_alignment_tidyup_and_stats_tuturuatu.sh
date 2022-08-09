#!/bin/bash -e 

#20 July 2022
#Tuturuatu alignment stats script adapted from https://github.com/janawold1/PhD_thesis/blob/main/Chapter_4/3_align_stats.md
#Jana Wold, adapted by Olivia Janes

sppdir=~/data/tuturuatu/

mergedbamdir=${sppdir}merged_bam_files/
                #directory that holds the aligned, sorted and merged bam files

mkdir -p ${sppdir}nodup_bam/ ${sppdir}nodup_bam_stats/
    #Making directories to hold the output bam files and the stats generated.

#Aligned bam files are sorted, mates fixed and PCR duplicates removed.
for bam in ${mergedbamdir}*_merged.bam
do
    base=$(basename ${bam} _merged.bam)
    echo "Now preparing to mark duplicates for ${base}..."
    #Sort alignment file based on the names of reads
    samtools sort -@ 8 -n -o ${sppdir}nodup_bam/${base}.nsorted.bam ${bam}
    #Fixmate to fill in mate coordinates and insert size fields
    samtools fixmate -@ 8 -r -m -c ${sppdir}nodup_bam/${base}.nsorted.bam \
        ${sppdir}nodup_bam/${base}.fixmate.bam
    #Sort based on chromosome number and coordinates
    samtools sort -@ 8 -o ${sppdir}nodup_bam/${base}.fixmate.sorted.bam \
        ${sppdir}nodup_bam/${base}.fixmate.bam
    #Remove duplicate alignments and print basic stats (-s flag)
    samtools markdup -@ 8 -r -s ${sppdir}nodup_bam/${base}.fixmate.sorted.bam \
        ${sppdir}nodup_bam/${base}_nodup.bam
    samtools index -@ 16 -b ${sppdir}nodup_bam/${base}_nodup.bam
    samtools stats ${bam} > ${sppdir}nodup_bam_stats/${base}_nodup.stats #Note: this may have to go in below for loop? Had issues with the nodup file not being found.
done
echo "Sorting, fixing, and duplicating has finished. Time to run stats, stat!"

#Running stats to observe the outcome of fixed bam files.
for bam in ${sppdir}nodup_bam/*_nodup.bam
    do
    base=$(basename ${bam} _nodup.bam)
    echo "Running Qualimap for ${base}..."
    qualimap bamqc \
        -bam ${bam} \
        -nw 10000 \
        -nt 16 -c \
        -outdir ${sppdir}nodup_bam_stats/${base}.graphmap \
        --java-mem-size=8G
    echo "Running calculating stats for ${base}..."
    mosdepth --threads 24 --fast-mode --by 50 ${sppdir}nodup_bam_stats/${base} ${bam}
done
echo "Stats have finished."

#Plotting mosdepth outputs.
python ~/data/general_scripts/plot-dist.py ${sppdir}nodup_bam_stats/*.global.dist.txt