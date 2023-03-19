#!/bin/bash -e 

#20 July 2022

#Tuturuatu alignment stats script adapted from https://github.com/janawold1/PhD_thesis/blob/main/Chapter_4/3_align_stats.md
#Jana Wold, adapted by Olivia Janes

sppdir=~/data/tuturuatu_all/

mkdir -p ${sppdir}nodup_bam/ ${sppdir}nodup_bam_stats/
    #Making directories to hold the output bam files and the stats generated.
mkdir -p ${sppdir}nodup_bam/intermediate_files/

mergedbamdir=${sppdir}merged_bam_files/
                #directory that holds the aligned, sorted and merged bam files
nodupbamdir=${sppdir}nodup_bam/

#Aligned bam files are sorted, mates fixed and PCR duplicates removed.
for bam in ${mergedbamdir}*_merged.bam
do
    base=$(basename ${bam} _merged.bam)
    echo "Now preparing to mark duplicates for ${base}..."
    #Sort alignment file based on the names of reads
    samtools sort -@ 16 -n -o ${nodupbamdir}intermediate_files/${base}_sorted.bam ${bam}
    #Fixmate to fill in mate coordinates and insert size fields
    samtools fixmate -@ 16 -r -m -c ${nodupbamdir}intermediate_files/${base}_sorted.bam \
        ${nodupbamdir}intermediate_files/${base}_fixmate.bam
    #Sort based on chromosome number and coordinates
    samtools sort -@ 16 -o ${nodupbamdir}intermediate_files/${base}_fixmate_sorted.bam \
        ${nodupbamdir}intermediate_files/${base}_fixmate.bam
    #Remove duplicate alignments and print basic stats (-s flag)
    samtools markdup -@ 16 -r -s ${nodupbamdir}intermediate_files/${base}_fixmate_sorted.bam \
        ${sppdir}nodup_bam/${base}_nodup.bam
    samtools index -@ 16 -b ${sppdir}nodup_bam/${base}_nodup.bam
done  
echo "Sorting, fixing, and removing duplicates has finished. Time to run stats, stat!"


#Running stats to observe the outcome of fixed bam files.
for bam in ${sppdir}nodup_bam/*_nodup.bam
    do
    base=$(basename ${bam} _nodup.bam)
    echo "Running samtools stats for ${base}"
    samtools stats ${bam} > ${sppdir}nodup_bam_stats/${base}_nodup.stats
    echo "Running Qualimap for ${base}..."
    qualimap bamqc \
        -bam ${bam} \
        -nw 10000 \
        -nt 16 -c \
        -outdir ${sppdir}nodup_bam_stats/${base}_nodup.graphmap \
        --java-mem-size=8G
    echo "Running calculating stats for ${base}..."
    mosdepth --threads 24 --fast-mode --by 50 ${sppdir}nodup_bam_stats/${base}_nodup ${bam}
done
echo "Stats have finished."

#Plotting mosdepth outputs.
python ~/data/general_scripts/plot-dist.py ${sppdir}nodup_bam_stats/*.global.dist.txt > ${sppdir}nodup_bam_stats/nodup_bam_plot_dist_tutu_all.html

echo "Done plotting with python. Sssscript hasss finished....."


###     Outputs     ###
# _nodup.stats                      from samtools stats
# _nodup.graphmap                   from qualimap
# _nodup.mosdepth.summary.txt       from mosdepth
# .per-base.bed.gz                  from mosdepth?
# .regions.bed.gz                   from mosdepth?
# scripts/dist.html                 from python plot-dist.py
