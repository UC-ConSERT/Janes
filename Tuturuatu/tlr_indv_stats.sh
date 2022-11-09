#!/bin/bash -e

#9 Nov 2022
#Olivia Janes
#Creating individual-level stats for TLR regions within duplicate removed (nodup.bam) bam files.
#This will be used to analyse the quality of reads within the TLR regions in tuturuatu individuals, to assess if poor overall
#   genome quality correlates to poor quality within the TLR regions.
#The output file can be run through multiqc for a visual report.

sppdir=~/data/tuturuatu/

nodupbamdir=${sppdir}nodup_bam/
        #directory that holds the merged bam files that have been sorted, fixed and had duplicates removed.
tlr_regions=${sppdir}bcf/tlr_regions.bed

mkdir -p ${nodupbamdir}tlr_indv_stats/

for file in ${nodupbamdir}*nodup.bam
do
    base=$(basename ${file} _nodup.bam)
    echo "Creating TLR region stats for ${base}..."
    samtools stats ${file} -@ 64 -t ${tlr_regions} >> ${nodupbamdir}tlr_indv_stats/tlr_indv_stats.stat
done

echo "Samtools stats has finished. Time to run multiqc!!!!!"
