#!/bin/bash -e 

#16 Mar 2023
## For CAPTIVE samples from 2019 (IKMB) and 2021 (LIC) that need to be merged before variant calling with the rest. ##
## Using 2_1_reheader.sh, the 2019 (IKMB) files have had their headers (SM tags) changed to match the 2021 (LIC) files.

#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu alignment from bwa_alignment_tara_iti_oj.sh


sppdir=~/data/tuturuatu_roh/

datadir=${sppdir}to_merge/
         #directory with bam files to be merged
mergedbamdir=${sppdir}merged_bam_files/
                #directory that holds the aligned, sorted and merged bam files
species="Tuturuatu"

## Files to merge are named CT01.bam and CT01mr.bam, for example. 
##      m = to merge, r = reheadered.


#Merging two samples over two lanes of the same individual (L001 & L002).
        ######### Must be edited to be sample specific ######
for file in ${datadir}*mr.bam
do
        base=$(basename $file mr.bam) 
        echo "Merging file $base"
        samtools merge -@ 32 ${mergedbamdir}${base}_merged.bam \
                ${datadir}${base}_merged.bam \
                ${datadir}${base}mr.bam
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

