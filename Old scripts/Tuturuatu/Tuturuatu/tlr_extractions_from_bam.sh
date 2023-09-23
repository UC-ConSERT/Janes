#!/bin/bash -e

#8 Nov 2022
#Olivia Janes
#Extracting TLR regions from duplicate removed (nodup.bam) bam files, for further analysis in the Integrative Genomics Viewer (IGV).
#This will be used to analyse the quality of reads within the TLR regions in a subset of tuturuatu individuals, to assess if poor overall
#   genome quality correlates to poor quality within the TLR regions.

sppdir=~/data/tuturuatu/

nodupbamdir=${sppdir}nodup_bam/
        #directory that holds the merged bam files that have been sorted, fixed and had duplicates removed.
tlr_regions=${sppdir}bcf/tlr_regions.bed

mkdir -p ${nodupbamdir}tlr_only_bam/

#Define files for the analysis. I have chosen 5 poor depth individuals, and 4 higher depth individuals.
filelist="CR12_nodup.bam CR11_nodup.bam CR02_nodup.bam CR18_nodup.bam CR10_nodup.bam \
            CR20_nodup.bam CR06_nodup.bam CT04_nodup.bam CR15_nodup.bam"

#Important to change file directory so that for loop can find files in filelist.
cd ${nodupbamdir}

#Extract only the regions defined in the TLR bed file from bam files into a new file and folder.
for file in ${filelist}
do
    base=$(basename ${file} _nodup.bam)
    echo "Extracting TLR region for ${base}..."
    samtools view -L ${tlr_regions} -@ 64 -o ${nodupbamdir}tlr_only_bam/${base}_nodup_tlr.bam ${file}
    echo "Indexing TLR bam file $base"
    samtools index -@ 16 -b ${nodupbamdir}tlr_only_bam/${base}_nodup_tlr.bam
done


echo "Script has finished."

#Upload TLR bam files, TLR bam index files, reference genome and ref genome index file to IGV at https://igv.org/app/.