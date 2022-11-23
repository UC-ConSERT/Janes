#!/bin/bash -e 

#12 June 2022
#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu alignment from bwa_alignment_tara_iti_oj.sh

sppdir=~/data/tuturuatu_trial_2/

mkdir -p ${sppdir}ref_genome/ ${sppdir}fq_gz_files/ ${sppdir}sam_files/ ${sppdir}bam_files/ \
        ${sppdir}processed_bam_files/ ${sppdir}bcf/ ${sppdir}chunks/ ${sppdir}to_merge/

ref=${sppdir}ref_genome/Maui_merged_assembly.fa
         #reference genome for alignment
         ##### Must be edited to be sample specific #####
datadir=${sppdir}fq_gz_files/
         #directory with trimmed fastq data
samdir=${sppdir}sam_files/
         #sam file directory
bamdir=${sppdir}bam_files/
         #bam file directory
processedbamdir=${sppdir}processed_bam_files/ 
        #a directory for processed bam files to protect the originals
mergedbamdir=${sppdir}merged_bam_files/
        #directory that holds the aligned, sorted and merged bam files
fq1=R1.fq.gz
        #Read 1 suffix
        ##### Must be edited to be sample specific #####
fq2=R2.fq.gz
        #Read 2 suffix
        ##### Must be edited to be sample specific #####
platform="Illumina"
species="Tuturuatu"

#follow 2_0_0_renaming.sh first

#first index the reference genome
#bwa index $ref  (already indexed in tuturuatu_trial_2)

#retrieving read group and instrument information.
echo "Retrieving read information"
for samp in ${datadir}*${fq1}
do

        #extract sample information from fastq name
        base=$(basename $samp _R1.fq.gz)
        infoline=$(zcat ${samp} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo $infoline | cut -d ':' -f2`
        flowcell=`echo $infoline | cut -d ':' -f3`
        lane=`echo $infoline | cut -d ':' -f4`
        index=`echo $infoline | cut -d ':' -f10`
        #the _apr/aug indicates the two files of the same indv to be merged later on.
        name0=$(echo ${base} | sed 's/_apr//g')
        name=$(echo ${name0} | sed 's/_aug//g')

        #now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${name}"

        echo "  $rgid   $rgpl   $rgpu   $rglb   $rgsm" >> ${bamdir}/read_groups.txt

        echo "Aligning reads for $base" #align paired sample reads with ref genome using bwa v 0.7.17
        bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 64 $ref $samp ${datadir}${base}*${fq2} > ${samdir}${base}.sam
        echo "Converting sam file to bam file for $base"
        samtools view -@ 16 -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam
done


#Sorting bam files
for file in ${bamdir}*.bam
do
        base=$(basename $file .bam)
        echo "Sorting file $base"
        samtools sort -@ 16 -o ${processedbamdir}${base}_aligned_sorted.bam ${bamdir}${base}.bam
        rm ${samdir}${base}.sam
done
echo "Sorting bam files is complete"

<<"COMMENTS"

#Some Apr files have no august to be merged with, so removing _apr from their name.
mv ${bamdir}A11_apr_aligned_sorted.bam ${bamdir}A11_aligned_sorted.bam
mv ${bamdir}F09_apr_aligned_sorted.bam ${bamdir}F09_aligned_sorted.bam
mv ${bamdir}F11_apr_aligned_sorted.bam ${bamdir}F11_aligned_sorted.bam
mv ${bamdir}G09_apr_aligned_sorted.bam ${bamdir}G09_aligned_sorted.bam
mv ${bamdir}G11_apr_aligned_sorted.bam ${bamdir}G11_aligned_sorted.bam
mv ${bamdir}H09_apr_aligned_sorted.bam ${bamdir}H09_aligned_sorted.bam


#Not all files are merged, but they all must finish with the same file name (_merged.bam), so separating the ones to merge.
mv *_apr_aligned_sorted.bam ${sppdir}to_merge/
mv *_aug_aligned_sorted.bam ${sppdir}to_merge/

#Merging two samples over two lanes of the same individual (_apr & _aug).
for file in ${sppdir}to_merge/*_apr_aligned_sorted.bam
do
        base=$(basename $file _apr_aligned_sorted.bam) 
        echo "Merging file $base"
        samtools merge -@ 32 ${mergedbamdir}${base}_merged.bam \
                ${sppdir}to_merge/${base}_apr_aligned_sorted.bam \
                ${sppdir}to_merge/${base}_aug_aligned_sorted.bam
done
echo "Merging is complete"

#Moving and renaming the files that don't need to be merged
for file in ${processedbamdir}*_aligned_sorted.bam
do
        base=$(basename $file _aligned_sorted.bam)
        mv ${file} ${base}_merged.bam
        mv ${base}_merged.bam ${mergedbamdir}
done

#Indexing the merged bam file
for file in ${mergedbamdir}*_merged.bam
do
        base=$(basename $file _merged.bam)
        echo "Indexing merged bam file $base"
        samtools index -@ 16 -b ${mergedbamdir}${base}_merged.bam
done
echo "Indexing merged bam files is complete"

COMMENTS

echo "Script is complete."