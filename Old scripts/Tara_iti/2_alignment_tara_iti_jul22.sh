#!/bin/bash -e 

#26 July 2022
#Olivia Janes adapted from Molly Magid and Jana Wold
#Tara iti and Aus Fairy Tern alignment from bwa_alignment_tara_iti_oj.sh using common tern genome
sppdir=~/data/tara_iti_jul22/

mkdir -p ${sppdir}ref_genome/ ${sppdir}fq_gz_files/ ${sppdir}sam_files/ ${sppdir}bam_files/ \
        ${sppdir}processed_bam_files/ ${sppdir}bcf/ ${sppdir}chunks/

ref=~/data/tara_iti_jul22/ref_genome/bSteHir1.pri.cur.20190820.fasta
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
species="Tara iti"

#first index the reference genome
bwa index $ref

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
        name=$(echo ${base} | sed 's/_lib[1-9]//g')  ####Must be edited to be sample specific

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
        samtools sort -@ 16 -o ${processedbamdir}${base}.aligned.sorted.bam ${bamdir}${base}.bam
        rm ${samdir}${base}.sam
done
echo "Sorting bam files is complete"


#Merging two samples over two lanes of the same individual (lib1 & lib2).
        ######### Must be edited to be sample specific ######
for file in ${processedbamdir}*_lib1.aligned.sorted.bam
do
        base=$(basename $file _lib1.aligned.sorted.bam) 
        echo "Merging file $base"
        samtools merge -@ 32 ${mergedbamdir}${base}_merged.bam \
                ${processedbamdir}${base}_lib1.aligned.sorted.bam \
                ${processedbamdir}${base}_lib2.aligned.sorted.bam
done
echo "merging is complete"

#Indexing the merged bam file
for file in ${mergedbamdir}*_merged.bam
do
        base=$(basename $file _merged.bam)
        echo "Indexing merged bam file $base"
        samtools index -@ 16 -b ${mergedbamdir}${base}_merged.bam
done
echo "Indexing merged bam files is complete"

