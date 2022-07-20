#!/bin/bash -e 

#12 June 2022
#Olivia Janes, adapted from Molly Magid
#Tuturuatu alignment from bwa_alignment_tara_iti_oj.sh
sppdir=~/data/tuturuatu/

mkdir -p ${sppdir}ref_genome_tutu/ ${sppdir}fq_gz_files/ ${sppdir}sam_files/ ${sppdir}bam_files/ \
        ${sppdir}processed_bam_files/ ${sppdir}bcf/ ${sppdir}chunks/

ref=~/data/tuturuatu/ref_genome_tutu/Maui_merged_assembly.fa
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
chunksdir=${sppdir}chunks/
        #a directory to hold the chunked bam files
bcf_file=${sppdir}bcf/ 
        #bcf file output
fq1=R1.fq.gz
        #Read 1 suffix
        ##### Must be edited to be sample specific #####
fq2=R2.fq.gz
        #Read 2 suffix
        ##### Must be edited to be sample specific #####
platform="Illumina"
species="Tuturuatu"

<<'COMMENTS'

#rename files to remove unnecessary text
        ##### Must be edited to be sample specific #####
for sample in ${datadir}*_L001_R1.fq.gz
do 
        echo $sample
        base=$(basename $sample _L001_R1.fq.gz)
        echo $base
        name=$(echo $base | sed 's/_S[0-9][0-9]/_S1/g')
        name=$(echo $base | sed 's/_S[0-9]//g')
        echo $name
        rename "s/${base}/${name}/g" ${datadir}/${base}* 
done

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
        name=$(echo ${base} | sed 's/_L00[1-9]//g')

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


#Aligning and indexing bam files
for file in ${bamdir}*.bam
do
         echo "Aligning and indexing file"
        base=$(basename $file .bam)
        samtools sort -@ 16 -o ${processedbamdir}${base}.aligned.sorted.bam ${bamdir}${base}.bam
        samtools index -@ 16 -b ${processedbamdir}${base}.aligned.sorted.bam
        rm ${samdir}${base}.sam
done

#Merging two samples over two lanes of the same individual (L001 & L002).
        ######### Must be edited to be sample specific ######
for file in ${processedbamdir}*_L001.aligned.sorted.bam
do
        base=$(basename $file _L001.aligned.sorted.bam) 
        echo "Merging file $base"
        samtools merge -@ 32 ${processedbamdir}${base}_merged.bam \
                ${processedbamdir}${base}_L001.aligned.sorted.bam \
                ${processedbamdir}${base}_L002.aligned.sorted.bam
done
echo "merging is complete"

#QC using qualimap and mosdepth (install)
#Mark duplicates (nodup using Jana Ch 4_3 aligned stats)
#QC using qualmap and mosdepth (install)

#chunk bam files for mpileup
ls ${processedbamdir}*.aligned.sorted.bam > ${processedbamdir}${species}bam_list.txt
perl ~/data/general_scripts/split_bamfiles_tasks.pl -b ${processedbamdir}${species}bam_list.txt -g $ref -n 13 -o ${chunksdir} | parallel -j 12 {}

#run mpileup on chunks of bam files
for ((i=1; i<=12; i++)); do
        bcftools mpileup -O b -f $ref -a AD,ADF,ADR,DP,SP -o ${bcf_file}${species}_${i}_raw.bcf  ${chunksdir}${i}/* &
done
wait
echo “mpileup is done running”

COMMENTS

#variant calling on bcf files
#for file in ${bcf_file}*.bcf
#do
#base=$(basename $file .bcf)
#bcftools call $file -mv -O b -f GQ -o ${bcf_file}${base}_VariantCalls.bcf &
#done
#wait
#echo “variant calling is complete”

#Molly doesn't have either of the below 2 lines in her script. it appears the list of bcf is created below on line 131
#>${bcf_file}list_of_bcf.txt   #do i need this part? I think this creates the list_of_bcf.txt file
#ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt #And this writes all of the variant calls names to it. However,
        #this may be done below after reheader. As in Molly's script. Final verdict: I THINK remove this part.


#prepare files for filtering with bgzip and indexing
#for file in ${bcf_file}*.vcf
#do
#base=$(basename $file _raw_VariantCalls.vcf)
#bcftools reheader -s ${bamdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
#wait
#>list_of_bcf.txt
#ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt
#done

#concatenate the chunked vcf files
#bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O v -o ${bcf_file}${species}_VariantCalls_concat.vcf --threads 16
#echo “bcf file is ready for filtering!”
