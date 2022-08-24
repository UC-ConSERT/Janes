#!/bin/bash -e 

################################
#Not sure if bcf/vcf formats are right - what format should files be in at each stage?


#22 Aug 2022
#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu variant calling and preparing files for filtering adapted from the end of bwa_alignment_tara_iti_oj.sh

sppdir=~/data/tuturuatu/

ref=${sppdir}ref_genome/Maui_merged_assembly.fa
         #reference genome for alignment
         ##### Must be edited to be sample specific #####
nodupbamdir=${sppdir}nodup_bam/
        #directory that holds the merged bam files that have been sorted, fixed and had duplicates removed.
chunksdir=${sppdir}chunks/
        #a directory to hold the chunked bam files
bcf_file=${sppdir}bcf/
        #bcf file output
species="Tuturuatu"


#chunk bam files for mpileup
ls ${nodupbamdir}*_nodup.bam > ${nodupbamdir}${species}_bam_list.txt
perl ~/data/general_scripts/split_bamfiles_tasks.pl \
        -b ${nodupbamdir}${species}_bam_list.txt \
        -g $ref -n 16 -o ${chunksdir} | parallel -j 16 {}

#run mpileup on chunks of bam files
for ((i=1; i<=16; i++)); do
        bcftools mpileup \
                --threads 16 \
                -f $ref \
                -a AD,ADF,ADR,DP,SP \
                -O b -o ${bcf_file}${species}_${i}_raw.bcf \
                ${chunksdir}${i}/* &
done
wait
echo “mpileup is done running. Beginning variant calling...”


#variant calling on bcf files
for file in ${bcf_file}*.bcf
do
    base=$(basename $file .bcf)
    bcftools call --threads 16 $file -mv -O v -f GQ -o ${bcf_file}${base}_VariantCalls.vcf &    
                                                ##Molly's Tara iti script has output as -O b / .bcf BUT github Molly and Jana is vcf#######
done
wait
echo “Variant calling is complete. Preparing files for filtering...”


#prepare files for filtering
for file in ${bcf_file}*.vcf
do
base=$(basename $file _VariantCalls.vcf)
#reheader each chunked bcf so it has the same sample names
bcftools reheader -s ${nodupbamdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf  
                        #####I think the output of this (*_reheader.bcf) are actually .vcf, based on checking file type with "file *_reheader.bcf"##############
wait
#put bcf files names into a list for concatenation
ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
done
echo "Preparing files complete. Concatenating chunked bcf files"

#concatenate the chunked vcf files
bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O b -o ${bcf_file}${species}_VariantCalls_concat.bcf --threads 16
echo “bcf file is ready for filtering!”