#!/bin/bash -e 

################################
#Need to change all merged/processed bam dir to the nodupbamdir
##################################


#12 June 2022
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
perl ~/data/general_scripts/split_bamfiles_tasks.pl -b ${nodupbamdir}${species}_bam_list.txt -g $ref -n 13 -o ${chunksdir} | parallel -j 12 {}

#run mpileup on chunks of bam files
for ((i=1; i<=12; i++)); do
        bcftools mpileup -O b -f $ref -a AD,ADF,ADR,DP,SP -o ${bcf_file}${species}_${i}_raw.bcf  ${chunksdir}${i}/* &
done
wait
echo “mpileup is done running”


#variant calling on bcf files
for file in ${bcf_file}*.bcf
do
    base=$(basename $file .bcf)
    bcftools call $file -mv -O b -f GQ -o ${bcf_file}${base}_VariantCalls.vcf &
done
wait
echo “variant calling is complete”


<<"COMMENTS"
#Molly doesn't have either of the below 2 lines in her script. it appears the list of bcf is created below on line 131
#>${bcf_file}list_of_bcf.txt   #do i need this part? I think this creates the list_of_bcf.txt file
#ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt #And this writes all of the variant calls names to it. However,
        #this may be done below after reheader. As in Molly's script. Final verdict: I THINK remove this part.


#prepare files for filtering with bgzip and indexing
for file in ${bcf_file}*.vcf
do
    base=$(basename $file _raw_VariantCalls.vcf)  #will it have raw_? Apparently not according to above output, but double check.
    bcftools reheader -s ${nodupbamdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
    wait
    >list_of_bcf.txt
    ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt #thinkng that the .bcf -> .vcf?
done

####THIS IS A REPEAT of above paragraph, just from Molly's github. 
#### I have already changed the bcftools call output from .bcf to .vcf, so may have to edit in above text (should match below)
#prepare files for filtering
for file in ${bcf_file}*.vcf
do
base=$(basename $file .vcf)
#reheader each chunked bcf so it has the same sample names
bcftools reheader -s ${bamdir}OFK_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
wait
#put bcf files names into a list for concatenation
ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
done

COMMENTS

##### What I have come up with as a combo from my paragraph and Molly's github.
#prepare files for filtering
for file in ${bcf_file}*.vcf
do
base=$(basename $file _VariantCalls.vcf)
#reheader each chunked bcf so it has the same sample names
bcftools reheader -s ${nodupbamdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
wait
#put bcf files names into a list for concatenation
ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
done

#concatenate the chunked vcf files
#bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O v -o ${bcf_file}${species}_VariantCalls_concat.vcf --threads 16
#echo “bcf file is ready for filtering!”