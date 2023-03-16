#!/bin/bash -e

#14 Mar 2023
#Molly Magid, adapted by Olivia Janes
## note to oj: find this script in Molly's github https://github.com/UC-ConSERT/Magid_et_al/blob/main/4_haplotype_creation.sh
##      or in manuscript (different)
## I'm pretty sure this is done one tlr at a time, so ive done tlr7 for a start


sppdir=~/data/tuturuatu_trial_2/
finaldir=${sppdir}bcf_final/


#phase each haplotype into a filtered vcf file with Beagle v. 5.2, beagle.21Apr21.304.jar
beagle gt=${finaldir}Tuturuatu_VariantCalls_final_variants.bcf.gz out=tuturuatu_phased.bcf.gz
##OJ if using this, change all below inputs to above phased output


#create consensus sequences with all ALT variants in a population using samtools faidx
    #this will create a file for comparison to ref genome sequence
samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 | bcftools consensus ${finaldir}tuturuatu_phased.bcf.gz -M N >> ${finaldir}TLR_variants.fasta

#creates TLR haplotypes for all individuals and output sequences into a MEGA file format to use in DNAsp
counter=1 #counter adds unique number labels to each individual for MEGA format
echo "#MEGA" > ${finaldir}tlr_haplotypes.meg

#for loop to output haplotypes for each individual bam file within the final vcf for one TLR region
for file in ${sppdir}nodup_bam/*_nodup.bam #loop through the population bam files
do
    base=$(basename $file _nodup.bam)
    echo "#"${counter}"-1" >> ${finaldir}tlr_haplotypes.meg
        #outputs sequence for phased haplotype 1 for one individual
    samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 | bcftools consensus ${finaldir}tuturuatu_phased.vcf.gz -s $base -M N -H 1pIu >> ${finaldir}tlr_haplotypes.meg
    echo "#"${counter}"-2" >> ${finaldir}tlr_haplotypes.meg
        #outputs sequence for phased haplotype 2 for one individual
    samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 | bcftools consensus ${finaldir}tuturuatu_phased.vcf.gz -s $base -M N -H 2pIu >> ${finaldir}tlr_haplotypes.meg
    counter=$(($counter+1))
done