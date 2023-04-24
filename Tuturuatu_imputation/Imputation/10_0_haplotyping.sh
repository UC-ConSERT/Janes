#!/bin/bash -e
set -e

#25 April 2023
#Olivia Janes, adapted from Molly Magid
## note to oj: find this script in Molly's github https://github.com/UC-ConSERT/Magid_et_al/blob/main/4_haplotype_creation.sh
##      or in manuscript (different)
## I'm pretty sure this is done one tlr at a time, so ive done tlr7 for a start


sppdir=~/data/tuturuatu_all_vcf/
bamdir=~/data/tuturuatu_all/nodup_bam/
    #Directory containing all individuals. To extract indv names from.
finaldir=${sppdir}impute/beagle_imputations/filtered/miss_filter_trial/
    #Directory holding the final phased, imputed and filtered vcf.
run=low_cov

#Setting up directories for use
mkdir -p ${sppdir}final_ouputs/fasta/ ${sppdir}final_ouputs/meg/
outdir=${sppdir}final_ouputs/


#Moving files out of the nodup_bam/not_var_calling/ folder (as sometimes files moved by previous scripts)
    if [ -n "$(ls -A ${bamdir}not_var_calling)" ]; then
        echo "Moving files from ${bamdir}not_var_calling to another folder"
        mv ${bamdir}not_var_calling/* ${bamdir}
    else
        echo "There are no files in not_var_calling folder"
    fi


for i in {0.9,0.8,0.5,0}
do
    echo ""; echo "Outputting files for ${i} missingness"
    
    #Defining output files
        fasta=${outdir}fasta/TLR_variants_${i}miss.fasta
        meg=${outdir}meg/tlr_haplotypes_${i}miss.meg

    #create consensus sequences with all ALT variants in a population using samtools faidx
        #this will create a file for comparison to ref genome sequence
        echo ""; echo "Outputting fasta file for ${i} missingness"
        samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 \
            | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -M N >> ${fasta}

    #creates TLR haplotypes for all individuals and output sequences into a MEGA file format to use in DNAsp
        counter=1 #counter adds unique number labels to each individual for MEGA format
        echo "#MEGA" > ${meg}

    #for loop to output haplotypes for each individual bam file within the final vcf for one TLR region
        echo ""; echo "Outputting meg files for ${i} missingness"
        for file in ${sppdir}nodup_bam/*_nodup.bam #loop through the population bam files
        do
            base=$(basename $file _nodup.bam)
            echo "#"${counter}"-1" >> ${meg}

                #outputs sequence for phased haplotype 1 for one individual
            samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 \
                | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -s $base -M N -H 1pIu >> ${meg}

            echo "#"${counter}"-2" >> ${meg}

                #outputs sequence for phased haplotype 2 for one individual
            samtools faidx ${sppdir}ref_genome/Maui_merged_assembly.fa jcf7180002696332:371554-374695 \
                | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -s $base -M N -H 2pIu >> ${meg}
            counter=$(($counter+1))
        done

done

echo ""; echo "Script is complete."