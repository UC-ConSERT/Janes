#!/bin/bash -e
set -e

#25 April 2023
#Olivia Janes, adapted from Molly Magid
## note to oj: find this script in Molly's github https://github.com/UC-ConSERT/Magid_et_al/blob/main/4_haplotype_creation.sh
##      or in manuscript (different)
## I'm pretty sure this is done one tlr at a time, so ive done tlr7 for a start

#Environment: samtools

sppdir=~/data/tuturuatu_all_vcf/
bamdir=~/data/tuturuatu_all/nodup_bam/
    #Directory containing all individuals. To extract indv names from.
ref_genome=~/data/tuturuatu_all/ref_genome/Maui_merged_assembly.fa
finaldir=${sppdir}impute/beagle_imputations/filtered/miss_filter_trial/
    #Directory holding the final phased, imputed and filtered vcf
run=low_cov

#Setting up directories for use
mkdir -p ${sppdir}final_outputs/fasta/ ${sppdir}final_outputs/meg/
outdir=${sppdir}final_outputs/

#Creating the TLR file list
tlrlist=${sppdir}final_outputs/tlr_list.txt
echo "jcf7180002669510:1-2191 TLR1A" > ${tlrlist}
echo "jcf7180002696225:79300-80387 TLR1B" >> ${tlrlist}
echo "jcf7180002687310:7857-10285 TLR2A" >> ${tlrlist}
echo "jcf7180002686685:232999-235155 TLR2B" >> ${tlrlist}
echo "jcf7180002693511:279461-279884 TLR3_1" >> ${tlrlist}
echo "jcf7180002693511:284293-286123 TLR3_2" >> ${tlrlist}
echo "jcf7180002693511:286306-286420 TLR3_3" >> ${tlrlist}
echo "jcf7180002688524:101290-103549 TLR4" >> ${tlrlist}
echo "jcf7180002688267:101472-103857 TLR5" >> ${tlrlist}
echo "jcf7180002696332:371554-374695 TLR7" >> ${tlrlist}
echo "jcf7180002694481:2000-3649 TLR21" >> ${tlrlist}
echo "jcf7180002693589:0-2103 TLR15_part" >> ${tlrlist}


#Moving files out of the nodup_bam/not_var_calling/ folder (as sometimes files moved by previous scripts)
    if [ -n "$(ls -A ${bamdir}not_var_calling)" ]; then
        echo "Moving files from ${bamdir}not_var_calling to another folder"
        mv ${bamdir}not_var_calling/* ${bamdir}
    else
        echo "There are no files in not_var_calling folder"
    fi

#Creating output files
while read tlr_pos tlr_name
    #Cycle through TLRs listed in ${tlrlist}
do
    echo ""; echo "Working on ${tlr_name}"
    
    
    for i in {0.9,0.8,0.5,0}
    do
        echo ""; echo "Outputting files for ${i} missingness"
        
        #Defining output files
            fasta=${outdir}fasta/TLR_variants_${i}miss_${tlr_name}.fasta
            meg=${outdir}meg/tlr_haplotypes_${i}miss_${tlr_name}.meg

        #create consensus sequences with all ALT variants in a population using samtools faidx
            #this will create a file for comparison to ref genome sequence
            echo ""; echo "Outputting fasta file for ${i} missingness"
            samtools faidx ${ref_genome} ${tlr_pos} \
                | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -M N >> ${fasta}

        #creates TLR haplotypes for all individuals and output sequences into a MEGA file format to use in DNAsp
            counter=1 #counter adds unique number labels to each individual for MEGA format
            echo "#MEGA" > ${meg}

        #for loop to output haplotypes for each individual bam file within the final vcf for one TLR region
            echo ""; echo "Outputting meg files for ${i} missingness"
            for file in ${bamdir}*_nodup.bam #loop through the population bam files
            do
                base=$(basename ${file} _nodup.bam)
                echo "#"${counter}"-1" >> ${meg}

                    #outputs sequence for phased haplotype 1 for one individual
                samtools faidx ${ref_genome} ${tlr_pos} \
                    | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -s ${file} -M N -H 1pIu >> ${meg}

                echo "#"${counter}"-2" >> ${meg}

                    #outputs sequence for phased haplotype 2 for one individual
                samtools faidx ${ref_genome} ${tlr_pos} \
                    | bcftools consensus ${finaldir}Tuturuatu_VariantCalls_${i}miss_${run}_final.vcf.gz -s ${file} -M N -H 2pIu >> ${meg}
                counter=$(($counter+1))
            done
    done
done < ${tlrlist}

echo ""; echo "Script is complete."