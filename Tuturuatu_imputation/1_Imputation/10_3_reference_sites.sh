#!/bin/bash -e
set -e

#25 April 2023
#Olivia Janes
#Extracting the TLR regions from the reference genome into a fasta file for input into Geneius

#Environment: N/A

sppdir=~/data/tuturuatu_all_vcf/
ref=~/data/tuturuatu_all/ref_genome/Maui_merged_assembly.fa

#Setting up directories for use
mkdir -p ${sppdir}final_outputs/ref_fasta
outdir=${sppdir}final_outputs/ref_fasta/

#Defining the TLR file list createc in 10_0_haplotyping or 10_2_tlr_snp_counts
tlrlist=${sppdir}final_outputs/tlr_list.txt

#Creating fasta files for each TLR
    while read tlr_pos tlr_name
        #Cycle through TLRs listed in ${tlrlist}
    do
        echo ""; echo "Working on ${tlr_name}"
        samtools faidx ${ref} ${tlr_pos} > ${outdir}${tlr_name}.fasta
    done < ${tlrlist}

echo ""; echo "Script is complete."