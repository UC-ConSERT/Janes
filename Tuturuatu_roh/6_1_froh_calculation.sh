#!/bin/bash -e 
set -e

#12 May 2023
#Olivia Janes
#Calculating genome size of contigs over 300kb and 1000kb for FROH.

#Environment: NA

sppdir=~/data/tuturuatu_roh/
    ## Must be edited to be run specific
parameters=5x_0.1miss

#Defining directories
    ref=${sppdir}ref_genome/Shore_Golden_Plover_Pseudogenome.fasta
    rohdir=${sppdir}roh/
    contigdir=${rohdir}ref_genome_contig_size/

#Calculating length of all contigs.
    bcftools view -h ${sppdir}bcf/Tuturuatu_VariantCalls_concat.vcf.gz | grep '^##contig' | cut -f 1-3 | sed 's/##contig=<ID=//' | sed 's/,.*length=/\t/' | sed 's/>//' | sort | uniq > ${contigdir}contig_size_all.txt

#Calculating length of contigs over 300kb and 1000kb.
    contigs_file="${contigdir}contig_size_all.txt"

    #For 300kb
        size_threshold_300=300000
        sum_300=0

        while IFS=$'\t' read -r contig size; do
            if (( size > size_threshold_300 )); then
                sum_300=$(( sum_300 + size ))
            fi
        done < "$contigs_file"

        echo "Sum of contig sizes over 300kb: $sum_300"

    #For 1000kb
        size_threshold_1000=1000000
        sum_1000=0

        while IFS=$'\t' read -r contig size; do
            if (( size > size_threshold_1000 )); then
                sum_1000=$(( sum_1000 + size ))
            fi
        done < "$contigs_file"
        
        echo "Sum of contig sizes over 1000kb: $sum_1000"

#Calculate individual FROH
    for i in {300,1000}
    do
        file=${rohdir}tuturuatu_roh_${parameters}_${i}kb_window.hom.indiv
        echo "Calculating FROH for individuals in tuturuatu_roh_${parameters}_${i}kb_window.hom.indiv"
        sum_var="sum_$i"
        sum="${!sum_var}"
        echo "$sum"
        echo "ID roh_kb roh_size contigs_total_size FROH" > ${contigdir}contig_size_${parameters}_${i}kb.txt
        awk -v sum="$sum" '{print $1, $5, $5 * 1000, sum, ($5 * 1000) / sum}' ${file} >> ${contigdir}contig_size_${parameters}_${i}kb.txt
    done

echo "Extract with:"; echo "rsync -rav rccuser:${contigdir} ./"
echo "FROH: Functionality Reached, Operation Halted. Final Results Obtained, Huzzah!"